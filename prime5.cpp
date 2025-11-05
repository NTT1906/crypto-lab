#include <array>
#include <vector>
#include <cstdio>

#define CIGINT_IMPLEMENTATION
#define CIGINT_STRIP_PREFIX
#define CIGINT_N (512 / 32)
#include "cigint.h"

// we use the same small r as your Python
static constexpr int R = 14;

struct Poly {
	std::array<Cigint, R> c;
};

// -------- basic Cigint helpers --------
static inline Cigint z0() { return CIGINT_ZERO(); }
static inline Cigint z1() { return cigint_from_u32(1u); }

#include <iostream>

// anti-overflow add mod
// u32 addmod(u32 a, u32 b, u32 m) {
//     if (a >= m) a %= m;
//     if (b >= m) b %= m;
//     u32 c = a + b;
//     if (c < a) c -= m;
//     if (c >= m) c -= m;
//     return c;
// }

void cigint_add_wrap(Cigint *lhs, const Cigint *rhs) {
	u64 carry = 0;
	for (size_t i = CIGINT_N; i-- > 0;) {
		carry = (u64) lhs->data[i] + (u64) rhs->data[i] + carry;
		lhs->data[i] = (u32) carry;
		carry >>= 32;
	}
}

// res = (a + b) mod m (all fixed width)
static inline void add_mod(Cigint& res, Cigint& b, const Cigint &m) {
	Cigint res_c = res;
	cigint_divmod_ref(&res, &m, NULL, &res);
	cigint_divmod_ref(&b, &m, NULL, &b);
	cigint_add_wrap(&res, &b);
	if (res < res_c) {
		cigint_sub_ref(&res, &m);
	}
	if (cigint_cmp(res, m) >= 0) {
		cigint_sub_ref(&res, &m);
	}
}

// static Cigint mul_mod_safe(Cigint a, Cigint b, const Cigint &m) {
// 	// cigint_divmod_ref(&a, &m, NULL, &a);
// 	// cigint_divmod_ref(&b, &m, NULL, &b);
// 	Cigint res = z0();
// 	u32 bits = cigint_highest_order(b);
// 	for (u32 i = 0; i < bits; ++i) {
// 		if (cigint_get_bit(b, i)) {
// 			add_mod(res, a, m);
// 		}
// 		add_mod(a, a, m);
// 	}
// 	return res;
// }

static inline u32 cigint_add_wrap_carry(Cigint *lhs, const Cigint *rhs) {
    u64 carry = 0;
    for (size_t i = CIGINT_N; i-- > 0;) {
        u64 t = (u64)lhs->data[i] + (u64)rhs->data[i] + carry;
        lhs->data[i] = (u32)t;
        carry = t >> 32;
    }
    return (u32)carry;
}

/* res = (res + b) mod m
 * PRE: 0 <= res < m, 0 <= b < m
 * COST: one limb add + at most one subtract
 */
static inline void add_mod_fast(Cigint &res, const Cigint &b, const Cigint &m) {
    u32 c = cigint_add_wrap_carry(&res, &b);
    if (c || cigint_cmp(res, m) >= 0) cigint_sub_ref(&res, &m);
}

/* x = 2x mod m
 * PRE: 0 <= x < m
 */
static inline void dbl_mod_fast(Cigint &x, const Cigint &m) {
    u32 carry = 0;
    for (size_t i = CIGINT_N; i-- > 0;) {
        u32 v = x.data[i];
        u32 nv = (v << 1) | carry;
        carry = v >> 31;
        x.data[i] = nv;
    }
    if (carry || cigint_cmp(x, m) >= 0) cigint_sub_ref(&x, &m);
}

/* ====== 4-bit window mul-mod (MSW-first), standard-only ======
 * r = (a * b) % m
 * - Reduce a,b once (so they stay < m)
 * - Precompute T[d] = d*a (mod m) for d in {0..15}
 * - Scan b nibble-by-nibble MSW→LSW
 * Complexity: O(N²) with small constants, ~4× fewer adds vs bit-by-bit
 */
static Cigint mul_mod_win4(Cigint a, Cigint b, const Cigint &m) {
    a = cigint_mod(a, m);
    b = cigint_mod(b, m);
    if (cigint_is0(a) || cigint_is0(b)) return CIGINT_ZERO();

    /* table T[0..15] */
    Cigint T[16];
    for (int i = 0; i < 16; ++i) T[i] = CIGINT_ZERO();
    T[1] = a;
    T[2] = a; dbl_mod_fast(T[2], m);            // 2a
    T[4] = T[2]; dbl_mod_fast(T[4], m);         // 4a
    T[8] = T[4]; dbl_mod_fast(T[8], m);         // 8a
    for (int d = 3; d < 16; ++d) {
        if (d == 4 || d == 8) continue;
        Cigint acc = CIGINT_ZERO();
        if (d & 1) add_mod_fast(acc, T[1], m);
        if (d & 2) add_mod_fast(acc, T[2], m);
        if (d & 4) add_mod_fast(acc, T[4], m);
        if (d & 8) add_mod_fast(acc, T[8], m);
        T[d] = acc;               // T[d] < m
    }

    Cigint res = CIGINT_ZERO();
    u32 bits = cigint_highest_order(b);        // 0..(32*CIGINT_N)
    if (bits == 0) return res;
    u32 nibbles = (bits + 3u) / 4u;
    u32 capacity = (u32)CIGINT_N * 8u;         // 8 nibbles per 32-bit limb
    u32 skip = capacity - nibbles;             // leading zero nibbles to skip

    bool first = true;
    u32 seen = 0;

    for (size_t i = 0; i < CIGINT_N; ++i) {
        u32 v = b.data[i];
        for (int k = 0; k < 8; ++k) {          // MSW nibble first
            if (seen++ < skip) continue;
            u32 d = (v >> (28 - 4*k)) & 0xFu;
            if (first) {
                if (d) res = T[d];             // res < m
                first = false;
            } else {
                // multiply by 16 via 4 doublings
                dbl_mod_fast(res, m);
                dbl_mod_fast(res, m);
                dbl_mod_fast(res, m);
                dbl_mod_fast(res, m);
                if (d) add_mod_fast(res, T[d], m);
            }
        }
    }
    if (first) return CIGINT_ZERO();           // b == 0
    return res;                                // < m
}

/* (optional) keep the old name but route to the fast one */
static inline Cigint mul_mod_safe(Cigint a, Cigint b, const Cigint &m) {
    return mul_mod_win4(a, b, m);
}


// poly multiplication mod (x^R - 1, m)
// C[k] = sum_{i+j≡k mod R} A[i]*B[j] (mod m)
static Poly poly_mul(const Poly &A, const Poly &B, const Cigint &m) {
	Poly C;
	for (int i = 0; i < R; ++i) C.c[i] = z0();

	for (int i = 0; i < R; ++i) {
		if (cigint_is0(A.c[i])) continue;
		for (int j = 0; j < R; ++j) {
			if (cigint_is0(B.c[j])) continue;
			int idx = (i + j) % R;
			Cigint prod = mul_mod_safe(A.c[i], B.c[j], m);
			add_mod(C.c[idx], prod, m);
		}
	}
	return C;
}

// compute n mod R using your own bigint division
static u32 cigint_mod_small(const Cigint &n, u32 mod) {
	Cigint q;
	u32 r = 0;
	cigint_sdivmod(n, mod, &q, &r);  // r = n % mod
	return r;
}

// check poly == 1 + x^(n mod R)
static bool poly_is_xn_plus_1(const Poly &P, const Cigint &n, const Cigint &m) {
	u32 k = cigint_mod_small(n, R); // exact n % R
	for (int i = 0; i < R; ++i) {
		if (i == 0 || i == (int)k) {
			if (cigint_cmp(P.c[i], z1()) != 0) return false;
		} else {
			if (!cigint_is0(P.c[i])) return false;
		}
	}
	return true;
}

// compute (1 + x)^n mod (x^R - 1, n)
static Poly poly_pow_1x(const Cigint &n, const Cigint &modn) {
	// base = 1 + x
	Poly base;
	for (int i = 0; i < R; ++i) base.c[i] = z0();
	base.c[0] = z1();
	base.c[1] = z1();

	// res = 1
	Poly res;
	for (int i = 0; i < R; ++i) res.c[i] = z0();
	res.c[0] = z1();

	// exponent = n (big)
	Cigint e = n;
	while (!cigint_is0(e)) {
		if (cigint_get_bit(e, 0)) {
			res = poly_mul(res, base, modn);
		}
		base = poly_mul(base, base, modn);
		e = cigint_shr(e, 1);  // e >>= 1
	}
	return res;
}

static bool aks_like_cigint(const Cigint &n) {
	// reject 0,1
	if (cigint_is0(n)) return false;
	if (cigint_cmp(n, z1()) == 0) return false;

	Poly p = poly_pow_1x(n, n);
	return poly_is_xn_plus_1(p, n, n);
}

#include <random>

void cigint_fill_more_random(Cigint *cig) {
	std::random_device rd;             // Random device
	std::mt19937 gen(rd());            // Mersenne Twister engine
	std::uniform_int_distribution<u32> dist(0, UINT32_MAX);  // Distribution for random u32
	size_t randomLimb = 1 + gen() % CIGINT_N;
	for (u32 &i : cig->data) i = 0;
	// Generate random values for cig->data
	for (size_t i = randomLimb; i-- > 0;) {
		cig->data[i] = dist(gen);  // Using the Mersenne Twister to get random values
	}
}

Cigint get_random_odd() {
	Cigint cig;
	cigint_fill_more_random(&cig);
	cigint_set_bit_ref(&cig, 0, 1);
	return cig;
}

static bool has_small_factor(const Cigint &n) {
	static const int SMALL_PRIMES[] = {
		2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97
		// extend to 1000 if you like
	};
	for (int p : SMALL_PRIMES) {
		u32 r = 0;
		cigint_sdivmod(n, (u32)p, NULL, &r);
		if (r == 0) return cigint_cmp(n, cigint_from_u32((u32)p)) != 0;
	}
	return false;
}

Cigint gen_prime() {
	Cigint cig;
	do {
		cig = get_random_odd();
	} while (has_small_factor(cig) || !aks_like_cigint(cig));
	return cig;
}

int main() {
	// Cigint n = cigint_from_dec("7669373798138436444425917457853162445620122597663674808454541043947311884377"); //up to 76digit
	// Cigint n = cigint_from_dec("771458402252315418489174122874340496560720991135040845033515625723932549727443357205751720461233037696609749135759573352834404574203828395431821751181407"); //up to 154 digit
	// Cigint n = cigint_from_dec("12598680682438752944055149498662865476737841742262510750565369004945449194176119763492007280023309651826310840991041734855800903377926420446403888376677867"); //up to 154 digit

	Cigint n = gen_prime();
	// Cigint a = cigint_from_dec("229911617100");
	// Cigint b = cigint_from_dec("1199165601554601993");
	// Cigint r = cigint_mul_2_mod(a, b, n);
	// std::string r_str = r.toDecStr();
	// printf("r = %s\n", r_str.c_str());
	printf("n = ");
	cigint_print10(n);
	bool ok = aks_like_cigint(n);
	printf(" -> %s\n", ok ? "maybe prime" : "composite-ish");
	return 0;
}
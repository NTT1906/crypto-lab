#include <array>
#include <cstdio>

#define CIGINT_IMPLEMENTATION
#define CIGINT_STRIP_PREFIX
#define CIGINT_N (512 / 32)
#include "../cigint.h"

// we use the same small r as your Python
static constexpr int R = 14;

struct Poly{
	std::array<Cigint, R> c;
};

static inline Cigint z0() { return CIGINT_ZERO(); }
static inline Cigint z1() { return cigint_from_u32(1u); }

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
        u32 nv = v << 1 | carry;
        carry = v >> 31;
        x.data[i] = nv;
    }
    if (carry || cigint_cmp(x, m) >= 0) cigint_sub_ref(&x, &m);
}

// compute n mod R using your own bigint division
static u32 cigint_mod_small(const Cigint &n, u32 m) {
	u32 r = 0;
	cigint_sdivmod(n, m, NULL, &r); // r = n % mod
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

/* a,b < m required in the trusted version (use safe version if unsure) */
static inline Cigint mul_mod_win4_trusted(Cigint a, const Cigint &b, const Cigint &m) {
	if (cigint_is0(a) || cigint_is0(b)) return CIGINT_ZERO();

	Cigint T[16]; for (auto & i : T) i = CIGINT_ZERO();
	T[1] = a;
	T[2] = a; dbl_mod_fast(T[2], m);
	T[4] = T[2]; dbl_mod_fast(T[4], m);
	T[8] = T[4]; dbl_mod_fast(T[8], m);
	for (int d = 3; d < 16; ++d) if (d != 4 && d != 8) {
		Cigint acc = CIGINT_ZERO();
		if (d & 1) add_mod_fast(acc, T[1], m);
		if (d & 2) add_mod_fast(acc, T[2], m);
		if (d & 4) add_mod_fast(acc, T[4], m);
		if (d & 8) add_mod_fast(acc, T[8], m);
		T[d] = acc;
	}

	Cigint res = CIGINT_ZERO();
	u32 bits = cigint_highest_order(b);
	if (!bits) return res;
	u32 total_nibbles = (bits + 3u) / 4u;
	u32 skip = 8u * (u32)CIGINT_N - total_nibbles;

	u32 seen = 0; bool first = true;
	for (unsigned int v : b.data) {
			for (int k = 0; k < 8; ++k) {
			if (seen++ < skip) continue;
			u32 d = (v >> (28 - 4*k)) & 0xFu;
			if (first) {
				if (d) res = T[d];
				first = false;
			} else {
				dbl_mod_fast(res, m); dbl_mod_fast(res, m);
				dbl_mod_fast(res, m); dbl_mod_fast(res, m);
				if (d) add_mod_fast(res, T[d], m);
			}
		}
	}
	return first ? CIGINT_ZERO() : res;
}

/* C[k] = sum_{i+j≡k (mod R)} A[i]*B[j] (mod m), with every coeff kept < m */
static inline int add_modR(int t) { return (t >= R) ? (t - R) : t; }

static Poly poly_mul_fast(const Poly &A, const Poly &B, const Cigint &m) {
	Poly C; for (int i = 0; i < R; ++i) C.c[i] = CIGINT_ZERO();

	for (int i = 0; i < R; ++i) {
		if (cigint_is0(A.c[i])) continue;
		for (int j = 0; j < R; ++j) {
			if (cigint_is0(B.c[j])) continue;
			const int idx = add_modR(i + j);
			Cigint prod = mul_mod_win4_trusted(A.c[i], B.c[j], m);
			add_mod_fast(C.c[idx], prod, m);
		}
	}
	return C;
}

#include <cmath>
using namespace std;

// === Karatsuba cyclic multiply for R=14 over Z_n ===
// Assumes: struct Poly { std::array<Cigint, 14> c; };  // R == 14
// Requires: cigint_cmp, cigint_add_ref, cigint_sub_ref, cigint_divmod_ref,
//           CIGINT_ZERO(), CIGINT_N, and your limb layout (MSW at data[0]).

static inline void add_mod_big(Cigint& x, const Cigint& y, const Cigint& modn) {
    uint64_t carry=0;
    for (size_t i=CIGINT_N; i-- > 0; ) {
        uint64_t t = (uint64_t)x.data[i] + y.data[i] + carry;
        x.data[i] = (uint32_t)t;
        carry = t >> 32;
    }
    if (carry || cigint_cmp(x, modn) >= 0) cigint_sub_ref(&x, &modn);
}

// --- replace sub_mod_big, add_mod_big (unchanged), mul_mod_big+linconv7_mod ---

static inline void sub_noborrow_ip(Cigint &x, const Cigint &y) {
    uint64_t borrow = 0;
    for (size_t i = CIGINT_N; i-- > 0;) {
        uint64_t xi = x.data[i];
        uint64_t yi = (uint64_t)y.data[i] + borrow;
        uint64_t d  = xi - yi;
        borrow = (xi < yi);
        x.data[i] = (uint32_t)d;
    }
}

static inline void add_nocarry_ip(Cigint &x, const Cigint &y) {
    uint64_t carry = 0;
    for (size_t i = CIGINT_N; i-- > 0;) {
        uint64_t t = (uint64_t)x.data[i] + y.data[i] + carry;
        x.data[i] = (uint32_t)t;
        carry = t >> 32;
    }
}

/* x = (x - y) mod modn; requires 0<=x<modn and 0<=y<modn on entry */
static inline void sub_mod_big(Cigint& x, const Cigint& y, const Cigint& modn) {
    if (cigint_cmp(x, y) >= 0) {
        // x := x - y
        sub_noborrow_ip(x, y);
    } else {
        // x := x + modn - y
        Cigint tmp = x;
        add_nocarry_ip(tmp, modn);   // tmp = x + modn  (no overflow since x<modn)
        sub_noborrow_ip(tmp, y);     // tmp -= y
        x = tmp;
    }
    // normalize: result is < 2*modn; bring it into [0, modn)
    if (cigint_cmp(x, modn) >= 0) cigint_sub_ref(&x, &modn);
}

/* use your trusted window-4 mul here for consistency/speed */
static inline Cigint mul_mod_big(Cigint a, const Cigint& b, const Cigint& modn) {
    return mul_mod_win4_trusted(a, b, modn);
}

/* linear convolution 7×7 over Z_modn (coeffs always raeduced) */
static void linconv7_mod(const Cigint A[7], const Cigint B[7], const Cigint& modn, Cigint out[13]) {
    for (int i = 0; i < 13; ++i) out[i] = CIGINT_ZERO();
    for (int i = 0; i < 7; ++i) {
        if (cigint_is0(A[i])) continue;
        for (int j = 0; j < 7; ++j) {
            if (cigint_is0(B[j])) continue;
            Cigint p = mul_mod_big(A[i], B[j], modn);
            add_mod_big(out[i + j], p, modn);
        }
    }
}


// Karatsuba cyclic multiply for R = 14: C(x) = A(x) * B(x) mod (x^14 - 1, n)
static Poly poly_mul_cyclic_kara14(const Poly& A, const Poly& B, const Cigint& modn) {
    // Split A,B into halves of 7
    Cigint A0[7], A1[7], B0[7], B1[7];
    for (int i = 0; i < 7; ++i) {
        A0[i] = A.c[i];
        A1[i] = A.c[i + 7];
        B0[i] = B.c[i];
        B1[i] = B.c[i + 7];
    }

    // Z0 = A0 * B0 (len 13), Z2 = A1 * B1 (len 13)
    Cigint Z0[13], Z2[13];
    linconv7_mod(A0, B0, modn, Z0);
    linconv7_mod(A1, B1, modn, Z2);

    // S = A0 + A1, T = B0 + B1  (mod n)
    Cigint S[7], T[7];
    for (int i = 0; i < 7; ++i) {
        S[i] = A0[i]; add_mod_big(S[i], A1[i], modn);
        T[i] = B0[i]; add_mod_big(T[i], B1[i], modn);
    }

    // Z1 = (A0+A1)*(B0+B1) (len 13)
    Cigint Z1[13]; linconv7_mod(S, T, modn, Z1);

    // M = Z1 - Z0 - Z2 (len 13), modulo n
    Cigint M[13];
    for (int i = 0; i < 13; ++i) {
        M[i] = Z1[i];
        sub_mod_big(M[i], Z0[i], modn);
        sub_mod_big(M[i], Z2[i], modn);
    }

    // Combine into cyclic length-14:
    // C[k] = Z0[k] + Z2[k] + M[(k-7) mod 14], but M is defined 0..12 -> use (t+7) mod 14 placement.
    Poly C; for (int k=0;k<14;++k) C.c[k] = CIGINT_ZERO();

    for (int t = 0; t < 13; ++t) { // Z0, Z2 land at indices t (no extra folding needed)
        add_mod_big(C.c[t], Z0[t], modn);
        add_mod_big(C.c[t], Z2[t], modn);
    }
    // M rotated by +7: add M[t] to C[(t+7) % 14]
    for (int t = 0; t < 13; ++t) {
        int k = t + 7; if (k >= 14) k -= 14;
        add_mod_big(C.c[k], M[t], modn);
    }
    // Note: maximum degree in Z* arrays is 12 (since 6+6), so no further folding needed.

    return C;
}

// Hook: use this instead of your previous poly_mul_fast in the exponentiation core.
static inline Poly POLY_MUL(const Poly& A, const Poly& B, const Cigint& modn) {
    return poly_mul_cyclic_kara14(A, B, modn);
}


/* (1+x)^n mod (x^R-1, n) with 4-bit window on polynomials */
static Poly poly_pow_1x_win4(const Cigint &n, const Cigint &modn) {
	// base = 1 + x
	Poly base; for (int i = 0; i < R; ++i) base.c[i] = CIGINT_ZERO();
	base.c[0] = z1(); base.c[1] = z1();

	// precompute P[1..15] = base^d
	Poly P[16];
	P[0] = {}; for (int i = 0; i < R; ++i) P[0].c[i] = CIGINT_ZERO();
	P[1] = base;
	for (int d = 2; d <= 15; ++d) P[d] = POLY_MUL(P[d - 1], base, modn);

	// res = 1
	Poly res; for (int i = 0; i < R; ++i) res.c[i] = CIGINT_ZERO();
	res.c[0] = z1();

	// scan exponent n by nibbles (MSW → LSW)
	u32 bits = cigint_highest_order(n);
	if (bits == 0) return res;
	u32 n_nibbles = (bits + 3u) / 4u;
	u32 skip = 8u * (u32)CIGINT_N - n_nibbles;

	bool first = true;
	u32 seen = 0;

	for (unsigned int v : n.data) {
			for (int k = 0; k < 8; ++k) {
			if (seen++ < skip) continue;
			u32 d = (v >> (28 - 4*k)) & 0xFu;

			if (first) {
				if (d) res = P[d];
				first = false;
			} else {
				// 4 squaring (res = res^(16))
				res = POLY_MUL(res, res, modn);
				res = POLY_MUL(res, res, modn);
				res = POLY_MUL(res, res, modn);
				res = POLY_MUL(res, res, modn);
				if (d) res = POLY_MUL(res, P[d], modn);
			}
		}
	}
	return first ? res : res;
}

static bool aks_like_cigint(const Cigint &n) {
	// reject 0,1
	if (cigint_is0(n)) return false;
	if (cigint_cmp(n, z1()) == 0) return false;

	Poly p = poly_pow_1x_win4(n, n);
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
		 2, 3, 5, 7,11,13,17,19,23,29,31,37,41,
		43,47,53,59,61,67,71,73,79,83,89,97
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
	// Cigint n = cigint_from_dec("68785536162660004414084120045689084892162254686597696425215272205120048508574098397551086658564179126105400178601774223071432354833665740594087341354341133906161408274754546505739193243119243225335933");
	// Cigint n = cigint_from_dec("12598680682438752944055149498662865476737841742262510750565369004945449194176119763492007280023309651826310840991041734855800903377926420446403888376677867"); //up to 154 digit
	Cigint n = cigint_from_dec("137"); //up to 154 digit

	// Cigint n = gen_prime();
	printf("n = ");
	cigint_print10(n);
	bool ok = aks_like_cigint(n);
	printf(" -> %s\n", ok ? "maybe prime" : "composite-ish");
	return 0;
}
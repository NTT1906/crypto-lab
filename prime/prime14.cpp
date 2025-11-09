#include <array>
#include <cstdio>

#define CIGINT_IMPLEMENTATION
#define CIGINT_STRIP_PREFIX
#define CIGINT_N (1024 / 32)
#include "../cigint.h"

static constexpr int R = 14;

struct Poly {
	std::array<Cigint, R> c;
};

static inline Cigint z0() { return CIGINT_ZERO(); }
static inline Cigint z1() { return cigint_from_u32(1u); }

static inline void add_mod_fast(Cigint &res, const Cigint &b, const Cigint &m) {
	u64 carry = 0;
	for (size_t i = CIGINT_N; i-- > 0;) {
		u64 t = (u64)res.data[i] + (u64)b.data[i] + carry;
		res.data[i] = (u32)t;
		carry = t >> SIZEOF_U32;
	}
	if (carry || cigint_cmp(res, m) >= 0) cigint_sub_ref(&res, &m);
}

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

static void sub_noborrow_ip(Cigint &x, const Cigint &y) {
	u64 borrow = 0;
	for (size_t i = CIGINT_N; i-- > 0;) {
		u64 xi = x.data[i];
		u64 yi = (u64)y.data[i] + borrow;
		u64 d  = xi - yi;
		borrow = (xi < yi);
		x.data[i] = (u32)d;
	}
}

static void add_nocarry_ip(Cigint &x, const Cigint &y) {
	u64 carry = 0;
	for (size_t i = CIGINT_N; i-- > 0;) {
		u64 t = (u64)x.data[i] + y.data[i] + carry;
		x.data[i] = (u32)t;
		carry = t >> 32;
	}
}

static inline void sub_mod_big(Cigint& x, const Cigint& y, const Cigint& modn) {
	if (cigint_cmp(x, y) >= 0) {
		sub_noborrow_ip(x, y);
	} else {
		Cigint tmp = x;
		add_nocarry_ip(tmp, modn);
		sub_noborrow_ip(tmp, y);
		x = tmp;
	}
	if (cigint_cmp(x, modn) >= 0) cigint_sub_ref(&x, &modn);
}

// From Step 1: Add Support for Full-Size Products
typedef uint32_t Limb;
#define CIGINT_2N (2 * CIGINT_N)

struct BigTwice {
	Limb data[CIGINT_2N];
};

BigTwice big_zero() {
	BigTwice res = {0};
	return res;
}

BigTwice full_mul(const Cigint &lhs, const Cigint &rhs) {
	BigTwice res = big_zero();
	for (size_t i = 0; i < CIGINT_N; ++i) {
		u64 carry = 0;
		for (size_t j = 0; j < CIGINT_N; ++j) {
			size_t idx = i + j;
			if (idx >= CIGINT_2N) break;
			u64 acc = (u64)lhs.data[i] * rhs.data[j] +
					  res.data[idx] + carry;
			res.data[idx] = (u32)acc;
			carry = acc >> 32;
			size_t k = idx + 1;
			while (carry && k < CIGINT_2N) {
				u64 sum = res.data[k] + carry;
				res.data[k] = (u32)sum;
				carry = sum >> 32;
				++k;
			}
		}
	}
	return res;
}

void big_add_ref(BigTwice *lhs, const BigTwice *rhs) {
	u64 carry = 0;
	for (size_t i = CIGINT_2N; i-- > 0;) {
		u64 sum = (u64)lhs->data[i] + rhs->data[i] + carry;
		lhs->data[i] = (u32)sum;
		carry = sum >> 32;
	}
}

BigTwice big_shr(const BigTwice &lhs, u32 amnt) {
	if (amnt == 0) return lhs;
	size_t limb_off = amnt / 32;
	u32 bits = amnt % 32;
	if (limb_off >= CIGINT_2N) return big_zero();
	BigTwice r = big_zero();
	for (size_t i = CIGINT_2N; i-- > limb_off; )
		r.data[i] = lhs.data[i - limb_off];
	if (bits) {
		u32 snapshot[CIGINT_2N];
		memcpy(snapshot, r.data, sizeof(snapshot));
		for (size_t i = 0; i < CIGINT_2N; ++i) {
			const u32 lo = snapshot[i] >> bits;
			const u32 hi = i > 0 ? snapshot[i - 1] << (32 - bits) : 0u;
			r.data[i] = hi | lo;
		}
	}
	return r;
}

Cigint cigint_mul_low(Cigint lhs, const Cigint &rhs) {
	Cigint tmp = CIGINT_ZERO();
	u64 carry = 0;
	for (size_t k = 0; k < CIGINT_N; ++k) {
		u64 acc = carry;
		for (size_t i = 0; i <= k; ++i) {
			size_t li = CIGINT_N - 1 - i;
			size_t ri = CIGINT_N - 1 - (k - i);
			acc += (u64)lhs.data[li] * rhs.data[ri];
		}
		size_t out_idx = CIGINT_N - 1 - k;
		tmp.data[out_idx] = (u32)acc;
		carry = acc >> 32;
	}
	size_t idx = 0;
	while (carry && idx < CIGINT_N) {
		u64 sum = (u64)tmp.data[idx] + carry;
		tmp.data[idx] = (u32)sum;
		carry = sum >> 32;
		++idx;
	}
	return tmp;
}

// From Step 2: Modular Inverse mod 2^{512}
// Cigint mod_inverse_2k(const Cigint &m) {
//     Cigint x = cigint_from_u32(1u);
//     Cigint two = cigint_from_u32(2u);
//     for (u32 j = 1; j < SIZEOF_CIGINT; j <<= 1) {  // Double j each time for efficiency
//         Cigint mx = cigint_mul_low(m, x);
//         Cigint neg_one = cigint_not(CIGINT_ZERO());  // -1 mod 2^{512}
//         Cigint diff = cigint_add(mx, neg_one);
//         Cigint corr = cigint_mul_low(diff, x);
//         x = cigint_sub(x, corr);
//         // Truncate to current bit level if needed, but low mul handles it
//     }
//     return x;
// }

Cigint mod_inverse_2k(const Cigint &m) {
	// Newton iteration for 1/m mod 2^512  (9 steps → 2^512)
	Cigint x = cigint_from_u32(1u);
	for (int i = 0; i < 9; ++i) {
		Cigint mx = cigint_mul_low(m, x);          // m·x mod 2^512
		Cigint two_minus = cigint_from_u32(2u);
		cigint_sub_ref(&two_minus, &mx);           // 2 – m·x  (wrap-around OK)
		x = cigint_mul_low(x, two_minus);          // x ← x·(2 – m·x)
	}
	// -----  we need –x  -----
	Cigint minv = CIGINT_ZERO();
	cigint_sub_ref(&minv, &x);                     // minv = –x  (two’s complement)
	return minv;
}

// From Step 3: Montgomery Multiplication
Cigint mont_mul(const Cigint &a, const Cigint &b, const Cigint &m, const Cigint &minv) {
	BigTwice p = full_mul(a, b);
	Cigint u;
	memcpy(u.data, p.data, sizeof(u.data)); // ← low 512 bits
	Cigint v = cigint_mul_low(u, minv);
	BigTwice vm = full_mul(v, m);
	BigTwice t = p;
	big_add_ref(&t, &vm);
	Cigint res;
	memcpy(res.data, &t.data[CIGINT_N], sizeof(res.data));  // ← high 512 bits
	if (cigint_cmp(res, m) >= 0) cigint_sub_ref(&res, &m);
	return res;
}

// Modified linconv7_mod (use mont_mul)
static void linconv7_mod(const Cigint A[7], const Cigint B[7], const Cigint& modn, const Cigint& minv, Cigint out[13]) {
	for (int i = 0; i < 13; ++i) out[i] = CIGINT_ZERO();
	for (int i = 0; i < 7; ++i) {
		if (cigint_is0(A[i])) continue;
		for (int j = 0; j < 7; ++j) {
			if (cigint_is0(B[j])) continue;
			Cigint p = mont_mul(A[i], B[j], modn, minv);
			add_mod_fast(out[i + j], p, modn);
		}
	}
}

// Modified square7_mod (use mont_mul)
static void square7_mod(const Cigint A[7], const Cigint &modn, const Cigint& minv, Cigint out[13]) {
	for (int t = 0; t < 13; ++t) out[t] = CIGINT_ZERO();

	Cigint sq[7];
	for (int i = 0; i < 7; ++i) {
		sq[i] = mont_mul(A[i], A[i], modn, minv);
	}
	for (int i = 0; i < 7; ++i) {
		add_mod_fast(out[2*i], sq[i], modn);
	}

	for (int i = 0; i < 7; ++i) {
		if (cigint_is0(A[i])) continue;
		for (int j = i + 1; j < 7; ++j) {
			if (cigint_is0(A[j])) continue;
			Cigint p = mont_mul(A[i], A[j], modn, minv);
			dbl_mod_fast(p, modn);
			add_mod_fast(out[i + j], p, modn);
		}
	}
}

// Modified poly_mul_cyclic_kara14 (pass minv to linconv7_mod)
static Poly poly_mul_cyclic_kara14(const Poly& A, const Poly& B, const Cigint& modn, const Cigint& minv) {
	const Cigint *A0 = A.c.data();
	const Cigint *A1 = A0 + 7;
	const Cigint *B0 = B.c.data();
	const Cigint *B1 = B0 + 7;

	Cigint S[7], T[7];
	for (int i = 0; i < 7; ++i) {
		S[i] = A0[i];
		add_mod_fast(S[i], A1[i], modn);
		T[i] = B0[i];
		add_mod_fast(T[i], B1[i], modn);
	}

	Cigint Z0[13], Z2[13], Z1[13];

	#pragma omp parallel sections
	{
		#pragma omp section
		linconv7_mod(A0, B0, modn, minv, Z0);
		#pragma omp section
		linconv7_mod(A1, B1, modn, minv, Z2);
		#pragma omp section
		linconv7_mod(S, T, modn, minv, Z1);
	}

	Cigint M[13];
	for (int i = 0; i < 13; ++i) {
		M[i] = Z1[i];
		sub_mod_big(M[i], Z0[i], modn);
		sub_mod_big(M[i], Z2[i], modn);
	}

	Poly C{};
	for (int t = 0; t < 13; ++t) {
		add_mod_fast(C.c[t], Z0[t], modn);
		add_mod_fast(C.c[t], Z2[t], modn);
	}
	for (int t = 0; t < 13; ++t) {
		add_mod_fast(C.c[(t + 7) % R], M[t], modn);
	}
	return C;
}

// Modified poly_square_cyclic_kara14 (pass minv to square7_mod)
static Poly poly_square_cyclic_kara14(const Poly& A, const Cigint& modn, const Cigint& minv) {
	const Cigint *A0 = A.c.data();
	const Cigint *A1 = A0 + 7;

	Cigint S[7];
	for (int i = 0; i < 7; ++i) {
		S[i] = A0[i];
		add_mod_fast(S[i], A1[i], modn);
	}

	Cigint Z0[13], Z2[13], Z1[13];
	#pragma omp parallel sections
	{
		#pragma omp section
		square7_mod(A0, modn, minv, Z0);
		#pragma omp section
		square7_mod(A1, modn, minv, Z2);
		#pragma omp section
		square7_mod(S, modn, minv, Z1);
	}

	Cigint M[13];
	for (int t = 0; t < 13; ++t) {
		M[t] = Z1[t];
		sub_mod_big(M[t], Z0[t], modn);
		sub_mod_big(M[t], Z2[t], modn);
	}

	Poly C = {};
	for (int t = 0; t < 13; ++t) {
		add_mod_fast(C.c[t], Z0[t], modn);
		add_mod_fast(C.c[t], Z2[t], modn);
	}
	for (int t = 0; t < 13; ++t) {
		add_mod_fast(C.c[(t + 7) % R], M[t], modn);
	}
	return C;
}

// Modified POLY_SQR and POLY_MUL (pass minv)
static inline Poly POLY_SQR(const Poly& A, const Cigint& modn, const Cigint& minv) {
	return poly_square_cyclic_kara14(A, modn, minv);
}

static inline Poly POLY_MUL(const Poly& A, const Poly& B, const Cigint& modn, const Cigint& minv) {
	return poly_mul_cyclic_kara14(A, B, modn, minv);
}

// Modified poly_pow_1x_win4 (work in Mont form, pass minv)
static Poly poly_pow_1x_win4(const Cigint &n, const Cigint &modn, const Cigint& minv, const Cigint& r2) {
	Poly base{};
	base.c[0] = z1(); base.c[1] = z1();

	// Convert base to Mont form
	Poly base_mont = base;
	for (int i = 0; i < R; ++i) {
		base_mont.c[i] = mont_mul(base.c[i], r2, modn, minv);
	}

	Poly P[16];
	P[1] = base_mont;
	for (int d = 2; d <= 15; ++d) P[d] = POLY_MUL(P[d - 1], base_mont, modn, minv);

	// res starts as Mont(1)
	Poly res{};
	res.c[0] = mont_mul(z1(), r2, modn, minv);

	// Compute bit length: floor(log2(n)) + 1
	u32 bit_len = 0;
	for (int i = CIGINT_N - 1; i >= 0; --i) {
		if (n.data[i] != 0) {
			bit_len = (CIGINT_N - 1 - i) * 32 + 32 - __builtin_clz(n.data[i]);
			break;
		}
	}
	if (bit_len == 0) return res;

	u32 total_nibbles = (bit_len + 3) / 4;
	u32 skip = 8 * CIGINT_N - total_nibbles;

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
				res = POLY_SQR(res, modn, minv);
				res = POLY_SQR(res, modn, minv);
				res = POLY_SQR(res, modn, minv);
				res = POLY_SQR(res, modn, minv);
				if (d) res = POLY_MUL(res, P[d], modn, minv);
			}
		}
	}
	return first ? res : res;
}

static bool poly_is_xn_plus_1(const Poly &P, const Cigint &n, const Cigint &m) {
	u32 k;
	cigint_sdivmod(n, R, NULL, &k);
	for (int i = 0; i < R; ++i) {
		bool should_be_one = (i == 0) || (i == (int)k && k != 0);
		if (cigint_cmp(P.c[i], should_be_one ? z1() : z0()) != 0) {
			return false;
		}
	}
	return true;
}

// Modified aks_like_cigint (add precomputations, convert back at end)
static bool aks_like_cigint(const Cigint &n) {
	if (cigint_is0(n)) return false;
	if (cigint_cmp(n, z1()) == 0) return false;

	// Precomputations for Montgomery
	Cigint minv = mod_inverse_2k(n);  // returns –n⁻¹
	Cigint r = z1();
	for (int i = 0; i < 512; ++i) dbl_mod_fast(r, n);
	Cigint r2 = r;
	for (int i = 0; i < 512; ++i) dbl_mod_fast(r2, n);  // 2^{1024} mod n

	Poly p = poly_pow_1x_win4(n, n, minv, r2);

	// Convert p back from Mont form
	for (int i = 0; i < R; ++i) {
		p.c[i] = mont_mul(p.c[i], z1(), n, minv);
	}

	return poly_is_xn_plus_1(p, n, n);
}

#include <random>
static void cigint_randomize(Cigint &x) {
	std::random_device rd; std::mt19937 gen(rd());
	std::uniform_int_distribution<u32> dist(0, UINT32_MAX);
	size_t limbs = 1 + gen() % (CIGINT_N / 2);
	for (u32 &i : x.data) i = 0;
	for (size_t i = limbs; i-- > 0;) x.data[i] = dist(gen);
}

static Cigint random_odd() {
	Cigint x; cigint_randomize(x);
	cigint_set_bit_ref(&x, 0, 1);
	return x;
}

static bool has_small_factor(const Cigint &n) {
	static const int SMALL_PRIMES[] = {
		 2, 3, 5, 7,11,13,17,19,23,29,31,37,41,
		43,47,53,59,61,67,71,73,79,83,89,97
	};
	for (int p : SMALL_PRIMES) {
		u32 r = 0;
		cigint_sdivmod(n, (u32)p, NULL, &r);
		if (r == 0) return cigint_cmp(n, cigint_from_u32((u32)p)) != 0;
	}
	return false;
}

static Cigint gen_prime() {
	Cigint x;
	do { x = random_odd(); }
	while (has_small_factor(x) || !aks_like_cigint(x));
	return x;
}

#include "../benchmark.h"
// observable sink to block DCE
static volatile u32 g_sink = 0;
inline void consume(const Cigint& x) {
	u32 acc = 0;
	for (unsigned int i : x.data) acc ^= i; g_sink ^= acc;
}

template<typename F>
void bench_func_ref(const std::string& name, F f, const Cigint& n) {
	volatile F fp = f;
	Cigint x;
	bench(name, [&]{
		Cigint x = n;
		(*fp)(x);
		consume(x);
	}, 10);
}

int main() {
	Cigint n = cigint_from_dec("13");
	// Cigint n = gen_prime();
	// bench_func_ref("aks", aks_like_cigint, n);
	printf("n = ");
	cigint_print10(n);
	bool ok = aks_like_cigint(n);
	printf("\nAKS   : %s\n", ok ? "maybe prime" : "composite-ish");
	return 0;
}
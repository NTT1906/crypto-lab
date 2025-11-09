#include <array>
#include <cstdio>

#define CIGINT_IMPLEMENTATION
#define CIGINT_STRIP_PREFIX
#define CIGINT_N (512 / 32)
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

static inline Cigint mul_mod_big(Cigint a, const Cigint &b, const Cigint &modn) {
	return mul_mod_win4_trusted(a, b, modn);
}

static inline void add_mod_big(Cigint &x, const Cigint &y, const Cigint &modn) {
	u64 carry = 0;
	for (size_t i = CIGINT_N; i-- > 0; ) {
		u64 t = (u64)x.data[i] + (u64)y.data[i] + carry;
		x.data[i] = (u32)t;
		carry = t >> 32;
	}
	if (carry || cigint_cmp(x, modn) >= 0) cigint_sub_ref(&x, &modn);
}

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

static Poly poly_mul_cyclic_kara14(const Poly& A, const Poly& B, const Cigint& modn) {
	Cigint A0[7], A1[7], B0[7], B1[7];
	for (int i = 0; i < 7; ++i) {
		A0[i] = A.c[i]; A1[i] = A.c[i + 7];
		B0[i] = B.c[i]; B1[i] = B.c[i + 7];
	}

	Cigint S[7], T[7];
	for (int i = 0; i < 7; ++i) {
		S[i] = A0[i]; add_mod_big(S[i], A1[i], modn);
		T[i] = B0[i]; add_mod_big(T[i], B1[i], modn);
	}

	Cigint Z0[13], Z2[13], Z1[13];

	#pragma omp parallel sections
	{
		#pragma omp section
		linconv7_mod(A0, B0, modn, Z0);
		#pragma omp section
		linconv7_mod(A1, B1, modn, Z2);
		#pragma omp section
		linconv7_mod(S, T, modn, Z1);
	}

	Cigint M[13];
	for (int i = 0; i < 13; ++i) {
		M[i] = Z1[i];
		sub_mod_big(M[i], Z0[i], modn);
		sub_mod_big(M[i], Z2[i], modn);
	}

	Poly C;
	for (int k = 0; k < R; ++k) C.c[k] = CIGINT_ZERO();

	for (int t = 0; t < 13; ++t) {
		add_mod_big(C.c[t], Z0[t], modn);
		add_mod_big(C.c[t], Z2[t], modn);
	}
	for (int t = 0; t < 13; ++t) {
		int k = t + 7; if (k >= R) k -= R;
		add_mod_big(C.c[k], M[t], modn);
	}
	return C;
}

static void square7_mod(const Cigint A[7], const Cigint &modn, Cigint out[13]) {
	for (int t = 0; t < 13; ++t) out[t] = CIGINT_ZERO();

	Cigint sq[7];
	for (int i = 0; i < 7; ++i) {
		sq[i] = mul_mod_big(A[i], A[i], modn);
	}
	for (int i = 0; i < 7; ++i) {
		add_mod_big(out[2*i], sq[i], modn);
	}

	for (int i = 0; i < 7; ++i) {
		if (cigint_is0(A[i])) continue;
		for (int j = i + 1; j < 7; ++j) {
			if (cigint_is0(A[j])) continue;
			Cigint p = mul_mod_big(A[i], A[j], modn);
			dbl_mod_fast(p, modn);
			add_mod_big(out[i + j], p, modn);
		}
	}
}

static Poly poly_square_cyclic_kara14(const Poly& A, const Cigint& modn) {
	Cigint A0[7], A1[7];
	for (int i = 0; i < 7; ++i) { A0[i] = A.c[i]; A1[i] = A.c[i + 7]; }
	Cigint S[7];
	for (int i = 0; i < 7; ++i) { S[i] = A0[i]; add_mod_big(S[i], A1[i], modn); }

	Cigint Z0[13], Z2[13], Z1[13];
	#pragma omp parallel sections
	{
		#pragma omp section
		square7_mod(A0, modn, Z0);
		#pragma omp section
		square7_mod(A1, modn, Z2);
		#pragma omp section
		square7_mod(S, modn, Z1);
	}

	Cigint M[13];
	for (int t = 0; t < 13; ++t) {
		M[t] = Z1[t];
		sub_mod_big(M[t], Z0[t], modn);
		sub_mod_big(M[t], Z2[t], modn);
	}

	Poly C; for (int k = 0; k < R; ++k) C.c[k] = CIGINT_ZERO();

	for (int t = 0; t < 13; ++t) {
		add_mod_big(C.c[t], Z0[t], modn);
		add_mod_big(C.c[t], Z2[t], modn);
	}
	for (int t = 0; t < 13; ++t) {
		int k = t + 7; if (k >= R) k -= R;
		add_mod_big(C.c[k], M[t], modn);
	}
	return C;
}

static inline Poly POLY_SQR(const Poly& A, const Cigint& modn) {
	return poly_square_cyclic_kara14(A, modn);
}

static inline Poly POLY_MUL(const Poly& A, const Poly& B, const Cigint& modn) {
	return poly_mul_cyclic_kara14(A, B, modn);
}

static Poly poly_pow_1x_win4(const Cigint &n, const Cigint &modn) {
	Poly base; for (int i = 0; i < R; ++i) base.c[i] = CIGINT_ZERO();
	base.c[0] = z1(); base.c[1] = z1();

	Poly P[16];
	P[0] = Poly{}; for (int i = 0; i < R; ++i) P[0].c[i] = CIGINT_ZERO();
	P[1] = base;
	for (int d = 2; d <= 15; ++d) P[d] = POLY_MUL(P[d - 1], base, modn);

	Poly res; for (int i = 0; i < R; ++i) res.c[i] = CIGINT_ZERO();
	res.c[0] = z1();

	u32 bits = cigint_highest_order(n);
	if (bits == 0) return res;
	u32 total_nibbles = (bits + 3u) / 4u;
	u32 skip = 8u * (u32)CIGINT_N - total_nibbles;

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
				res = POLY_SQR(res, modn);
				res = POLY_SQR(res, modn);
				res = POLY_SQR(res, modn);
				res = POLY_SQR(res, modn);
				if (d) res = POLY_MUL(res, P[d], modn);
			}
		}
	}
	return first ? res : res;
}

static bool poly_is_xn_plus_1(const Poly &P, const Cigint &n, const Cigint &m) {
	u32 k;
	cigint_sdivmod(n, R, NULL, &k);
	for (int i = 0; i < R; ++i) {
		if (i == 0 || i == (int)k) {
			if (cigint_cmp(P.c[i], z1()) != 0) return false;
		} else {
			if (!cigint_is0(P.c[i])) return false;
		}
	}
	return true;
}

static bool aks_like_cigint(const Cigint &n) {
	if (cigint_is0(n)) return false;
	if (cigint_cmp(n, z1()) == 0) return false;

	Poly p = poly_pow_1x_win4(n, n);
	return poly_is_xn_plus_1(p, n, n);
}

#include <random>
static void cigint_randomize(Cigint &x) {
	std::random_device rd; std::mt19937 gen(rd());
	std::uniform_int_distribution<u32> dist(0, UINT32_MAX);
	size_t limbs = 1 + gen() % CIGINT_N;
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

int main() {
	Cigint n = cigint_from_dec("12598680682438752944055149498662865476737841742262510750565369004945449194176119763492007280023309651826310840991041734855800903377926420446403888376677867");
	// Cigint n = gen_prime();
	printf("n = ");
	cigint_print10(n);
	bool ok = aks_like_cigint(n);
	printf(" -> %s\n", ok ? "maybe prime" : "composite-ish");
	return 0;
}
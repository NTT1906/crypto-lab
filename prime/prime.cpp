#include <array>
#include <vector>
#include <cstdio>

#define CIGINT_IMPLEMENTATION
#define CIGINT_STRIP_PREFIX
#define CIGINT_N (1024 / 32)
#include "../cigint.h"

// we use the same small r as your Python
static constexpr int R = 14;

struct Poly {
	std::array<Cigint, R> c;
};

// -------- basic Cigint helpers --------
static inline Cigint z0() { return CIGINT_ZERO(); }
static inline Cigint z1() { return cigint_from_u32(1u); }

// res = (a + b) mod m (all fixed width)
static inline void add_mod(Cigint& res, const Cigint& b, const Cigint &m) {
	Cigint res_c = res;
	cigint_add_ref(&res, &b);
	if (res < res_c) {
		cigint_printf("res = %Cd\n", res);
		assert(false && "add_mod overflow");
	}
	if (cigint_cmp(res, m) >= 0) {
		cigint_sub_ref(&res, &m);
	}
}

void cigint_mul_ref_2_mod(Cigint *lhs, Cigint *rhs, const Cigint *m) {
	cigint_divmod_ref(lhs, m, NULL, lhs);
	cigint_divmod_ref(rhs, m, NULL, rhs);
	__uint128_t carry = 0;
	for (size_t k = 0; k < CIGINT_N; ++k) {
		__uint128_t acc = carry;
		for (size_t i = 0; i <= k; ++i) {
			acc += (__uint128_t)lhs->data[CIGINT_N - 1 - i]
				 * (__uint128_t)rhs->data[CIGINT_N - 1 - (k - i)];
		}
		carry = acc >> 32;
		lhs->data[CIGINT_N - 1 - k] = (u32)acc;
	}
	cigint_divmod_ref(lhs, m, NULL, lhs);
}

Cigint cigint_mul_2_mod(Cigint lhs, Cigint rhs, Cigint m) {
	cigint_mul_ref_2_mod(&lhs, &rhs, &m);
	return lhs;
}

// poly multiplication mod (x^R - 1, m)
// C[k] = sum_{i+j≡k mod R} A[i]*B[j] (mod m)
void poly_mul_mod(Poly *A, Poly *B, Cigint *m) {
	Poly C;
	for (int i = 0; i < R; ++i) C.c[i] = z0();
	for (int i = 0; i < R; ++i) {
		if (cigint_is0(A->c[i])) continue;
		for (int j = 0; j < R; ++j) {
			if (cigint_is0(B->c[j])) continue;
			int idx = (i + j) % R;
			Cigint r = cigint_mul_2_mod(A->c[i], B->c[j], *m);
			// cigint_divmod_ref(&C.c[idx], m, NULL, &C.c[idx]);
			// cigint_divmod_ref(&r, m, NULL, &r);
			add_mod(C.c[idx], r, *m);
		}
	}
	*A = C;
}

// compute n mod R using your own bigint division
static u32 cigint_mod_small(const Cigint &n, u32 m) {
	u32 r = 0;
	cigint_sdivmod(n, m, NULL, &r);  // r = n % mod
	return r;
}

// check poly == 1 + x^(n mod R)
static bool poly_is_xn_plus_1(const Poly &P, const Cigint &n, const Cigint &m) {
	if (cigint_cmp(P.c[0], z1()) != 0) return false;
	u32 k = cigint_mod_small(n, R);   // exact n % R
	if (cigint_cmp(P.c[k], z1()) != 0) return false;
	for (int i = 0; i < R; ++i) {
		if (!cigint_is0(P.c[i])) return false;
	}
	return true;
}

void pow_mod_2(Cigint *x, Cigint *e, Cigint *m) {
	Cigint ans = 1;
	u32 hob = cigint_highest_order(*e);
	for (size_t i = 0; i < hob; ++i) {
		if (cigint_get_bit_ref(e, i)) {
			cigint_mul_ref_2_mod(&ans, x, m);
		}
		cigint_mul_ref_2_mod(x, x, m);
	}
	*x = ans;
}

// compute (1 + x)^n mod (x^R - 1, n)
static Poly poly_pow_1x(Cigint &n, Cigint &modn) {
	// base = 1 + x
	Poly base;
	for (int i = 0; i < R; ++i) base.c[i] = z0();
	base.c[0] = z1();
	base.c[1] = z1();

	// res = 1
	Poly res;
	for (int i = 0; i < R; ++i) res.c[i] = z0();
	res.c[0] = z1();

	u32 hob = cigint_highest_order(n);
	for (size_t i = 0; i < hob; ++i) {
		if (cigint_get_bit(n, i)) {
			poly_mul_mod(&res, &base, &modn);
		}
		Poly tmp = base;
		poly_mul_mod(&base, &tmp, &modn);
	}
	return res;
}

static bool aks_like_cigint(Cigint &n) {
	// reject 0,1
	if (cigint_is0(n)) return false;
	if (cigint_cmp(n, z1()) == 0) return false;

	Poly p = poly_pow_1x(n, n);
	return poly_is_xn_plus_1(p, n, n);
}

int main() {
	// your prime
	// Cigint n = cigint_from_dec("7669373798138436444425917457853162445620122597663674808454541043947311884377"); //up to 76digit
	// Cigint n = cigint_from_dec("771458402252315418489174122874340496560720991135040845033515625723932549727443357205751720461233037696609749135759573352834404574203828395431821751181407"); //up to 154 digit
	// Cigint n = cigint_from_dec("12598680682438752944055149498662865476737841742262510750565369004945449194176119763492007280023309651826310840991041734855800903377926420446403888376677867"); //up to 154 digit
	Cigint n = cigint_from_dec("5"); //up to 154 digit
	printf("n = ");
	cigint_print10(n);
	bool ok = aks_like_cigint(n);
	printf(" -> %s\n", ok ? "maybe prime" : "composite-ish");
	return 0;
}

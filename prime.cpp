// #include <chrono>
#include <iostream>
#define BI_BIT 512
#include "bigint.h"

constexpr u32 POLY_R = 14;
struct Poly : std::array<bui, POLY_R>{};

// void printBuiA(const bui *poly, int n) {
// 	printf("{");
// 	for (int i = 0; i < n; ++i) {
// 		printf("%s, ", bui_to_dec(poly[i]).c_str());
// 	}
// 	printf("}\n");
// }

// static void poly_mul_mod_ip(Poly &A, const Poly &B, const bui& m) {
// 	bool a_skips[POLY_R] = {false};
// 	bool b_skips[POLY_R] = {false};
// 	for (int i = 0; i < POLY_R; i++) {
// 		a_skips[i] = bui_is0(A[i]);
// 		b_skips[i] = bui_is0(B[i]);
// 	}
// 	Poly C = {};
// 	for (int i = 0; i < POLY_R; ++i) {
// 		if (a_skips[i]) continue;
// 		for (int j = 0; j < POLY_R; ++j) {
// 			if (b_skips[j]) continue;
// 			bui p = A[i];
// 			mul_mod_ip(p, B[j], m);
// 			add_mod_ip(C[(i + j) % POLY_R], p, m);
// 		}
// 	}
// 	A = C;
// }

// static void poly_sqr_mod_ip(Poly &A, const bui& m) {
// 	bool a_skips[POLY_R] = {false};
// 	for (int i = 0; i < POLY_R; i++) {
// 		a_skips[i] = bui_is0(A[i]);
// 	}
// 	Poly C = {};
// 	for (int i = 0; i < POLY_R; ++i) {
// 		if (a_skips[i]) continue;
// 		for (int j = 0; j < POLY_R; ++j) {
// 			if (a_skips[j]) continue;
// 			bui p = A[i];
// 			mul_mod_ip(p, A[j], m);
// 			add_mod_ip(C[(i + j) % POLY_R], p, m);
// 		}
// 	}
// 	A = C;
// }
//
// Poly poly_pow_1x(const bui &n) {
// 	Poly base{};
// 	base[0] = bui1(); base[1] = bui1();
// 	Poly res{};
// 	res[0] = bui1();
// 	u32 hb = highest_bit(n);
// 	for (u32 i = 0; i < hb; ++i) {
// 		if (get_bit(n, i)) {
// 			poly_mul_mod_ip(res, base, n);
// 		}
// 		poly_sqr_mod_ip(base, n);
// 	}
// 	return res;
// }

// a = (a + b) % m
inline void add_true_mod_ip(bui &a, bui b, const bui &m) {
	a = mod_native(a, m);
	b = mod_native(b, m);
	add_ip_n_imp(a.data(), b.data(), BI_N);
	a = mod_native(a, m);
}

static void poly_mul_mod_mont_ip(Poly &A, const Poly &B, const MontgomeryReducer &mr) {
	bool a_skips[POLY_R] = {false};
	bool b_skips[POLY_R] = {false};
	for (int i = 0; i < POLY_R; ++i) {
		a_skips[i] = bui_is0(A[i]);
		b_skips[i] = bui_is0(B[i]);
	}
	Poly C{};
	for (int i = 0; i < POLY_R; ++i) {
		if (a_skips[i]) continue;
		for (int j = 0; j < POLY_R; ++j) {
			if (b_skips[j]) continue;
			bui p = mr.multiply(A[i], B[j]); // Montgomery product (result in Mont. form)
			add_true_mod_ip(C[(i + j) % POLY_R], p, mr.modulus); // addition mod m (same domain)
		}
	}
	A = C;
}

static void poly_sqr_mod_mont_ip(Poly &A, const MontgomeryReducer &mr) {
	bool a_skips[POLY_R] = {false};
	for (int i = 0; i < POLY_R; ++i) a_skips[i] = bui_is0(A[i]);
	Poly C{};
	for (int i = 0; i < POLY_R; ++i) {
		if (a_skips[i]) continue;
		for (int j = 0; j < POLY_R; ++j) {
			if (a_skips[j]) continue;
			bui p = mr.multiply(A[i], A[j]); // montgomery multiply
			add_true_mod_ip(C[(i + j) % POLY_R], p, mr.modulus);
		}
	}
	A = C;
}

Poly poly_pow_1x_mont(const bui &n) {
	MontgomeryReducer mr(n);
	Poly base{}; base[0] = mr.convertedOne; base[1] = mr.convertedOne;
	Poly res{};  res[0] = mr.convertedOne;

	u32 hb = highest_bit(n);
	for (u32 i = 0; i < hb; ++i) {
		if (get_bit(n, i)) {
			poly_mul_mod_mont_ip(res, base, mr);
		}
		poly_sqr_mod_mont_ip(base, mr);
	}
	// convert result coefficients back to standard form
	for (int i = 0; i < POLY_R; ++i) {
		if (!bui_is0(res[i])) res[i] = mr.convertOut(res[i]);
	}
	return res;
}

// TODO: gate-away small number?
static bool aks_like_prime(const bui &n) {
	if (!get_bit(n, 0)) return false;
	Poly p = poly_pow_1x_mont(n);
	bui b1 = bui1();
	if (cmp(p[0], b1) != 0) return false;
	u32 k;
	bui q;
	u32divmod(n, POLY_R, q,  k);
	if (cmp(p[k], b1) != 0) return false;
	for (u32 i = 1; i < POLY_R; ++i) {
		if (i == k) continue;
		if (!bui_is0(p[i])) return false;
	}
	return true;
}

// static bool has_small_factor(const bui &n) {
// 	static const int SMALL_PRIMES[] = {
// 		2, 3, 5, 7,11,13,17,19,23,29,31,37,41,
// 	   43,47,53,59,61,67,71,73,79,83,89,97
//    };
// 	for (int p : SMALL_PRIMES) {
// 		u32 r = 0; bui tmp;
// 		u32divmod(n, (u32)p, tmp, r);
// 		if (r == 0) return cmp(n, bui_from_u32((u32)p)) != 0;
// 	}
// 	return false;
// }

// static bui gen_prime() {
// 	bui x;
// 	do {
// 		x = random_odd();
// 		// printf("Testing: %s\n", bui_to_dec(x).c_str());
// 	}
// 	while (has_small_factor(x) || !aks_like_prime(x));
// 	return x;
// }

int main(int argc, char* argv[]) {
	if (argc != 3) {
		return 1;
	}

	if (!freopen(argv[1], "r", stdin)) return 1;
	if (!freopen(argv[2], "w", stdout)) return 1;

	std::ios::sync_with_stdio(false);
	std::cin.tie(nullptr);

	bui p = read_bui_le();
	// printf("%s\n", bui_to_dec(p).c_str());
	std::cout << (aks_like_prime(p) ? "1" : "0") << '\n';
	return 0;
}

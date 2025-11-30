// #include <chrono>
#include <iostream>
#define BI_BIT 512
#include "bigint.h"

constexpr u32 POLY_R = 14;
struct Poly : std::array<bui, POLY_R>{};

void printBuiA(const bui *poly, int n) {
	printf("{");
	for (int i = 0; i < n; ++i) {
		printf("%s, ", bui_to_dec(poly[i]).c_str());
	}
	printf("}\n");
}

static void poly_mul_mod_ip(Poly &A, const Poly &B, const bui& m) {
	bool a_skips[POLY_R] = {false};
	bool b_skips[POLY_R] = {false};
	for (int i = 0; i < POLY_R; i++) {
		a_skips[i] = bui_is0(A[i]);
		b_skips[i] = bui_is0(B[i]);
	}
	Poly C = {};
	for (int i = 0; i < POLY_R; ++i) {
		if (a_skips[i]) continue;
		for (int j = 0; j < POLY_R; ++j) {
			if (b_skips[j]) continue;
			bui p = A[i];
			mul_mod_ip(p, B[j], m);
			add_mod_ip(C[(i + j) % POLY_R], p, m);
		}
	}
	A = C;
}

static void poly_sqr_mod_ip(Poly &A, const bui& m) {
	bool a_skips[POLY_R] = {false};
	for (int i = 0; i < POLY_R; i++) {
		a_skips[i] = bui_is0(A[i]);
	}
	Poly C = {};
	for (int i = 0; i < POLY_R; ++i) {
		if (a_skips[i]) continue;
		for (int j = 0; j < POLY_R; ++j) {
			if (a_skips[j]) continue;
			bui p = A[i];
			mul_mod_ip(p, A[j], m);
			add_mod_ip(C[(i + j) % POLY_R], p, m);
		}
	}
	A = C;
}

Poly poly_pow_1x(const bui &n) {
	Poly base{};
	base[0] = bui1(); base[1] = bui1();
	Poly res{};
	res[0] = bui1();
	u32 hb = highest_bit(n);
	for (u32 i = 0; i < hb; ++i) {
		if (get_bit(n, i)) {
			poly_mul_mod_ip(res, base, n);
			printf("%5d: R1: ", i);
			printBuiA(res.data(), res.size());
		}
		poly_sqr_mod_ip(base, n);
		printf("%5d: B1: ", i);
		printBuiA(base.data(), base.size());
	}
	return res;
}

// a = (a + b) % m
inline void add_true_mod_ip(bui &a, bui b, const bui &m) {
	if (add_ip_n_imp(a.data(), b.data(), BI_N) || cmp(a, m) >= 0) {
		sub_ip(a, m);
	}
}

static void poly_mul_mod_mont_ip(Poly &A, const Poly &B, const MontgomeryReducer &mr) {
	bool a_skips[POLY_R] = {false};
	bool b_skips[POLY_R] = {false};
	for (int i = 0; i < POLY_R; ++i) {
		a_skips[i] = bui_is0(A[i]);
		b_skips[i] = bui_is0(B[i]);
	}
	Poly C{};
	// Poly cv{};
	// Poly cv2{};
	// Poly cC{};
	// for (int i = 0; i < POLY_R; ++i) {
	// 	if (!bui_is0(A[i])) cv[i] = mr.convertOut(A[i]);
	// 	if (!bui_is0(B[i])) cv2[i] = mr.convertOut(B[i]);
	// }
	for (int i = 0; i < POLY_R; ++i) {
		if (a_skips[i]) continue;
		for (int j = 0; j < POLY_R; ++j) {
			if (b_skips[j]) continue;
			bui p = mr.multiply(A[i], B[j]);
			add_true_mod_ip(C[(i + j) % POLY_R], p, mr.modulus);
			// bui tp = cv[i];
			// mul_mod_ip(tp, cv2[j], mr.modulus);
			// if (cmp(mr.convertOut(p), tp) != 0) {
				// printf("NO MUL: %5d - %5d\n", i, j);
				// exit(1);
			// }
			// add_mod_ip(cC[(i + j) % POLY_R], tp, mr.modulus);
			// if (cmp(mr.convertOut(C[(i + j) % POLY_R]), cC[(i + j) % POLY_R]) != 0) {
				// printf("NOE MUL: %5d - %5d\n", i, j);
				// exit(1);
			// }
		}
	}
	A = C;
}

static void poly_sqr_mod_mont_ip(Poly &A, const MontgomeryReducer &mr) {
	bool a_skips[POLY_R] = {false};
	for (int i = 0; i < POLY_R; ++i) a_skips[i] = bui_is0(A[i]);
	Poly C{};
	// Poly cv{};
	// Poly cC{};
	// for (int i = 0; i < POLY_R; ++i) {
		// if (!bui_is0(A[i])) cv[i] = mr.convertOut(A[i]);
	// }
	for (int i = 0; i < POLY_R; ++i) {
		if (a_skips[i]) continue;
		for (int j = 0; j < POLY_R; ++j) {
			if (a_skips[j]) continue;
			bui p = mr.multiply(A[i], A[j]);
			add_true_mod_ip(C[(i + j) % POLY_R], p, mr.modulus);
			// bui tp = cv[i];
			// mul_mod_ip(tp, cv[j], mr.modulus);
			// if (cmp(mr.convertOut(p), tp) != 0) {
				// printf("NO SQR: %5d - %5d\n", i, j);
			// }
			// add_mod_ip(cC[(i + j) % POLY_R], tp, mr.modulus);
			// if (cmp(mr.convertOut(C[(i + j) % POLY_R]), cC[(i + j) % POLY_R]) != 0) {
				// printf("NOE SQR: %5d - %5d\n", i, j);
			// }
		}
	}
	A = C;
}

void printOg(Poly p, const MontgomeryReducer& mr) {
	for (int i = 0; i < POLY_R; ++i) {
		if (!bui_is0(p[i])) p[i] = mr.convertOut(p[i]);
	}
	printBuiA(p.data(), p.size());
}

Poly poly_pow_1x_mont(const bui &n) {
	MontgomeryReducer mr(n);
	Poly base{}; base[0] = mr.convertedOne; base[1] = mr.convertedOne;
	Poly res{};  res[0] = mr.convertedOne;

	u32 hb = highest_bit(n);
	for (u32 i = 0; i < hb; ++i) {
		if (get_bit(n, i)) {
			poly_mul_mod_mont_ip(res, base, mr);
			// printf("%5d: R2: ", i);
			// printOg(res, mr);
		}
		poly_sqr_mod_mont_ip(base, mr);
		// printf("%5d: B2: ", i);
		// printOg(base, mr);
	}
	// printBuiA(res.data(), res.size());
	// convert result coefficients back to standard form
	for (int i = 0; i < POLY_R; ++i) {
		if (!bui_is0(res[i])) res[i] = mr.convertOut(res[i]);
	}
	// printBuiA(res.data(), res.size());
	return res;
}

// TODO: gate-away small number?
static bool aks_like_prime(const bui &n) {
	if (!get_bit(n, 0)) return false;
	Poly p = poly_pow_1x_mont(n);
	// p = poly_pow_1x(n);
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
	// if (argc != 3) {
		// return 1;
	// }

	// if (!freopen(argv[1], "r", stdin)) return 1;
	// if (!freopen(argv[2], "w", stdout)) return 1;
	// if (!freopen("out.txt", "w", stdout)) return 1;

	// std::ios::sync_with_stdio(false);
	// std::cin.tie(nullptr);

	// bui p = read_bui_le();
	bui p = bui_from_dec("9862580434556848933093118044369795906452209005604134993142891065799068045921485909427627718142455707644541651618163328127809698482899632857003280134349623");
	printf("%s\n", bui_to_dec(p).c_str());
	std::cout << (aks_like_prime(p) ? "1" : "0") << '\n';
	return 0;
}
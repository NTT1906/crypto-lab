#define BI_BIT 1024
#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "bigint.h"

constexpr u32 POLY_R = 14;
struct Poly : std::array<bui, POLY_R>{};

// static void linconv7_mod(const bui A[7], const bui B[7], const bui& modn, bui out[13]) {
// 	bool a_skips[7] = {false};
// 	bool b_skips[7] = {false};
// 	for (int i = 0; i < 7; i++) {
// 		a_skips[i] = bui_is0(A[i]);
// 		b_skips[i] = bui_is0(B[i]);
// 	}
// 	for (int i = 0; i < 13; ++i) out[i] = {};
// 	for (int i = 0; i < 7; ++i) {
// 		if (a_skips[i]) continue;
// 		for (int j = 0; j < 7; ++j) {
// 			if (b_skips[j]) continue;
// 			bui p = A[i];
// 			mul_mod_ip(p, B[j], modn);
// 			add_mod_ip(out[i + j], p, modn);
// 		}
// 	}
// }

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
		}
		poly_sqr_mod_ip(base, n);
	}
	return res;
}

// a = (a + b) % m
inline void add_true_mod_ip(bui &a, bui b, const bui &m) {
	a = mod_native(a, m);
	b = mod_native(b, m);
	add_ip_n_imp(a.data(), b.data(), BI_N);
	a = mod_native(a, m);
}

static void poly_mul_mod_mont_ip(Poly &A, const Poly &B, MontgomeryReducer &mr) {
	bool a_skips[POLY_R] = {false};
	bool b_skips[POLY_R] = {false};
	for (int i = 0; i < POLY_R; ++i) {
		a_skips[i] = bui_is0(A[i]);
		b_skips[i] = bui_is0(B[i]);
	}
	Poly C = {};
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

static void poly_sqr_mod_mont_ip(Poly &A, MontgomeryReducer &mr) {
	bool a_skips[POLY_R] = {false};
	for (int i = 0; i < POLY_R; ++i) a_skips[i] = bui_is0(A[i]);
	Poly C = {};
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
void printBuiA(bui *buis, int n);

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

	// printf("Pre: ");
	// printBuiA(res.data(), POLY_R);

	// Convert result coefficients back to standard form
	for (int i = 0; i < POLY_R; ++i) {
		if (!bui_is0(res[i])) res[i] = mr.convertOut(res[i]);
	}
	return res;
}


void printBuiA(bui *buis, int n) {
	printf("{");
	for (int i = 0; i < n; ++i) {
		printf("%s, ", bui_to_dec(buis[i]).c_str());
	}
	printf("}\n");
}

static bool aks_like_prime(const bui &n) {
	if (!get_bit(n, 0)) return false;
	Poly p = poly_pow_1x_mont(n);
	// printBuiA(p.data(), POLY_R);
	// p = poly_pow_1x(n);
	// printBuiA(p.data(), POLY_R);
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

using namespace std;

string normalize_hex_le_to_be(string s) {
	string hex;
	for (char c : s) {
		if (!isspace((unsigned char)c)) hex.push_back(c);
	}
	if (hex.empty()) return string("0");
	reverse(hex.begin(), hex.end());
	return hex;
}

bui read_bui_le() {
	string line;
	getline(cin, line);
	string be_hex = normalize_hex_le_to_be(line);
	return bui_from_hex(be_hex);
}

int main(int argc, char* argv[]) {
	if (argc != 3) {
		return 1;
	}

	if (!freopen(argv[1], "r", stdin)) return 1;
	// if (!freopen(argv[2], "w", stdout)) return 1;

	// std::ios::sync_with_stdio(false);
	// std::cin.tie(nullptr);

	bui p = read_bui_le();
	printf("%s\n", bui_to_dec(p).c_str());
	std::cout << (aks_like_prime(p) ? "prime" : "composite") << '\n';

	// Poly a = {bui1(), bui1()};
	// printBuiA(a.data(), POLY_R);
	// poly_sqr_mod_ip(a, bui_from_u32(100007));
	// printBuiA(a.data(), POLY_R);
	//
	// // Benchmark optimized version
	// Poly b;
	// int N = 3;
	// bui n = bui_from_dec("6274904334290417405341624571932224150456224549917673444239237760272785701939526927698156030175801211624849856326839256526253153336777911614668501375751381");
	bui n = bui_from_dec("91");
	printf("%s\n", bui_to_dec(n).c_str());
	if (aks_like_prime(n)) {
		printf("PRIME\n");
	} else {
		printf("COMPOSITE\n");
	}
	// // bui n = bui_from_dec("3");
	// auto start = std::chrono::high_resolution_clock::now();
	// b = poly_pow_1x(n);
	// auto end = std::chrono::high_resolution_clock::now();
	// auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
	// // Print results
	// std::cout << "Dur: " << duration.count() / N << " ns per run\n";
	// printBuiA(b.data(), b.size());
	// start = std::chrono::high_resolution_clock::now();
	// b = poly_pow_1x_mont(n);
	// end = std::chrono::high_resolution_clock::now();
	// duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
	// // Print results
	// std::cout << "Dur: " << duration.count() / N << " ns per run\n";
	// printBuiA(b.data(), b.size());
	//
	// for (int i = 1; i <= 7; ++i) {
	// 	printf("%2d: ", i);
	// 	Poly poly = poly_pow_1x(bui_from_u32(i));
	// 	printBuiA(poly.data(), poly.size());
	// }
	return 0;
}

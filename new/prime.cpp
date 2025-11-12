#define BI_BIT 512
#include <iostream>

#include "bigint.h"

constexpr u32 POLY_R = 14;
struct Poly : std::array<bui, POLY_R>{};

// inline bui mod(const bul &x, const bui &m) {
// 	bui r = bul_high(x);
// 	bui low = bul_low(x);
//
// 	for (u32 i = 0; i < BI_N; ++i) {
// 		// shift r left by 32b
// 		for (int j = BI_N - 1; j > 0; --j) r[j] = r[j-1];
// 		r[BI_N - 1] = low[i];
// 		if (cmp(r, m) >= 0) {
// 			sub_ip(r, m);
// 		}
// 	}
// 	while (cmp(r, m) >= 0) {
// 		sub_ip(r, m);
// 	}
// 	return r;
// }

Poly poly_pow_1x(const bui &n) {
	Poly base{};
	base[0] = bui1(); base[1] = bui1();
	Poly res{};
	res[0] = bui1();
	return res;
}

// bool is_prime(const bui &x) {
// 	if (cmp(x, bui1()) <= 0) return false;
//
// }

int main() {
	bui a = bui_from_dec("98765432123456789");
	bui b = bui_from_dec("123456789");
	printf("A = %s\n", bui_to_dec(a).c_str());
	printf("B = %s\n", bui_to_dec(b).c_str());
	bui q, r;

	// divmod_fast(a, b, &q, &r);
	divmod_limb(a, b, q, r);

	printf("Q = %s\n", bui_to_dec(q).c_str());
	printf("R = %s\n", bui_to_dec(r).c_str());
}

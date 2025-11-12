#include <iostream>
#include <iomanip>
#include <chrono>  // For benchmarking

// #define BI_N (16384 / 32)
// #define BI_N (32768 / 32)
// #define BI_N (512 / 32)
// fail with N odd
#include "bigint.h"

int main() {
	// bui a = bui_from_u32(2);
	// bui b = bui_from_u32(3);
	// bui a = bui_from_dec("340282366920938463463374607461832982528");
	// bui b = bui_from_dec("340282366920938463463374607461832981000");
	bui a = bui_from_dec("115792089237316195423570985008687907853269984665640564039457584007913129639935");
	bui b = bui_from_dec("115792089237316195423570985008687907853269984665640564039457584007913129639900");

	std::cout << "Initial Values:" << '\n';
	std::cout << "a = " << bui_to_dec(a) << '\n';
	std::cout << "b = " << bui_to_dec(b) << '\n';

	// Benchmarking the three multiplication methods:

	// 1. Test Multiplication: mul(a, b);
	auto start = std::chrono::high_resolution_clock::now();
	bul p = mul(a, b);
	auto end = std::chrono::high_resolution_clock::now();
	auto mul_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
	bui low = bul_low(p);
	std::cout << "\nSchoolbook:" << '\n';
	std::cout << "a * b = " << bui_to_dec(low) << '\n';
	std::cout << mul_duration.count() << " ns" << '\n';

	// 3. Test Multiplication: karatsuba(a, b, BI_N);
	start = std::chrono::high_resolution_clock::now();
	bul p3 = karatsuba_be_top(a, b);
	end = std::chrono::high_resolution_clock::now();
	mul_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
	low = bul_low(p3);
	std::cout << "\nKaratsuba:" << '\n';
	std::cout << "a * b = " << bui_to_dec(low) << '\n';
	std::cout << mul_duration.count() << " ns" << '\n';

	return 0;
}

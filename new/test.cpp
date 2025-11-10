#include <iostream>
#include <iomanip>
#include <chrono>  // For benchmarking

#define BI_N (1024 / 32)
#include "bigint.h"

int main() {
	// bui a = bui_from_u32(17);
	// bui b = bui_from_u32(42);
	bui a = bui_from_dec("115792089237316195423570985008687907853269984665640564039457584007913129639936");
	bui b = bui_from_dec("115792089237316195423570985008687907853269984665640564039457584007913129639936");

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
	std::cout << "\nAfter Multiplication (mul):" << '\n';
	std::cout << "a * b = " << bui_to_dec(low) << '\n';
	std::cout << "Time taken by mul: " << mul_duration.count() << " ns" << '\n';

	// 2. Test Multiplication: mul_low(a, b);
	start = std::chrono::high_resolution_clock::now();
	bui p2 = mul_low(a, b);
	end = std::chrono::high_resolution_clock::now();
	mul_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
	std::cout << "\nAfter Multiplication (mul_low):" << '\n';
	std::cout << "a * b = " << bui_to_dec(p2) << '\n';
	std::cout << "Time taken by mul_low: " << mul_duration.count() << " ns" << '\n';

	// 3. Test Multiplication: karatsuba(a, b, BI_N);
	start = std::chrono::high_resolution_clock::now();
	bul p3 = karatsuba_test(a, b);
	end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> karatsuba_duration = end - start;
	low = bul_low(p3);
	std::cout << "\nAfter Multiplication (karatsuba):" << '\n';
	std::cout << "a * b = " << bui_to_dec(low) << '\n';
	std::cout << "Time taken by karatsuba: " << karatsuba_duration.count() << " seconds" << '\n';

	return 0;
}

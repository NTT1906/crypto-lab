#define CIGINT_IMPLEMENTATION
#define CIGINT_STRIP_PREFIX
#define CIGINT_N (512 / 32)
#define CIGINT_IMPLEMENTATION
#include <chrono>
#include <cmath>
#include <random>
#include <iostream>
#include "cigint.h"

using u32 = uint32_t;
static inline uint64_t now_ns() {
	using namespace std::chrono;
	return std::chrono::duration_cast<std::chrono::nanoseconds>(
		std::chrono::steady_clock::now().time_since_epoch()).count();
}

// Cigint cigint_isqrt(Cigint a) {
// 	/* 0,1 fast-path */
// 	if (cigint_is0(a)) return a;
// 	Cigint one = (Cigint){0}; one.data[CIGINT_N - 1] = 1;
// 	if (cigint_cmp(a, one) <= 0) return a;
//
// 	uint bits = cigint_highest_order(a);              /* 1.. */
// 	uint up   = (bits + 1u) / 2u;                     /* ceil(bits/2) */
// 	Cigint low = {0};
// 	Cigint high = cigint_pow2(up);                    /* 2^up */
// 	high = cigint_sub(high, one);                     /* 2^up - 1 */
//
// 	while (cigint_cmp(low, high) < 0) {
// 		/* mid = (low + high + 1) >> 1 */
// 		Cigint mid = cigint_add(low, high);
// 		mid = cigint_add(mid, one);
// 		mid = cigint_shr(mid, 1);
//
// 		Cigint mid2 = cigint_mul(mid, mid);
// 		if (cigint_cmp(mid2, a) <= 0) {
// 			low = mid;
// 		} else {
// 			high = cigint_sub(mid, one);
// 		}
// 	}
// 	return low;
// }

Cigint cigint_isqrt(Cigint a) {
	if (cigint_is0(a)) return a;

	Cigint res = (Cigint){0};            // “trial root × 2”
	u32 bits = cigint_highest_order(a);
	u32 k = bits ? (bits - 1u) : 0u;    // top bit index
	if (k & 1u) k--;                     // make even → largest power of 4 ≤ a
	Cigint bit = cigint_pow2(k);         // starts at 4^⌊(bits-1)/2⌋

	while (!cigint_is0(bit)) {
		Cigint trial = cigint_add(res, bit);
		if (cigint_cmp(a, trial) >= 0) {
			a = cigint_sub(a, trial);
			res = cigint_add(cigint_shr(res, 1), bit);
		} else {
			res = cigint_shr(res, 1);
		}
		bit = cigint_shr(bit, 2);
	}
	return res;
}

Cigint cigint_isqrt3(Cigint a) {
	if (cigint_is0(a)) return a;

	Cigint rem = a;                      // running remainder
	Cigint res = (Cigint){0};            // “trial root × 2”

	u32 bits = cigint_highest_order(a);
	u32 k = bits ? (bits - 1u) : 0u;    // top bit index
	if (k & 1u) k--;                     // make even → largest power of 4 ≤ a
	Cigint bit = cigint_pow2(k);         // starts at 4^⌊(bits-1)/2⌋

	while (!cigint_is0(bit)) {
		Cigint trial = cigint_add(res, bit);
		if (cigint_cmp(rem, trial) >= 0) {
			rem = cigint_sub(rem, trial);
			res = cigint_add(cigint_shr(res, 1), bit);
		} else {
			res = cigint_shr(res, 1);
		}
		bit = cigint_shr(bit, 2);
	}
	return res;
}

// Cigint cigint_isqrt2(Cigint a) {
// 	return 1 << cigint_highest_order(a) / 2;
// }

Cigint cigint_isqrt2(Cigint a) {
	if (cigint_is0(a)) return a;

	// Initial rough guess: 2^(highest_order/2)
	u32 bits = cigint_highest_order(a);
	Cigint x = cigint_pow2(bits / 2);

	// Refine with Newton–Raphson iteration:
	// x_{n+1} = (x + a/x) / 2
	for (int i = 0; i < 10; i++) {
		Cigint q = cigint_div(a, x);
		Cigint sum = cigint_add(x, q);
		Cigint newx = cigint_shr(sum, 1);
		if (cigint_cmp(newx, x) == 0) break;
		x = newx;
	}

	// Make sure we don't overshoot
	Cigint x2 = cigint_mul(x, x);
	if (cigint_cmp(x2, a) > 0)
		x = cigint_sub(x, 1u);

	return x;
}

int main() {
	std::mt19937_64 rng(12345);
	auto rand_cig = [&]{
		Cigint x = {0};
		for (size_t i=0;i<CIGINT_N;++i) x.data[i] = (u32)rng();
		x.data[CIGINT_N-1] |= 1u;  // ensure non-zero LSW
		return x;
	};

	const int ITERS = 20000;
	uint64_t t1=0, t2=0, t3=0;

	// warmup + correctness
	// for (int i=0;i<1000;++i){
	for (int i=0;i<10;++i){
		auto x = rand_cig();
		auto r1 = cigint_isqrt(x);
		auto r2 = cigint_isqrt2(x);
		auto r3 = cigint_isqrt3(x);
		// Optional: check consistency
		assert(cigint_cmp(r1,r3)==0);
		assert(cigint_cmp(r1,r2)==0);
	}

	// benchmarking
	for (int i=0;i<ITERS;++i){
		auto x = rand_cig();

		auto ts = now_ns();
		cigint_isqrt(x); t1 += now_ns()-ts;
		ts = now_ns(); cigint_isqrt2(x); t2 += now_ns()-ts;
		ts = now_ns(); cigint_isqrt3(x); t3 += now_ns()-ts;
	}

	std::cout << "avg ns/op (lower is better)\n";
	std::cout << "cigint_isqrt:  " << (t1/ITERS) << "\n";
	std::cout << "cigint_isqrt2: " << (t2/ITERS) << "\n";
	std::cout << "cigint_isqrt3: " << (t3/ITERS) << "\n";

	return 0;
}

// PS D:\code\clion\rsa\pj01> .\bench.exe
// avg ns/op (lower is better)
// cigint_isqrt:  5617300
// cigint_isqrt2: 306424
// cigint_isqrt3: 33874
#include <chrono>
#include <iostream>
#include <string>
#define CIGINT_IMPLEMENTATION
#define CIGINT_STRIP_PREFIX
#define CIGINT_N (512 / 32)
#include <iomanip>

#include "cigint.h"

using namespace std::chrono;

constexpr int ITER = 1000000;

void cigint_fill_random(Cigint* cig) {
	for (size_t i = 0; i < CIGINT_N; i++) {
		cig->data[i] = rand() % 0xFFFFFFFF;
	}
}

// Reusable benchmark template for returning-by-value functions
template<typename Func>
void bench_func(const std::string& name, Func f, const Cigint& a, const Cigint& b) {
	Cigint r = a;
	auto start = high_resolution_clock::now();
	for (int i = 0; i < ITER; i++)
		r = f(a, b);
	auto end = high_resolution_clock::now();

	auto duration = duration_cast<nanoseconds>(end - start).count();
	std::cout << std::left << std::setw(25) << name
			  << " took " << std::right << std::setw(10) << duration << " ns\n";
}

// Benchmark for in-place reference functions
template<typename Func>
void bench_func_ref(const std::string& name, Func f, const Cigint& a, const Cigint& b) {
	Cigint x = a, y = b;
	auto start = high_resolution_clock::now();
	for (int i = 0; i < ITER; i++) {
		Cigint x_copy = x;
		Cigint y_copy = y;
		f(&x_copy, &y_copy);
	}
	auto end = high_resolution_clock::now();

	auto duration = duration_cast<nanoseconds>(end - start).count();
	std::cout << std::left << std::setw(25) << name
			  << " took " << std::right << std::setw(10) << duration << " ns\n";
}

Cigint cigint_add_fast_1(Cigint a, Cigint b) {
	Cigint r = CIGINT_ZERO();
	uint64_t carry = 0;
	for (size_t i = CIGINT_N; i-- > 0;) {
		uint64_t sum = (uint64_t)a.data[i] + (uint64_t)b.data[i] + carry;
		r.data[i] = (uint)sum;
		carry = sum >> 32;
	}
	return r;
}

Cigint cigint_add_fast_2(Cigint a, Cigint b) {
	Cigint r = CIGINT_ZERO();
	uint32_t carry = 0;
	for (size_t i = CIGINT_N; i-- > 0;) {
		uint32_t sum = a.data[i] + b.data[i] + carry;
		r.data[i] = sum;
		carry = sum < a.data[i] || sum < b.data[i];
	}
	return r;
}

// In-place (Reference) Functions
void cigint_add_fast_1_r(Cigint *a, Cigint *b) {
	uint64_t carry = 0;
	for (size_t i = CIGINT_N; i-- > 0;) {
		uint64_t sum = (uint64_t)a->data[i] + (uint64_t)b->data[i] + carry;
		a->data[i] = (uint)sum;
		carry = sum >> 32;
	}
}

void cigint_add_fast_1_a(Cigint *a, Cigint *b) {
	uint64_t sum = 0;
	for (size_t i = CIGINT_N; i-- > 0;) {
		sum = (uint64_t)a->data[i] + (uint64_t)b->data[i] + (sum >> 32);
		a->data[i] = (uint)sum;
	}
}

void cigint_add_fast_1_b(Cigint *a, Cigint *b) {
	for (size_t i = CIGINT_N; i-- > 0;) {
		a->data[i] = (uint) (uint64_t)a->data[i] + (uint64_t)b->data[i] + ((uint64_t) a->data[i + 1] >> 32);
	}
}

void cigint_add_fast_2_r(Cigint *a, Cigint *b) {
	uint32_t carry = 0;
	for (size_t i = CIGINT_N; i-- > 0;) {
		uint32_t sum = a->data[i] + b->data[i] + carry;
		carry = sum < a->data[i] || sum < b->data[i];
		a->data[i] = sum;
	}
}

void cigint_add_fast_2_a(Cigint *a, Cigint *b) {
	uint32_t sum = -1;
	for (size_t i = CIGINT_N; i-- > 0;) {
		sum = a->data[i] + b->data[i] + (sum < a->data[i] || sum < b->data[i]);
		a->data[i] = sum;
	}

}

void cigint_add_fast_2_b(Cigint *a, Cigint *b) {
	uint32_t carry = 0;
	for (size_t i = CIGINT_N; i-- > 0;) {
		a->data[i] += b->data[i] + carry;
		carry = a->data[i] < b->data[i] || (carry && a->data[i] < carry);  // Carry if sum overflows
	}
}

// void cigint_add_fast_2_c(Cigint *a, Cigint *b) {
//     uint32_t carry = 0;
//     for (size_t i = CIGINT_N; i-- > 0;) {
//         uint32_t prev = a->data[i];
//         a->data[i] += b->data[i] + carry;
//         carry = (a->data[i] < prev) || (carry && a->data[i] == prev);
//     }
// }

void cigint_add_fast_2_c(Cigint *a, Cigint *b) {
	for (size_t i = CIGINT_N; i-- > 0;) {
		a->data[i] = (a->data[i] += b->data[i]) + (a->data[i] < b->data[i]);
	}
}

void cigint_add_fast_2_d(Cigint *a, Cigint *b) {
	uint32_t carry = 0;
	for (size_t i = CIGINT_N; i-- > 0;) {
		uint32_t sum = a->data[i] + b->data[i] + carry;
		carry = sum < a->data[i] || sum < b->data[i];
		a->data[i] = sum;
	}
}

Cigint cigint_add_fast_2_d_c(Cigint a, Cigint b) {
	cigint_add_fast_2_d(&a, &b);
	return a;
}

Cigint cigint_add_fast_2_d_o(Cigint a, Cigint b) {
	uint32_t carry = 0;
	for (size_t i = CIGINT_N; i-- > 0;) {
		uint32_t sum = a.data[i] + b.data[i] + carry;
		carry = sum < a.data[i] || sum < b.data[i];
		a.data[i] = sum;
	}
	return a;
}

#include <immintrin.h>

static inline void cigint_add_adc_ip(Cigint* a, const Cigint* b) {
	unsigned char c = 0;
	for (size_t i = CIGINT_N; i-- > 0;) {
		c = _addcarry_u32(c, a->data[i], b->data[i], &a->data[i]);
	}
}

static inline Cigint cigint_add_adc(Cigint x, const Cigint& y) {
	cigint_add_adc_ip(&x, &y);
	return x;
}

void benchmark() {
	Cigint a, b;
	cigint_fill_random(&a);
	cigint_fill_random(&b);
	Cigint r = a;
	cprintf("a = %Cd\n", a);
	cprintf("b = %Cd\n", b);
	cigint_add_fast_1_a(&r, &b);
	cprintf("+ = %Cd\n", r);

	a = 1;
	b = 0xFFFFFFFF;
	r = a;
	r = cigint_add_fast_1(r, b);
	cprintf("1  = %Cd\n", r);
	r = a;
	r = cigint_add_fast_2(r, b);
	cprintf("2  = %Cd\n", r);
	r = a;
	cigint_add_fast_1_a(&r, &b);
	cprintf("1a = %Cd\n", r);
	r = a;
	cigint_add_fast_2_r(&r, &b);
	cprintf("2r = %Cd\n", r);
	r = a;
	cigint_add_fast_2_d(&r, &b);
	cprintf("2d = %Cd\n", r);
	r = a;
	cigint_add_fast_2_b(&r, &b);
	cprintf("2b = %Cd\n", r);

	bench_func("Original", cigint_add, a, b);
	bench_func("Fast_1 (uint64_t)", cigint_add_fast_1, a, b);
	bench_func("Fast_2 (uint32_t carry)", cigint_add_fast_2, a, b);
	bench_func_ref("Fast_1_r", cigint_add_fast_1_r, a, b);
	bench_func_ref("Fast_1_a", cigint_add_fast_1_a, a, b);
	bench_func_ref("Fast_2_r", cigint_add_fast_2_r, a, b);
	bench_func_ref("Fast_2_a", cigint_add_fast_2_a, a, b);
	bench_func_ref("Fast_2_b", cigint_add_fast_2_b, a, b);
	bench_func_ref("Fast_2_c", cigint_add_fast_2_c, a, b);
	bench_func_ref("Fast_2_d", cigint_add_fast_2_d, a, b);
	bench_func("Fast_2_d_c", cigint_add_fast_2_d_c, a, b);
	bench_func("Fast_2_d_o", cigint_add_fast_2_d_o, a, b);
	bench_func_ref("add_adc_ip", cigint_add_adc_ip, a, b);
	bench_func("add_adc", cigint_add_adc, a, b);
}

int main() {
	srand((unsigned int)time(NULL));  // Initialize random number generator
	benchmark();
	return 0;
}
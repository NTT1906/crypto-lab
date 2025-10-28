#define CIGINT_IMPLEMENTATION
#define CIGINT_STRIP_PREFIX
#define CIGINT_N (512 / 32)
#include <chrono>
#include <random>
#include <iostream>
#include "cigint.h"

using namespace std::chrono;

// Fast Addition Approaches

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
		// Update a->data[i] directly and calculate carry in one line
		a->data[i] += b->data[i] + carry;
		carry = (a->data[i] < b->data[i]) || (carry && a->data[i] < carry);  // Carry if sum overflows
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

// Fill Cigint structure with random values
void cigint_fill_random(Cigint* cig) {
	for (size_t i = 0; i < CIGINT_N; i++) {
		cig->data[i] = rand() % 0xFFFFFFFF;
	}
}

// Benchmarking function with high resolution timing (ns)
void benchmark() {
	Cigint a, b;
	cigint_fill_random(&a);
	cigint_fill_random(&b);
	Cigint r = a;
	// cprintf("a = %Cd\n", a);
	// cprintf("b = %Cd\n", b);
	// cigint_add_fast_1_a(&r, &b);
	// cprintf("+ = %Cd\n", r);
	// r = a;
	// cigint_add_fast_2_r(&r, &b);
	// cprintf("+ = %Cd\n", r);
	// r = a;
	// r = cigint_add_fast_2(r, b);
	// cprintf("+ = %Cd\n", r);
	// r = a;
	// cigint_add_fast_2_b(&r, &b);
	// cprintf("+ = %Cd\n", r);
	// r = a;
	// cigint_add_fast_2_c(&r, &b);
	// cprintf("+ = %Cd\n", r);

	a = 1;
	b = 0xFFFFFFFF;
	r = a;
	r = cigint_add_fast_1(r, b);
	cprintf("r = %Cd\n", r);
	r = a;
	r = cigint_add_fast_2(r, b);
	cprintf("r = %Cd\n", r);
	r = a;
	cigint_add_fast_1_a(&r, &b);
	cprintf("r = %Cd\n", r);
	r = a;
	cigint_add_fast_2_r(&r, &b);
	cprintf("r = %Cd\n", r);
	r = a;
	cigint_add_fast_2_d(&r, &b);
	cprintf("r = %Cd\n", r);

	// return;

	using namespace std::chrono;
	const int ITER = 1000000;

	// Benchmark First Approach (uint64_t method)
	auto start_time_1 = high_resolution_clock::now();
	for (int i = 0; i < ITER; i++) {
		Cigint r = cigint_add_fast_1(a, b);
	}
	auto end_time_1 = high_resolution_clock::now();
	auto duration_1 = duration_cast<nanoseconds>(end_time_1 - start_time_1);

	// Benchmark Second Approach (uint32_t carry propagation)
	auto start_time_2 = high_resolution_clock::now();
	for (int i = 0; i < ITER; i++) {
		Cigint r = cigint_add_fast_2(a, b);
	}
	auto end_time_2 = high_resolution_clock::now();
	auto duration_2 = duration_cast<nanoseconds>(end_time_2 - start_time_2);

	// Benchmark Third Approach (Bitwise addition)
	auto start_time_3 = high_resolution_clock::now();
	for (int i = 0; i < ITER; i++) {
		Cigint r = cigint_add(a, b);
	}
	auto end_time_3 = high_resolution_clock::now();
	auto duration_3 = duration_cast<nanoseconds>(end_time_3 - start_time_3);

	// Benchmark First Approach REF (uint64_t method) by reference
	auto start_time_1_r = high_resolution_clock::now();
	for (int i = 0; i < ITER; i++) {
		cigint_add_fast_1_r(&a, &b);
	}
	auto end_time_1_r = high_resolution_clock::now();
	auto duration_1_r = duration_cast<nanoseconds>(end_time_1_r - start_time_1_r);

	// Benchmark First Approach A (uint64_t method) by reference
	auto start_time_1_a = high_resolution_clock::now();
	for (int i = 0; i < ITER; i++) {
		cigint_add_fast_1_a(&a, &b);
	}
	auto end_time_1_a = high_resolution_clock::now();
	auto duration_1_a = duration_cast<nanoseconds>(end_time_1_a - start_time_1_a);

	// Benchmark First Approach B (uint64_t method) by reference
	auto start_time_1_b = high_resolution_clock::now();
	for (int i = 0; i < ITER; i++) {
		cigint_add_fast_1_b(&a, &b);
	}
	auto end_time_1_b = high_resolution_clock::now();
	auto duration_1_b = duration_cast<nanoseconds>(end_time_1_b - start_time_1_b);

	// Benchmark Second Approach REF (uint32_t carry propagation) by reference
	auto start_time_2_r = high_resolution_clock::now();
	for (int i = 0; i < ITER; i++) {
		cigint_add_fast_2_r(&a, &b);
	}
	auto end_time_2_r = high_resolution_clock::now();
	auto duration_2_r = duration_cast<nanoseconds>(end_time_2_r - start_time_2_r);

	// Benchmark Second Approach A (uint32_t carry propagation) by reference
	auto start_time_2_a = high_resolution_clock::now();
	for (int i = 0; i < ITER; i++) {
		cigint_add_fast_2_a(&a, &b);
	}
	auto end_time_2_a = high_resolution_clock::now();
	auto duration_2_a = duration_cast<nanoseconds>(end_time_2_a - start_time_2_a);

	// Benchmark Second Approach B (uint32_t carry propagation) by reference
	auto start_time_2_b = high_resolution_clock::now();
	for (int i = 0; i < ITER; i++) {
		cigint_add_fast_2_b(&a, &b);
	}
	auto end_time_2_b = high_resolution_clock::now();
	auto duration_2_b = duration_cast<nanoseconds>(end_time_2_b - start_time_2_b);

	// Benchmark Third Approach REF (uint32_t carry propagation) by reference
	auto start_time_2_c = high_resolution_clock::now();
	for (int i = 0; i < ITER; i++) {
		cigint_add_fast_2_c(&a, &b);
	}
	auto end_time_2_c = high_resolution_clock::now();
	auto duration_2_c = duration_cast<nanoseconds>(end_time_2_c - start_time_2_c);


	auto start_time_2_d = high_resolution_clock::now();
	for (int i = 0; i < ITER; i++) {
		cigint_add_fast_2_d(&a, &b);
	}
	auto end_time_2_d = high_resolution_clock::now();
	auto duration_2_d = duration_cast<nanoseconds>(end_time_2_d - start_time_2_d);

	auto start_time_2_d_c = high_resolution_clock::now();
	for (int i = 0; i < ITER; i++) {
		cigint_add_fast_2_d_c(a, b);
	}
	auto end_time_2_d_c = high_resolution_clock::now();
	auto duration_2_d_c = duration_cast<nanoseconds>(end_time_2_d_c - start_time_2_d_c);

	auto start_time_2_d_o = high_resolution_clock::now();
	for (int i = 0; i < ITER; i++) {
		cigint_add_fast_2_d_o(a, b);
	}
	auto end_time_2_d_o = high_resolution_clock::now();
	auto duration_2_d_o = duration_cast<nanoseconds>(end_time_2_d_o - start_time_2_d_o);

	// Print results
	std::cout << "First Approach (uint64_t method) took                : " << duration_1.count() << " ns\n";
	std::cout << "Second Approach (uint32_t carry propagation) took    : " << duration_2.count() << " ns\n";
	std::cout << "First Approach Ref (uint64_t method) took            : " << duration_1_r.count() << " ns\n";
	std::cout << "First Approach A   (uint64_t method) took            : " << duration_1_a.count() << " ns\n";
	std::cout << "First Approach B   (uint64_t method) took            : " << duration_1_b.count() << " ns\n";
	std::cout << "Second Approach Ref (uint32_t carry propagation) took: " << duration_2_r.count() << " ns\n";
	std::cout << "Second Approach A   (uint32_t carry propagation) took: " << duration_2_a.count() << " ns\n";
	std::cout << "Second Approach B   (uint32_t carry propagation) took: " << duration_2_b.count() << " ns\n";
	std::cout << "Second Approach C   (uint32_t carry propagation) took: " << duration_2_c.count() << " ns\n";
	// std::cout << "Third Approach (Bitwise addition) took               : " << duration_3.count() << " ns\n";

	std::cout << "Second Approach D                 took               : " << duration_2_d.count() << " ns\n";
	std::cout << "Second Approach Dc                took               : " << duration_2_d_c.count() << " ns\n";
	std::cout << "Second Approach Do                took               : " << duration_2_d_o.count() << " ns\n";
}

int main() {
	srand((unsigned int)time(NULL));  // Initialize random number generator
	benchmark();
	return 0;
}
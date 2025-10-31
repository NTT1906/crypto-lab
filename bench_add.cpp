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

#include <atomic>   // for atomic_signal_fence
using namespace std::chrono;

// observable sink to block DCE
static volatile u32 g_sink = 0;
inline void consume(const Cigint& x) {
	u32 acc = 0;
	for (size_t i = 0; i < CIGINT_N; ++i) acc ^= x.data[i];
	g_sink ^= acc;
}

// by-value functions (any of: (Cigint, Cigint) or (Cigint, const Cigint&))
template<typename F>
void bench_func(const std::string& name, F f, const Cigint& a, const Cigint& b) {
	volatile F fp = f;                   // blocks inlining/const-folding
	Cigint r = a;

	std::atomic_signal_fence(std::memory_order_seq_cst);
	auto t0 = steady_clock::now();
	for (int i = 0; i < ITER; ++i) {
		r = (*fp)(a, b);                 // call through the pointer
		consume(r);                      // make every iter observable
	}
	auto t1 = steady_clock::now();
	std::atomic_signal_fence(std::memory_order_seq_cst);

	auto ns = duration_cast<nanoseconds>(t1 - t0).count();
	std::cout << std::left << std::setw(25) << name
			  << " took " << std::right << std::setw(15) << ns << " ns\n";
}

// in-place functions (any of: (Cigint*, Cigint*) or (Cigint*, const Cigint*))
template<typename F>
void bench_func_ref(const std::string& name, F f, const Cigint& a, const Cigint& b) {
	volatile F fp = f;
	Cigint x = a, y = b;

	std::atomic_signal_fence(std::memory_order_seq_cst);
	auto t0 = steady_clock::now();
	for (int i = 0; i < ITER; ++i) {
		Cigint x_copy = x;               // ensure we do work
		(*fp)(&x_copy, &y);
		consume(x_copy);
	}
	auto t1 = steady_clock::now();
	std::atomic_signal_fence(std::memory_order_seq_cst);

	auto ns = duration_cast<nanoseconds>(t1 - t0).count();
	std::cout << std::left << std::setw(25) << name
			  << " took " << std::right << std::setw(15) << ns << " ns\n";
}

// For functions with signature like:
// void func(const Cigint* lhs, const Cigint* rhs, Cigint* res);
template<typename F>
void bench_func_refex(const std::string& name, F f,
					  const Cigint& a, const Cigint& b, Cigint& r)
{
	volatile F fp = f; // prevent optimization

	std::atomic_signal_fence(std::memory_order_seq_cst);
	auto t0 = steady_clock::now();
	for (int i = 0; i < ITER; ++i) {
		Cigint lhs = a;   // ensure same input every time
		Cigint rhs = b;
		Cigint res = r;   // temporary result buffer
		(*fp)(&lhs, &rhs, &res);
		consume(res);     // prevent DCE (dead code elimination)
	}
	auto t1 = steady_clock::now();
	std::atomic_signal_fence(std::memory_order_seq_cst);

	auto ns = duration_cast<nanoseconds>(t1 - t0).count();
	std::cout << std::left << std::setw(25) << name
			  << " took " << std::right << std::setw(15) << ns << " ns\n";
}

Cigint cigint_add_fast_1(Cigint a, Cigint b) {
	Cigint r = CIGINT_ZERO();
	uint64_t carry = 0;
	for (size_t i = CIGINT_N; i-- > 0;) {
		uint64_t sum = (uint64_t)a.data[i] + (uint64_t)b.data[i] + carry;
		r.data[i] = (u32)sum;
		carry = sum >> 32;
	}
	return r;
}

Cigint cigint_add_fast_2(Cigint a, Cigint b) {
	Cigint r = CIGINT_ZERO();
	u32 carry = 0;
	for (size_t i = CIGINT_N; i-- > 0;) {
		u32 sum = a.data[i] + b.data[i] + carry;
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
		a->data[i] = (u32)sum;
		carry = sum >> SIZEOF_U32;
	}
}

void cigint_add_fast_1_a(Cigint *a, Cigint *b) {
	uint64_t sum = 0;
	for (size_t i = CIGINT_N; i-- > 0;) {
		sum = (uint64_t)a->data[i] + (uint64_t)b->data[i] + (sum >> SIZEOF_U32);
		a->data[i] = (u32)sum;
	}
}

void cigint_add_fast_1_b(Cigint *a, Cigint *b) {
	for (size_t i = CIGINT_N; i-- > 0;) {
		a->data[i] = (u32) (uint64_t)a->data[i] + (uint64_t)b->data[i] + ((uint64_t) a->data[i + 1] >> 32);
	}
}

void cigint_add_fast_2_r(Cigint *a, Cigint *b) {
	u32 carry = 0;
	for (size_t i = CIGINT_N; i-- > 0;) {
		u32 sum = a->data[i] + b->data[i] + carry;
		carry = sum < a->data[i] || sum < b->data[i];
		a->data[i] = sum;
	}
}

void cigint_add_fast_2_a(Cigint *a, Cigint *b) {
	u32 sum = -1;
	for (size_t i = CIGINT_N; i-- > 0;) {
		sum = a->data[i] + b->data[i] + (sum < a->data[i] || sum < b->data[i]);
		a->data[i] = sum;
	}
}

void cigint_add_fast_2_a_i(Cigint *a, Cigint *b) {
	u32 sum = -1;
	u32 carry = 0;
	for (size_t i = CIGINT_N; i-- > 0;) {
		sum = a->data[i] + b->data[i] + carry;
		carry = sum < a->data[i] | sum < b->data[i];
		a->data[i] = sum;
	}
}

void cigint_add_fast_2_b(Cigint *a, Cigint *b) {
	u32 carry = 0;
	for (size_t i = CIGINT_N; i-- > 0;) {
		a->data[i] += b->data[i] + carry;
		carry = a->data[i] < b->data[i] || (carry && a->data[i] < carry);
	}
}

void cigint_add_fast_2_d(Cigint *a, Cigint *b) {
	u32 carry = 0;
	for (size_t i = CIGINT_N; i-- > 0;) {
		u32 sum = a->data[i] + b->data[i] + carry;
		carry = sum < a->data[i] || sum < b->data[i];
		a->data[i] = sum;
	}
}

void cigint_add_fast_2_b_i(Cigint* a, Cigint* b) {
	u32 carry = 0;
	for (size_t i = CIGINT_N; i-- > 0;) {
		u32 x = a->data[i];
		u32 y = b->data[i];

		u32 s  = x + y;          // first add
		u32 c1 = (s < x);        // overflow from x+y?
		u32 t  = s + carry;      // add incoming carry
		u32 c2 = (t < s);        // overflow from +carry?

		a->data[i] = t;
		carry = c1 | c2;              // next carry
	}
}

void cigint_add_fast_2_b_iii(Cigint* a, Cigint* b) {
	u32 carry = 0;
	for (size_t i = CIGINT_N; i-- > 0;) {
		u32 x = a->data[i];
		u32 y = b->data[i];
		u32 s  = x + y;
		u32 t  = s + carry;
		a->data[i] = t;
		carry = s < x | t < s;
	}
}

void cigint_add_fast_2_b_iv(Cigint* a, Cigint* b) {
	u32 carry = 0;
	for (size_t i = CIGINT_N; i-- > 0;) {
		u32 s  = a->data[i] + b->data[i];
		u32 t  = s + carry;
		carry = s < a->data[i] | t < s;
		a->data[i] = t;
	}
}

void cigint_add_fast_2_b_ii(Cigint* a, Cigint* b) {
	u32 carry = 0;
	for (size_t i = CIGINT_N; i-- > 0;) {
		u32 x = a->data[i], y = b->data[i];
		u32 t = x + y + carry;       // full sum
		a->data[i] = t;
		// carry-out of x + y + carry (2's complement identity)
		carry = (x & y | ((x | y) & ~t)) >> 31;
	}
}

void cigint_add_fast_2_b_ptr(Cigint* a, Cigint* b) {
	u32 carry = 0;
	u32*       pa = a->data + CIGINT_N;
	const u32* pb = b->data + CIGINT_N;
	do {
		--pa; --pb;
		u32 x = *pa, y = *pb;
		u32 s = x + y;          u32 c1 = (s < x);
		u32 t = s + carry;      u32 c2 = (t < s);
		*pa = t;
		carry = c1 | c2;
	} while (pa != a->data);
}

void cigint_add_fast_2_b_ptr2(Cigint* a, Cigint* b) {
	u32 carry = 0;
	u32*       pa = a->data + CIGINT_N;
	const u32* pb = b->data + CIGINT_N;
	do {
		--pa; --pb;
		u32 x = *pa, y = *pb;
		u32 s = x + y;
		u32 t = s + carry;
		*pa = t;
		carry = s < x | t < s;
	} while (pa != a->data);
}

void cigint_add_fast_2_c(Cigint *a, Cigint *b) {
    u32 carry = 0;
    for (size_t i = CIGINT_N; i-- > 0;) {
        u32 prev = a->data[i];
        a->data[i] += b->data[i] + carry;
        carry = a->data[i] < prev || (carry && a->data[i] == prev);
    }
}

// void cigint_add_fast_2_c(Cigint *a, Cigint *b) {
// 	for (size_t i = CIGINT_N; i-- > 0;) {
// 		a->data[i] = (a->data[i] += b->data[i]) + (a->data[i] < b->data[i]);
// 	}
// }

Cigint cigint_add_fast_2_d_c(Cigint a, Cigint b) {
	cigint_add_fast_2_d(&a, &b);
	return a;
}

Cigint cigint_add_fast_2_d_o(Cigint a, Cigint b) {
	u32 carry = 0;
	for (size_t i = CIGINT_N; i-- > 0;) {
		u32 sum = a.data[i] + b.data[i] + carry;
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
	// cprintf("a = %Cd\n", a);
	// cprintf("b = %Cd\n", b);
	// cigint_add_fast_1_a(&r, &b);
	// cprintf("+ = %Cd\n", r);
	//
	// // a = 1;
	// a = 0xFFFFFFFF;
	// b = 0xFFFFFFFF;
	// r = a;
	// r = cigint_add_fast_1(r, b);
	// cprintf("1  = %Cd\n", r);
	// r = a;
	// r = cigint_add_fast_2(r, b);
	// cprintf("2  = %Cd\n", r);
	// r = a;
	// cigint_add_fast_1_a(&r, &b);
	// cprintf("1a = %Cd\n", r);
	// r = a;
	// cigint_add_fast_1_r(&r, &b);
	// cprintf("1r = %Cd\n", r);
	// r = a;
	// cigint_add_fast_2_r(&r, &b);
	// cprintf("2r = %Cd\n", r);
	// // r = a;
	// // cigint_add_fast_2_a(&r, &b);
	// // cprintf("2a = %Cd\n", r);
	// r = a;
	// cigint_add_fast_2_a_i(&r, &b);
	// cprintf("2ai= %Cd\n", r);
	// r = a;
	// cigint_add_fast_2_b(&r, &b);
	// cprintf("2b = %Cd\n", r);
	// r = a;
	// cigint_add_fast_2_b_i(&r, &b);
	// cprintf("2bi= %Cd\n", r);
	// r = a;
	// cigint_add_fast_2_b_iii(&r, &b);
	// cprintf("2b3= %Cd\n", r);
	// r = a;
	// cigint_add_fast_2_b_iv(&r, &b);
	// cprintf("2b4= %Cd\n", r);
	// r = a;
	// cigint_add_fast_2_b_ii(&r, &b);
	// cprintf("2b2= %Cd\n", r);
	// r = a;
	// cigint_add_fast_2_b_ptr(&r, &b);
	// cprintf("2b2= %Cd\n", r);
	// r = a;
	// cigint_add_fast_2_b_ptr2(&r, &b);
	// cprintf("2b2= %Cd\n", r);
	// r = a;
	// cigint_add_fast_2_d(&r, &b);
	// cprintf("2d = %Cd\n", r);
	// r = a;
	// cigint_add_fast_2_c(&r, &b);
	// cprintf("2c = %Cd\n", r);

	// a = 1;
	// b = 2;
	// r = cigint_mul(a, b);
	// cprintf("mul = %Cd\n", r);
	// r = a;
	// cigint_mul_ref(&r, &b);
	// cprintf("mul_ref = %Cd\n", r);
	// r = a;
	// a = 0;
	// b = 0xFFFFFFFF;
	// r = cigint_mul(a, b);
	// cprintf("mul = %Cd\n", r);
	// r = a;
	// cigint_mul_ref(&r, &b);
	// cprintf("mul_ref = %Cd\n", r);

	cigint_fill_random(&a);
	cigint_fill_random(&b);
	cprintf("a = %Cd\n", a);
	cprintf("b = %Cd\n", b);

	// bench_func("Original", cigint_add, a, b);
	// bench_func("Fast_1 (u64)", cigint_add_fast_1, a, b);
	// bench_func("Fast_2 (u32)", cigint_add_fast_2, a, b);
	// bench_func_ref("Fast_1_r", cigint_add_fast_1_r, a, b);
	// bench_func_ref("Fast_1_a", cigint_add_fast_1_a, a, b);
	// bench_func_ref("Fast_2_r", cigint_add_fast_2_r, a, b);
	// // bench_func_ref("Fast_2_a", cigint_add_fast_2_a, a, b); S
	// bench_func_ref("Fast_2_a_i", cigint_add_fast_2_a_i, a, b);
	// bench_func_ref("Fast_2_b", cigint_add_fast_2_b, a, b);
	// bench_func_ref("Fast_2_b_i", cigint_add_fast_2_b_i, a, b);
	// bench_func_ref("Fast_2_b_ii", cigint_add_fast_2_b_ii, a, b);
	// bench_func_ref("Fast_2_b_iii", cigint_add_fast_2_b_iii, a, b);
	// bench_func_ref("Fast_2_b_iv", cigint_add_fast_2_b_iv, a, b);
	// bench_func_ref("Fast_2_b_ptr", cigint_add_fast_2_b_ptr, a, b);
	// bench_func_ref("Fast_2_b_ptr2", cigint_add_fast_2_b_ptr, a, b);
	// bench_func_ref("Fast_2_c", cigint_add_fast_2_c, a, b);
	// bench_func_ref("Fast_2_d", cigint_add_fast_2_d, a, b);
	// bench_func("Fast_2_d_c", cigint_add_fast_2_d_c, a, b);
	// bench_func("Fast_2_d_o", cigint_add_fast_2_d_o, a, b);
	// bench_func_ref("add_adc_ip", cigint_add_adc_ip, a, b);
	// bench_func("add_adc", cigint_add_adc, a, b);
	// bench_func("sub1", cigint_sub, a, b);
	// bench_func_ref("sub_ref", cigint_sub_ref, a, b);
	bench_func_ref("add_ref", cigint_add_ref, a, b);
	bench_func_ref("mul_ref", cigint_mul_ref, a, b);
	bench_func_refex("mul_refex", cigint_mul_refex, a, b, r);
	// bench_func("mul", cigint_mul, a, b);
}

int main() {
	srand((unsigned int)time(NULL));  // Initialize random number generator
	benchmark();
	return 0;
}
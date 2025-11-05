#include <chrono>
#include <iostream>
#include <string>
#define CIGINT_IMPLEMENTATION
#define CIGINT_STRIP_PREFIX
#define CIGINT_N (512 / 32)
#include <iomanip>

#include "cigint.h"

using namespace std::chrono;

void cigint_fill_random(Cigint *cig) {
	for (size_t i = 0; i < CIGINT_N; i++) {
		u32 r = (u32)rand();
		r = (r << 16) ^ (u32)rand();
		cig->data[i] = r;
	}
}

void cigint_fill_more_random(Cigint *cig) {
	size_t randomLimb = 1 + rand() % CIGINT_N;
	for (unsigned int & i : cig->data) i = 0;
	for (size_t i = 0; i < randomLimb; i++) cig->data[i] =  (u32)rand() << 16 ^ (u32)rand();
}

#include "benchmark.h"
// observable sink to block DCE
static volatile u32 g_sink = 0;
inline void consume(const Cigint& x) {
	u32 acc = 0;
	for (size_t i = 0; i < CIGINT_N; ++i) acc ^= x.data[i]; g_sink ^= acc;
}

// ------------------------------------------------------------
// 1) by-value functions: Cigint f(Cigint, Cigint)
// ------------------------------------------------------------
template<typename F>
void bench_func(const std::string& name, F f, const Cigint& a, const Cigint& b) {
    volatile F fp = f; // block inlining
	bench(name, [&]{
        Cigint r = (*fp)(a, b);
		consume(r);
	});
}

// ------------------------------------------------------------
// 2) in-place, RHS is effectively const: void f(Cigint*, const Cigint*)
//    (your old bench_func_ref)
// ------------------------------------------------------------
template<typename F>
void bench_func_ref(const std::string& name, F f, const Cigint& a, const Cigint& b) {
    volatile F fp = f;
	bench(name, [&]{
        Cigint x = a;
		Cigint y = b;
        (*fp)(&x, &y);
		consume(x);
	});
}

// ------------------------------------------------------------
// 4) 3-arg: void f(const Cigint*, const Cigint*, Cigint*)
//    (your bench_func_refex)
// ------------------------------------------------------------
template<typename F>
void bench_func_refex(const std::string& name, F f, const Cigint& a, const Cigint& b, const Cigint& r0) {
	volatile F fp = f;
	bench(name, [&]{
		Cigint lhs = a;
		Cigint rhs = b;
		Cigint res = r0;
		(*fp)(&lhs, &rhs, &res);
		consume(res);
	});
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

// #include <immintrin.h>
// static inline void cigint_add_adc_ip(Cigint* a, const Cigint* b) {
// 	unsigned char c = 0;
// 	for (size_t i = CIGINT_N; i-- > 0;) {
// 		c = _addcarry_u32(c, a->data[i], b->data[i], &a->data[i]);
// 	}
// }
//
// static inline Cigint cigint_add_adc(Cigint x, const Cigint& y) {
// 	cigint_add_adc_ip(&x, &y);
// 	return x;
// }

void cigint_add_ref_old(Cigint *lhs, Cigint *rhs) {
	while (!cigint_is0_ref(rhs)) {
		Cigint carry = cigint_and(*lhs, *rhs);
		cigint_xor_ref(lhs, rhs);
		*rhs = cigint_shl(carry, 1);
	}
}

Cigint cigint_add_old(Cigint lhs, Cigint rhs) {
	while (!cigint_is0(rhs)) {
		Cigint carry = cigint_and(lhs, rhs);
		lhs = cigint_xor(lhs, rhs);
		rhs = cigint_shl(carry, 1);
	}
	return lhs;
}

void cigint_mul_ref_4(Cigint *lhs, const Cigint *rhs) {
	u32 carry = 0;

	for (size_t k = 0; k < CIGINT_N; ++k) {
		u64 acc_lo = carry;  /* low 32 of current column */
		u32 acc_hi = 0;      /* high 32 of current column */

		for (size_t i = 0; i <= k; ++i) {
			u32 a = lhs->data[CIGINT_N - 1 - i];
			u32 b = rhs->data[CIGINT_N - 1 - (k - i)];
			u64 p = (u64)a * (u64)b;
			u32 p_lo = (u32)p;
			u32 p_hi = (u32)(p >> 32);

			/* add low */
			u32 old = acc_lo;
			acc_lo = old + p_lo;
			u32 c1 = (acc_lo < old);

			/* add high + carry from low */
			acc_hi = acc_hi + p_hi + c1;
		}

		lhs->data[CIGINT_N - 1 - k] = acc_lo;
		carry = acc_hi;
	}
}

void mul_mode_school(Cigint *lhs, const Cigint *rhs, const Cigint *m) {
    Cigint x = cigint_mod(*lhs, *m);      // a mod m
    Cigint p = CIGINT_ZERO();
    u32 nbits = cigint_highest_order(*rhs);
    for (u32 i = 0; i < nbits; ++i) {
        if (cigint_get_bit_ref(rhs, i)) {
            cigint_add_ref(&p, &x);
            p = cigint_mod(p, *m);
        }
        x = cigint_shl(x, 1);
        x = cigint_mod(x, *m);
    }
    *lhs = p;
}

void mul_mode_school_test(Cigint *lhs, const Cigint *rhs) {
	const Cigint m = cigint_from_u32(23127);
	mul_mode_school(lhs, rhs, &m);
}

u32 mul_mod_u32_school(const Cigint &a, const Cigint &b, u32 m32) {
	// a_mod = a % m32; // rem only (no big-int needed afterward)
	u32 a_mod;
	cigint_sdivmod(a, m32, NULL, &a_mod);
	u32 acc = 0;
	u32 nbits = cigint_highest_order(b);
	for (u32 i = 0; i < nbits; ++i) {
		if (cigint_get_bit_ref(&b, i)) acc = (acc + a_mod) % m32;
		a_mod = (2u * a_mod) % m32;
	}
	return acc;
}

void mul_mode_u32_school_test(Cigint *lhs, const Cigint *rhs) {
	const u32 m32 = 23127;
	u32 r = mul_mod_u32_school(*lhs, *rhs, m32);
	*lhs = cigint_from_u32(r);
}

u32 mul_mod_u32_school_2(const Cigint &a, const Cigint &b, u32 m32) {
	u32 a_mod;
	u32 b_mod;
	cigint_sdivmod(a, m32, NULL, &a_mod);
	cigint_sdivmod(b, m32, NULL, &b_mod);
	return ((u64)a_mod * (u64)b_mod) % m32;
}

void mul_mode_u32_school_2_test(Cigint *lhs, const Cigint *rhs) {
	const u32 m32 = 23127;
	u32 r = mul_mod_u32_school_2(*lhs, *rhs, m32);
	*lhs = cigint_from_u32(r);
}

void cigint_mul_mod_wordwise(Cigint *lhs, const Cigint *rhs, const Cigint *m) {
	Cigint res = CIGINT_ZERO();

	for (size_t k = 0; k < CIGINT_N; ++k) {
		u64 acc_lo = 0;
		u32 acc_hi = 0;

		for (size_t i = 0; i <= k; ++i) {
			u32 a = lhs->data[CIGINT_N - 1 - i];
			u32 b = rhs->data[CIGINT_N - 1 - (k - i)];
			u64 p = (u64)a * (u64)b;

			u32 p_lo = (u32)p;
			u32 p_hi = (u32)(p >> 32);

			u32 old_lo = acc_lo;
			acc_lo += p_lo;
			u32 c1 = (acc_lo < old_lo);

			acc_hi += p_hi + c1;
		}

		res.data[CIGINT_N - 1 - k] = (u32)acc_lo;

		if (k + 1 < CIGINT_N) {
			res.data[CIGINT_N - 2 - k] += acc_hi;
		}
	}

	res = cigint_mod(res, *m);
	*lhs = res;
}

void cigint_mul_ref_2_mod(Cigint *lhs, Cigint *rhs, const Cigint *m) {
	cigint_divmod_ref(lhs, m, NULL, lhs);
	cigint_divmod_ref(rhs, m, NULL, rhs);
	__uint128_t carry = 0;
	for (size_t k = 0; k < CIGINT_N; ++k) {
		__uint128_t acc = carry;
		for (size_t i = 0; i <= k; ++i) {
			acc += (__uint128_t)lhs->data[CIGINT_N - 1 - i]
				 * (__uint128_t)rhs->data[CIGINT_N - 1 - (k - i)];
		}
		carry = acc >> 32;
		lhs->data[CIGINT_N - 1 - k] = (u32)acc;
	}
	cigint_divmod_ref(lhs, m, NULL, lhs);
}

void mul_mode_wordwise_test(Cigint *lhs, Cigint *rhs) {
	Cigint m = cigint_from_u32(23127);
	cigint_mul_ref_2_mod(lhs, rhs, &m);
}

void divmod_ref_test(Cigint *lhs, Cigint *rhs) {
	cigint_divmod_ref(lhs, rhs, rhs, lhs);
}

void pow_mod(Cigint *x, Cigint *e, Cigint *m) {
	Cigint ans = 1;
	while (!cigint_is0_ref(e)) {
		if (cigint_get_bit_ref(e, 0)) {
			cigint_mul_ref_2_mod(&ans, x, m);
		}
		cigint_mul_ref_2_mod(x, x, m);
		*e = cigint_shr(*e, 1);
	}
	*x = ans;
}

void pow_mod_2(Cigint *x, Cigint *e, Cigint *m) {
	Cigint ans = 1;
	u32 hob = cigint_highest_order(*e);
	for (size_t i = 0; i < hob; ++i) {
		if (cigint_get_bit_ref(e, i)) {
			cigint_mul_ref_2_mod(&ans, x, m);
		}
		cigint_mul_ref_2_mod(x, x, m);
	}
	*x = ans;
}

void pow_mod_test(Cigint *lhs, Cigint *rhs) {
	Cigint m = 23127;
	pow_mod(lhs, rhs, &m);
}

void pow_mod_2_test(Cigint *lhs, Cigint *rhs) {
	Cigint m = 23127;
	pow_mod_2(lhs, rhs, &m);
}

void cigint_add_wrap(Cigint *lhs, const Cigint *rhs) {
	u64 carry = 0;
	for (size_t i = CIGINT_N; i-- > 0;) {
		carry = (u64) lhs->data[i] + (u64) rhs->data[i] + carry;
		lhs->data[i] = (u32) carry;
		carry >>= 32;
	}
}

void benchmark() {
	Cigint a, b;
	cigint_fill_more_random(&a);
	cigint_fill_more_random(&b);
	std::cout << "a = " << a << "\n";
	std::cout << "b = " << b << "\n";
	bench_func_ref("add_ref", cigint_add_ref, a, b);
	bench_func_ref("add_wrap", cigint_add_wrap, a, b);
	// bench_func("add", cigint_add, a, b);
	// bench_func_ref("mul_ref", cigint_mul_ref, a, b);
	// bench_func_refex("mul_refex", cigint_mul_refex, a, b, r);
	// bench_func_ref("mul_ref2", cigint_mul_ref2, a, b);
	bench_func_ref("mul_2_ref2_b", cigint_mul_ref2_b, a, b);
	// bench_func_ref("mul_ref3_b", cigint_mul_ref3_b, a, b);
	// bench_func_ref("mul_ref4", cigint_mul_ref_4, a, b);
	bench_func_ref("mul_mode_school", mul_mode_school_test, a, b);
	bench_func_ref("mul_mode_u32_school", mul_mode_u32_school_test, a, b);
	bench_func_ref("mul_mode_wordwise", mul_mode_wordwise_test, a, b);
	bench_func_ref("mul_mode_u32_school_2", mul_mode_u32_school_2_test, a, b);
	bench_func_ref("divmod", divmod_ref_test, a, b);
	bench_func_ref("pow_mod", pow_mod_test, a, b);
	bench_func_ref("pow_mod_2", pow_mod_2_test, a, b);
	auto x = a;
	auto y = b;
	mul_mode_school_test(&x, &y);
	cigint_printf("rc = a * b % 23127 = %Cd\n", x);
	x = a;
	y = b;
	mul_mode_u32_school_test(&x, &y);
	cigint_printf("ru = a * b % 23127 = %Cd\n", x);
	x = a;
	y = b;
	mul_mode_u32_school_2_test(&x, &y);
	cigint_printf("ru2= a * b % 23127 = %Cd\n", x);
	x = a;
	y = b;
	mul_mode_wordwise_test(&x, &y);
	cigint_printf("rw = a * b % 23127 = %Cd\n", x);
	x = a;
	y = b;
	pow_mod_test(&x, &y);
	cigint_printf("p  = a ^ b % 23127 = %Cd\n", x);
	x = a;
	y = b;
	pow_mod_2_test(&x, &y);
	cigint_printf("p2 = a ^ b % 23127 = %Cd\n", x);
}

int main() {
	srand((unsigned int)time(NULL));  // Initialize random number generator
	benchmark();
	return 0;
}

// a = 6285088557015142153805641923869937749706537137630097017862647915330197022639331637859852227293232937634720834295549987992776817383678384221041455562941383
// b = 3917538045356230381857991288974190634479233992622172085610044456706146309726126356257005088871671263022323328662412525459899758020885878247371443781912954
// add_ref                   | first:        140 ns | total_loop:       939610 ns | mean_loop:       9.40 ns | total_body:      2714170 ns | mean_body:      27.14 ns
// add_ref                   | first:         54 ns | total_loop:      1010836 ns | mean_loop:      10.11 ns | total_body:      2878229 ns | mean_body:      28.78 ns
// add_ref_old               | first:        287 ns | total_loop:      6693175 ns | mean_loop:      66.93 ns | total_body:      7280860 ns | mean_body:      72.81 ns
// add_old                   | first:        102 ns | total_loop:      4053913 ns | mean_loop:      40.54 ns | total_body:      6301981 ns | mean_body:      63.02 ns
// add                       | first:        807 ns | total_loop:      1257885 ns | mean_loop:      12.58 ns | total_body:      3485199 ns | mean_body:      34.85 ns
// mul_ref                   | first:        783 ns | total_loop:      4550324 ns | mean_loop:      45.50 ns | total_body:      6882402 ns | mean_body:      68.82 ns
// mul_refex                 | first:        579 ns | total_loop:      4226172 ns | mean_loop:      42.26 ns | total_body:      5553441 ns | mean_body:      55.53 ns
// mul_ref2                  | first:        896 ns | total_loop:      8820169 ns | mean_loop:      88.20 ns | total_body:      9800960 ns | mean_body:      98.01 ns
// mul_2_ref2_b              | first:        505 ns | total_loop:      6378370 ns | mean_loop:      63.78 ns | total_body:      9341184 ns | mean_body:      93.41 ns
// mul_ref3                  | first:        337 ns | total_loop:      7810444 ns | mean_loop:      78.10 ns | total_body:     12243382 ns | mean_body:     122.43 ns
// mul_ref_old               | first:       7485 ns | total_loop:    485269560 ns | mean_loop:    4852.70 ns | total_body:    446230310 ns | mean_body:    4462.30 ns
// mul_old                   | first:       3345 ns | total_loop:    353873425 ns | mean_loop:    3538.73 ns | total_body:    362419415 ns | mean_body:    3624.19 ns


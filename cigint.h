#ifndef CIGINT_H
#define CIGINT_H

#ifdef __cplusplus
#include <cassert>
#include <cstdio>
#include <cstdarg>
#include <cstring>
#include <cstdint>
#include <string>
#else
#include <assert.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <stdint.h>
#endif

#ifndef CIGINT_N
#define CIGINT_N 8
#endif

#ifdef __cplusplus
#define sassert static_assert
#else
#define sassert _Static_assert
#endif
sassert(CIGINT_N > 1, "CIGINT_N > 1");

typedef uint32_t u32;
typedef uint64_t u64;
#define SIZEOF_U32 (8 * sizeof(u32))
#define SIZEOF_U64 (8 * sizeof(u64))
#define SIZEOF_CIGINT (CIGINT_N * SIZEOF_U32)

typedef struct Cigint {
	u32 data[CIGINT_N];

#ifdef __cplusplus
	Cigint() = default;
	Cigint(const Cigint &rhs) : data{} {
		memcpy(this->data, rhs.data, sizeof(rhs.data));
	}
	Cigint &operator=(const Cigint &rhs) {
		if (this != &rhs) {
			memcpy(this->data, rhs.data, sizeof(rhs.data));
		}
		return *this;
	}
	Cigint(u32 rhs) : data{} {
		this->data[CIGINT_N - 1] = rhs;
	}
	[[nodiscard]] std::string toDecStr() const;
	[[nodiscard]] std::string toHexStr(bool uppercase = false) const;
	[[nodiscard]] std::string toBinStr() const;
#endif
} Cigint;

#ifdef __cplusplus
#define CIGINT_ZERO() Cigint{}
#else
#define CIGINT_ZERO() (Cigint){0}
#endif

#ifdef __cplusplus
#define FREF(type) type&
#define CFREF(type) const FREF(type)
#else
#define FREF(type) type
#define CFREF(type) FREF(type)
#endif

u32 cigint_get_bit_ref(const Cigint *a, u32 pos);
u32 cigint_get_bit(CFREF(Cigint) a, u32 pos);
Cigint *cigint_set_bit_ref(Cigint *a, u32 pos, u32 val);
Cigint cigint_set_bit(Cigint a, u32 pos, u32 val);
u32 cigint_print10(CFREF(Cigint) a);
u32 cigint_print2(CFREF(Cigint) a);
u32 cigint_print16(CFREF(Cigint) a);
u32 cigint_print16_upper(CFREF(Cigint) a);
int cigint_is0(CFREF(Cigint) a);
int cigint_cmp(CFREF(Cigint) lhs, CFREF(Cigint) rhs);
u32 cigint_highest_order(CFREF(Cigint) num);

void cigint_and_ref(Cigint *lhs, const Cigint *rhs);
Cigint cigint_and(Cigint lhs, CFREF(Cigint) rhs);
void cigint_not_ref(Cigint *a);
Cigint cigint_not(Cigint a);
void cigint_xor_ref(Cigint *lhs, const Cigint *rhs);
Cigint cigint_xor(Cigint lhs, CFREF(Cigint) rhs);
void cigint_or_ref(Cigint *lhs, const Cigint *rhs);
Cigint cigint_or(Cigint lhs, CFREF(Cigint) rhs);

Cigint cigint_shl(CFREF(Cigint) lhs, u32 amnt);
Cigint cigint_shr(CFREF(Cigint) lhs, u32 amnt);
Cigint cigint_pow2(u32 amnt);

Cigint cigint_from_u32(u32 x);

void cigint_add_ref(Cigint *lhs, const Cigint *rhs);
Cigint cigint_add(Cigint lhs, CFREF(Cigint) rhs);
void cigint_sub_ref(Cigint *lhs, const Cigint *rhs);
Cigint cigint_sub(Cigint lhs, CFREF(Cigint) rhs);
void cigint_mul_ref(Cigint *lhs, const Cigint *rhs);
Cigint cigint_mul(Cigint lhs, CFREF(Cigint) rhs);
void cigint_sqr_ref(Cigint *base);
void cigint_pow_ref(Cigint *base, u32 exp);
Cigint cigint_pow(Cigint base, u32 exp);
void cigint_divmod_ref(const Cigint *lhs, const Cigint *rhs, Cigint *q, Cigint *r);
void cigint_divmod(CFREF(Cigint) lhs, CFREF(Cigint) rhs, Cigint *q, Cigint *r);
Cigint cigint_div(CFREF(Cigint) lhs, CFREF(Cigint) rhs);
Cigint cigint_mod(CFREF(Cigint) lhs, CFREF(Cigint) rhs);
Cigint cigint_sqr(Cigint base);
// u32 cigint_printf(const char *fmt, ...);

Cigint cigint_from_dec(const char *s);
Cigint cigint_from_hex(const char *s);
Cigint cigint_from_bin(const char *s);
size_t cigint_to_dec(const Cigint *x, char *buf, size_t buf_size);
size_t cigint_to_hex(const Cigint *x, char *buf, size_t buf_size, int uppercase);
size_t cigint_to_bin(const Cigint *x, char *buf, size_t buf_size);

#ifdef CIGINT_STRIP_PREFIX
	#define from_u32 cigint_from_u32
	#define get_bit_ref cigint_get_bit_ref
	#define get_bit cigint_get_bit
	#define set_bit cigint_set_bit

	/* not to be confused with C++ macros */
	#define cand cigint_and
	#define cor cigint_or
	#define cxor cigint_xor
	#define cnot cigint_not

	#define shl cigint_shl
	#define shr cigint_shr
	#define highest_order cigint_highest_order
	#define pow2 cigint_pow2
	#define cmp cigint_cmp
	#define is0 cigint_is0

	#define add cigint_add
	#define addr cigint_add_ref
	#define sub cigint_sub
	#define subr cigint_sub_ref
	#define mul cigint_mul
	#define mulr cigint_mul_ref

	#define cdivmod cigint_divmod
	#define cdivmodr cigint_divmod_ref
	#define cdiv cigint_div
	#define cmod cigint_mod

	#define cpow cigint_pow
	#define c cigint_pow_ref
	#define csqr cigint_sqr

	#define print2 cigint_print2
	#define print10 cigint_print10
	#define print16 cigint_print16
	#define print16U cigint_print16_upper
	#define cprintf cigint_printf
#endif

#ifdef CIGINT_IMPLEMENTATION
u32 u1_get_bit(u32 num, u32 pos) { return (num >> pos) & (u32)1; }

u32 u1_set_bit(u32 num, u32 pos, u32 val) {
	if (pos >= SIZEOF_U32) return num;
	u32 mask = (u32)1 << pos;
	return (num & ~mask) | ((val & 1u) ? mask : 0u);
}

// reverse bits by divide-and-conquer
u32 u1_bit_reverse(u32 num) {
	num = ((num & 0x55555555U) << 1) | ((num >> 1) & 0x55555555U);
	num = ((num & 0x33333333U) << 2) | ((num >> 2) & 0x33333333U);
	num = ((num & 0x0F0F0F0FU) << 4) | ((num >> 4) & 0x0F0F0F0FU);
	num = ((num & 0x00FF00FFU) << 8) | ((num >> 8) & 0x00FF00FFU);
	return (num << 16) | (num >> 16);
}

u32 u1_highest_order(u32 num) {
	u32 pos = 0;
	while (num > 0) {
		++pos;
		num >>= 1;
	}
	return pos;
}

u32 u1_get_last_nbits(u32 num, u32 amnt) {
	if (amnt >= SIZEOF_U32) {
		return 0;
	}
	return num & ((1ul << amnt) - 1);
}

u32 cigint_get_bit_ref(const Cigint *a, u32 pos) {
	assert(pos < CIGINT_N * SIZEOF_U32);
	size_t data_index = CIGINT_N - 1 - pos / SIZEOF_U32;
	return u1_get_bit(a->data[data_index], pos % SIZEOF_U32);
}

u32 cigint_get_bit(CFREF(Cigint) a, u32 pos) {
	return cigint_get_bit_ref(&a, pos);
}

Cigint *cigint_set_bit_ref(Cigint *a, u32 pos, u32 val) {
	assert(pos < CIGINT_N * SIZEOF_U32);
	size_t data_index = CIGINT_N - 1 - pos / SIZEOF_U32;
	a->data[data_index] = u1_set_bit(a->data[data_index], pos % SIZEOF_U32, val);
	return a;
}

Cigint cigint_set_bit(Cigint a, u32 pos, u32 val) {
	return *cigint_set_bit_ref(&a, pos, val);
}

void cigint_bit_reverse_ref(Cigint *a) {
	u32 *start = a->data;
	u32 *end = a->data + CIGINT_N - 1;
	while (start < end) {
		u32 tmp = u1_bit_reverse(*start);
		*start++ = u1_bit_reverse(*end);
		*end-- = tmp;
	}
	if (CIGINT_N % 2 != 0) { // reverse the middle one if odd n
		*start = u1_bit_reverse(*start);
	}
}

void cigint_bit_reverse_n_ref(Cigint *a, u32 n) {
	for (u32 i = 0; i < n / 2; ++i) {
		u32 j = n - 1 - i;
		u32 bit_i = cigint_get_bit_ref(a, i);
		u32 bit_j = cigint_get_bit_ref(a, j);
		if (bit_i != bit_j) {
			cigint_set_bit_ref(a, i, bit_j);
			cigint_set_bit_ref(a, j, bit_i);
		}
	}
}

Cigint cigint_bit_reverse_n(FREF(Cigint) a, u32 n) {
	cigint_bit_reverse_n_ref(&a, n);
	return a;
}

Cigint cigint_bit_reverse_high(FREF(Cigint) a) {
	cigint_bit_reverse_n_ref(&a, cigint_highest_order(a));
	return a;
}

// TODO: Write doc
void cigint_and_ref(Cigint *lhs, const Cigint *rhs) {
	for (size_t i = 0; i < CIGINT_N; ++i) {
		lhs->data[i] &= rhs->data[i];
	}
}

// TODO: Write doc
Cigint cigint_and(Cigint lhs, CFREF(Cigint) rhs) {
	cigint_and_ref(&lhs , &rhs);
	return lhs;
}

/**
 * @brief Performs a bitwise OR operation (lhs |= rhs) between two Cigint values.
 * @param lhs Pointer to the left-hand operand; stores the result.
 * @param rhs Pointer to the right-hand operand.
 */
void cigint_or_ref(Cigint *lhs, const Cigint *rhs) {
	for (size_t i = 0; i < CIGINT_N; ++i) {
		lhs->data[i] |= rhs->data[i];
	}
}

/**
 * @brief Returns the bitwise OR (r = lhs | rhs) of two Cigint values.
 * This is a copy wrapper of ::cigint_or()
 * @param lhs Left-hand operand (copied before modification).
 * @param rhs Right-hand operand.
 */
Cigint cigint_or(Cigint lhs, CFREF(Cigint) rhs) {
	cigint_or_ref(&lhs, &rhs);
	return lhs;
}

// TODO: Write doc
void cigint_not_ref(Cigint *a) {
	for (size_t i = 0; i < CIGINT_N; ++i)
		a->data[i] = ~a->data[i];
}

// TODO: Write doc
Cigint cigint_not(Cigint a) {
	cigint_not_ref(&a);
	return a;
}

/**
 * @brief Performs a bitwise XOR operation (lhs ^= rhs) between two Cigint values.
 * @param lhs Pointer to the left-hand operand; stores the result.
 * @param rhs Pointer to the right-hand operand.
 */
void cigint_xor_ref(Cigint *lhs, const Cigint *rhs) {
	for (size_t i = 0; i < CIGINT_N; ++i) {
		lhs->data[i] ^= rhs->data[i];
	}
}

/**
 * @brief Returns the bitwise XOR (r = lhs ^ rhs) of two Cigint values.
 * This is a copy wrapper of ::cigint_xor()
 * @param lhs Left-hand operand (copied before modification).
 * @param rhs Right-hand operand.
 */
Cigint cigint_xor(Cigint lhs, CFREF(Cigint) rhs) {
	cigint_xor_ref(&lhs, &rhs);
	return lhs;
}

/* Two-phase, MSW-first left shift: limb move plus intra-word stitch. */
Cigint cigint_shl(CFREF(Cigint) lhs, u32 amnt) {
	if (amnt == 0) return lhs;

	const u32 limbs = amnt / SIZEOF_U32;
	const u32 bits = amnt % SIZEOF_U32;

	if (limbs >= CIGINT_N) return CIGINT_ZERO(); // everything shifted out
	Cigint r = (Cigint){0};
	// Phase 1: limb-only move (toward MSW)
	for (size_t i = 0; i + limbs < CIGINT_N; ++i)
		r.data[i] = lhs.data[i + limbs];
	// Phase 2: intra-word stitch (only if bits != 0)
	if (bits) {
		u32 snapshot[CIGINT_N];
		// for (size_t i = 0; i < CIGINT_N; ++i) snapshot[i] = r.data[i];
		memcpy(snapshot, r.data, sizeof snapshot);
		for (size_t i = 0; i < CIGINT_N; ++i) {
			const u32 hi = snapshot[i] << bits;
			const u32 lo = i + 1 < CIGINT_N ? snapshot[i + 1] >> (SIZEOF_U32 - bits) : 0u;
			r.data[i] = hi | lo;
		}
	}
	return r;
}

/* Two-phase, MSW-first right shift: limb move plus intra-word stitch */
Cigint cigint_shr(CFREF(Cigint) lhs, u32 amnt) {
	if (amnt == 0) return lhs;
	const size_t limb_off = amnt / SIZEOF_U32;
	const u32 bits = amnt % SIZEOF_U32;
	if (limb_off >= CIGINT_N) return CIGINT_ZERO(); // everything shifted out
	Cigint r = CIGINT_ZERO();
	// Phase 1: limb-only move (toward LSW)
	for (size_t i = CIGINT_N; i-- > limb_off; )
		r.data[i] = lhs.data[i - limb_off];
	// Phase 2: intra-word stitch (only if bits != 0)
	if (bits) {
		u32 snapshot[CIGINT_N];
		memcpy(snapshot, r.data, sizeof snapshot);
		for (size_t i = 0; i < CIGINT_N; ++i) {
			const u32 lo = snapshot[i] >> bits;
			const u32 hi = i > 0 ? snapshot[i - 1] << (SIZEOF_U32 - bits) : 0u;
			r.data[i] = hi | lo;
		}
	}
	return r;
}

u32 cigint_highest_order_ref(const Cigint *num) {
	for (size_t i = 0; i < CIGINT_N; ++i) {
		if (num->data[i] > 0) {
			return u1_highest_order(num->data[i]) + (CIGINT_N - i - 1) * SIZEOF_U32;
		}
	}
	return 0;
}

u32 cigint_highest_order(CFREF(Cigint) num) {
	return cigint_highest_order_ref(&num);
}

u32 cigint_highest_word_ref(const Cigint *num) {
	for (size_t i = 0; i < CIGINT_N; ++i) {
		if (num->data[i] > 0) return CIGINT_N - i - 1;
	}
	return 0;
}

u32 cigint_highest_word(CFREF(Cigint) num) {
	return cigint_highest_word_ref(&num);
}

Cigint cigint_pow2(u32 amnt) {
	assert(amnt < CIGINT_N * SIZEOF_U32);
	Cigint r = CIGINT_ZERO();
	size_t limb = CIGINT_N - 1 - (amnt / SIZEOF_U32);
	r.data[limb] = (u32)(1u << (amnt % SIZEOF_U32));
	return r;
}

int cigint_cmp_ref(const Cigint *lhs, const Cigint *rhs) {
	for (size_t i = 0; i < CIGINT_N; ++i) {
		if (lhs->data[i] != rhs->data[i]) {
			return lhs->data[i] > rhs->data[i] ? 1 : -1;
		}
	}
	return 0;
}

int cigint_cmp(CFREF(Cigint) lhs, CFREF(Cigint) rhs) {
	return cigint_cmp_ref(&lhs, &rhs);
}

int cigint_is0_ref(const Cigint *a) {
	for (size_t i = 0; i < CIGINT_N; ++i)
		if (a->data[i] != 0) return 0;
	return 1;
}

int cigint_is0(CFREF(Cigint) a) {
	return cigint_is0_ref(&a);
}

void cigint_add_ref(Cigint *lhs, const Cigint *rhs) {
	u64 sum = 0;
	for (size_t i = CIGINT_N; i-- > 0;) {
		sum = (u64) lhs->data[i] + (u64) rhs->data[i] + (sum >> SIZEOF_U32);
		lhs->data[i] = (u32) sum;
	}
}

Cigint cigint_add(Cigint lhs, CFREF(Cigint) rhs) {
	cigint_add_ref(&lhs, &rhs);
	return lhs;
}

void cigint_sub_ref(Cigint *lhs, const Cigint *rhs) {
	u64 borrow = 0;
	for (size_t i = CIGINT_N; i-- > 0;) {
		u64 a = lhs->data[i];
		u64 t = (u64)rhs->data[i] + borrow;
		lhs->data[i] = (uint32_t)(a - t);
		borrow = a < t;
	}
}

Cigint cigint_sub(Cigint lhs, CFREF(Cigint) rhs) {
	cigint_sub_ref(&lhs, &rhs);
	return lhs;
}

// void cigint_mul_ref_old(Cigint *lhs, Cigint *rhs) {
// 	Cigint res = CIGINT_ZERO();
// 	while (!cigint_is0_ref(rhs)) {
// 		if (u1_get_bit(rhs->data[CIGINT_N - 1], 0) == 1) {
// 			cigint_add_ref(&res, lhs);
// 		}
// 		*lhs = cigint_shl(*lhs, 1);
// 		*rhs = cigint_shr(*rhs, 1);
// 	}
// 	*lhs = res;
// }

// Cigint cigint_mul_old(Cigint lhs, Cigint rhs) {
// 	Cigint res = CIGINT_ZERO();
// 	while (!cigint_is0(rhs)) {
// 		if (u1_get_bit(rhs.data[CIGINT_N - 1], 0) == 1) {
// 			cigint_add(res, lhs);
// 		}
// 		lhs = cigint_shl(lhs, 1);
// 		rhs = cigint_shr(rhs, 1);
// 	}
// 	return lhs;
// }

void cigint_mul_ref(Cigint *lhs, const Cigint *rhs) {
	Cigint tmp = CIGINT_ZERO();
	u64 carry = 0;

	for (size_t k = 0; k < CIGINT_N; ++k) {
		u64 acc = carry;
		for (size_t i = 0; i <= k; ++i) {
			// printf("acc-b: %llu\n", acc);
			acc += (u64) lhs->data[CIGINT_N - 1 - i] * (u64) rhs->data[CIGINT_N - 1 - (k - i)];
			// printf("acc-a: %llu\n", acc);
		}
		tmp.data[CIGINT_N - 1 - k] = (u32) acc;
		carry = acc >> SIZEOF_U32;
		// printf("Carry: %llu\n", carry);
	}

	if (carry != 0) {
		// printf("Carry: %llu\n", carry);
		// fprintf(stderr, "Cigint overflow in multiplication: result truncated\n");
		// assert(0 && "cigint_mul_ref Cigint overflow – increase CIGINT_N or use dynamic bigint");
	}

	*lhs = tmp;
}

void cigint_mul_ref2(Cigint *lhs, const Cigint *rhs) {
	Cigint tmp = CIGINT_ZERO();
	__uint128_t carry = 0;

	for (size_t k = 0; k < CIGINT_N; ++k) {
		__uint128_t acc = carry;
		for (size_t i = 0; i <= k; ++i) {
			acc += (__uint128_t)lhs->data[CIGINT_N - 1 - i]
				 * (__uint128_t)rhs->data[CIGINT_N - 1 - (k - i)];
		}
		tmp.data[CIGINT_N - 1 - k] = (u32)acc;
		carry = acc >> 32;
	}
	*lhs = tmp;
}

void cigint_mul_ref2_b(Cigint *lhs, const Cigint *rhs) {
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
}

void cigint_mul_ref3(Cigint *lhs, const Cigint *rhs) {
	Cigint tmp = CIGINT_ZERO();
	u64 carry = 0;

	// Loop from least significant position (rightmost index) to most (leftmost)
	// But since data[0] is MSB, we reverse index order properly.
	for (ssize_t k = CIGINT_N - 1; k >= 0; --k) {
		u64 acc = carry;

		// inner loop handles all digit pairs whose product contributes to position k
		for (ssize_t i = CIGINT_N - 1; i >= k; --i) {
			ssize_t j = (CIGINT_N - 1) - (i - k); // ensures i + j = k + (N-1)
			if (j < 0 || j >= (ssize_t) CIGINT_N) continue;
			acc += (u64)lhs->data[i] * (u64)rhs->data[j];
		}

		tmp.data[k] = (u32)acc;
		carry = acc >> 32;
	}

	if (carry != 0) {
		printf("Carry: %lu\n", carry);
		fprintf(stderr, "Cigint overflow in multiplication: result truncated\n");
		assert(0 && "cigint_mul_ref3: Cigint overflow – increase CIGINT_N or use dynamic bigint");
	}

	*lhs = tmp;
}

void cigint_mul_ref3_b(Cigint *lhs, const Cigint *rhs) {
	u64 carry = 0;
	// Loop from the least significant position (rightmost index) to most (leftmost)
	// But since data[0] is MSB, we reverse index order properly.
	for (ssize_t k = CIGINT_N - 1; k >= 0; --k) {
		u64 acc = carry;
		// inner loop handles all digit pairs whose product contributes to position k
		for (ssize_t i = CIGINT_N - 1; i >= k; --i) {
			ssize_t j = (CIGINT_N - 1) - (i - k); // ensures i + j = k + (N-1)
			if (j < 0 || j >= (ssize_t) CIGINT_N) continue;
			acc += (u64)lhs->data[i] * (u64)rhs->data[j];
		}
		lhs->data[k] = (u32)acc;
		carry = acc >> 32;
	}

	if (carry != 0) {
		// fprintf(stderr, "Cigint overflow in multiplication: result truncated\n");
		// assert(0 && "Cigint overflow – increase CIGINT_N or use dynamic bigint");
	}
}

void cigint_mul_refex(const Cigint *lhs, const Cigint *rhs, Cigint *res) {
	u64 carry = 0;
	for (size_t k = 0; k < CIGINT_N; ++k) {
		u64 acc = carry;
		for (size_t i = 0; i <= k; ++i) {
			acc += (u64) lhs->data[CIGINT_N - 1 - i] * (u64) rhs->data[CIGINT_N - 1 - (k - i)];
		}
		res->data[CIGINT_N - 1 - k] = (uint32_t) acc;
		carry = acc >> SIZEOF_U32;
	}
}

Cigint cigint_mul(Cigint lhs, CFREF(Cigint) rhs) {
	cigint_mul_ref(&lhs, &rhs);
	return lhs;
}

void cigint_sqr_ref(Cigint *base) {
	Cigint op = *base;
	cigint_mul_ref(base, &op);
}

Cigint cigint_sqr(Cigint base) {
	cigint_sqr_ref(&base);
	return base;
}

void cigint_pow_ref(Cigint *base, u32 exp) {
	Cigint acc = CIGINT_ZERO();
	acc.data[CIGINT_N - 1] = 1;
	Cigint x = *base;
	while (exp) {
		if (exp % 2 == 1)
			acc = cigint_mul(acc, x);
		x = cigint_mul(x, x);
		exp >>= 1u;
		// if (!exp) break;
	}
	*base = acc;
}

Cigint cigint_pow(Cigint base, u32 exp) {
	cigint_pow_ref(&base, exp);
	return base;
}

Cigint cigint_pow_v2(Cigint *base, u32 exp) {
	Cigint res = CIGINT_ZERO();
	res.data[CIGINT_N - 1] = 1;
	for (size_t i = 0; i < exp; ++i) {
		cigint_mul_ref(&res, base);
	}
	return res;
}

Cigint cigint_pow_v3(Cigint base, u32 exp) {
	if (exp == 0) return cigint_from_u32(1);
	Cigint res = cigint_pow_v3(base, exp / 2);
	Cigint tmp = cigint_mul(res, res);
	// cigint_printf("exp :%d\n", exp);
	// cigint_printf("res :%Cd\n", res);
	// cigint_printf("tmp :%Cd\n", tmp);
	if (exp % 2 == 0) {
		return tmp;
	}
	return cigint_mul(base, tmp);
}

inline Cigint cigint_from_u32(u32 x) {
	Cigint tmp = CIGINT_ZERO();
	tmp.data[CIGINT_N - 1] = x;
	return tmp;
}

/* bitwise restoring long division
 * q,r can be NULL independent
 */
void cigint_divmod_ref(const Cigint *lhs, const Cigint *rhs, Cigint *q, Cigint *r) {
	assert(!cigint_is0_ref(rhs));
	int cmp = cigint_cmp_ref(lhs, rhs);
	if (cmp < 0) {
		if (q) *q = CIGINT_ZERO();
		if (r) *r = *lhs;
		return;
	}
	if (cmp == 0) {
		if (q) {
			Cigint one = CIGINT_ZERO();
			one.data[CIGINT_N - 1] = 1;
			*q = one;
		}
		if (r) *r = CIGINT_ZERO();
		return;
	}

	Cigint quotient = CIGINT_ZERO();
	Cigint rem = CIGINT_ZERO();

	int top = (int) cigint_highest_order_ref(lhs) - 1;
	for (int bit = top; bit >= 0; --bit) {
		u32 carry = cigint_get_bit_ref(lhs, (u32) bit);
		for (size_t i = CIGINT_N; i-- > 0; ) {
			u32 v = rem.data[i];
			rem.data[i] = (v << 1) | carry;
			carry = v >> 31;
		}

		if (cigint_cmp_ref(&rem, rhs) >= 0) {
			cigint_sub_ref(&rem, rhs);
			quotient = cigint_set_bit(quotient, (u32) bit, 1u);
		}
	}

	if (q) *q = quotient;
	if (r) *r = rem;
}

void cigint_divmod(CFREF(Cigint) lhs, CFREF(Cigint) rhs, Cigint *q, Cigint *r) {
	cigint_divmod_ref(&lhs, &rhs, q, r);
}

/* quotient only */
Cigint cigint_div(CFREF(Cigint) lhs, CFREF(Cigint) rhs) {
	Cigint q;
	cigint_divmod_ref(&lhs, &rhs, &q, NULL);
	return q;
}

/* remainder only */
Cigint cigint_mod(CFREF(Cigint) lhs, CFREF(Cigint) rhs) {
	Cigint r;
	cigint_divmod_ref(&lhs, &rhs, NULL, &r);
	return r;
}

/* single-limb divisor, fixed version */
void cigint_sdivmod_ref(const Cigint *lhs, u32 rhs, Cigint *q, u32 *r) {
	assert(rhs != 0);
	u64 rem = 0;
	Cigint quo = CIGINT_ZERO();

	u32 hord = cigint_highest_order_ref(lhs);
	// walk bits from MSB -> LSB
	for (int bit = (int) hord - 1; bit >= 0; --bit) {
		rem = (rem << 1) | cigint_get_bit_ref(lhs, (u32) bit);
		if (rem >= rhs) {
			rem -= rhs;
			quo = cigint_set_bit(quo, (u32)bit, 1u);
		}
	}

	if (q) *q = quo;
	if (r) *r = rem;
}

void cigint_sdivmod(CFREF(Cigint) lhs, const u32 rhs, Cigint *q, u32 *r) {
	cigint_sdivmod_ref(&lhs, rhs, q, r);
}

u32 cigint_print2(CFREF(Cigint) a) {
	u32 counter = printf("0b");
	u32 ho = cigint_highest_order(a);
	if (ho == 0) {
		counter += (u32)putchar('0');
		return counter;
	}
	for (u32 bit = ho; bit-- > 0;) {
		int d = (int) cigint_get_bit_ref(&a, bit);
		counter += printf("%d", d);
	}
	return counter;
}

// u32 cigint_print2(CFREF(Cigint) a) {
// 	u32 counter = printf("0b"), old_counter = counter;
// 	int bit_index = cigint_highest_order(a) - 1;
//
// 	while (bit_index >= 0) {
// 		/* TODO: use %2 */
// 		int digit = cigint_get_bit_ref(&a, bit_index);
// 		counter += printf("%d", digit);
// 		--bit_index;
// 	}
// 	if (counter == old_counter) {
// 		counter += putchar('0');
// 	}
// 	return counter;
// }

// Each Cigint has CIGINT_N * 32 bits.
// Each 10^8 chunk holds ~26.6 bits (log2(10^8) ≈ 26.6). So max chunks = ceil(total_bits / 26.6)
#define CIGINT_PRINT10_CHUNKS ((CIGINT_N * SIZEOF_U32 + 26) / 27)

u32 cigint_print10(CFREF(Cigint) a) {
	if (cigint_is0(a)) {
		putchar('0');
		return 1;
	}
	u32 remainders[CIGINT_PRINT10_CHUNKS];
	size_t count = 0;
	Cigint q, n = a;
	while (!cigint_is0(n)) {
		u32 r;
		cigint_sdivmod(n, 100000000U, &q, &r);
		remainders[count++] = r;
		n = q;
	}

	u32 counter = 0;
	counter += printf("%u", remainders[count - 1]);
	for (size_t i = count - 1; i-- > 0;) {
		counter += printf("%08u", remainders[i]);
	}
	return counter;
}

u32 cigint_print16_impl(CFREF(Cigint) a, int upper) {
	u32 counter = 0;
	counter += printf("0x");
	size_t i = 0;
	while (i < CIGINT_N && a.data[i] == 0u) ++i;
	if (i == CIGINT_N) {
		counter += putchar('0');
		return counter;
	}
	counter += printf(upper ? "%X" : "%x", a.data[i]);
	for (size_t j = i + 1; j < CIGINT_N; ++j) {
		counter += printf(upper ? "%08X" : "%08x", a.data[j]);
	}
	return counter;
}

u32 cigint_print16(CFREF(Cigint) a)       { return cigint_print16_impl(a, 0); }
u32 cigint_print16_upper(CFREF(Cigint) a) { return cigint_print16_impl(a, 1); }

#ifdef __cplusplus
inline u32 cigint_print_spec(char spec, const Cigint &x) {
	switch (spec) {
		case 'b': return cigint_print2(x);
		case 'd': return cigint_print10(x);
		case 'x': return cigint_print16(x);
		case 'X': return cigint_print16_upper(x);
		default:  return 0;
	}
}

inline u32 cigint_printf(const char *fmt) {
	u32 counter = 0;
	while (*fmt) {
		if (*fmt == '%' && *(fmt + 1) == '%') {
			putchar('%');
			++counter; fmt += 2;
		} else {
			putchar(*fmt++);
			++counter;
		}
	}
	return counter;
}

template <typename... Rest>
inline typename std::enable_if<(sizeof...(Rest) >= 0), u32>::type
cigint_printf(const char *fmt, const Cigint &arg, Rest&&... rest) {
	u32 counter = 0;

	while (*fmt) {
		if (*fmt == '%') {
			++fmt;
			// "%%"
			if (*fmt == '%') {
				putchar('%');
				++counter;
				++fmt;
				continue;
			}
			// "%C?"
			if (*fmt == 'C') {
				char spec = *(fmt + 1);
				if (spec == 'b' || spec == 'd' || spec == 'x' || spec == 'X') {
					fmt += 2;
					counter += cigint_print_spec(spec, arg);
					return counter + cigint_printf(fmt, std::forward<Rest>(rest)...);
				}
				putchar('%'); ++counter;
				putchar('C'); ++counter;
				continue;
			}
			putchar('%'); ++counter;
			if (*fmt) {
				putchar(*fmt); ++counter;
				++fmt;
			}
		} else {
			putchar(*fmt++);
			++counter;
		}
	}
	// we reached the end of fmt but still had an arg — just ignore the arg
	return counter;
}
#else
u32 cigint_printf(const char *fmt, ...) {
	u32 counter = 0;
	va_list lst;
	va_start(lst, fmt);
	while (*fmt != '\0') {
		switch (*fmt) {
			case '%': {
				fmt++;
				if (*fmt == '%') {
					putchar('%');
					counter++;
				} else if (*fmt == 'C') {
					char next = *(fmt + 1);
					if (next == 'b' || next == 'd' || next == 'x' || next == 'X') {
						fmt++;
						Cigint num = va_arg(lst, Cigint);
						if (*fmt == 'b')      counter += cigint_print2(num);
						else if (*fmt == 'd') counter += cigint_print10(num);
						else if (*fmt == 'x') counter += cigint_print16(num);
						else                  counter += cigint_print16_upper(num); /* 'X' */
					}
				} else if (*fmt == 'c') {
					int ch = va_arg(lst, int);
					counter += putchar(ch);
				} else if (*fmt == 'd' || *fmt == 'i') {
					int num = va_arg(lst, int);
					counter += printf("%d", num);
				} else if (*fmt == 's') {
					char *str = (char*) va_arg(lst, char*);
					counter += printf("%s", str);
				}
				break;
			}
			default: {
				counter += putchar(*fmt);
			}
		}
		fmt++;
	}

	va_end(lst);
	return counter;
}
#endif

int is_space_c(unsigned char c) {
	return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\f' || c=='\v';
}

int hex_val(unsigned char c) {
	if (c >= '0' && c <= '9') return c - '0';
	if (c >= 'a' && c <= 'f') return 10 + (c - 'a');
	if (c >= 'A' && c <= 'F') return 10 + (c - 'A');
	return -1;
}

Cigint cigint_from_dec(const char *s) {
	assert(s && "cigint_from_dec: null pointer");
	Cigint out = CIGINT_ZERO();

	// skip leading spaces and optional '+'
	while (is_space_c(*s)) ++s;
	if (*s == '+') ++s;
	assert(*s != '-' && "cigint_from_dec: negative not supported");

	// skip leading zeros, underscores, spaces
	const char *p = s;
	int any = 0;
	while (*p == '0' || *p == '_' || is_space_c(*p)) {
		if (*p == '0') any = 1;
		++p;
	}
	Cigint n10 = cigint_from_u32(10u);
	Cigint tmp = CIGINT_ZERO();
	for (; *p; ++p) {
		char c = *p;
		if (c == '_' || is_space_c(c)) continue;
		if (c < '0' || c > '9') break;
		any = 1;
		cigint_mul_ref(&out, &n10);
		tmp.data[CIGINT_N - 1] = c - '0';
		cigint_add_ref(&out, &tmp);
	}

	// trailing junk? allow only spaces/underscores
	while (*p && (is_space_c(*p) || *p == '_')) ++p;
	assert(any && "cigint_from_dec: no digits");
	return out;
}

Cigint cigint_from_hex(const char *s) {
	assert(s && "cigint_from_hex: null pointer");
	Cigint out = CIGINT_ZERO();

	while (is_space_c(*s)) ++s;
	if (*s == '+') ++s;
	assert(*s != '-' && "cigint_from_hex: negative not supported");

	if (s[0] == '0' && (s[1] == 'x' || s[1] == 'X')) s += 2;

	const char *p = s;
	int any = 0;

	// skip leading zeros/sep
	while (*p == '0' || *p == '_' || is_space_c(*p)) {
		if (*p == '0') any = 1;
		++p;
	}

	Cigint n16 = cigint_from_u32(16u);
	Cigint tmp = CIGINT_ZERO();
	for (; *p; ++p) {
		char c = *p;
		if (c == '_' || is_space_c(c)) continue;
		int d = hex_val(c);
		if (d < 0) break;
		any = 1;
		cigint_mul_ref(&out, &n16);
		tmp.data[CIGINT_N - 1] = d;
		cigint_add_ref(&out, &tmp);
	}

	while (*p && (is_space_c(*p) || *p == '_')) ++p;
	assert(any && "cigint_from_hex: no digits");
	return out;
}

Cigint cigint_from_bin(const char *s) {
	assert(s && "cigint_from_bin: null pointer");
	Cigint out = CIGINT_ZERO();
	while (is_space_c(*s)) ++s;
	if (s[0] == '0' && (s[1] == 'b' || s[1] == 'B')) s += 2;
	while (*s == '0' || *s == '_' || is_space_c(*s)) ++s;
	size_t k = 0;
	for (const char *p = s; *p; ++p) {
		if (*p == '0' || *p == '1') {
			cigint_set_bit_ref(&out, k++, *p == '1');
		} else if (*p != '_' && !is_space_c(*p)) {
			break; // stop at the 1st invalid chr
		}
	}
	assert(k && "cigint_from_bin: no digits");
	cigint_bit_reverse_n_ref(&out, k);
	return out;
}

Cigint cigint_from_str(const char *s) {
	assert(s && "cigint_from_str: null pointer");
	const char *p = s;
	while (is_space_c(*p)) ++p;
	if (*p == '+') ++p;
	if (*p == '-') {
		assert(0 && "cigint_from_str: negative not supported");
	}

	if (p[0] == '0' && (p[1] == 'x' || p[1] == 'X')) return cigint_from_hex(p);
	if (p[0] == '0' && (p[1] == 'b' || p[1] == 'B')) return cigint_from_bin(p);
	return cigint_from_dec(p);
}

/* --------------------------------------------------------------------- */ /* HEX → string */
/* returns the number of chars written (not counting '\0') */
/* --------------------------------------------------------------------- */
size_t cigint_to_hex(const Cigint *x, char *buf, size_t buf_size, int uppercase) {
	if (buf_size == 0) return 0;
	size_t pos = 0; /* skip leading zero limbs (MSW first) */
	size_t i = 0;
	while (i < CIGINT_N && x->data[i] == 0u) ++i;
	if (i == CIGINT_N) { /* x == 0 */
		if (buf_size > 1) {
			buf[0] = '0';
			buf[1] = '\0';
			return 1;
		}
		buf[0] = '\0';
		return 0;
	} /* first limb: no zero-padding */
	{
		int written = snprintf(buf + pos, buf_size - pos, uppercase ? "%X" : "%x", x->data[i]);
		if (written <= 0 || (size_t)written >= buf_size - pos) { /* truncated */
			buf[buf_size - 1] = '\0';
			return buf_size - 1;
		}
		pos += (size_t)written;
	} /* remaining limbs: always 8 hex digits */
	for (size_t j = i + 1; j < CIGINT_N; ++j) {
		int written = snprintf(buf + pos, buf_size - pos, uppercase ? "%08X" : "%08x", x->data[j]);
		if (written <= 0 || (size_t)written >= buf_size - pos) {
			buf[buf_size - 1] = '\0';
			return buf_size - 1;
		}
		pos += (size_t)written;
	}
	buf[pos] = '\0';
	return pos;
}

/* ---------------------------------------------------------------------
 * BINARY → string
 * returns the number of chars written (not counting '\0')
 * --------------------------------------------------------------------- */
size_t cigint_to_bin(const Cigint *x, char *buf, size_t buf_size) {
	if (buf_size == 0) return 0;
	u32 ho = cigint_highest_order_ref(x); /* number of bits, 0.. */
	if (ho == 0) { /* x == 0 */
		if (buf_size > 1) {
			buf[0] = '0';
			buf[1] = '\0';
			return 1;
		}
		buf[0] = '\0';
		return 0;
	}
	size_t pos = 0;
	for (u32 bit = ho; bit-- > 0;) {
		if (pos + 1 >= buf_size) {
			buf[buf_size - 1] = '\0';
			return buf_size - 1;
		}
		buf[pos++] = (char)('0' + cigint_get_bit_ref(x, bit));
	}
	buf[pos] = '\0';
	return pos;
}

/* --------------------------------------------------------------------- */ /* DECIMAL → string */
/* uses the same 1e8-chunk trick as your cigint_print10() */ /* returns the number of chars written (not counting '\0') */
/* --------------------------------------------------------------------- */
size_t cigint_to_dec(const Cigint *x, char *buf, size_t buf_size) {
	if (buf_size == 0) return 0;
	if (cigint_is0_ref(x)) {
		if (buf_size > 1) {
			buf[0] = '0';
			buf[1] = '\0';
			return 1;
		}
		buf[0] = '\0';
		return 0;
	} /* max chunks like in header */
	enum { BASE = 100000000U }; /* 1e8 */ /* enough for any CIGINT_N; same formula as in cigint.h */
	u32 rems[(CIGINT_N * 32 + 26) / 27];
	size_t count = 0;
	Cigint n = *x;
	Cigint q;
	while (!cigint_is0_ref(&n)) {
		u32 r;
		cigint_sdivmod_ref(&n, BASE, &q, &r);
		rems[count++] = r;
		n = q;
	}
	size_t pos = 0; /* print most significant chunk without leading zeros */
	{
		char tmp[16];
		int written = snprintf(tmp, sizeof(tmp), "%u", rems[count - 1]);
		if (written <= 0) {
			buf[0] = '\0';
			return 0;
		}
		if ((size_t)written >= buf_size) { /* truncated */
			memcpy(buf, tmp, buf_size - 1);
			buf[buf_size - 1] = '\0';
			return buf_size - 1;
		}
		memcpy(buf, tmp, (size_t)written);
		pos += (size_t)written;
	} /* remaining chunks with leading zeros (always 8 digits) */
	for (size_t i = count - 1; i-- > 0;) {
		if (pos + 8 >= buf_size) { /* we would overflow, truncate */
			size_t remain = buf_size - 1 - pos; /* write as much as we can from 8-digit block */
			char tmp[9];
			snprintf(tmp, sizeof(tmp), "%08u", rems[i]);
			memcpy(buf + pos, tmp, remain);
			pos += remain;
			buf[pos] = '\0';
			return pos;
		} /* safe to write 8 digits */
		snprintf(buf + pos, buf_size - pos, "%08u", rems[i]);
		pos += 8;
	}
	buf[pos] = '\0';
	return pos;
}
#endif

#ifdef __cplusplus
inline Cigint operator+(const Cigint &lhs, const Cigint &rhs) {
	return cigint_add(lhs, rhs);
}

inline const Cigint& operator+=(Cigint &lhs, const Cigint &rhs) {
	cigint_add_ref(&lhs, &rhs);
	return lhs;
}

/* TODO: find an ideal way when calculate with u32 */
inline Cigint operator+(const Cigint &lhs, u32 rhs) {
	Cigint tmp = cigint_from_u32(rhs);
	return cigint_add(lhs, tmp);
}

inline const Cigint& operator+=(Cigint &lhs, u32 rhs) {
	Cigint tmp = cigint_from_u32(rhs);
	cigint_add_ref(&lhs, &tmp);
	return lhs;
}

inline Cigint operator-(const Cigint &lhs, const Cigint &rhs) {
	return cigint_sub(lhs, rhs);
}

inline const Cigint &operator-=(Cigint &lhs, const Cigint &rhs) {
	cigint_sub_ref(&lhs, &rhs);
	return lhs;
}

inline Cigint operator*(const Cigint &lhs, const Cigint &rhs) {
	return cigint_mul(lhs, rhs);
}

inline const Cigint &operator*=(Cigint &lhs, const Cigint &rhs) {
	cigint_mul_ref(&lhs, &rhs);
	return lhs;
}

inline Cigint operator/(const Cigint &lhs, const Cigint &rhs) {
	return cigint_div(lhs, rhs);
}

inline const Cigint& operator/=(Cigint &lhs, const Cigint &rhs) {
	lhs = cigint_div(lhs, rhs);
	return lhs;
}

inline bool operator==(const Cigint &lhs, const Cigint &rhs) {
	return cigint_cmp_ref(&lhs, &rhs) == 0;
}

inline bool operator!=(const Cigint &lhs, const Cigint &rhs) {
	return cigint_cmp_ref(&lhs, &rhs) != 0;
}

inline bool operator>(const Cigint &lhs, const Cigint &rhs) {
	return cigint_cmp_ref(&lhs, &rhs) > 0;
}

inline bool operator<(const Cigint &lhs, const Cigint &rhs) {
	return cigint_cmp_ref(&lhs, &rhs) < 0;
}

inline bool operator>=(const Cigint &lhs, const Cigint &rhs) {
	return cigint_cmp_ref(&lhs, &rhs) >= 0;
}

inline bool operator<=(const Cigint &lhs, const Cigint &rhs) {
	return cigint_cmp_ref(&lhs, &rhs) <= 0;
}

/* pick a safe buffer size for decimal:
 * max decimal digits ≈ CIGINT_N * 32 * log10(2) ≈ CIGINT_N * 10
 * so 16 * CIGINT_N is very safe for 512-bit; add a bit more.
 */
// inline
std::string Cigint::toDecStr() const {
	char buf[16 * CIGINT_N + 64];
	cigint_to_dec(this, buf, sizeof(buf));
	return buf;
}

// inline
std::string Cigint::toHexStr(bool uppercase) const {
	/* each limb → up to 8 hex chars; +1 for '\0' */
	char buf[8 * CIGINT_N + 8];
	cigint_to_hex(this, buf, sizeof(buf), uppercase ? 1 : 0);
	return buf;
}

// inline
std::string Cigint::toBinStr() const {
	/* worst case: CIGINT_N*32 bits + 1 */
	char buf[CIGINT_N * 32 + 2];
	cigint_to_bin(this, buf, sizeof(buf));
	return buf;
}

inline std::ostream& operator<<(std::ostream& os, const Cigint& x) {
	return os << x.toDecStr();
}

#endif
#endif
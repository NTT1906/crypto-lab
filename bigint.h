#ifndef _BIGINT_H_
#define _BIGINT_H_
#include <algorithm>
#include <iomanip>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <future>
#include <string>
#include <sstream>
#include <random>
#include <cctype>
#include <cstring>
#include <iostream>
#include <vector>

typedef uint32_t u32;
typedef uint64_t u64;

// #if !defined(DEBUG) && !defined(NDEBUG)
// #define NDEBUG
// #endif

#define SU32 sizeof(u32)
#define SBU32 32

#ifndef BI_BIT
#define BI_BIT 512
#endif
#ifndef BI_N
#define BI_N (BI_BIT / 32)
#endif

#if defined(_MSC_VER)
#define ALWAYS_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#define ALWAYS_INLINE inline __attribute__((always_inline))
#else
#define ALWAYS_INLINE inline
#endif

#if !defined(_MSC_VER) && defined(__cpp_constexpr)
#define OP_CONSTEXPR constexpr
#else
#define OP_CONSTEXPR
#endif

// big endian: data[0] = MSW
// eg: assign 1 to bui: a[BI_N - 1] = 1;
// eg: assign 0x12345678'9ABCDEF0'11223344'55667788 to bui
// a[BI_N - 1] = 0x55667788u;
// a[BI_N - 2] = 0x11223344u;
// a[BI_N - 3] = 0x9ABCDEF0u;
// a[BI_N - 4] = 0x12345678u;
struct bui : std::array<u32, BI_N> {
	using base_type = std::array<u32, BI_N>;
	OP_CONSTEXPR bui() : base_type{} {}
	OP_CONSTEXPR static bui zero() { return {}; }
	OP_CONSTEXPR static bui one() {
		bui r{}; r[BI_N - 1] = 1;
		return r;
	}
	OP_CONSTEXPR static bui from_u32(const u32 x) {
		bui r{}; r[BI_N - 1] = x;
		return r;
	}
};
struct bul : std::array<u32, BI_N * 2> {
	using base_type = std::array<u32, BI_N * 2>;
	OP_CONSTEXPR bul() : base_type{} {}
	OP_CONSTEXPR static bul zero() { return {}; }
	OP_CONSTEXPR static bul one() {
		bul r{}; r[BI_N * 2 - 1] = 1;
		return r;
	}
	OP_CONSTEXPR static bul from_u32(const u32 x) {
		bul r{}; r[BI_N * 2 - 1] = x;
		return r;
	}
};
struct MontgomeryReducer;

std::string bui_to_dec(const bui& x);
std::string bui_to_hex(const bui& a, bool split = false);
bui bui_from_dec(const std::string& s);
bui bui_from_hex(const std::string& s);
bul bui_to_bul(const bui& x);

int cmp(const bui &a, const bui &b);
void add_ip(bui& a, const bui& b);
void add_ip(bul& a, const bul& b);
void sub_ip(bui& a, const bui& b);
void add_mod_ip(bui& a, const bui& b, const bui& m);
void sub_mod_ip(bui& a, const bui& b, const bui& m);
bui mod_native(bui x, const bui& m);
bui mod_native(bul x, const bui& m);
void mul_mod_ip(bui &a, bui b, const bui &m);
void mul_ref(const bui &a, const bui &b, bul &r);
bui bui_pow2(u32 k);
bul bul_pow2(u32 k);
bui bui_binary_flood1(u32 k);
bul bul_binary_flood1(u32 k);
void dbl_ip(bui &x);
u32 u32_divmod_bul(const bul &a, u32 d, bul &q);

ALWAYS_INLINE bui bui0() { return bui::zero(); }
ALWAYS_INLINE bui bui1() { return bui::one(); }
ALWAYS_INLINE bui bui_from_u32(const u32 x) { return bui::from_u32(x); }
ALWAYS_INLINE bul bul0() { return bul::zero(); }
ALWAYS_INLINE bul bul1() { return bul::one(); }
ALWAYS_INLINE bul bul_from_u32(const u32 x) { return bul::from_u32(x); }

inline u32 get_bit(const u32 num, const u32 pos) { return num >> pos & 1; }

inline u32 set_bit(const u32 num, const u32 pos, const u32 val) {
	if (pos >= 32) return num;
	u32 mask = (u32)1 << pos;
	return (num & ~mask) | (val & 1u ? mask : 0u);
}

inline u32 get_bit(const bui &a, const u32 pos) {
	assert(pos < BI_N * SBU32);
	u32 k = BI_N - 1 - pos / SBU32;
	return get_bit(a[k], pos % SBU32);
}

// set in-place
inline void set_bit_ip(bui &a, const u32 pos, const u32 val) {
	assert(pos < BI_N * SBU32 && "Cannot set bit outside the scope of the big integer");
	u32 k = BI_N - 1 - pos / SBU32;
	a[k] = set_bit(a[k], pos % 32, val);
}

inline void set_bit_ip(bul &a, const u32 pos, const u32 val) {
	assert(pos < BI_N * 2 * SBU32 && "Cannot set bit outside the scope of the big integer");
	u32 k = BI_N * 2 - 1 - pos / SBU32;
	a[k] = set_bit(a[k], pos % 32, val);
}

inline bui set_bit(bui a, const u32 pos, const u32 val) {
	set_bit_ip(a, pos, val);
	return a;
}

inline u32 highest_bit(u32 x) {
#if defined(__GNUC__) || defined(__clang__)
	if (x == 0) return 0;
	return SBU32 - __builtin_clz(x);
#else
	u32 pos = 0;
	while (x > 0) {
		++pos;
		x >>= 1;
	}
	return pos;
#endif
}

inline u32 highest_bit(const bui &x) {
	for (size_t i = 0; i < BI_N; ++i) {
		if (x[i] != 0)
			return highest_bit(x[i]) + (BI_N - i - 1) * SBU32;
	}
	return 0; // all limbs zero
}

inline u32 highest_bit(const bul &x) {
	for (size_t i = 0; i < BI_N * 2; ++i) {
		if (x[i] != 0)
			return highest_bit(x[i]) + (BI_N * 2 - i - 1) * SBU32;
	}
	return 0; // all limbs zero
}

inline void bitwise_and_ip(bui &a, const bui &b) {
	for (u32 i = BI_N; i-- > 0;) {
		a[i] &= b[i];
	}
}

// find highest (MSB) limb
inline u32 highest_limb(const bui &x) {
	for (size_t i = 0; i < BI_N; ++i)
		if (x[i] > 0) return BI_N - i - 1;
	return 0;
}

// find highest (MSB) limb
inline u32 highest_limb(const bul &x) {
	for (size_t i = 0; i < BI_N * 2; ++i) {
		if (x[i] > 0)
			return (BI_N * 2 - 1) - i;
	}
	return 0;
}

// Big long: shift left by l whole limbs (each limb is 32 bits) in big‑endian representation.
// Storage is [x[0] = MSW, ..., x[2*BI_N-1] = LSW].
// eg: n = 5, l = 1
//   before: index  0    1    2    3    4
//           value  a0   a1   a2   a3   a4
//   after:         a1   a2   a3   a4   0   // multiplied by 2^(32*l)
inline void shift_limb_left(bul &x, const u32 l) {
	if (l == 0) return;
	if (l >= BI_N * 2) {
		std::fill(x.begin(), x.end(), 0);
		return;
	}
	std::copy(x.begin() + l, x.end(), x.begin());
	std::fill(x.end() - l, x.end(), 0);
}

// shift left in-place (x *= 2^k)
ALWAYS_INLINE void shift_left_ip_imp(u32 *x, const u32 n, const u32 k) {
	if (k == 0) return;
	const u32 limbs = k / SBU32;
	if (limbs >= n) {
		memset(x, 0, n * SU32);
		return;
	}
	const u32 bits = k % SBU32;
	// limb-only move (toward MSW)
	if (limbs) {
		memmove(x, x + limbs, (n - limbs) * SU32);
		memset(x + n - limbs, 0, limbs * SU32);
	}
	// intra-word stitch (only if bits != 0)
	if (bits) {
		u32 c = 0, i = n;
		while (i-- > 0) {
			u32 tmp = x[i];
			x[i] = tmp << bits | c;
			c = tmp >> (SBU32 - bits);
		}
	}
}

inline void shift_left_ip(bui &x, const u32 k) {
	shift_left_ip_imp(x.data(), BI_N, k);
}

inline void shift_left_ip(bul &x, const u32 k) {
	shift_left_ip_imp(x.data(), BI_N * 2, k);
}

// shift left (r = x * 2^k)
inline bui shift_left(bui x, const u32 k) {
	assert(k < BI_N - 1 && "Cannot shift left by big amount (k > BIN_N - 1)");
	if (k == 0) return x;
	u32 limbs = k / SBU32;
	if (limbs >= BI_N) return {};
	u32 bits = k % SBU32;
	bui r{};
	// limb-only move (toward MSW)
	std::copy(x.begin() + limbs, x.end(), r.begin());
	// intra-word stitch (only if bits != 0)
	if (bits) {
		u32 c = 0, i = BI_N;
		while (i-- > 0) {
			u32 tmp = r[i];
			r[i] = tmp << bits | c;
			c = tmp >> (32 - bits);
		}
	}
	return r;
}

// Experiment: shift left expand from bui to bul (r = x * 2^k)
inline bul shift_left_expand(bui x, const u32 k) {
	assert(k < BI_N * 2 - 1 && "Cannot shift left by big amount (k > 2xBIN_N - 1)");
	if (k == 0) return bui_to_bul(x);
	u32 limbs = k / SBU32;
	if (limbs >= BI_N * 2) return {};
	u32 bits = k % SBU32;
	bul r{};
	// limb-only move (toward MSW)
	std::copy_backward(x.begin() + (limbs > BI_N) * (limbs - BI_N), x.end(), r.begin() + BI_N * 2 - limbs);
	// intra-word stitch (only if bits != 0)
	if (bits) {
		u32 c = 0, i = BI_N * 2;
		while (i-- > 0) {
			u32 tmp = r[i];
			r[i] = tmp << bits | c;
			c = tmp >> (32 - bits);
		}
	}
	return r;
}

// shift left mod (r = x * 2^k mod m)
inline bui shift_left_mod(bui x, const u32 k, const bui& m) {
	assert(k < BI_N * 2 && "Cannot shift left by big amount (k > 2xBIN_N - 1)");
	bul p2 = bul_pow2(k);
	bui p2m = mod_native(p2, m);
	x = mod_native(x, m);
	mul_ref(x, p2m, p2);
	p2m = mod_native(p2, m);
	return p2m;
}

// shift right in-place (x /= 2^k)
ALWAYS_INLINE void shift_right_ip_imp(u32 *x, const u32 n, const u32 k) {
	if (k == 0) return;
	const u32 limbs = k / SBU32;
	if (limbs >= n) {
		memset(x, 0, n * SU32);
		return;
	}
	const u32 bits = k % SBU32;
	// limb-only move (toward MSW)
	if (limbs) {
		memmove(x + limbs, x, (n - limbs) * SU32);
		memset(x, 0, limbs * SU32);
	}
	// intra-word stitch (only if bits != 0)
	if (bits) {
		u32 carry = 0;
		for (u32 i = 0; i < n; ++i) {
			u32 v = x[i];
			u32 new_val = v >> bits | carry;
			carry = v << (SBU32 - bits);
			x[i] = new_val;
		}
	}
}

// Big int: shift right in-place (x /= 2^k)
inline void shift_right_ip(bui& x, const u32 k) {
	shift_right_ip_imp(x.data(), BI_N, k);
}

// Big long: shift right in-place (x /= 2^k)
inline void shift_right_ip(bul& x, const u32 k) {
	shift_right_ip_imp(x.data(), BI_N * 2, k);
}

// Checking if input bigint is zero
inline bool bu_is0(const u32 *x, u32 n) {
	while (n-- > 0)
		if (x[n] != 0) return false;
	return true;
}

// Checking if input bui is zero
inline bool bui_is0(const bui& x) {
	for (const u32 val : x)
		if (val != 0) return false;
	return true;
}

// Return low-part of bul as bui
inline bui bul_low(const bul& x) {
	bui r{};
	std::copy(x.begin() + BI_N, x.end(), r.begin());
	return r;
}

// Return high-part of bul as bui
inline bui bul_high(const bul& x) {
	bui r{};
	std::copy_n(x.begin(), BI_N, r.begin());
	return r;
}

// Return new bul with low-part being input bui x
inline bul bui_to_bul(const bui& x) {
	bul r{};
	std::copy(x.begin(), x.end(), r.begin() + BI_N);
	return r;
}

// Compare between two bui
inline int cmp(const bui &a, const bui &b) {
	for (u32 i = 0; i < BI_N; ++i) {
		if (a[i] != b[i])
			return a[i] > b[i] ? 1 : -1;
	}
	return 0;
}

// Compare between two bul
inline int cmp(const bul &a, const bul &b) {
	for (u32 i = 0; i < BI_N * 2; ++i) {
		if (a[i] != b[i])
			return a[i] > b[i] ? 1 : -1;
	}
	return 0;
}

// Compare between bul and bui
inline int cmp(const bul& a, const bui& b) {
	for (int i = 0; i < BI_N; ++i) {
		if (a[i] != 0) return 1;
	}
	for (int i = 0; i < BI_N; ++i) {
		u32 al = a[BI_N + i], bl = b[i];
		if (al != bl) return al > bl ? 1 : -1;
	}
	return 0;
}

inline void randomize_ip(bui &x) {
	std::random_device rd; std::mt19937 gen(rd());
	std::uniform_int_distribution<u32> dist(0, UINT32_MAX);
	size_t limbs = 1 + gen() % BI_N;
	for (u32 &i : x) i = 0;
	for (size_t i = limbs; i-- > 0;) x[i] = dist(gen);
}

inline bui random_odd() {
	bui x; randomize_ip(x);
	set_bit_ip(x, 0, 1);
	return x;
}

ALWAYS_INLINE u32 add_ip_n_imp(u32* a, const u32* b, u32 n) {
	u32 c = 0;
	while (n-- > 0) {
		u64 s = (u64)a[n] + b[n] + c;
		a[n] = (u32)s;
		c = s >> 32;
	}
	return c;
}

inline void add_ip_n(u32* a, const u32* b, const u32 n) {
	add_ip_n_imp(a, b, n);
}

// a += b;
inline void add_ip(bui& a, const bui& b) {
	add_ip_n_imp(a.data(), b.data(), BI_N);
}

// a += b
inline void add_ip(bul& a, const bul& b) {
	add_ip_n_imp(a.data(), b.data(), BI_N * 2);
}

// r.size = 2n
inline void add_n(const u32* a, const u32* b, u32* r, const u32 n) {
	std::copy_n(a, n, r);
	add_ip_n_imp(r, b, n);
}

// r = a + b
inline bui add(bui a, const bui& b) {
	add_ip(a, b);
	return a;
}

// a = (a + b) % m
inline void add_mod_ip(bui &a, const bui &b, const bui &m) {
	if (add_ip_n_imp(a.data(), b.data(), BI_N) || cmp(a, m) >= 0) {
		sub_ip(a, m);
	}
}

ALWAYS_INLINE u32 sub_ip_n_imp(u32* a, const u32* b, u32 n) {
	u32 borrow = 0;
	while (n-- > 0) {
		u64 d = (u64)a[n] - b[n] - borrow;
		a[n] = (u32)d;
		borrow = d >> 32 & 1; // borrow occurs if 32nd bit is 1
	}
	return borrow;
}

// a -= b; // assume a > b
inline void sub_ip(bui& a, const bui& b) {
	sub_ip_n_imp(a.data(), b.data(), BI_N);
}

inline void sub_ip(bul& a, const bul& b) {
	sub_ip_n_imp(a.data(), b.data(), BI_N * 2);
}

// a -= b; // assume a > b
inline void sub_n(const u32* a, const u32* b, u32* r, u32 n) {
	std::copy_n(a, n, r);
	sub_ip_n_imp(r, b, n);
 }

inline bui sub(bui a, const bui& b) {
	sub_ip(a, b);
	return a;
}

inline void sub_mod_ip(bui& a, const bui& b, const bui& m) {
	if (cmp(a, b) >= 0) {
		sub_ip(a, b);
	} else {
		bui t = m, bb = b;
		sub_ip(bb, a); // bb = b - a
		sub_ip(t, bb); // t  = m - (b - a)
		a = t;
	}
}

ALWAYS_INLINE void mul_imp(const u32* a, const u32* b, u32* r, const u32 n) {
	std::fill_n(r, 2 * n, 0);
	for (u32 i = n; i-- > 0;) {
		if (!a[i]) continue;
		u32 c = 0;
		for (u32 j = n; j-- > 0;) {
			u64 p = (u64)a[i] * b[j] + r[i + j + 1] + c;
			r[i + j + 1] = (u32)p;
			c = p >> 32;
		}
		u32 k = i;
		while (c) {
			u64 s = (u64)r[k] + c;
			r[k] = (u32)s;
			c = s >> 32;
			if (k == 0) break;
			k--;
		}
	}
}

inline void mul_ref(const bui &a, const bui &b, bul &r) {
	mul_imp(a.data(), b.data(), r.data(), BI_N);
}

inline void mul_ip(bui &a, const bui &b) {
	bul r{};
	mul_ref(a, b, r);
	a = bul_low(r);
}

inline bul mul(const bui& a, const bui& b) {
	bul r{};
	mul_ref(a, b, r);
	return r;
}

inline bui mul_low(const bui& a, const bui& b) {
	bul r = mul(a, b);
	return bul_low(r);
}

inline bui mul_low_fast(const bui& a, const bui& b) {
	bui r{};
	for (u32 i = 0; i < BI_N; ++i) {
		if (!a[BI_N - 1 - i]) continue;
		u32 c = 0;
		for (u32 j = 0; j < BI_N; ++j) {
			if (i + j >= BI_N) continue;
			u64 p = (u64)a[BI_N - 1 - i] * b[BI_N - 1 - j] + r[BI_N - 1 - (i + j)] + c;
			r[BI_N - 1 - (i + j)] = (u32)p;
			c = p >> 32;
		}
	}
	return r;
}

inline void mul_mod_ip(bui &a, bui b, const bui &m) {
	a = mod_native(a, m);
	b = mod_native(b, m);
	bul r;
	mul_ref(a, b, r);
	a = mod_native(r, m);
}

inline bui mod_native(bui x, const bui& m) {
	long long shift = (long long) highest_bit(x) - highest_bit(m);
	if (shift < 0) return x;

	for (; shift >= 0; --shift) {
		bui tmp = m;
		shift_left_ip(tmp, shift);
		if (cmp(x, tmp) >= 0)
			sub_ip(x, tmp);
	}
	return x;
}

inline bui mod_native(bul x, const bui& m) {
	long long shift = (long long) highest_bit(x) - highest_bit(m);
	if (shift < 0) return bul_low(x);

	for (; shift >= 0; --shift) {
		bul tmp = bui_to_bul(m);
		shift_left_ip(tmp, shift);
		if (cmp(x, tmp) >= 0)
			sub_ip(x, tmp);
	}
	return bul_low(x);
}

inline void mod_native_ip(bui& x, const bui& m) {
	long long shift = (long long) highest_bit(x) - highest_bit(m);
	if (shift < 0) return;
	for (; shift >= 0; --shift) {
		bui tmp = m;
		shift_left_ip(tmp, shift);
		if (cmp(x, tmp) >= 0)
			sub_ip(x, tmp);
	}
}

inline void mod_native_ip(bul& x, const bui& m) {
	long long shift = (long long) highest_bit(x) - highest_bit(m);
	if (shift < 0) return;

	for (; shift >= 0; --shift) {
		bul tmp = bui_to_bul(m);
		shift_left_ip(tmp, shift);
		if (cmp(x, tmp) >= 0)
			sub_ip(x, tmp);
	}
}

inline void divmod(const bui& a, const bui& b, bui &q, bui &r) {
	q = {};
	r = a;
	long long shift = (long long) highest_bit(a) - highest_bit(b);
	if (shift < 0) return;
	for (; shift >= 0; --shift) {
		bui tmp = b;
		shift_left_ip(tmp, shift);
		if (cmp(r, tmp) >= 0) {
			sub_ip(r, tmp);
			set_bit_ip(q, shift, 1);
		}
	}
}

inline void divmod(const bul& a, const bui& b, bui &q, bul &r) {
	q = {};
	r = a;
	long long shift = highest_bit(a) - highest_bit(b);
	if (shift < 0) return;
	bul bb = bui_to_bul(b);
	for (; shift >= 0; --shift) {
		bul tmp = bb;
		shift_left_ip(tmp, shift);
		if (cmp(r, tmp) >= 0) {
			sub_ip(r, tmp);
			set_bit_ip(q, shift, 1);
		}
	}
}

inline void u32divmod(const bui& a, u32 b, bui& q, u32& r) {
	q = {}, r = 0;
	for (int i = 0; i < BI_N; ++i) {
		u64 dividend = (u64)r << 32 | a[i];
		q[i] = (u32)(dividend / b);
		r = (u32)(dividend % b);
	}
}

// Big int: return 2^k
inline bui bui_pow2(const u32 k) {
	assert(k < BI_BIT && "Input size must be in data range!");
	bui r{};
	set_bit_ip(r, k, 1);
	return r;
}

// Big long: return 2^k
inline bul bul_pow2(const u32 k) {
	assert(k < BI_BIT * 2);
	bul r{};
	set_bit_ip(r, k, 1);
	return r;
}

// Return 2^k - 1, k smaller than BI_BIN
inline bui bui_binary_flood1(const u32 k) {
	bui r{};
	u32 l = k / SBU32;
	u32 b = k % SBU32;
	if (l) std::fill_n(r.data() + BI_N - l, l, 0xffffffff);
	if (b) r[BI_N - 1 - l] = (1 << b) - 1;
	return r;
}

// Return 2^k - 1, k smaller than 2xBI_BIN
inline bul bul_binary_flood1(const u32 k) {
	bul r{};
	u32 l = k / SBU32;
	u32 b = k % SBU32;
	if (l) std::fill_n(r.data() + BI_N * 2 - l, l, 0xffffffff);
	if (b) r[BI_N * 2 - 1 - l] = (1 << b) - 1;
	return r;
}

// Return pow_mod (x^e % m)
inline bui pow_mod(bui x, const bui& e, const bui &m) {
	bui r = bui1();
	u32 hb = highest_bit(e);
	for (u32 i = 0; i < hb; ++i) {
		if (get_bit(e, i)) {
			mul_mod_ip(r, x, m);
		}
		mul_mod_ip(x, x, m);
	}
	return r;
}

// Knuth Algorithm D
// Donald E. Knuth, The Art of Computer Programming, Volume 2: Seminumerical Algorithms
// Section: 4.3.1, Algorithm D (Division of large integers).
// https://skanthak.hier-im-netz.de/division.html
// inline void divmod_knuth(const bui& a, const bui& b, bui& quot, bui& rem) {
//     assert(!bui_is0(b));
//     int cm = cmp(a, b);
//     if (cm < 0) {
//        quot = {};
//        rem = a;
//        return;
//     }
//     if (cm == 0) { // a == b
//        quot = bui1();
//        rem = {};
//        return;
//     }
//
//     // normalize
//     bul r = bui_to_bul(a);
//     bui d = b;
//     u32 d_lead_pow = highest_limb(b);
//     u32 d_msw_idx = BI_N - 1 - d_lead_pow;
//     u32 d0 = d[d_msw_idx];
//     const u32 norm_shift = (d0 == 0) ? 0 : 32 - highest_bit(d0);
//
//     if (norm_shift > 0) {
//         shift_left_ip(d, norm_shift);
//         shift_left_ip(r, norm_shift);
//     }
//     // recalculate divisor info after normalization
//     d_lead_pow = highest_limb(d);
//     d_msw_idx = BI_N - 1 - d_lead_pow;
//     d0 = d[d_msw_idx];
//     const u32 d1 = (d_msw_idx + 1 < BI_N) ? d[d_msw_idx + 1] : 0;
//     const u32 n = d_lead_pow + 1; // number of limbs in divisor
//
//     // init
//     quot = {};
//     u32 r_lead_pow = highest_limb(r); // use bul highest_limb
//     const int m = (int)r_lead_pow - (int)d_lead_pow;
//
//     // *Fix*: Handle n=1 (short division) by calling u32divmod
//     // Knuth's algorithm (with the d1 correction) requires n > 1.
//     // Your test case is n=1, which was breaking the logic.
//     if (n == 1) {
//         u32 r_denorm;
//         // We must divide the *full* normalized 'r' (bul), not just the low part
//         // Since u32divmod only takes a bui, we have to simulate it on the bul
//         bul q_bul = {}; // quotient (bul)
//         u32 r_temp = 0; // remainder (u32)
//         u32 d_val = d[BI_N-1]; // The single normalized divisor limb
//
//         for (int i = 0; i < BI_N * 2; ++i) {
//            u64 dividend = (u64)r_temp << 32 | r[i];
//            q_bul[i] = (u32)(dividend / d_val);
//            r_temp = (u32)(dividend % d_val);
//         }
//         // The quotient is in the low half of q_bul
//         std::copy(q_bul.begin() + BI_N, q_bul.end(), quot.begin());
//
//         // Denormalize the remainder
//         // The *un-normalized* remainder is r_temp >> norm_shift
//         // But we must also account for the bits from r that were shifted out
//         if (norm_shift > 0) {
//             shift_right_ip(r, norm_shift);
//             rem = bul_low(r);
//             // The last `norm_shift` bits of the original `a` are in `rem`
//             // and the `r_temp` from the division needs to be shifted back.
//             // This is tricky. Let's just re-calculate the remainder from scratch.
//             bui q_times_b = mul_low(quot, b);
//             rem = sub(a, q_times_b);
//         } else {
//             rem = bui_from_u32(r_temp);
//         }
//         return; // Exit after handling short division
//     }
//
//     // j = Knuth's j (power of current quotient digit), from m down to 0
//     for (int j = m; j >= 0; --j) {
//         // calculate q_hat (guess)
//         // *Fix*: This is the correct, non-sliding index for u_j
//         u32 r_idx = (BI_N * 2 - 1) - (j + n);
//
//         u32 u_jn = r[r_idx];
//         u32 u_jn1 = (r_idx + 1 < BI_N * 2) ? r[r_idx + 1] : 0;
//         u32 u_jn2 = (r_idx + 2 < BI_N * 2) ? r[r_idx + 2] : 0;
//
//         u64 r_top = (u64)u_jn << SBU32 | u_jn1;
//         u64 qhat, rhat;
//
//         if (u_jn >= d0) {
//             qhat = 0xFFFFFFFFULL;
//             rhat = ((u64)(u_jn - d0) << SBU32) + u_jn1 + d0;
//         } else {
//             qhat = r_top / d0;
//             rhat = r_top % d0;
//         }
//
//         // knuth's correction step
//         while (qhat >= (1ULL << 32) || (n > 1 && qhat * d1 > (rhat << 32) + u_jn2)) {
//             qhat--;
//             rhat += d0;
//             if (rhat >= (1ULL << 32)) break;
//         }
//
//         // multiply and subtract (r -= qhat * d * b^j)
//         u64 borrow = 0;
//         u32 d_lsw_idx = BI_N - 1;
//
//         // *Fix*: The window to subtract from starts at r_idx
//         for (u32 i = 0; i < n; ++i) {
//             u32 r_i = r_idx + n - i; // LSW of window is r_idx + n
//             u32 d_i = d_lsw_idx - i;
//             // No OOB check needed now, r_idx is static for the j-loop
//
//             u64 prod = qhat * d[d_i];
//             u64 diff = (u64)r[r_i] - (prod & 0xFFFFFFFF) - borrow;
//             r[r_i] = (u32)diff;
//             borrow = (prod >> 32) - (diff >> 32);
//         }
//
//         u64 final_diff = (u64)r[r_idx] - borrow;
//         r[r_idx] = (u32)final_diff;
//
//         // store quotient digit
//         if (j < BI_N) {
//             u32 q_idx = BI_N - 1 - j;
//             quot[q_idx] = (u32)qhat;
//         }
//
//         // add back (if guess was too high)
//         if (final_diff >> 32 & 1) {
//             if (j < BI_N) {
//                 u32 q_idx = BI_N - 1 - j;
//                 --quot[q_idx];
//             }
//
//             u64 carry = 0;
//             for (u32 i = 0; i < n; ++i) {
//                 u32 r_i = r_idx + n - i;
//                 u32 d_i = d_lsw_idx - i;
//
//                 u64 sum = (u64)r[r_i] + d[d_i] + carry;
//                 r[r_i] = (u32)sum;
//                 carry = sum >> 32;
//             }
//             r[r_idx] = (u32)((u64)r[r_idx] + carry);
//         }
//     }
//
//     // denormalize, remainder is in the low half of 'r'
//     if (norm_shift > 0) {
//         shift_right_ip(r, norm_shift);
//     }
//     rem = bul_low(r); // copy low half of 'r' into 'rem'
// }

ALWAYS_INLINE u32 dbl_ip_imp(bui &x) {
	u32 c = 0, i = BI_N;
	while (i-- > 0) {
		u32 v = x[i];
		u32 nv = v << 1 | c;
		c = v >> 31;
		x[i] = nv;
	}
	return c;
}

// x = 2x (= x << 1)
inline void dbl_ip(bui &x) {
	dbl_ip_imp(x);
}

// x = (2x) % m
static void dbl_mod_ip(bui &x, const bui &m) {
	if (dbl_ip_imp(x) || cmp(x, m) >= 0)
		sub_ip(x, m);
}

inline int hex_val(unsigned char c) {
	if (c >= '0' && c <= '9') return c - '0';
	if (c >= 'a' && c <= 'f') return 10 + (c - 'a');
	if (c >= 'A' && c <= 'F') return 10 + (c - 'A');
	return -1;
}

inline std::string bui_to_dec(const bui& x) {
	if (bui_is0(x)) return "0";
	std::vector<u32> parts;
	bui n = x, q{};
	u32 r = 0;
	while (!bui_is0(n)) {
		u32 BASE = 1000000000u;
		u32divmod(n, BASE, q, r);
		parts.push_back(r);
		n = q;
	}
	std::ostringstream oss;
	oss << parts.back();
	for (int i = (int)parts.size() - 2; i >= 0; --i)
		oss << std::setw(9) << std::setfill('0') << parts[i];
	return oss.str();
}

// Convert bul to decimal string using base 1e9 chunks.
inline std::string bul_to_dec(const bul& x) {
	if (bu_is0(x.data(), BI_N * 2)) return "0";
	std::vector<u32> parts;
	bul n = x, q{};
	while (!bu_is0(n.data(), BI_N * 2)) {
		u32 BASE = 1000000000u;
		u32 r = u32_divmod_bul(n, BASE, q);
		parts.push_back(r);
		n = q;
	}
	std::ostringstream oss;
	oss << parts.back();
	for (int i = (int)parts.size() - 2; i >= 0; --i) {
		oss << std::setw(9) << std::setfill('0') << parts[i];
	}
	return oss.str();
}

inline std::string bui_to_hex(const bui &a, bool split) {
	std::ostringstream o;
	o << std::hex << std::setfill('0');
	for (u32 i = 0; i < BI_N; ++i) {
		o << std::setw(8) << a[i];
		if (split) o << ' ';
	}
	return o.str();
}

static std::string normalize_hex_le_to_be(const std::string& s) {
	std::string hex;
	for (char c : s) {
		if (!isspace((unsigned char)c)) hex.push_back(c);
	}
	if (hex.empty()) return "0";
	reverse(hex.begin(), hex.end());
	return hex;
}

inline bui read_bui_le() {
	std::string line;
	getline(std::cin, line);
	std::string be_hex = normalize_hex_le_to_be(line);
	return bui_from_hex(be_hex);
}

// Divide a double-width big-int (bul, MSW at index 0) by a 32-bit divisor.
// q := a / d (quotient), returns remainder r = a % d.
// Requires: d != 0
inline u32 u32_divmod_bul(const bul &a, u32 d, bul &q) {
	u64 rem = 0;
	for (u32 i = 0; i < BI_N * 2; ++i) q[i] = 0;

	for (u32 i = 0; i < BI_N * 2; ++i) {
		rem = (rem << 32) | (u64)a[i]; // bring down next limb
		// quotient limb fits in 32 bits because rem < d * 2^32 here
		u32 qi = (u32)(rem / d);
		q[i] = qi;
		rem = rem - (u64)qi * (u64)d; // rem = rem % d
	}
	return (u32)rem;
}

// Big int: return bui from dec string
inline bui bui_from_dec(const std::string& s) {
	assert(!s.empty() && "bui_from_dec: empty string");
	u32 i = 0;
	// skip leading spaces and optional '+'
	while (isspace(s[i])) ++i;
	if (s[i] == '+') ++i;
	assert(s[i] != '-' && "bui_from_dec: negative not supported");
	// skip leading zeros, underscores, spaces
	while (s[i] == '0' || s[i] == '_') ++i;
	bool any_digit = false;
	// process each digit in the decimal string
	bui n10 = bui_from_u32(10u);
	bui out{}, tmp{};
	for (; i < s.size(); ++i) {
		char c = s[i];
		if (c == '_' || isspace(c)) continue;
		if (c < '0' || c > '9') break;
		any_digit = true;
		mul_ip(out, n10);
		tmp[BI_N - 1] = c - '0';
		add_ip(out, tmp);
	}
	assert(any_digit && "bui_from_dec: no digits found");
	return out;
}

// Big int: return bui from hex string
inline bui bui_from_hex(const std::string& s) {
	assert(!s.empty() && "bui_from_hex: empty string");
	u32 i = 0;
	// skip leading spaces
	while (i < s.size() && isspace(s[i])) ++i;
	// optional "0x" or "0X" prefix
	if (i + 1 < s.size() && s[i] == '0' && (s[i+1] == 'x' || s[i+1] == 'X')) i += 2;

	bool any_digit = false;
	bui out{};
	bui tmp{};
	bui n16 = bui_from_u32(16u);

	for (; i < s.size(); ++i) {
		char c = s[i];
		if (c == '_' || isspace(c)) continue;
		int val = hex_val(c);
		if (val < 0) break;
		any_digit = true;
		mul_ip(out, n16);
		tmp[BI_N - 1] = (u32)val;
		add_ip(out, tmp);
	}
	assert(any_digit && "bui_from_hex: no digits found");
	return out;
}

// ALWAYS_INLINE void split_bui(const bui &x, bui &high, bui &low, u32 n) {
// 	std::copy_n(x.begin(), n, high.begin() + (BI_N - n));
// 	std::copy(x.begin() + n, x.end(), low.begin() + n);
// }
//
// inline bul karatsuba(const bui &a, const bui &b, const u32 n) {
// 	bul r{};
// 	if (n <= 16) return mul(a, b);
// 	u32 half = n / 2;
//
// 	bui a1{}, a0{}, b1{}, b0{};
// 	split_bui(a, a1, a0, half);
// 	split_bui(b, b1, b0, half);
//
// 	bul z2 = karatsuba(a1, b1, half);
// 	bul z0 = karatsuba(a0, b0, half);
//
// 	bui a_sum = add(a1, a0);
// 	bui b_sum = add(b1, b0);
// 	bul z1 = karatsuba(a_sum, b_sum, half);
// 	sub_ip(z1, z2);
// 	sub_ip(z1, z0);
//
// 	shift_limb_left(z2, 2 * half);
// 	shift_limb_left(z1, half);
//
// 	add_ip(r, z0);
// 	add_ip(r, z1);
// 	add_ip(r, z2);
// 	return r;
// }
//
// size_t KARATSUBA_CUTOFF = 2;
// // size_t KARATSUBA_CUTOFF = 4;
// // size_t KARATSUBA_CUTOFF = 8;
// // size_t KARATSUBA_CUTOFF = 16;
// // size_t KARATSUBA_CUTOFF = 32; // tune this experimentally
// inline void karatsuba_be_rec_old(const u32* a, const u32* b, u32* r, const u32 n, u32* scratch) {
//     if (n <= KARATSUBA_CUTOFF) {
//         mul_imp(a, b, r, n);
//         return;
//     }
//
//     u32 half = n / 2;
//     const u32* a1 = a;
//     const u32* a0 = a + half;
//     const u32* b1 = b;
//     const u32* b0 = b + half;
//
//     u32* z2 = r;       // high part
//     u32* z0 = r + n;   // low part
//     u32* z1 = scratch; // middle temp (2*half)
//
// 	const u32 maxlen = std::max(half, half);
//
//     u32* tmp_a = z1 + 2 * maxlen;
//     u32* tmp_b = tmp_a + maxlen;
//     u32* tmp_scratch = tmp_b + maxlen;
//
//     karatsuba_be_rec_old(a0, b0, z0, half, tmp_scratch); // z0 = a0 * b0
//     karatsuba_be_rec_old(a1, b1, z2, half, tmp_scratch); // z2 = a1 * b1
//
// 	add_n(a1, a0, tmp_a, half); // tmp_a = a1 + a0
// 	add_n(b1, b0, tmp_b, half); // tmp_b = b1 + b0
//     karatsuba_be_rec_old(tmp_a, tmp_b, z1, half, tmp_scratch); // z1 = (a1 + a0) * (b1 + b0)
//
//     // z1 = z1 - z2 - z0
// 	sub_n(z1, z2, z1, 2 * half);
// 	sub_n(z1, z0, z1, 2 * half);
//     // combine: r = z2 << (2*half*32) + z1 << (half*32) + z0
// 	add_n(r + half, z1, r + half, 2 * half);
// }
//
// inline void karatsuba_be_rec(const u32* a, const u32* b, u32* r, u32 n, u32* scratch) {
// 	if (bu_is0(a, n) || bu_is0(b, n)) {
// 		std::fill_n(r, 2 * n, 0);
// 		return;
// 	}
// 	if (n <= KARATSUBA_CUTOFF) {
// 		mul_imp(a, b, r, n);
// 		return;
// 	}
//
// 	const u32 half  = n / 2;
// 	const u32 other = n - half;       // may be half+1 if n is odd
//
// 	// Big-endian split
// 	const u32* a1 = a;           // high half
// 	const u32* a0 = a + half;    // low  half
// 	const u32* b1 = b;
// 	const u32* b0 = b + half;
//
// 	// workspace layout
// 	u32* z0 = scratch;                 // size 2*other
// 	u32* z1 = z0 + 2 * other;          // size 2*other
// 	u32* z2 = z1 + 2 * other;          // size 2*half
// 	u32* tmp_a = z2 + 2 * half;
// 	u32* tmp_b = tmp_a + other;
// 	u32* subscratch = tmp_b + other;
//
// 	// z0 = a0 * b0
// 	karatsuba_be_rec(a0, b0, z0, other, subscratch);
//
// 	// z2 = a1 * b1
// 	karatsuba_be_rec(a1, b1, z2, half, subscratch);
//
// 	// tmp_a = a0 + a1 (aligned to low indices)
// 	std::fill_n(tmp_a, other, 0);
// 	std::fill_n(tmp_b, other, 0);
// 	std::copy(a0 + (other - half), a0 + other, tmp_a + (other - half));
// 	for (u32 i = 0; i < half; ++i)
// 		tmp_a[i] += a1[i];
// 	std::copy(b0 + (other - half), b0 + other, tmp_b + (other - half));
// 	for (u32 i = 0; i < half; ++i)
// 		tmp_b[i] += b1[i];
//
// 	// z1 = (a0+a1)*(b0+b1)
// 	karatsuba_be_rec(tmp_a, tmp_b, z1, other, subscratch);
//
// 	// z1 = z1 - z2 - z0
// 	sub_n(z1 + (2 * other - 2 * half), z2, z1 + (2 * other - 2 * half), 2 * half);
// 	sub_n(z1, z0, z1, 2 * other);
//
// 	// clear result
// 	std::fill_n(r, 2 * n, 0);
//
// 	// combine (big-endian)
// 	// copy z0 → low end
// 	std::copy(z0 + 2 * other - n, z0 + 2 * other, r + 2 * n - 2 * other);
//
// 	// add z1 shifted by (other limbs)
// 	add_n(r + (n - other), z1 + 2 * other - n, r + (n - other), 2 * other);
//
// 	// add z2 shifted by (2*other limbs)
// 	add_n(r, z2 + 2 * half - n, r, 2 * half);
// }
//
// inline u32 next_pow2(u32 x) {
// 	if (x == 0) return 1;
// 	x--;
// 	x |= x >> 1;
// 	x |= x >> 2;
// 	x |= x >> 4;
// 	x |= x >> 8;
// 	x |= x >> 16;
// 	x++;
// 	return x;
// }
//
// inline bul karatsuba_be_top(const bui& a, const bui& b) {
//     size_t n = BI_N;
// 	// u32 n = std::max(highest_limb(a), highest_limb(b));
// 	// n = next_pow2(n) * 2;
//     bul r{};
//     std::array<u32, 6 * BI_N> scratch{};
//     // std::array<u32, 8 * BI_N> scratch{};
//     karatsuba_be_rec_old(a.data(), b.data(), r.data(), n, scratch.data());
//     return r;
// }
//
// // for compatibility with your test code
// inline bul karatsu_test(const bui& a, const bui& b) {
//     return karatsuba_be_top(a, b);
// }

// Extended Euclidean algorithm
inline bool mod_inverse_old(bui a, const bui &m, bui &inv_out) {
	// invalid modulus or zero
	if (bui_is0(m)) return false;
	if (cmp(a, m) >= 0) a = mod_native(a, m);
	if (bui_is0(a)) return false; // zero has no inverse
	bui r0 = m, r1 = a, t0{}, t1 = bui1();
	while (!bui_is0(r1)) {
		// q = r0 / r1, rem = r0 % r1
		bui q, rem;
		divmod(r0, r1, q, rem);
		// r0, r1 = r1, rem
		r0 = r1, r1 = rem;
		// t_new = (t0 - q * t1) mod m
		// compute q * t1 -> bul, then reduce modulo m to get r_qt (bui)
		bul prod{};
		mul_ref(q, t1, prod);  // prod = q * t1 (2N words)
		auto qtm_rem = mod_native(prod, m); // qtm_rem = (prod) % m

		// t_new = t0 - qtm_rem mod m
		bui tnew = t0;
		if (cmp(tnew, qtm_rem) >= 0) {
			sub_ip(tnew, qtm_rem);
		} else {
			// tnew = (t0 - qtm_rem) mod m = m - (qtm_rem - t0)
			tnew = m;
			sub_ip(qtm_rem, t0);
			sub_ip(tnew, qtm_rem);
		}
		t0 = t1;
		t1 = tnew;
	}

	// r0 = gcd(a, m) so if gcd != 1 -> no inverse
	if (cmp(r0, bui1()) != 0) return false;
	inv_out = t0;
	return true;
}

// Extended Euclidean algorithm: find a^{-1} mod m
inline bool mod_inverse(const bui& a_in, const bui& m, bui& inv_out) {
	if (bui_is0(m)) return false; // invalid modulus
	bui a = a_in;
	if (cmp(a, m) >= 0) a = mod_native(a, m);      // reduce a mod m
	if (bui_is0(a)) return false; // 0 has no inverse

	bui r0 = m, r1 = a;
	bui t0{}, t1 = bui1();

	while (!bui_is0(r1)) {
		bui q, rem;
		divmod(r0, r1, q, rem); // r0 = q*r1 + rem
		r0 = r1;
		r1 = rem;

		// t_new = (t0 - q * t1) mod m
		bui qt = t1;
		mul_mod_ip(qt, q, m); // qt = (q * t1) mod m

		bui tnew;
		if (cmp(t0, qt) >= 0) {
			tnew = t0;
			sub_ip(tnew, qt); // tnew = t0 - qt
		} else {
			tnew = m;
			sub_ip(qt, t0);   // qt = qt - t0
			sub_ip(tnew, qt); // tnew = m - (qt - t0)
		}

		t0 = t1;
		t1 = tnew;
	}

	if (cmp(r0, bui1()) != 0) return false;        // gcd != 1 -> no inverse
	inv_out = t0;                                  // already in [0, m)
	return true;
}

/**
 * @brief Montgomery modular arithmetic helper for fixed-size big integers.
 *
 * This struct precomputes all constants needed to perform fast modular
 * multiplication and exponentiation modulo an odd big integer using the
 * Montgomery reduction algorithm in the bui/bul representation.
 *
 * Ref:
 * [1] https://en.wikipedia.org/wiki/Montgomery_modular_multiplication
 * [2] https://cp-algorithms.com/algebra/montgomery_multiplication.html
 * [3] https://en.algorithmica.org/hpc/number-theory/montgomery/
 * [4] MVP: https://www.nayuki.io/page/montgomery-reduction-algorithm
 */
struct MontgomeryReducer {
	bui modulus;        // must be odd >= 3
	bul reducer{};      // power of 2
	bui mask{};         // reducer - 1
	u32 reducerBits;    // log2(reducer)
	bui reciprocal{};   // reducer^-1 mod modulus
	bui factor{};       // (reducer * reciprocal - 1) / modulus
	bui convertedOne{}; // convertIn(1) aka reducer mod modulus
	static bui modInverse(const bui& a, const bui& m);

	MontgomeryReducer(const bui& modulus) : modulus(modulus) {
		assert(get_bit(modulus, 0) && cmp(modulus, bui1()) == 1);
		reducerBits = (highest_bit(modulus) / SBU32 + 1) * SBU32;
		if (reducerBits > BI_BIT) reducerBits = BI_BIT;
		reducer = bul_pow2(reducerBits);
		mask = bui_binary_flood1(reducerBits);
		convertedOne = mod_native(reducer, modulus);
		mod_inverse_old(convertedOne, modulus, reciprocal); // reducer^-1 mod modulus

		auto tmp = bui_to_bul(reciprocal);
		shift_left_ip(tmp, reducerBits);
		sub_ip(tmp, bul1());
		bul rem;
		divmod(tmp, modulus, factor, rem);
		// std::cout << "modulus      = " << bui_to_dec(modulus)      << "\n";
		// std::cout << "reducer      = " << bul_to_dec(reducer)      << "\n";
		// std::cout << "mask         = " << bui_to_dec(mask)         << "\n";
		// std::cout << "reducerBits  = " << reducerBits              << "\n";
		// std::cout << "reciprocal   = " << bui_to_dec(reciprocal)   << "\n";
		// std::cout << "factor       = " << bui_to_dec(factor)       << "\n";
		// std::cout << "convertedOne = " << bui_to_dec(convertedOne) << "\n";
	}

	// convert a standard integer into Montgomery form
	bui convertIn(const bui& x) const {
		return shift_left_mod(x, reducerBits, modulus);
	}

	// convert a Montgomery form integer back to standard form
	bui convertOut(bui x) const {
		mul_mod_ip(x, reciprocal, modulus);
		return x;
	}

	// Multiply two Montgomery-form numbers
	bui multiply(const bui& x, const bui& y) const {
		assert(cmp(x, modulus) < 0 && cmp(y, modulus) < 0);
		bul product = mul(x, y);
		bui t_low = bul_low(product);
		bitwise_and_ip(t_low, mask);
		t_low = mul_low_fast(t_low, factor);
		bitwise_and_ip(t_low, mask);
		auto tmp2 = mul(t_low, modulus);
		u32 c = add_ip_n_imp(product.data(), tmp2.data(), BI_N * 2);
		shift_right_ip(product, reducerBits);
		if (c) {
			bul carry = bul_pow2(BI_BIT * 2 - reducerBits);
			add_ip(product, carry);
		}
		if (cmp(product, modulus) >= 0) {
			sub_ip(product, bui_to_bul(modulus));
		}
		return bul_low(product);
	}

	// Montgomery exponentiation: x^e (e standard, x and result in Montgomery form)
	bui pow(bui x, const bui& e) const {
		bui r = convertedOne;
		u32 hb = highest_bit(e) + 1;
		for (u32 i = 0; i < hb; ++i) {
			if (get_bit(e, i)) {
				r = multiply(r, x);
			}
			x = multiply(x, x);
		}
		return r;
	}
};

// m > 1 and m is odd
inline bool is_valid_modulus(const bui &m) {
	return cmp(m, bui1()) > 0 && get_bit(m, 0);
}

// Montgomery power (faster than naive version for big num), m must be odd
inline bui mr_pow_mod(bui x, const bui& e, const bui& m) {
	if (!is_valid_modulus(m)) {
		return pow_mod(x, e, m);
	}
	MontgomeryReducer mr(m);
	x = mr.convertIn(x);
	bui r = mr.pow(x, e);
	return mr.convertOut(r);
}

#endif
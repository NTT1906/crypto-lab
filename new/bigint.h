#ifndef _BIGINT_H_
#define _BIGINT_H_
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <string>

typedef uint32_t u32;
typedef uint64_t u64;

#ifndef BI_N
#define BI_N (512 / 32)
#endif

#if defined(__GNUC__) || defined(__clang__)
#define GNU_BUILTIN
#endif

#if defined(_MSC_VER)
#define ALWAYS_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#define ALWAYS_INLINE inline __attribute__((always_inline))
#else
#define ALWAYS_INLINE inline
#endif

// big endian: data[0] = MSW
struct bui : std::array<u32, BI_N> {};
struct bul : std::array<u32, BI_N * 2> {};

constexpr bui bui0() { return {}; }

constexpr bui bui1() {
	bui r = {}; r[BI_N - 1] = 1;
	return r;
}

constexpr bui bui_from_u32(const u32 x) {
	bui r = {}; r[BI_N - 1] = x;
	return r;
}

constexpr bul bul0() { return {}; }

constexpr bul bul1() {
	bul r = {}; r[BI_N * 2 - 1] = 1;
	return r;
}

inline bul bul_from_u32(const u32 x) {
	bul r = {}; r[BI_N * 2 - 1] = x;
	return r;
}

inline bool bui_is0(const bui& x) {
	for (const u32 val : x)
		if (val != 0) return false;
	return true;
}

inline bui bul_low(bul& x) {
	bui r{};
	std::copy(x.begin() + BI_N, x.end(), r.begin());
	return r;
}

inline bui bul_high(const bul& x) {
	bui r{};
	std::copy_n(r.begin(), BI_N, r.begin());
	return r;
}

inline bul bui_to_bul(const bui& input) {
	bul r{};
	std::copy(input.begin(), input.end(), r.begin() + BI_N);
	return r;
}

std::string bui_to_dec(const bui& x);
bui bui_from_dec(const std::string& s);

int cmp(const bui &a, const bui &b);
void add_ip(bui& a, const bui& b);
void add_ip(bul& a, const bul& b);
void sub_ip(bui& a, const bui& b);

inline int cmp(const bui &a, const bui &b) {
	for (u32 i = 0; i < BI_N; ++i) {
		if (a[i] != b[i])
			return a[i] > b[i] ? 1 : -1;
	}
	return 0;
}

// a += b;
inline void add_ip(bui& a, const bui& b) {
	u32 c = 0, i = BI_N;
	while (i-- > 0) {
#ifdef GNU_BUILTIN
		c = __builtin_add_overflow(a[i], b[i] + c, &a[i]);
#else
		u64 s = (u64)a[i] + b[i] + c;
		a[i] = (u32)s;
		c = s >> 32;
#endif
	}
}

// a += b
inline void add_ip(bul& a, const bul& b) {
	u32 c = 0, i = BI_N * 2;
	while (i-- > 0) {
#ifdef GNU_BUILTIN
		c = __builtin_add_overflow(a[i], b[i] + c, &a[i]);
#else
		u64 s = (u64)a[i] + b[i] + c;
		a[i] = (u32)s;
		c = s >> 32;
#endif
	}
}

// r = a + b
inline bui add(bui a, const bui& b) {
	add_ip(a, b);
	return a;
}

// a = (a + b) % m
inline void add_mod_ip(bui &a, const bui &b, const bui &m) {
	u32 c = 0;
	for (u32 i = BI_N; i-- > 0;) {
#ifdef GNU_BUILTIN
		c = __builtin_add_overflow(a[i], b[i] + c, &a[i]);
#else
		u64 s = (u64)a[i] + (u64)b[i] + c;
		a[i] = (u32)s;
		c = s >> 32;
#endif
	}
	if (c || cmp(a, m) >= 0) {
		sub_ip(a, m);
	}
}

// a -= b; // assume a > b
#ifdef GNU_BUILTIN
inline void sub_ip(bui& a, const bui& b) {
	u32 borrow = 0, i = BI_N;
	while (i-- > 0)
		borrow = __builtin_sub_overflow(a[i], b[i] + borrow, &a[i]);
}
#else
inline void sub_ip(bui& a, const bui& b) {
	u32 borrow = 0, i = BI_N;
	while (i-- > 0) {
		u64 d = (u64)a[i] - b[i] - borrow;
		a[i] = static_cast<u32>(d);
		borrow = d > a[i] ? 1 : 0; // borrow occurs if there is underflow
	}
}
#endif

inline bui sub(bui a, const bui& b) {
	sub_ip(a, b);
	return a;
}

inline void mul_ref(const bui &a, const bui &b, bul &r) {
	for (u32 i = BI_N; i-- > 0;) {
		u32 c = 0;
		for (u32 j = BI_N; j-- > 0;) {
			u64 p = (u64)a[i] * b[j] + r[i + j + 1] + c;
			r[i + j + 1] = (u32)p;
			c = p >> 32;
		}
		r[i] += c;
	}
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

inline void bui_u32divmod(const bui& a, u32 b, bui& q, u32& r) {
	q = {};
	r = 0;
	for (int i = 0; i < BI_N; ++i) {
		u64 dividend = (u64)r << 32 | a[i];
		q[i] = (u32)(dividend / b);
		r = (u32)(dividend % b);
	}
}

inline bool is_space_c(char c) {
	return c == ' ' || c == '\t';
}

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

std::string bui_to_dec(const bui& x) {
	std::string result;
	if (bui_is0(x)) return "0";
	u32 rems[(BI_N * 32 + 26) / 27];
	size_t count = 0;
	bui n = x;
	bui q;
	while (!bui_is0(n)) {
		u32 r;
		bui_u32divmod(n, 100000000U, q, r);
		rems[count++] = r;
		n = q;
	}
	// first chunk is printed without leading zeros
	result += std::to_string(rems[count - 1]);
	// remaining chunks, padded with leading zeros
	for (size_t i = count - 2; i < count; --i) {
		result += std::string(8 - std::to_string(rems[i]).size(), '0') + std::to_string(rems[i]);
	}
	return result;
}

bui bui_from_dec(const std::string& s) {
	assert(!s.empty() && "bui_from_dec: empty string");
	size_t i = 0;
	// skip leading spaces and optional '+'
	while (is_space_c(s[i])) ++i;
	if (s[i] == '+') ++i;
	assert(s[i] != '-' && "bui_from_dec: negative not supported");
	// skip leading zeros, underscores, spaces
	while (s[i] == '0' || s[i] == '_') ++i;
	bool any_digit = false;
	// process each digit in the decimal string
	constexpr bui n10 = bui_from_u32(10u);
	bui out{};
	bui tmp{};
	for (; i < s.size(); ++i) {
		char c = s[i];
		if (c == '_' || is_space_c(c)) continue;
		if (c < '0' || c > '9') break;  // Stop if non-digit is encountered
		any_digit = true;
		mul_ip(out, n10);
		tmp[BI_N - 1] = c - '0';
		add_ip(out, tmp);
	}
	assert(any_digit && "bui_from_dec: no digits found");
	return out;
}

inline u32 highest_limb(bui &x) {
	for (size_t i = 0; i < BI_N; ++i)
		if (x[i] > 0) return BI_N - i - 1;
	return 0;
}

inline void shift_left_limb(bul &x, u32 l) {
	if (l == 0) return;
	std::copy(x.begin() + l, x.end(), x.begin());
	std::fill(x.end() - l, x.end(), 0);
}




inline bul karatsuba(const bui& a, const bui& b, u32 size) {
	bul result = {};  // The final result will be a bul.

	// Base case: simple multiplication for small numbers (when BI_N is small enough)
	if (size <= 2) {
		return mul(a, b);  // Use your existing multiplication function for small numbers
	}

	int half = BI_N / 2;

	// Split a into a1 and a0 (top and bottom halves)
	bui a1 = {}, a0 = {};
	std::copy_n(a.begin(), half, a1.begin() + BI_N - half);
	std::copy(a.begin() + half, a.end(), a0.begin() + BI_N - half);


	// Split b into b1 and b0 (top and bottom halves)
	bui b1 = {}, b0 = {};
	std::copy_n(b.begin(), half, b1.begin() + BI_N - half);
	std::copy(b.begin() + half, b.end(), b0.begin() + BI_N - half);

	// Compute z2 = a1 * b1 (this is a 'bui' result)
	bul tmp = karatsuba(a1, b1, half);
	bui z2 = bul_low(tmp);

	// Compute z0 = a0 * b0 (this is also a 'bui' result)
	tmp = karatsuba(a0, b0, half);
	bui z0 = bul_low(tmp);

	// Compute z1 = (a1 + a0) * (b1 + b0)
	bui a1_plus_a0 = add(a1, a0);
	bui b1_plus_b0 = add(b1, b0);
	tmp = karatsuba(a1_plus_a0, b1_plus_b0, half);
	bui z1 = bul_low(tmp);

	// Adjust the results with shifts to combine them into the final result
	// z2 needs to be shifted by (2 * half * 32) bits
	bui h = bui_from_u32(1 << (2 * half * 32));
	bul temp1 = mul(z2, h); // Shift z2 left by 2 * half * 32

	// z0 needs to be shifted by (half * 32) bits
	bul temp2 = mul(z0, h);  // Shift z0 left by half * 32

	// Now combine: result = temp1 + (z1 - z2 - z0) + temp2
	sub_ip(z1, z2); // z1 = z1 - z2
	sub_ip(z1, z0); // z1 = z1 - z0
	add_ip(result, temp1);  // result += temp1 (z2 * 10^(2m))
	add_ip(result, bui_to_bul(z1));     // result += z1 (intermediate sum)
	add_ip(result, temp2);  // result += temp2 (z0 * 10^m)
	return result;
}

inline bul karatsuba_test(bui& a, bui& b) {
	u32 n = std::max(highest_limb(a), highest_limb(b));
	return karatsuba(a, b, n);
}

#endif
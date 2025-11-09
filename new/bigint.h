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
	bul r = {}; r[BI_N - 1] = 1;
	return r;
}

inline bul bul_from_u32(const u32 x) {
	bul r = {}; r[BI_N - 1] = x;
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
	bul result = {};
	std::copy(input.begin(), input.end(), result.begin());
	return result;
}

// a += b;
inline void add_ip(bui& a, const bui& b) {
	u32 carry = 0;
	for (int i = BI_N - 1; i >= 0; --i) {  // Iterate from LSW to MSW
		u64 sum = static_cast<u64>(a[i]) + b[i] + carry;
		a[i] = static_cast<u32>(sum);
		carry = static_cast<u32>(sum >> 32);  // Propagate the carry to the next iteration
	}
}

inline void add_ip(bul& a, const bul& b) {
	u32 carry = 0;
	for (int i = BI_N * 2 - 1; i >= 0; --i) {  // Iterate from LSW to MSW
		u64 sum = static_cast<u64>(a[i]) + b[i] + carry;
		a[i] = static_cast<u32>(sum);
		carry = static_cast<u32>(sum >> 32);  // Propagate the carry to the next iteration
	}
}

inline bui add(bui a, const bui& b) {
	add_ip(a, b);
	return a;
}

// a -= b; // assume a > b
inline void sub_ip(bui& a, const bui& b) {
	u32 borrow = 0;
	for (int i = BI_N - 1; i >= 0; --i) {  // Iterate from LSW to MSW
		u64 diff = static_cast<u64>(a[i]) - b[i] - borrow;
		a[i] = static_cast<u32>(diff);
		borrow = (diff > a[i]) ? 1 : 0;  // Borrow occurs if there is underflow
	}
}

inline bui sub(bui a, const bui& b) {
	sub_ip(a, b);
	return a;
}

inline bul mul(const bui& a, const bui& b) {
	bul result = {};
	for (int i = BI_N - 1; i >= 0; --i) {
		u32 carry = 0;
		for (int j = BI_N - 1; j >= 0; --j) {
			u64 product = static_cast<u64>(a[i]) * b[j] + result[i + j + 1] + carry;
			result[i + j + 1] = static_cast<u32>(product);
			carry = product >> 32;
		}
		result[i] += carry;
	}
	return result;
}

inline bui mul_low(const bui& a, const bui& b) {
	bul r = mul(a, b);
	return bul_low(r);
}

inline bool is_space_c(char c) {
	return c == ' ' || c == '\t';
}

inline void bui_u32divmod(const bui& a, u32 b, bui& q, u32& r) {
	q = bui{};
	r = 0;

	for (int i = 0; i < BI_N; ++i) {
		u64 dividend = (static_cast<u64>(r) << 32) | a[i];
		q[i] = static_cast<u32>(dividend / b);
		r = static_cast<u32>(dividend % b);
	}
}

std::string bui_to_dec(const bui& x) {
	std::string result;
	if (bui_is0(x)) return "0";
	u32 rems[(BI_N * 32 + 26) / 27];
	size_t count = 0;
	bui n = x;
	bui q;

	// Divide by BASE until the quotient is zero
	while (!bui_is0(n)) {
		u32 r;
		bui_u32divmod(n, 100000000U, q, r);
		rems[count++] = r;
		n = q;
	}

	// First chunk is printed without leading zeros
	result += std::to_string(rems[count - 1]);

	// Remaining chunks, padded with leading zeros
	for (size_t i = count - 2; i < count; --i) {
		result += std::string(8 - std::to_string(rems[i]).size(), '0') + std::to_string(rems[i]);
	}

	return result;
}

bui bui_from_dec(const std::string& s) {
	assert(!s.empty() && "bui_from_dec: empty string");

	bui out = bui0();
	size_t i = 0;

	// Skip leading spaces and optional '+' sign
	while (i < s.size() && is_space_c(s[i])) ++i;
	if (i < s.size() && s[i] == '+') ++i;
	assert(i < s.size() && "bui_from_dec: invalid string format");

	// Skip leading zeros and underscores
	while (i < s.size() && (s[i] == '0' || s[i] == '_')) ++i;

	bool any_digit = false;

	// Process each digit in the decimal string
	for (; i < s.size(); ++i) {
		char c = s[i];
		if (c == '_' || is_space_c(c)) continue;  // Ignore underscores and spaces
		if (c < '0' || c > '9') break;  // Stop if non-digit is encountered
		any_digit = true;

		// Convert the character to its numeric value
		u32 digit = c - '0';

		// Shift the current number by one place (multiply by 10)
		u32 carry = 0;
		for (int j = BI_N - 1; j >= 0; --j) {
			u64 temp = static_cast<u64>(out[j]) * 10 + carry;
			out[j] = static_cast<u32>(temp);
			carry = static_cast<u32>(temp >> 32);  // Propagate the carry
		}

		// Add the current digit (without carry)
		out[BI_N - 1] += digit;
	}

	assert(any_digit && "bui_from_dec: no digits found");
	return out;
}

inline u32 highest_limb(bui &x) {
	for (size_t i = 0; i < BI_N; ++i)
		if (x[i] > 0) return BI_N - i - 1;
	return 0;
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
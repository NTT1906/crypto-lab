#ifndef _BIGINT_H_
#define _BIGINT_H_
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <string>

typedef uint32_t u32;
typedef uint64_t u64;

#ifdef BI_BIT
#define BI_N (BI_BIT / 32)
#endif
#ifndef BI_N
#define BI_N (512 / 32)
#endif

#if defined(_MSC_VER)
#define ALWAYS_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#define ALWAYS_INLINE inline __attribute__((always_inline))
#else
#define ALWAYS_INLINE inline
#endif

// big endian: data[0] = MSW
// bui a = 1 -> a = [0,0,...,1] -> a[BI_N]=1
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

inline u32 get_bit(u32 num, u32 pos) { return num >> pos & 1; }

inline u32 set_bit(u32 num, u32 pos, u32 val) {
	if (pos >= 32) return num;
	u32 mask = (u32)1 << pos;
	return (num & ~mask) | ((val & 1u) ? mask : 0u);
}

inline u32 get_bit(const bui &a, u32 pos) {
	assert(pos < BI_N * 32);
	u32 k = BI_N - 1 - pos / 32;
	return get_bit(a[k], pos % 32);
}

// set in-place
inline void set_bit_ip(bui &a, u32 pos, u32 val) {
	assert(pos < BI_N * 32);
	u32 k = BI_N - 1 - pos / 32;
	a[k] = set_bit(a[k], pos % 32, val);
}

inline bui set_bit(bui a, u32 pos, u32 val) {
	set_bit_ip(a, pos, val);
	return a;
}

inline u32 highest_bit(u32 x) {
#if defined(__GNUC__) || defined(__clang__)
	if (x == 0) return 0;
	return 32 - __builtin_clz(x);
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
		if (x[i] > 0)
			return highest_bit(x[i]) + (BI_N - i - 1) * 32;
	}
	return 0; // all limbs zero
}

inline u32 highest_limb(const bui &x) {
	for (size_t i = 0; i < BI_N; ++i)
		if (x[i] > 0) return BI_N - i - 1;
	return 0;
}

inline void shift_limb_left(bul &x, u32 l) {
	if (l == 0) return;
	if (l >= BI_N * 2) {
		std::fill(x.begin(), x.end(), 0);
		return;
	}
	std::copy(x.begin() + l, x.end(), x.begin());
	std::fill(x.end() - l, x.end(), 0);
}

inline bool bu_is0(const u32 *x, u32 n) {
	assert(n <= BI_N);
	while (n-- > 0)
		if (x[n] != 0) return false;
	return true;
}

inline bool bui_is0(const bui& x) {
	for (const u32 val : x)
		if (val != 0) return false;
	return true;
}

inline bui bul_low(const bul& x) {
	bui r{};
	std::copy(x.begin() + BI_N, x.end(), r.begin());
	return r;
}

inline bui bul_high(const bul& x) {
	bui r{};
	std::copy_n(x.begin(), BI_N, r.begin());
	return r;
}

inline bul bui_to_bul(const bui& x) {
	bul r{};
	std::copy(x.begin(), x.end(), r.begin() + BI_N);
	return r;
}

std::string bui_to_dec(const bui& x);
bui bui_from_dec(const std::string& s);

int cmp(const bui &a, const bui &b);
void add_ip(bui& a, const bui& b);
void add_ip(bul& a, const bul& b);
void sub_ip(bui& a, const bui& b);

void dbl_ip(bui &x);

inline int cmp(const bui &a, const bui &b) {
	for (u32 i = 0; i < BI_N; ++i) {
		if (a[i] != b[i])
			return a[i] > b[i] ? 1 : -1;
	}
	return 0;
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

inline void add_ip_n(u32* a, const u32* b, u32 n) {
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

ALWAYS_INLINE void mul_imp(const u32* a, const u32* b, u32* r, const u32 n) {
	std::fill_n(r, 2 * n, 0);
	for (u32 i = n; i-- > 0;) {
		u32 c = 0;
		for (u32 j = n; j-- > 0;) {
			u64 p = (u64)a[i] * b[j] + r[i + j + 1] + c;
			r[i + j + 1] = (u32)p;
			c = p >> 32;
		}
		r[i] += c;
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

/* bitwise restoring long division
 * q,r can be NULL independent
 */
void divmod_ref(const bui &a, const bui &b, bui *q, bui *r) {
	assert(!bui_is0(b));
	int cm = cmp(a, b);
	if (cm < 0) {
		if (q) *q = {};
		if (r) *r = a;
		return;
	}
	if (cm == 0) {
		if (q) *q = bui1();
		if (r) *r = {};
		return;
	}
	bui quotient{};
	bui rem{};
	u32 n = highest_bit(a);
	while (n-- > 0) {
		u32 c = get_bit(a, n), i = BI_N;
		while (i-- > 0) {
			u32 v = rem[i];
			rem[i] = v << 1 | c;
			c = v >> 31;
		}
		if (cmp(rem, b) >= 0) {
			sub_ip(rem, b);
			quotient = set_bit(quotient, n, 1u);
		}
	}
	if (q) *q = quotient;
	if (r) *r = rem;
}

ALWAYS_INLINE void divmod_fast(const bui &a, const bui &b, bui *q, bui *r) {
	assert(!bui_is0(b));

	int cm = cmp(a, b);
	if (cm < 0) {
		if (q) *q = {};
		if (r) *r = a;
		return;
	}
	if (cm == 0) {
		if (q) *q = bui1();
		if (r) *r = {};
		return;
	}

	bui rem{};
	bui quotient{};

	const u32 a_bits = highest_bit(a);
	const u32 b_bits = highest_bit(b);
	const u32 shift = a_bits - b_bits;

	// fast path: normalize divisor
	bui divisor = b;
	if (shift > 0) {
		for (int i = 0; i < shift; ++i) dbl_ip(divisor);
	}

	for (u32 i = shift + 1; i-- > 0;) {
		// rem <<= 1; rem[0] += current bit from a
		u32 c = get_bit(a, i + b_bits - 1);
		u32 carry = c;
		for (u32 j = BI_N; j-- > 0;) {
			u32 v = rem[j];
			rem[j] = (v << 1) | carry;
			carry = v >> 31;
		}
		if (cmp(rem, divisor) >= 0) {
			sub_ip(rem, divisor);
			set_bit_ip(quotient, i + b_bits - 1, 1);
		}
		// shift divisor right by 1 (divide by 2)
		for (int j = 0; j < BI_N; ++j) {
			u32 v = divisor[j];
			u32 next = (j + 1 < BI_N) ? divisor[j + 1] : 0;
			divisor[j] = (v >> 1) | (next << 31);
		}
	}

	if (q) *q = quotient;
	if (r) *r = rem;
}

// Binary long division using limb-by-limb (Knuth's Algorithm D)
// Operates on 32-bit limbs, O(n*m) where n = limbs in dividend, m = limbs in divisor
inline void divmod_limb_bruh(const bui& num, const bui& den, bui& quot, bui& rem) {
    assert(!bui_is0(den));
    quot = bui0();
    rem = num;

    if (cmp(num, den) < 0) return;

    // Find leading limb of divisor
    u32 d_lead_limb = highest_limb(den);
    u32 d_lead = den[d_lead_limb];
    u32 n_lead_limb = highest_limb(num);

    // Normalize: make divisor's leading bit >= 16 (helps estimation)
    u32 shift = 32 - highest_bit(d_lead);
    bui d = den;
    bui r = num;
    if (shift > 0) {
        for (u32 i = 0; i < BI_N; ++i) {
            u64 v = (u64)d[i] << shift;
            if (i > 0) v |= d[i-1] >> (32 - shift);
            d[i] = (u32)v;
        }
        for (u32 i = 0; i < BI_N; ++i) {
            u64 v = (u64)r[i] << shift;
            if (i > 0) v |= r[i-1] >> (32 - shift);
            r[i] = (u32)v;
        }
    }
    u32 d0 = d[d_lead_limb];  // normalized leading limb

    // Main loop: process each limb of dividend from MSB
    for (int j = (int)n_lead_limb - (int)d_lead_limb; j >= 0; --j) {
        u32 idx = d_lead_limb + j;

        // Estimate quotient digit: qhat = floor( (r[idx]*2^32 + r[idx-1]) / d0 )
        u64 r_high = (u64)r[idx] << 32 | (idx > 0 ? r[idx-1] : 0);
        u64 qhat = r_high / d0;
        if (qhat >= (1ULL << 32)) qhat = (1ULL << 32) - 1;

        // Multiply and subtract: r -= qhat * d (shifted by j limbs)
        u64 borrow = 0;
        for (u32 k = 0; k < BI_N; ++k) {
            u32 d_val = (k <= d_lead_limb + j) ? d[k - j] : 0;
            u64 prod = (u64)qhat * d_val + borrow;
            u64 diff = (u64)r[k] - (u32)prod;
            r[k] = (u32)diff;
            borrow = (prod >> 32) + (diff >> 32);
        }
        // Add back if borrow occurred (qhat was too large)
        if (borrow) {
            qhat--;
            u32 carry = 0;
            for (u32 k = j; k < BI_N; ++k) {
                u64 sum = (u64)r[k] + (k - j < BI_N ? d[k - j] : 0) + carry;
                r[k] = (u32)sum;
                carry = sum >> 32;
            }
        }

        // Store quotient digit
        if (j + d_lead_limb < BI_N)
            quot[j + d_lead_limb] = (u32)qhat;
    }

    // Denormalize remainder
    if (shift > 0) {
        u32 carry = 0;
        for (u32 i = BI_N; i-- > 0;) {
            u32 v = r[i] >> shift | carry;
            carry = r[i] << (32 - shift);
            rem[i] = v;
        }
    } else {
        rem = r;
    }
}

inline void shift_left_ip(bui &x, u32 amnt) {
	if (amnt == 0) return;
	const u32 limbs = amnt / 32;
	if (limbs >= BI_N) { x.fill(0); return; }
	const u32 bits = amnt % 32;
	// limb-only move (toward MSW)
	if (limbs) {
		std::copy_backward(x.begin() + limbs, x.end(), x.begin() + BI_N - limbs);
		std::fill(x.end() - limbs, x.end(), 0);
	}
	// intra-word stitch (only if bits != 0)
	if (bits) {
		u32 c = 0, i = BI_N;
		while (i-- > 0) {
			u32 tmp = x[i];
			x[i] = tmp << bits | c;
			c = tmp >> (32 - bits);
		}
	}
}

bui shift_left(const bui &x, u32 amnt) {
	if (amnt == 0) return x;
	u32 limbs = amnt / 32;
	if (limbs >= BI_N) return {};
	u32 bits = amnt % 32;
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

/* -------------------------------------------------------------
   Big-endian in-place right shift (bits) on a bui
   ------------------------------------------------------------- */
inline void shift_right_ip(bui& x, u32 amnt) {
	if (amnt == 0) return;
	if (amnt >= 32 * BI_N) {
		x.fill(0);
		return;
	}
	u32 limbs = amnt / 32, bits = amnt % 32;
	if (limbs) {
		std::copy_backward(x.begin(), x.end() - limbs, x.end());
		std::fill_n(x.begin(), limbs, 0);
	}
	if (bits) {
		u32 c = 0;
		for (u32 i = BI_N; i-- > 0;) {
			u32 v = x[i];;
			x[i] = x[i] >> bits | c;
			c = v >> (32 - bits);
		}
	}
}

// ------------------------------------------------------------
//  Limb-by-limb long division (Knuth Algorithm D) – **big-endian**
//  Works directly on bui (data[0] = MSB, data[BI_N-1] = LSB)
// ------------------------------------------------------------
inline void divmod_limb(const bui& a, const bui& b, bui& quot, bui& rem) {
    assert(!bui_is0(b));
	int cm = cmp(a, b);
	if (cm < 0) {
		quot = {};
		rem = a;
		return;
	}
	if (cm == 0) {
		quot = bui1();
		rem = {};
		return;
	}

	quot = bui0();
	rem  = a;

    // ----- find leading limb indices (big-endian) -----
	u32 d_lead = highest_limb(b);
	u32 n_lead = highest_limb(a);
    u32 d0 = b[BI_N - d_lead - 1];

    // ----- normalize divisor (make MSB >= 2^31) -----
    const u32 norm_shift = d0 == 0 ? 0 : 32 - highest_bit(d0);
    bui d = b;
    bui r = a;

	if (norm_shift) {
		shift_left_ip(d, norm_shift);
		shift_left_ip(r, norm_shift);
		d_lead = highest_limb(d);
		d0 = d[BI_N - d_lead - 1];
	}

	u32 limb_diff = (n_lead > d_lead) ? (n_lead - d_lead) : 0;
	for (u32 j = 0; j <= limb_diff; ++j) {
		u32 r_msw_idx = BI_N - d_lead - 1 - j; // current MSW of remainder

		u64 r_high = (u64)r[r_msw_idx] << 32;
		if (r_msw_idx + 1 < BI_N) r_high |= r[r_msw_idx + 1];

		u64 qhat = r_high / d0;
		if (qhat >= (1ULL << 32)) qhat = (1ULL << 32) - 1;

		// subtract qhat * d shifted left by j limbs
		u64 borrow = 0;
		for (u32 k = 0; k < BI_N; ++k) {
			u32 d_idx = k + j;
			u32 dval = (d_idx < BI_N) ? d[d_idx] : 0;
			u64 prod = qhat * dval + borrow;
			borrow = prod >> 32;
			u32 p_low = (u32)prod;
			if (k < BI_N) {
				u64 diff = (u64)r[k] - p_low;
				r[k] = (u32)diff;
				if (diff > (u64)r[k]) ++borrow;
			}
		}

		if (borrow) {
			--qhat;
			u32 carry = 0;
			for (u32 k = j; k < BI_N; ++k) {
				u32 d_idx = k - j;
				u32 dval = (d_idx < BI_N) ? d[d_idx] : 0;
				u64 sum = (u64)r[k] + dval + carry;
				r[k] = (u32)sum;
				carry = sum >> 32;
			}
		}

		u32 q_idx = BI_N - d_lead - 1 - j;
		if (q_idx < BI_N) quot[q_idx] = (u32)qhat;
	}

    // ----- denormalize remainder (shift right) -----
	if (norm_shift) shift_right_ip(r, norm_shift);
    rem = r;
}

// ------------------------------------------------------------
//  Convenience wrappers (big-endian safe)
// ------------------------------------------------------------
inline bui div_limb(const bui& a, const bui& b) {
    bui q, r;
    divmod_limb(a, b, q, r);
    return q;
}

inline bui mod_limb(const bui& a, const bui& b) {
    bui q, r;
    divmod_limb(a, b, q, r);
    return r;
}

inline void u32divmod(const bui& a, u32 b, bui& q, u32& r) {
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
		u32divmod(n, 100000000U, q, r);
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
		if (c < '0' || c > '9') break;
		any_digit = true;
		mul_ip(out, n10);
		tmp[BI_N - 1] = c - '0';
		add_ip(out, tmp);
	}
	assert(any_digit && "bui_from_dec: no digits found");
	return out;
}

bui bui_from_hex(const std::string& s) {
	assert(!s.empty() && "bui_from_hex: empty string");
	size_t i = 0;
	// skip leading spaces
	while (i < s.size() && is_space_c(s[i])) ++i;
	// optional "0x" or "0X" prefix
	if (i + 1 < s.size() && s[i] == '0' && (s[i+1] == 'x' || s[i+1] == 'X')) i += 2;

	bool any_digit = false;
	bui out{};
	bui tmp{};
	constexpr bui n16 = bui_from_u32(16u);

	for (; i < s.size(); ++i) {
		char c = s[i];
		if (c == '_' || is_space_c(c)) continue;
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

ALWAYS_INLINE void split_bui(const bui &x, bui &high, bui &low, u32 n) {
	std::copy_n(x.begin(), n, high.begin() + (BI_N - n));
	std::copy(x.begin() + n, x.end(), low.begin() + n);
}

inline bul karatsuba(const bui &a, const bui &b, const u32 n) {
	bul r{};
	if (n <= 16) return mul(a, b);
	u32 half = n / 2;

	bui a1{}, a0{}, b1{}, b0{};
	split_bui(a, a1, a0, half);
	split_bui(b, b1, b0, half);

	bul z2 = karatsuba(a1, b1, half);
	bul z0 = karatsuba(a0, b0, half);

	bui a_sum = add(a1, a0);
	bui b_sum = add(b1, b0);
	bul z1 = karatsuba(a_sum, b_sum, half);
	sub_ip(z1, z2);
	sub_ip(z1, z0);

	shift_limb_left(z2, 2 * half);
	shift_limb_left(z1, half);

	add_ip(r, z0);
	add_ip(r, z1);
	add_ip(r, z2);
	return r;
}

constexpr size_t KARATSUBA_CUTOFF = 2;
// constexpr size_t KARATSUBA_CUTOFF = 4;
// constexpr size_t KARATSUBA_CUTOFF = 8;
// constexpr size_t KARATSUBA_CUTOFF = 16;
// constexpr size_t KARATSUBA_CUTOFF = 32; // tune this experimentally
inline void karatsuba_be_rec_old(const u32* a, const u32* b, u32* r, const u32 n, u32* scratch) {
    if (n <= KARATSUBA_CUTOFF) {
        mul_imp(a, b, r, n);
        return;
    }

    u32 half = n / 2;
    const u32* a1 = a;
    const u32* a0 = a + half;
    const u32* b1 = b;
    const u32* b0 = b + half;

    u32* z2 = r;       // high part
    u32* z0 = r + n;   // low part
    u32* z1 = scratch; // middle temp (2*half)

	const u32 maxlen = std::max(half, half);

    u32* tmp_a = z1 + 2 * maxlen;
    u32* tmp_b = tmp_a + maxlen;
    u32* tmp_scratch = tmp_b + maxlen;

    karatsuba_be_rec_old(a0, b0, z0, half, tmp_scratch); // z0 = a0 * b0
    karatsuba_be_rec_old(a1, b1, z2, half, tmp_scratch); // z2 = a1 * b1

	add_n(a1, a0, tmp_a, half); // tmp_a = a1 + a0
	add_n(b1, b0, tmp_b, half); // tmp_b = b1 + b0
    karatsuba_be_rec_old(tmp_a, tmp_b, z1, half, tmp_scratch); // z1 = (a1 + a0) * (b1 + b0)

    // z1 = z1 - z2 - z0
	sub_n(z1, z2, z1, 2 * half);
	sub_n(z1, z0, z1, 2 * half);
    // combine: r = z2 << (2*half*32) + z1 << (half*32) + z0
	add_n(r + half, z1, r + half, 2 * half);
}

inline void karatsuba_be_rec(const u32* a, const u32* b, u32* r, u32 n, u32* scratch) {
	if (bu_is0(a, n) || bu_is0(b, n)) {
		std::fill_n(r, 2 * n, 0);
		return;
	}
	if (n <= KARATSUBA_CUTOFF) {
		mul_imp(a, b, r, n);
		return;
	}

	const u32 half  = n / 2;
	const u32 other = n - half;       // may be half+1 if n is odd

	// Big-endian split
	const u32* a1 = a;           // high half
	const u32* a0 = a + half;    // low  half
	const u32* b1 = b;
	const u32* b0 = b + half;

	// workspace layout
	u32* z0 = scratch;                 // size 2*other
	u32* z1 = z0 + 2 * other;          // size 2*other
	u32* z2 = z1 + 2 * other;          // size 2*half
	u32* tmpa = z2 + 2 * half;
	u32* tmpb = tmpa + other;
	u32* subscratch = tmpb + other;

	// z0 = a0 * b0
	karatsuba_be_rec(a0, b0, z0, other, subscratch);

	// z2 = a1 * b1
	karatsuba_be_rec(a1, b1, z2, half, subscratch);

	// tmpa = a0 + a1 (aligned to low indices)
	std::fill_n(tmpa, other, 0);
	std::fill_n(tmpb, other, 0);
	std::copy(a0 + (other - half), a0 + other, tmpa + (other - half));
	for (u32 i = 0; i < half; ++i)
		tmpa[i] += a1[i];
	std::copy(b0 + (other - half), b0 + other, tmpb + (other - half));
	for (u32 i = 0; i < half; ++i)
		tmpb[i] += b1[i];

	// z1 = (a0+a1)*(b0+b1)
	karatsuba_be_rec(tmpa, tmpb, z1, other, subscratch);

	// z1 = z1 - z2 - z0
	sub_n(z1 + (2 * other - 2 * half), z2, z1 + (2 * other - 2 * half), 2 * half);
	sub_n(z1, z0, z1, 2 * other);

	// clear result
	std::fill_n(r, 2 * n, 0);

	// combine (big-endian)
	// copy z0 → low end
	std::copy(z0 + 2 * other - n, z0 + 2 * other, r + 2 * n - 2 * other);

	// add z1 shifted by (other limbs)
	add_n(r + (n - other), z1 + 2 * other - n, r + (n - other), 2 * other);

	// add z2 shifted by (2*other limbs)
	add_n(r, z2 + 2 * half - n, r, 2 * half);
}

inline u32 next_pow2(u32 x) {
	if (x == 0) return 1;
	x--;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	x++;
	return x;
}


inline bul karatsuba_be_top(const bui& a, const bui& b) {
    constexpr size_t n = BI_N;
	// u32 n = std::max(highest_limb(a), highest_limb(b));
	// n = next_pow2(n) * 2;
    bul r{};
    std::array<u32, 6 * BI_N> scratch{};
    // std::array<u32, 8 * BI_N> scratch{};
    karatsuba_be_rec_old(a.data(), b.data(), r.data(), n, scratch.data());
    return r;
}

// for compatibility with your test code
inline bul karatsu_test(const bui& a, const bui& b) {
    return karatsuba_be_top(a, b);
}

#endif
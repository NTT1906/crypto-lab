#define BI_BIT 512
#include <chrono>
#include <iostream>

#include "bigint.h"

constexpr u32 POLY_R = 14;
struct Poly : std::array<bui, POLY_R>{};

inline bui mod_reduce_fast(bui x, const bui& m) {
	int shift = highest_bit(x) - highest_bit(m);
	if (shift < 0) return x;

	for (; shift >= 0; --shift) {
		bui tmp = m;
		shift_left_ip(tmp, shift);
		if (cmp(x, tmp) >= 0)
			sub_ip(x, tmp);
	}
	return x;
}

inline bul mod_reduce_fast(bul x, const bui& m) {
	int shift = highest_bit(x) - highest_bit(m);
	if (shift < 0) return x;

	for (; shift >= 0; --shift) {
		bul tmp = bui_to_bul(m);
		shift_left_ip(tmp, shift);
		if (cmp(x, tmp) >= 0)
			sub_ip(x, tmp);
	}
	return x;
}

inline void divmod_reduce_fast(const bui& a, const bui& b, bui &q, bui &r) {
	q = {};
	r = a;
	int shift = highest_bit(a) - highest_bit(b);
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

// Montgomery helper: compute m_inv = -m^{-1} mod 2^32
u32 montgomery_inverse(u32 m_lsw) {
	u32 m_inv = 1;
	for (int i = 0; i < 5; ++i)
		m_inv *= 2 - m_lsw * m_inv;
	return -m_inv;
}

// Montgomery multiplication
// res = a * b * R^-1 mod m
bui mont_mul(const bui &a, const bui &b, const bui &m, u32 m_inv) {
	bul t;
	mul_ref(a, b, t);
	// Montgomery reduction
	for (int i = BI_N * 2 - 1; i >= BI_N; --i) {
		u32 u = (u32)((u64)t[i] * m_inv);
		u64 carry = 0;
		for (int j = BI_N - 1; j >= 0; --j) {
			size_t k = i - (BI_N - 1 - j);
			u64 sum = (u64)t[k] + (u64)u * m[j] + carry;
			t[k] = (u32)sum;
			carry = sum >> 32;
		}
		for (int k = i - BI_N; carry && k >= 0; --k) {
			u64 sum = (u64)t[k] + carry;
			t[k] = (u32)sum;
			carry = sum >> 32;
		}
	}

	// result = least significant BI_N words (end of array)
	bui res{};
	for (int i = 0; i < BI_N; ++i)
		res[BI_N - 1 - i] = t[BI_N * 2 - 1 - i];

	if (cmp(res, m) >= 0)
		sub_ip(res, m);

	return res;
}

// Convert to Montgomery form: x_bar = x * R mod m
// bui to_mont(const bui &x, const bui &m, u32 m_inv) {
// 	bui R{};
// 	R[BI_N - 1] = 1; // R = 2^(32*BI_N)
// 	return mont_mul(x, R, m, m_inv);
// }

bui to_mont(const bui &x, const bui &m, u32 m_inv) {
	// Compute R^2 mod m
	bui R2{};
	R2[BI_N - 1] = 1;
	for (int i = 0; i < 2 * BI_N * 32; ++i) {
		add_ip(R2, R2);
		if (cmp(R2, m) >= 0)
			sub_ip(R2, m);
	}
	return mont_mul(x, R2, m, m_inv);
}

// Modular exponentiation using Montgomery
bui mont_pow(bui base, const bui &exp, const bui &m) {
	u32 m_inv = montgomery_inverse(m[BI_N - 1]);

	bui one{};
	one[BI_N - 1] = 1;

	// compute R^2 mod m
	bui R2{};
	R2[BI_N - 1] = 1;
	for (int i = 0; i < 2 * BI_N * 32; ++i) {
		add_ip(R2, R2);
		if (cmp(R2, m) >= 0)
			sub_ip(R2, m);
	}

	bui base_m = mont_mul(base, R2, m, m_inv);
	bui res_m  = mont_mul(one, R2, m, m_inv);

	for (int i = 0; i < BI_N; ++i) {
		for (int bit = 31; bit >= 0; --bit) {
			res_m = mont_mul(res_m, res_m, m, m_inv);
			if ((exp[i] >> bit) & 1)
				res_m = mont_mul(res_m, base_m, m, m_inv);
		}
	}

	return mont_mul(res_m, one, m, m_inv); // back from Montgomery form
}

Poly poly_pow_1x(const bui &n) {
	Poly base{};
	base[0] = bui1(); base[1] = bui1();
	Poly res{};
	res[0] = bui1();
	return res;
}

void mul_mod_ip(bui &a, bui b, bui &m) {
	a = mod_reduce_fast(a, m);
	b = mod_reduce_fast(b, m);
	bul r;
	mul_ref(a, b, r);
	r = mod_reduce_fast(r, m);
	a = bul_low(r);
}

bui pow_mod(bui x, const bui& e, bui &m) {
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

struct MontgomeryReducer {
	bui modulus;      // must be odd >= 3
	bui reducer;      // power of 2
	int reducerBits;  // log2(reducer)
	bui reciprocal;   // reducer^-1 mod modulus
	bui mask;         // reducer - 1
	bui factor;       // (reducer * reciprocal - 1) / modulus
	bui convertedOne; // convertIn(1)
	static bui modInverse(const bui& a, const bui& m);

	MontgomeryReducer(const bui& mod) : modulus(mod) {
		assert(get_bit(modulus, 0) && cmp(modulus, bui1()) == 1);
		// compute reducer as a power of 2 bigger than modulus
		reducerBits = (highest_bit(modulus) / 8 + 1) * 8;  // multiple of 8
		reducer = shift_left(bui1(), reducerBits);
		mask = sub(reducer, bui1());                         // mask = reducer - 1
		// assert(gcd(reducer, modulus) == bui1()); m must be a prime thingy
		// other precomputations
		reciprocal = modInverse(reducer, modulus);         // reducer^-1 mod modulus
		factor = (reducer * reciprocal - oneBui()) / modulus;
		auto tmp = mul(reducer, reciprocal);
		sub_ip(tmp, bul1());
		bui tmp2;
		divmod_reduce_fast(tmp, modulus, factor, tmp2);
		convertedOne = reducer % modulus;
	}
};

bui MontgomeryReducer::modInverse(const bui& a, const bui& m) {
	bui t{}, newt = bui1();
	bui r = m, newr = a;

	while (!bui_is0(newr)) {
		bui quotient = r / newr;

		bui temp = t;
		t = newt;
		newt = temp - quotient * newt;

		temp = r;
		r = newr;
		newr = temp - quotient * newr;
	}

	if (r != oneBui()) throw std::runtime_error("Inverse does not exist");

	if (t < bui{}) t = t + m;  // make positive
	return t;
}

int main() {
	bui u = bui_from_dec("5");
	bui v = bui_from_dec("13");
	bui m = bui_from_dec("17");
	u32 m_inv = montgomery_inverse(m[BI_N - 1]);
	auto c = mont_mul(u, v, m, m_inv);
	printf("u = %s\n", bui_to_dec(u).c_str());
	printf("v = %s\n", bui_to_dec(v).c_str());
	printf("m = %s\n", bui_to_dec(m).c_str());
	printf("c = %s\n", bui_to_dec(c).c_str());

	// auto c1 = u;
	// mul_mod_ip(c1, v, m);
	auto c1 = pow_mod(u, v, m);
	printf("c1= %s\n", bui_to_dec(c1).c_str());

	return 0;
	bui a = bui_from_dec("115792089237316195423570985008687907853269984665640564039457584007913129639");
	bui b = bui_from_dec("123456789123456789");
	// bui b = bui_from_dec("2");
	printf("A = %s\n", bui_to_dec(a).c_str());
	printf("B = %s\n", bui_to_dec(b).c_str());
	bui r;
	bui q;
	divmod_knuth(a, b, q, r);
	auto t0 = mul(b, q);
	auto r3 = bui_to_bul(a);
	sub_ip(r3, t0);
	printf("Q = %s\n", bui_to_dec(q).c_str());
	printf("R = %s\n", bui_to_dec(r).c_str());
	printf("R3= %s\n", bui_to_dec(bul_low(r3)).c_str());
	bui q2, r2;
	divmod_reduce_fast(a, b, q2, r2);
	printf("Q2= %s\n", bui_to_dec(q2).c_str());
	printf("R2= %s\n", bui_to_dec(r2).c_str());
	r2 = mod_reduce_fast(a, b);
	printf("R2= %s\n", bui_to_dec(r2).c_str());
}

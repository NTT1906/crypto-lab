#define BI_BIT 512
#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "bigint.h"

constexpr u32 POLY_R = 14;
struct Poly : std::array<bui, POLY_R>{};

inline bui mod(bui x, const bui& m) {
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

inline bul mod(bul x, const bui& m) {
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

inline void divmod(const bui& a, const bui& b, bui &q, bui &r) {
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

inline void divmod(const bul& a, const bui& b, bui &q, bul &r) {
	q = {};
	r = a;
	int shift = highest_bit(a) - highest_bit(b);
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

Poly poly_pow_1x(const bui &n) {
	Poly base{};
	base[0] = bui1(); base[1] = bui1();
	Poly res{};
	res[0] = bui1();
	return res;
}

void mul_mod_ip(bui &a, bui b, const bui &m) {
	a = mod(a, m);
	b = mod(b, m);
	bul r;
	mul_ref(a, b, r);
	r = mod(r, m);
	a = bul_low(r);
}

bui pow_mod(bui x, const bui& e, const bui &m) {
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

bui bitwise_and(bui a, bui b) {
	bui r;
	for (u32 i = BI_N; i-- > 0;) {
		r[i] = a[i] & b[i];
	}
	return r;
}

bool modinv(bui a, const bui &m, bui &inv_out) {
	// invalid modulus or zero
	if (bui_is0(m)) return false;
	if (cmp(a, m) >= 0) a = mod(a, m);
	if (bui_is0(a)) return false; // zero has no inverse

	// Initialize: r0 = m, r1 = a; t0 = 0, t1 = 1
	bui r0 = m, r1 = a;
	bui t0{};
	bui t1 = bui1();

	while (!bui_is0(r1)) {
		// q = r0 / r1, rem = r0 % r1
		bui q, rem;
		divmod(r0, r1, q, rem);
		// r0, r1 = r1, rem
		r0 = r1;
		r1 = rem;

		// t_new = (t0 - q * t1) mod m
		// compute q * t1 -> bul, then reduce modulo m to get r_qt (bui)
		bul prod{};
		mul_ref(q, t1, prod);  // prod = q * t1 (2N words)
		auto qtm_rem = bul_low(mod(prod, m));// qtm_rem = (prod) % m

		// t_new = t0 - qtm_rem  (in modulo m arithmetic)
		bui tnew = t0;
		if (cmp(tnew, qtm_rem) >= 0) {
			sub_ip(tnew, qtm_rem);
		} else {
			// tnew = (t0 - qtm_rem) mod m = m - (qtm_rem - t0)
			bui tmp = qtm_rem;
			sub_ip(tmp, t0);   // tmp = qtm_rem - t0
			tnew = m;
			sub_ip(tnew, tmp); // tnew = m - tmp
		}

		// advance t's
		t0 = t1;
		t1 = tnew;
	}

	// r0 = gcd(a, m). If gcd != 1 -> no inverse.
	if (cmp(r0, bui1()) != 0) return false;

	// t0 is the inverse, ensure it's reduced < m
	if (cmp(t0, m) >= 0) {
		t0 = mod(t0, m);
	}
	inv_out = t0;
	return true;
}

bui shift_left_mod(bui x, int shift, const bui& m) {
	bui p2 = bui_pow2(shift);
	mul_mod_ip(x, p2, m);
	return x;
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

	MontgomeryReducer(const bui& modulus) : modulus(modulus) {
		assert(get_bit(modulus, 0) && cmp(modulus, bui1()) == 1);
		// compute reducer as a power of 2 bigger than modulus
		reducerBits = (highest_bit(modulus) / 8 + 1) * 8;  // multiple of 8
		reducer = shift_left(bui1(), reducerBits);
		mask = sub(reducer, bui1());                         // mask = reducer - 1
		// assert(gcd(reducer, modulus) == bui1()); m must be a prime thingy
		// other precomputations
		modinv(reducer, modulus, reciprocal);         // reducer^-1 mod modulus
		{
			auto tmp = mul(reducer, reciprocal);
			sub_ip(tmp, bul1());
			bul tmp2;
			divmod(tmp, modulus, factor, tmp2);
		}
		convertedOne = mod(reducer, modulus);
	}

	// Convert a standard integer into Montgomery form
	bui convertIn(bui x) const {
		// TODO: shift overflow problem
		return shift_left_mod(x, reducerBits, modulus);
		// shift_left_ip(x, reducerBits);
		// return mod(x, modulus);
	}

	// Convert a Montgomery form integer back to standard form
	bui convertOut(bui x) const {
		mul_mod_ip(x, reciprocal, modulus);
		return x;
	}

	// Multiply two Montgomery-form numbers
	bui multiply(const bui& x, const bui& y) const {
		assert(cmp(x, modulus) < 0 && cmp(y, modulus) < 0);
		// printf("x      = %s\n", bui_to_dec(x).c_str());
		// printf("y      = %s\n", bui_to_dec(y).c_str());
		bul product = mul(x, y);
		bui t_low = bul_low(product);
		// printf("1: p   = %s\n", bul_to_dec(product).c_str());
		t_low = bitwise_and(t_low, mask);
		// printf("2.1: t1= %s\n", bui_to_dec(t_low).c_str());
		t_low = mul_low(t_low, factor);
		// printf("2.2: t2= %s\n", bui_to_dec(t_low).c_str());
		t_low = bitwise_and(t_low, mask);
		// printf("2: temp= %s\n", bui_to_dec(t_low).c_str());
		auto tmp2 = mul(t_low, modulus);
		// printf("3.1: r1= %s\n", bul_to_dec(tmp2).c_str());
		add_ip(product, tmp2);
		// printf("3.2: r2= %s\n", bul_to_dec(product).c_str());
		shift_right_ip(product, reducerBits);
		// printf("3: redu= %s\n", bul_to_dec(product).c_str());
		if (cmp(product, modulus) >= 0) {
			sub_ip(product, bui_to_bul(modulus));
		}
		// printf("4: resu= %s\n", bul_to_dec(product).c_str());
		// if (cmp(product, modulus) >= 0) {
			// sub_ip(product, bui_to_bul(modulus));
			// printf("NO1\n");
		// }
		// if (cmp(product, modulus) >= 0) {
			// sub_ip(product, bui_to_bul(modulus));
			// printf("NO2\n");
		// }
		return bul_low(product);
	}

	// Montgomery exponentiation: x^e (e standard, x and result in Montgomery form)
	bui pow(bui x, const bui& e) const {
		bui r = convertedOne;
		u32 hb = highest_bit(e) + 1;
		// printf("Highest bit of e: {%d}\n", hb);
		for (u32 i = 0; i < hb; ++i) {
			if (get_bit(e, i)) {
				r = multiply(r, x);
			}
			x = multiply(x, x);
		}
		return r;
	}
};

std::string bui_to_hex(const bui &a) {
	std::ostringstream o;
	o << std::hex << std::setfill('0');
	for (u32 i = 0; i < BI_N; ++i) {
		o << std::setw(8) << a[i] << ' ';
	}
	return o.str();
}

int main() {
	bui u = bui_from_dec("260428835329122752520818469321216072583938198616075453742527759001901820374664228839496959095353854544158481165265966459");
	// bui u = bui_from_dec("260428835329122752520818469321216072583938198616075453742527");
	bui v = bui_from_dec("1312312317639123213");
	// bui m = bui_from_dec("108579795932485217312615519053");
	bui m = bui_from_dec("823887783191267813656599693818502133610549771176609410328824491902309472167445766968176579098197424208002485918297593219");

	// bui u = bui_from_dec("123456789");
	// bui v = bui_from_dec("6713");
	// bui m = bui_from_dec("896947");

	// auto sm = shift_mod(u, 400, m);
	// printf("sm = %s\n", bui_to_dec(sm).c_str());
	// return 0;

	printf("u = %s\n", bui_to_dec(u).c_str());
	printf("u = %s\n", bui_to_hex(u).c_str());
	printf("v = %s\n", bui_to_dec(v).c_str());
	printf("v = %s\n", bui_to_hex(v).c_str());
	printf("m = %s\n", bui_to_dec(m).c_str());


	auto start_mont = std::chrono::high_resolution_clock::now();
	MontgomeryReducer mr(m);
	// printf("reducerBits= %d\n", mr.reducerBits);
	// printf("reciprocal = %s\n", bui_to_dec(mr.reciprocal).c_str());
	// printf("mask       = %s\n", bui_to_dec(mr.mask).c_str());
	// printf("reducer    = %s\n", bui_to_dec(mr.reducer).c_str());
	// printf("factor     = %s\n", bui_to_dec(mr.factor).c_str());
	// printf("c1         = %s\n", bui_to_dec(mr.convertedOne).c_str());

	auto cu = mr.convertIn(u);
	// auto ou = mr.convertOut(cu);
	// printf("cu= %s\n", bui_to_dec(cu).c_str());
	// printf("ou= %s\n", bui_to_dec(ou).c_str());
	// printf("o1= %s\n", bui_to_dec(mr.convertedOne).c_str());
	// ---------- Time Montgomery exponentiation ----------
	auto c = mr.pow(cu, v);
	auto cc = mr.convertOut(c);
	auto end_mont = std::chrono::high_resolution_clock::now();
	auto duration_mont = std::chrono::duration_cast<std::chrono::nanoseconds>(end_mont - start_mont).count();

	printf("Montgomery pow: c    = %s\n", bui_to_dec(c).c_str());
	printf("Converted out: cc    = %s\n", bui_to_dec(cc).c_str());

	// ---------- Time standard modular exponentiation ----------
	auto start_std = std::chrono::high_resolution_clock::now();
	auto c1 = pow_mod(u, v, m);
	auto end_std = std::chrono::high_resolution_clock::now();
	auto duration_std = std::chrono::duration_cast<std::chrono::nanoseconds>(end_std - start_std).count();

	printf("Standard pow_mod: c1 = %s\n", bui_to_dec(c1).c_str());
	printf("--------------------------------\n");
	printf("Montgomery pow time  : %lld ns\n", duration_mont);
	printf("Standard pow_mod time: %lld ns\n", duration_std);

	return 0;
}

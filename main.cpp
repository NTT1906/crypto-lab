#define CIGINT_IMPLEMENTATION
#define CIGINT_STRIP_PREFIX
#define CIGINT_N (512 / 32)
#include "cigint.h"
#include <iostream>

/* Optional short aliases if you use CIGINT_STRIP_PREFIX like the rest */
#if defined(CIGINT_STRIP_PREFIX)
  #define from_dec  cigint_from_dec
  #define from_hex  cigint_from_hex
  #define from_bin  cigint_from_bin
  #define from_str  cigint_from_str
#endif

/* floor(sqrt(a)) via binary search: O(log n * M(n)) */
Cigint cigint_isqrt(CFREF(Cigint) a) {
	/* 0,1 fast-path */
	if (cigint_is0(a)) return a;
	Cigint one = (Cigint){0}; one.data[CIGINT_N - 1] = 1;
	if (cigint_cmp(a, one) <= 0) return a;

	u32 bits = cigint_highest_order(a);              /* 1.. */
	u32 up   = (bits + 1u) / 2u;                     /* ceil(bits/2) */
	Cigint low = {0};
	Cigint high = cigint_pow2(up);                    /* 2^up */
	high = cigint_sub(high, one);                     /* 2^up - 1 */

	while (cigint_cmp(low, high) < 0) {
		/* mid = (low + high + 1) >> 1 */
		Cigint mid = cigint_add(low, high);
		mid = cigint_add(mid, one);
		mid = cigint_shr(mid, 1);

		Cigint mid2 = cigint_mul(mid, mid);
		if (cigint_cmp(mid2, a) <= 0) {
			low = mid;
		} else {
			high = cigint_sub(mid, one);
		}
	}
	return low;
}

Cigint cigint_isqrt3(Cigint& a) {
	if (cigint_is0(a)) return a;

	Cigint rem = a;                      // running remainder
	Cigint res = (Cigint){0};            // "trial root × 2"

	u32 bits = cigint_highest_order(a);
	u32 k = bits ? (bits - 1u) : 0u;    // top bit index
	if (k & 1u) k--;                     // make even → largest power of 4 ≤ a
	Cigint bit = cigint_pow2(k);         // starts at 4^⌊(bits-1)/2⌋

	while (!cigint_is0(bit)) {
		Cigint trial = cigint_add(res, bit);
		if (cigint_cmp(rem, trial) >= 0) {
			rem = cigint_sub(rem, trial);
			res = cigint_add(cigint_shr(res, 1), bit);
		} else {
			res = cigint_shr(res, 1);
		}
		bit = cigint_shr(bit, 2);
	}
	return res;
}

Cigint cigint_isqrt2(const Cigint& a) {
	return 1 << cigint_highest_order(a) / 2;
}

Cigint cigint_isqrt4(Cigint a) {
	if (cigint_is0(a)) return a;

	Cigint res = CIGINT_ZERO();            // “trial root × 2"
	u32 bits = cigint_highest_order(a);
	u32 k = bits ? (bits - 1u) : 0u;    // top bit index
	if (k & 1u) k--;                     // make even → largest power of 4 ≤ a
	Cigint bit = cigint_pow2(k);         // starts at 4^⌊(bits-1)/2⌋

	while (!cigint_is0(bit)) {
		Cigint trial = cigint_add(res, bit);
		if (cigint_cmp(a, trial) >= 0) {
			a = cigint_sub(a, trial);
			res = cigint_add(cigint_shr(res, 1), bit);
		} else {
			res = cigint_shr(res, 1);
		}
		bit = cigint_shr(bit, 2);
	}
	return res;
}

Cigint euclidean_gcd(Cigint a, Cigint b) {
	while (!cigint_is0(b)) {
		Cigint tmp = b;
		b = cigint_mod(a, b);
		a = tmp;
	}
	return a;
}

bool isPrime(const Cigint &n) {
	Cigint p = {2};
	if (cigint_cmp(n, p) < 0) return false;
	while (cigint_cmp(p * p, n) <= 0) {
		if (cigint_mod(n, p) == 0) return false;
		p = cigint_add(p, 1);
	}
	return true;
}

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wnon-pod-varargs"
int main() {
	Cigint tmp = {0};
	tmp = cigint_from_dec("000001000A11");
	cprintf("TMP: %Cd\n", tmp);
	tmp = cigint_from_str("12345678901234567890");
	cprintf("TMP: %Cd\n", tmp);
	tmp = cigint_from_str("0xdead_beef_012345");
	cprintf("TMP: %Cd\n", tmp);
	cprintf("HEX: %Cx\n", tmp);
	cprintf("HEX: %CX\n", tmp);
	tmp = cigint_from_str("0b1010_1111_0011");
	cprintf("Dec: %Cd\n", tmp);
	cprintf("Bin: %Cb\n", tmp);
	auto r = cigint_bit_reverse_high(tmp);
	cprintf("Dec:_%Cd\n", r);
	cprintf("Bin:_%Cb\n", r);
	cprintf("Dec: %Cd\n", tmp);
	cprintf("Bin: %Cb\n", tmp);
	// cigint_bit_reverse_inplace(&tmp);
	cprintf("Dec: %Cd\n", tmp);
	cprintf("Bin: %Cb\n", tmp);
	return 0;
	while (cigint_cmp(tmp, 0) > 0) {
		tmp = cigint_shr(tmp, 1);
		cprintf("Dec: %Cd\n", tmp);
		cprintf("Bin: %Cb\n", tmp);
	}
	tmp = 1;
	tmp = cigint_shl(tmp, 1);
	cprintf("Dec: %Cd\n", tmp);
	cprintf("Bin: %Cb\n", tmp);
	tmp -= 1;
	cprintf("Dec: %Cd\n", tmp);
	cprintf("Bin: %Cb\n", tmp);
	tmp -= 1;
	cprintf("Dec: %Cd\n", tmp);
	cprintf("Bin: %Cb\n", tmp);
	tmp -= 1;
	cprintf("Dec: %Cd\n", tmp);
	cprintf("Bin: %Cb\n", tmp);
	// tmp = cigint_from_str("  42_000  ");
	// cprintf("TMP: %Cd\n", tmp);
	// tmp.data[CIGINT_N - 1] = (1ul << 31);
	// tmp = 2;
	// Cigint res1 = cigint_isqrt(tmp);
	// Cigint res2 = cigint_isqrt2(tmp);
	// Cigint res3 = cigint_isqrt3(tmp);
	// Cigint res4 = cigint_isqrt4(tmp);
	// cprintf("ISQRT1 of %Cd: %Cd\n", tmp, res1);
	// cprintf("ISQRT2 of %Cd: %Cd\n", tmp, res2);
	// cprintf("ISQRT2 of %Cd: %Cd\n", tmp, res3);
	// cprintf("ISQRT2 of %Cd: %Cd\n", tmp, res4);
	// cprintf("%Cb\n", tmp);
	// cprintf("%Cb\n", res1);
	// cprintf("%Cb\n", res2);
	// cprintf("%Cb\n", res3);
	// cprintf("%Cb\n", res4);
	return 0;
}
#pragma clang diagnostic pop
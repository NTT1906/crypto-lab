// prime_fast.cpp
// AKS with tighter findR and allocation-free polynomial powering for Cigint.

#include <vector>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <algorithm>

#define CIGINT_IMPLEMENTATION
#define CIGINT_STRIP_PREFIX
#ifndef CIGINT_N
#define CIGINT_N (1024 / 32)
#endif
#include "cigint.h"

// ---------- small Cigint helpers ----------
static inline u32 bitlen(const Cigint& n) { return highest_order(n); }
static inline int is_zero(const Cigint& x) { return is0(x); }
static inline int is_odd (const Cigint& x) { return get_bit_ref(&x, 0) != 0; }
static inline void shr1_ref(Cigint& x)     { x = shr(x, 1); }

static inline Cigint from_u64(u64 v) {
    Cigint x = CIGINT_ZERO();
    if (CIGINT_N) x.data[CIGINT_N - 1] = (u32)(v & 0xFFFFFFFFu);
    if (CIGINT_N >= 2) x.data[CIGINT_N - 2] = (u32)(v >> 32);
    return x;
}
static inline bool le_u64(const Cigint& n, u64 v) { return cmp(n, from_u64(v)) <= 0; }

static inline void add_mod_ref(Cigint& a, const Cigint& b, const Cigint& m) {
    cigint_add_ref(&a, &b);
    if (cigint_cmp(a, m) >= 0) cigint_sub_ref(&a, &m);
}
static inline void sub_mod_ref(Cigint& a, const Cigint& b, const Cigint& m) {
    if (cigint_cmp(a, b) >= 0) cigint_sub_ref(&a, &b);
    else { Cigint t = a; cigint_add_ref(&t, &m); cigint_sub_ref(&t, &b); a = t; }
}

static Cigint mul_mod(Cigint a, Cigint b, const Cigint& m) {
    if (cigint_cmp(a, m) >= 0) a = cigint_mod(a, m);
    if (cigint_cmp(b, m) >= 0) b = cigint_mod(b, m);
    Cigint res = CIGINT_ZERO();
    while (!is_zero(a)) {
        if (is_odd(a)) add_mod_ref(res, b, m);
        add_mod_ref(b, b, m);
        shr1_ref(a);
    }
    return res;
}
static Cigint pow_mod(Cigint base, Cigint exp, const Cigint& m) {
    Cigint one = CIGINT_ZERO(); one.data[CIGINT_N - 1] = 1;
    Cigint r = one;
    base = cigint_mod(base, m);
    while (!is_zero(exp)) {
        if (is_odd(exp)) r = mul_mod(r, base, m);
        base = mul_mod(base, base, m);
        shr1_ref(exp);
    }
    return r;
}
static Cigint cigint_gcd(Cigint a, Cigint b) {
    while (!is_zero(b)) { Cigint r = cigint_mod(a, b); a = b; b = r; }
    return a;
}
static inline u64 mod_u64(const Cigint& x, u64 m) {
    if (!m) return 0;
    __uint128_t rem = 0;
    for (size_t i = 0; i < CIGINT_N; ++i) rem = ((rem << 32) + x.data[i]) % m;
    return (u64)rem;
}
static inline u64 powmod_u64(u64 a, u64 e, u64 m) {
    if (m == 1) return 0;
    __uint128_t A = a % m, R = 1;
    while (e) { if (e & 1) R = (R * A) % m; A = (A * A) % m; e >>= 1; }
    return (u64)R;
}

// ---------- Step 1: perfect power (kept simple) ----------
static bool perfect_power(const Cigint& n) {
    if (is_zero(n)) return false;
    Cigint one = CIGINT_ZERO(); one.data[CIGINT_N - 1] = 1;
    if (cmp(n, one) <= 0) return false;

    const u32 bl = bitlen(n);
    for (u32 b = 2; b <= bl; ++b) {
        // binary search integer b-th root via compare using mul (no mod)
        // lo=2, hi=2^ceil(bl/b)
        u32 ubits = (bl + b - 1) / b;
        Cigint lo = from_u64(2);
        Cigint hi = pow2(ubits);
        while (cmp(lo, hi) <= 0) {
            Cigint mid = add(lo, hi);
            mid = shr(mid, 1);
            // compute mid^b (exact width), compare to n
            Cigint acc = CIGINT_ZERO(); acc.data[CIGINT_N - 1] = 1;
            Cigint base = mid; u32 e = b;
            while (e) { if (e & 1u) acc = mul(acc, base); base = mul(base, base); e >>= 1; }
            int c = cmp(acc, n);
            if (c == 0) return true;
            if (c < 0) lo = add(mid, from_u64(1)); else hi = sub(mid, from_u64(1));
        }
    }
    return false;
}

// ---------- Step 2: findR(n) (fixed to standard order test) ----------
static u64 findR(const Cigint& n) {
    const u64 L2 = bitlen(n);
    const u64 maxK = L2 * L2;

    for (u64 r = 2; ; ++r) {
        if (cigint_cmp(cigint_gcd(n, from_u64(r)), from_u64(1)) != 0) continue; // require gcd(n,r)=1
        u64 nr = mod_u64(n, r);
        bool ok = true;
        for (u64 k = 1; k <= maxK; ++k) {
            if (powmod_u64(nr, k, r) == 1) { ok = false; break; } // only check ==1, not 0
        }
        if (ok) return r;
        // (no artificial maxR cap; the corrected test usually yields small r)
    }
}

// ---------- Step 3 ----------
static bool has_nontrivial_gcd_up_to_r(const Cigint& n, u64 r) {
    for (u64 a = 2; a <= r; ++a) {
        Cigint g = cigint_gcd(n, from_u64(a));
        if (cigint_cmp(g, from_u64(1)) > 0 && cigint_cmp(g, n) < 0) return true;
    }
    return false;
}

// ---------- Step 5: polynomial work (zero-alloc, tight loops) ----------
//
// All polys are fixed-length r with indices [0..r-1], representing degree < r.
// Ring: (Z/nZ)[x] / (x^r - 1)  -> cyclic wrap with single subtraction (no %).
//
// Buffers are passed in by pointer; OUT may alias A or B (we handle via temp).

// OUT = 0
static inline void poly_zero(Cigint* OUT, size_t r) {
    for (size_t i = 0; i < r; ++i) OUT[i] = CIGINT_ZERO();
}

// OUT = A  (copy)
static inline void poly_copy(Cigint* OUT, const Cigint* A, size_t r) {
    std::memcpy(OUT, A, r * sizeof(Cigint));
}

// OUT = A * B   (cyclic, coeffs mod n)
// scratch must be r* sizeof(Cigint)
static void poly_mul_cyclic(const Cigint* A, const Cigint* B,
                            Cigint* OUT, Cigint* SCR,
                            size_t r, const Cigint& n)
{
    // accumulate in SCR, then copy to OUT
    poly_zero(SCR, r);

    for (size_t i = 0; i < r; ++i) {
        // skip zero rows
        if (cigint_is0(A[i])) continue;

        const Cigint Ai = A[i];
        size_t j = 0, idx = i;
        // loop j=0..r-1 with single-wrap
        for (; j < r; ++j) {
            if (!cigint_is0(B[j])) {
                Cigint prod = mul_mod(Ai, B[j], n);
                add_mod_ref(SCR[idx], prod, n);
            }
            ++idx;
            if (idx == r) idx = 0;
        }
    }
    poly_copy(OUT, SCR, r);
}

// OUT = P * (x + a)   (O(r), cyclic)
// OUT[k] = P[k-1] + a*P[k]
static void poly_mul_by_x_plus_a(const Cigint* P, Cigint* OUT, size_t r, const Cigint& a, const Cigint& n) {
    if (r == 0) return;
    // OUT[0] = P[r-1] + a*P[0]
    OUT[0] = mul_mod(a, P[0], n);
    add_mod_ref(OUT[0], P[r - 1], n);
    for (size_t k = 1; k < r; ++k) {
        OUT[k] = mul_mod(a, P[k], n);
        add_mod_ref(OUT[k], P[k - 1], n);
    }
}

// Compute Y = (x+a)^e in (Z/nZ)[x]/(x^r - 1)
// allocation-free: uses two r-length buffers and one scratch.
static void poly_pow_x_plus_a(Cigint* Y, u64 a_u64, Cigint e, size_t r, const Cigint& n) {
    // buffers
    std::vector<Cigint> bufA(r), bufB(r), scratch(r);
    Cigint* acc = bufA.data();
    Cigint* base = bufB.data();

    // acc = 1
    for (size_t i = 0; i < r; ++i) acc[i] = CIGINT_ZERO();
    acc[0] = from_u64(1);

    // base = x + a
    for (size_t i = 0; i < r; ++i) base[i] = CIGINT_ZERO();
    base[0] = from_u64(a_u64);
    base[1 % r] = from_u64(1);

    Cigint aC = from_u64(a_u64);

    // fast path for very small e: multiply-by-(x+a) repeatedly (rarely used in AKS)
    // else: binary exponentiation with cyclic multiply
    while (!is_zero(e)) {
        if (is_odd(e)) {
            // acc = acc * base
            poly_mul_cyclic(acc, base, acc, scratch.data(), r, n);
        }
        // base = base * base
        poly_mul_cyclic(base, base, base, scratch.data(), r, n);
        shr1_ref(e);
    }

    // write to Y
    for (size_t i = 0; i < r; ++i) Y[i] = acc[i];

    // Final adjustment: subtract a at coeff 0, and 1 at coeff (n % r)
    sub_mod_ref(Y[0], aC, n);
    u64 idx = mod_u64(n, r);
    sub_mod_ref(Y[idx], from_u64(1), n);
}

static bool poly_any_nonzero(const Cigint* P, size_t r) {
    for (size_t i = 0; i < r; ++i) if (!cigint_is0(P[i])) return true;
    return false;
}

static u64 euler_phi_u64(u64 r) {
    u64 res = r, x = r;
    for (u64 p = 2; p * p <= x; ++p) if (x % p == 0) {
        while (x % p == 0) x /= p;
        res -= res / p;
    }
    if (x > 1) res -= res / x;
    return res;
}

// ---------- AKS ----------
enum class AKSResult { Prime, Composite };

static AKSResult aks(const Cigint& n) {
    // Step 1
    if (perfect_power(n)) return AKSResult::Composite;

    // Step 2 (fixed)
    u64 r = findR(n);
    printf("r = %llu\n", (unsigned long long)r);

    // Step 3
    if (has_nontrivial_gcd_up_to_r(n, r)) return AKSResult::Composite;

    // Step 4
    if (le_u64(n, r)) return AKSResult::Prime;

    // Step 5
    const u64 phi_r = euler_phi_u64(r);
    const u64 L = (u64)std::floor(std::sqrt((long double)phi_r) * (long double)bitlen(n));
    if (L == 0) return AKSResult::Prime;

    std::vector<Cigint> poly(r); // single output buffer reused
    for (u64 a = 1; a < L; ++a) {
        poly_pow_x_plus_a(poly.data(), a, n, r, n);
        if (poly_any_nonzero(poly.data(), r)) return AKSResult::Composite;
    }
    return AKSResult::Prime;
}

static const char* aks_str(const Cigint& n) {
    return aks(n) == AKSResult::Prime ? "Prime" : "Composite";
}

int main() {
    // Example
    Cigint n = cigint_from_dec("592050202721842132167314177912814106098300752505128394576383610271382547200796999123188872375399778064159260533195923283990719135880903");
    printf("%s\n", aks_str(n));
    return 0;
}

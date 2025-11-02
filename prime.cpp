#include <array>
#include <vector>
#include <cstdio>

#define CIGINT_IMPLEMENTATION
#define CIGINT_STRIP_PREFIX
#define CIGINT_N (512 / 32)
#include "cigint.h"

// we use the same small r as your Python
static constexpr int R = 14;

struct Poly {
    std::array<Cigint, R> c;
};

// -------- basic Cigint helpers --------
static inline Cigint z0() { return CIGINT_ZERO(); }
static inline Cigint z1() { return cigint_from_u32(1u); }

// res = (a + b) mod m (all fixed width)
static inline void add_mod(Cigint& res, const Cigint& b, const Cigint &m) {
    cigint_add_ref(&res, &b);
    if (cigint_cmp(res, m) >= 0) {
        cigint_sub_ref(&res, &m);
    }
}

// safely modularly multiply: (a * b) mod m
// does NOT rely on a 2N-limb product; stays in N limbs the whole time
static Cigint mul_mod_safe(Cigint a, Cigint b, const Cigint &m) {
    // reduce inputs first
    a = cigint_mod(a, m);
    b = cigint_mod(b, m);

    Cigint res = z0();

    // how many bits to scan in b?
    u32 bits = cigint_highest_order(b);  // 0 if b==0
    for (u32 i = 0; i < bits; ++i) {
        if (cigint_get_bit(b, i)) {
            add_mod(res, a, m);
        }
        // a = (a + a) mod m
        a = cigint_add(a, a);
        if (cigint_cmp(a, m) >= 0) {
            a = cigint_sub(a, m);
        }
    }
    return res;
}

// poly multiplication mod (x^R - 1, m)
// C[k] = sum_{i+j≡k mod R} A[i]*B[j] (mod m)
static Poly poly_mul(const Poly &A, const Poly &B, const Cigint &m) {
    Poly C;
    for (int i = 0; i < R; ++i) C.c[i] = z0();

    for (int i = 0; i < R; ++i) {
        if (cigint_is0(A.c[i])) continue;
        for (int j = 0; j < R; ++j) {
            if (cigint_is0(B.c[j])) continue;
            int idx = (i + j) % R;
            Cigint prod = mul_mod_safe(A.c[i], B.c[j], m);
            add_mod(C.c[idx], prod, m);
        }
    }
    return C;
}

// compute n mod R using your own bigint division
static u32 cigint_mod_small(const Cigint &n, u32 mod) {
    Cigint q;
    u32 r = 0;
    cigint_sdivmod(n, mod, &q, &r);  // r = n % mod
    return r;
}

// check poly == 1 + x^(n mod R)
static bool poly_is_xn_plus_1(const Poly &P, const Cigint &n, const Cigint &m) {
    u32 k = cigint_mod_small(n, R);   // exact n % R
    for (int i = 0; i < R; ++i) {
        if (i == 0 || i == (int)k) {
            if (cigint_cmp(P.c[i], z1()) != 0) return false;
        } else {
            if (!cigint_is0(P.c[i])) return false;
        }
    }
    return true;
}

// compute (1 + x)^n mod (x^R - 1, n)
static Poly poly_pow_1x(const Cigint &n, const Cigint &modn) {
    // base = 1 + x
    Poly base;
    for (int i = 0; i < R; ++i) base.c[i] = z0();
    base.c[0] = z1();
    base.c[1] = z1();

    // res = 1
    Poly res;
    for (int i = 0; i < R; ++i) res.c[i] = z0();
    res.c[0] = z1();

    // exponent = n (big)
    Cigint e = n;
    while (!cigint_is0(e)) {
        if (cigint_get_bit(e, 0)) {
            res = poly_mul(res, base, modn);
        }
        base = poly_mul(base, base, modn);
        e = cigint_shr(e, 1);  // e >>= 1
    }
    return res;
}

static bool aks_like_cigint(const Cigint &n) {
    // reject 0,1
    if (cigint_is0(n)) return false;
    if (cigint_cmp(n, z1()) == 0) return false;

    Poly p = poly_pow_1x(n, n);
    return poly_is_xn_plus_1(p, n, n);
}

int main() {
    // your prime
    // Cigint n = cigint_from_dec("7669373798138436444425917457853162445620122597663674808454541043947311884377"); //up to 76digit
    Cigint n = cigint_from_dec("1897964784648536816215132903105500768398446447110886849902932658628213742067511093681881833444371252205590663430767842679721896956834629961605331435987959"); //up to 154 digit
    printf("n = ");
    cigint_print10(n);
    bool ok = aks_like_cigint(n);
    printf(" -> %s\n", ok ? "maybe prime" : "composite-ish");
    return 0;
}

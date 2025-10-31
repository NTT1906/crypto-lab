// aks_512.c
// AKS primality test using the 512-bit Cigint implementation
// requires: CIGINT_IMPLEMENTATION, CIGINT_STRIP_PREFIX, CIGINT_N = 16 (512/32)
// compile together with your cigint.h
//
// NOTE: this is the straight, slow, textbook-style AKS for your fixed-width bigint.
// It uses bigint-safe (a*b) mod n via shift-add, and polynomial arithmetic
// (mod x^r - 1, n) with a fixed max r. It's meant to be correct, not fast.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define AKS_MAX_R 2048    /* must be >= r we find; raise if needed */

#define CIGINT_IMPLEMENTATION
#define CIGINT_STRIP_PREFIX
#define CIGINT_N (512 / 32)
#include "cigint.h"

/* ---------- small helpers on Cigint ---------- */

static inline Cigint cigint_from_u32(u32 x) {
    Cigint t = CIGINT_ZERO();
    t.data[CIGINT_N - 1] = x;
    return t;
}

static inline int cigint_is_one(const Cigint *a) {
    for (size_t i = 0; i + 1 < CIGINT_N; ++i)
        if (a->data[i] != 0) return 0;
    return a->data[CIGINT_N - 1] == 1;
}

static inline int cigint_eq_u32(const Cigint *a, u32 v) {
    for (size_t i = 0; i + 1 < CIGINT_N; ++i)
        if (a->data[i] != 0) return 0;
    return a->data[CIGINT_N - 1] == v;
}

static inline int cigint_le_u32(const Cigint *a, u32 v) {
    /* a <= v ? */
    for (size_t i = 0; i + 1 < CIGINT_N; ++i) {
        if (a->data[i] != 0) return 0;
    }
    return a->data[CIGINT_N - 1] <= v;
}

static inline u32 cigint_mod_u32(const Cigint *a, u32 m) {
    Cigint q;
    u32 r = 0;
    cigint_sdivmod(*a, m, &q, &r);
    return r;
}

/* gcd(a,b) in big-int domain */
static Cigint cigint_gcd(Cigint a, Cigint b) {
    while (!cigint_is0(b)) {
        Cigint r = cigint_mod(a, b);
        a = b;
        b = r;
    }
    return a;
}

/* (a + b) % mod, with 0 <= a,b < mod, so at most one sub */
static inline void cigint_add_mod_inplace(Cigint *a, const Cigint *b, const Cigint *mod) {
    cigint_add_ref(a, b);
    if (cigint_cmp(*a, *mod) >= 0) {
        cigint_sub_ref(a, mod);
    }
}

/* x = 2*x mod mod, assuming x < mod */
static inline void cigint_double_mod_inplace(Cigint *x, const Cigint *mod) {
    cigint_add_ref(x, x);
    if (cigint_cmp(*x, *mod) >= 0) cigint_sub_ref(x, mod);
}

/* (a * b) % mod, bigint-safe, no 1024-bit intermediate:
 * Russian peasant / double-and-add modulo mod.
 * cost ~ (#bits of b) * (1 mul-by-2 + maybe 1 add) = OK for 512-bit.
 */
static Cigint cigint_mul_mod(const Cigint *a, const Cigint *b, const Cigint *mod) {
    Cigint res = CIGINT_ZERO();
    Cigint x = *a;
    Cigint y = *b;

    while (!cigint_is0(y)) {
        if (y.data[CIGINT_N - 1] & 1u) {
            cigint_add_mod_inplace(&res, &x, mod);
        }
        cigint_double_mod_inplace(&x, mod);
        y = cigint_shr(y, 1);
    }
    return res;
}

/* ---------- polynomials mod (x^r - 1, n) ---------- */

typedef struct {
    /* coeffs[0] + coeffs[1] x + ... + coeffs[r-1] x^{r-1}, all mod n */
    Cigint c[AKS_MAX_R];
} aks_poly;

/* p = 0 */
static inline void poly_zero(aks_poly *p, u32 r) {
    for (u32 i = 0; i < r; ++i) p->c[i] = CIGINT_ZERO();
}

/* p = 1 (constant) */
static inline void poly_one(aks_poly *p, u32 r) {
    poly_zero(p, r);
    p->c[0] = cigint_from_u32(1);
}

/* p = x + a (mod n) */
static inline void poly_x_plus_a(aks_poly *p, u32 r, u32 a, const Cigint *n) {
    poly_zero(p, r);
    /* coeff 0 = a mod n */
    Cigint ca = cigint_from_u32(a);
    if (cigint_cmp(ca, *n) >= 0) cigint_sub_ref(&ca, n);
    p->c[0] = ca;
    /* coeff 1 = 1 */
    p->c[1] = cigint_from_u32(1);
}

/* res = p * q mod (x^r - 1, n) */
static void poly_mul_mod(aks_poly *res,
                         const aks_poly *p,
                         const aks_poly *q,
                         u32 r,
                         const Cigint *n)
{
    aks_poly out;
    poly_zero(&out, r);

    for (u32 i = 0; i < r; ++i) {
        if (cigint_is0(p->c[i])) continue;
        for (u32 j = 0; j < r; ++j) {
            if (cigint_is0(q->c[j])) continue;
            u32 k = i + j;
            if (k >= r) k -= r;
            Cigint prod = cigint_mul_mod(&p->c[i], &q->c[j], n);
            cigint_add_mod_inplace(&out.c[k], &prod, n);
        }
    }
    *res = out;
}

/* res = base^exp mod (x^r - 1, n) */
static void poly_pow_mod(aks_poly *res,
                         const aks_poly *base,
                         const Cigint *exp,
                         u32 r,
                         const Cigint *n)
{
    aks_poly result;
    poly_one(&result, r);

    aks_poly b = *base;
    Cigint e = *exp;

    while (!cigint_is0(e)) {
        if (e.data[CIGINT_N - 1] & 1u) {
            aks_poly tmp;
            poly_mul_mod(&tmp, &result, &b, r, n);
            result = tmp;
        }
        e = cigint_shr(e, 1);
        if (!cigint_is0(e)) {
            aks_poly tmp;
            poly_mul_mod(&tmp, &b, &b, r, n);
            b = tmp;
        }
    }
    *res = result;
}

/* rhs = x^ (n mod r) + a  (mod n) */
static void poly_rhs_make(aks_poly *rhs,
                          const Cigint *n,
                          u32 r,
                          u32 a,
                          const Cigint *mod_n)
{
    poly_zero(rhs, r);
    /* x^{n mod r} */
    u32 nm = cigint_mod_u32(n, r);
    rhs->c[nm] = cigint_from_u32(1);
    /* + a */
    Cigint ca = cigint_from_u32(a);
    if (cigint_cmp(ca, *mod_n) >= 0) cigint_sub_ref(&ca, mod_n);
    rhs->c[0] = ca;
}

/* returns 1 if p == q (coeff-wise), else 0 */
static int poly_eq(const aks_poly *p, const aks_poly *q, u32 r) {
    for (u32 i = 0; i < r; ++i) {
        if (cigint_cmp(p->c[i], q->c[i]) != 0) return 0;
    }
    return 1;
}

/* ---------- AKS core ---------- */

/* step 1: (very light) check: if n is 0,1 -> composite; also small trial division */
static int aks_basic_checks(const Cigint *n) {
    if (cigint_is0(*n)) return 0;
    if (cigint_is_one(n)) return 0;

    /* small primes 2,3,5,7,... */
    static const u32 small_primes[] = {
        2,3,5,7,11,13,17,19,23,29,31,37,41,43,47
    };
    for (size_t i = 0; i < sizeof(small_primes)/sizeof(small_primes[0]); ++i) {
        u32 p = small_primes[i];
        if (cigint_eq_u32(n, p)) return 1;
        u32 r = cigint_mod_u32(n, p);
        if (r == 0) return 0;
    }
    return -1; /* unknown */
}

/* find the smallest r s.t. ord_r(n) > (log2 n)^2 */
static u32 aks_find_r(const Cigint *n) {
    u32 logn = cigint_highest_order(*n);
    u32 maxk = logn * logn;

    for (u32 r = 2; ; ++r) {
        if (r >= AKS_MAX_R) return r; /* cap, caller should check */

        int ok = 1;
        u32 n_mod_r = cigint_mod_u32(n, r);
        if (n_mod_r == 0) { ok = 0; goto cont; }

        {
            /* compute ord_r(n) */
            u32 x = n_mod_r;
            u32 ord = 1;
            u64 cur = x % r;
            while (cur != 1 && ord <= maxk) {
                cur = (cur * x) % r;
                ord++;
            }
            if (cur != 1) ok = 0;
            if (ok && ord <= maxk) ok = 0;
        }

    cont:
        if (ok) return r;
    }
}

/* gcd(n, a) for 2 <= a <= r; if 1 < gcd < n -> composite */
static int aks_gcd_step(const Cigint *n, u32 r) {
    for (u32 a = 2; a <= r; ++a) {
        Cigint ga = cigint_gcd(*n, cigint_from_u32(a));
        if (!cigint_is_one(&ga) && !cigint_eq_u32(&ga, 0) && !cigint_cmp(*n, ga)) {
            /* ga == n means a == n, ignore */
        }
        if (!cigint_is_one(&ga) && cigint_cmp(ga, *n) != 0) {
            /* 1 < gcd(n,a) < n */
            return 0;
        }
    }
    return 1;
}

/* main entry: 1 = prime, 0 = composite, -1 = "r too big / not handled" */
int aks_is_prime_512(const Cigint *n) {
    int bc = aks_basic_checks(n);
    if (bc == 1) return 1;
    if (bc == 0) return 0;

    /* step 2: find r */
    u32 r = aks_find_r(n);
    if (r >= AKS_MAX_R) {
        /* need to enlarge AKS_MAX_R */
        return -1;
    }

    /* step 3: check gcd(n, a) for a=2..r */
    if (!aks_gcd_step(n, r)) return 0;

    /* step 4: if n <= r, n is prime */
    if (cigint_le_u32(n, r)) return 1;

    /* step 5: polynomial check */
    u32 logn = cigint_highest_order(*n);
    /* phi(r) >= sqrt(r) heuristic: use sqrt(r) * logn */
    double maxA_d = floor(sqrt((double)r) * (double)logn);
    u32 maxA = (u32)maxA_d;
    if (maxA < 1) maxA = 1;

    for (u32 a = 1; a <= maxA; ++a) {
        aks_poly base, lhs, rhs;
        poly_x_plus_a(&base, r, a, n);
        poly_pow_mod(&lhs, &base, n, r, n);
        poly_rhs_make(&rhs, n, r, a, n);
        if (!poly_eq(&lhs, &rhs, r)) {
            return 0;
        }
    }
    return 1;
}

int main() {
     Cigint n = CIGINT_ZERO();
     // n = some 512-bit number, for demo just set small:
     n.data[CIGINT_N - 1] = 101; // 101 is prime
     int r = aks_is_prime_512(&n);
     printf("result = %d\n", r);
     return 0;
}

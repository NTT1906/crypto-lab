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

#include <iostream>

// anti-overflow add mod
// u32 addmod(u32 a, u32 b, u32 m) {
//     if (a >= m) a %= m;
//     if (b >= m) b %= m;
//     u32 c = a + b;
//     if (c < a) c -= m;
//     if (c >= m) c -= m;
//     return c;
// }

void cigint_add_wrap(Cigint *lhs, const Cigint *rhs) {
    u64 carry = 0;
    for (size_t i = CIGINT_N; i-- > 0;) {
        carry = (u64) lhs->data[i] + (u64) rhs->data[i] + carry;
        lhs->data[i] = (u32) carry;
        carry >>= 32;
    }
}

// res = (a + b) mod m (all fixed width)
static inline void add_mod(Cigint& res, Cigint& b, const Cigint &m) {
	Cigint res_c = res;
    cigint_divmod_ref(&res, &m, NULL, &res);
    cigint_divmod_ref(&b, &m, NULL, &b);
    cigint_add_wrap(&res, &b);
    if (res < res_c) {
        cigint_sub_ref(&res, &m);
        // cigint_printf("res = %Cd\n", res);
        // std::cout << (std::string("add_mod overflow\nres =") + res.toDecStr() + "\nresc=" + res_c.toDecStr()).c_str() << "\n";
        // exit(1);
        // assert(false && "add_mod overflow");
    }
    if (cigint_cmp(res, m) >= 0) {
        cigint_sub_ref(&res, &m);
    }
}

static Cigint mul_mod_safe(Cigint a, Cigint b, const Cigint &m) {
    a = cigint_mod(a, m);
    b = cigint_mod(b, m);
    Cigint res = z0();
    u32 bits = cigint_highest_order(b);
    for (u32 i = 0; i < bits; ++i) {
        if (cigint_get_bit(b, i)) {
            add_mod(res, a, m);
        }
        add_mod(a, a, m);
    }
    return res;
}

static inline u32 clz32(u32 x) {
#if defined(__GNUC__) || defined(__clang__)
    return x ? (u32)__builtin_clz(x) : 32u;
#else
    /* portable fallback */
    u32 n = 0; if (!x) return 32;
    if ((x & 0xFFFF0000u) == 0) { n += 16; x <<= 16; }
    if ((x & 0xFF000000u) == 0) { n +=  8; x <<=  8; }
    if ((x & 0xF0000000u) == 0) { n +=  4; x <<=  4; }
    if ((x & 0xC0000000u) == 0) { n +=  2; x <<=  2; }
    if ((x & 0x80000000u) == 0) { n +=  1; }
    return n;
#endif
}

/* ---- MSW<->LSW helpers (we keep scratch in LSW-first) ---- */
static inline void msw_to_lsw(const u32 *msw, u32 *lsw, size_t n) {
    for (size_t i = 0; i < n; ++i) lsw[i] = msw[n - 1 - i];
}
static inline void lsw_to_msw(const u32 *lsw, u32 *msw, size_t n) {
    for (size_t i = 0; i < n; ++i) msw[n - 1 - i] = lsw[i];
}

/* ---- left/right shift by s bits (0..31) on LSW-first arrays ---- */
static void lshift_bits(u32 *x, size_t n, u32 s) {
    if (s == 0) return;
    u32 carry = 0;
    for (size_t i = 0; i < n; ++i) {
        u32 v = x[i];
        x[i] = (v << s) | carry;
        carry = (s == 0) ? 0 : (v >> (32 - s));
    }
}
static u32 lshift_bits_with_out(u32 *x, size_t n, u32 s) {
    if (s == 0) return 0;
    u32 carry = 0;
    for (size_t i = 0; i < n; ++i) {
        u32 v = x[i];
        x[i] = (v << s) | carry;
        carry = (v >> (32 - s));
    }
    return carry; /* overflow word */
}
static void rshift_bits(u32 *x, size_t n, u32 s) {
    if (s == 0) return;
    u32 carry = 0;
    for (size_t i = n; i-- > 0;) {
        u32 v = x[i];
        x[i] = (v >> s) | (carry << (32 - s));
        carry = v & ((1u << s) - 1u);
    }
}

/* ---- schoolbook 2N-limb multiply: U = A * B (LSW-first) ---- */
static void mul_full_2n(const u32 *A, const u32 *B, u32 *U, size_t n) {
    memset(U, 0, sizeof(u32) * (2 * n));
    for (size_t i = 0; i < n; ++i) {
        u64 carry = 0;
        u64 ai = (u64)A[i];
        for (size_t j = 0; j < n; ++j) {
            u64 t = (u64)U[i + j] + ai * (u64)B[j] + carry;
            U[i + j] = (u32)t;
            carry = t >> 32;
        }
        U[i + n] += (u32)carry; /* carry fits in 32 bits */
    }
}

/* ---- compare two LSW-first arrays of length n ---- */
static int cmp_n(const u32 *x, const u32 *y, size_t n) {
    for (size_t i = n; i-- > 0;) {
        if (x[i] < y[i]) return -1;
        if (x[i] > y[i]) return 1;
    }
    return 0;
}

/* ---- x = x - y*q, y length n, q is 32-bit; return borrow (0/1) ---- */
static u32 submul_1(u32 *x, const u32 *y, size_t n, u32 q) {
    u64 borrow = 0;
    u64 qq = (u64)q;
    for (size_t i = 0; i < n; ++i) {
        u64 yiq = (u64)y[i] * qq;
        u64 cur = (u64)x[i];
        u64 t = cur - (yiq & 0xFFFFFFFFull) - borrow;
        x[i] = (u32)t;
        borrow = ((yiq >> 32) + ((t >> 32) & 1ull) /* 1 if underflow */);
    }
    /* Propagate borrow into the next word */
    u64 cur = (u64)x[n];
    u64 t = cur - borrow;
    x[n] = (u32)t;
    return (u32)((t >> 32) & 1ull); /* 1 if underflow after last subtract */
}

/* ---- x = x + y; y length n; return carry ---- */
static u32 add_n(u32 *x, const u32 *y, size_t n) {
    u64 carry = 0;
    for (size_t i = 0; i < n; ++i) {
        u64 t = (u64)x[i] + (u64)y[i] + carry;
        x[i] = (u32)t;
        carry = t >> 32;
    }
    return (u32)carry;
}

/* ---- remainder R = U mod V (V length n, U length 2n), LSW-first.
 *      Knuth Algorithm D (32-bit base). Only returns remainder.
 * -------------------------------------------------------------- */
static void mod_knuth_rem(const u32 *U2n, const u32 *Vn, u32 *R, size_t n) {
    /* handle V == 0? (undefined) */
    int all_zero = 1;
    for (size_t i = 0; i < n; ++i) if (Vn[i]) { all_zero = 0; break; }
    assert(!all_zero && "mod by zero");

    /* copy inputs because we'll normalize/modify */
    u32 *U = (u32*)alloca(sizeof(u32) * (2 * n + 1));
    u32 *V = (u32*)alloca(sizeof(u32) * n);
    memcpy(U, U2n, sizeof(u32) * (2 * n));
    U[2 * n] = 0;
    memcpy(V, Vn, sizeof(u32) * n);

    /* normalization shift s so that top limb of V has highest bit set */
    u32 vtop = V[n - 1];
    u32 s = (vtop == 0) ? 0u : clz32(vtop);
    if (s) {
        u32 carryU = lshift_bits_with_out(U, 2 * n, s);
        U[2 * n] = carryU;
        lshift_bits(V, n, s);
    }

    u32 v1 = V[n - 1];
    u32 v0 = (n >= 2) ? V[n - 2] : 0;

    /* Main loop: j = (2n - n) .. 0 => j = n .. 0 */
    for (size_t jj = n + 1; jj-- > 1; ) {
        size_t j = jj - 1;               /* j in [0..n] */
        u64 u2 = (u64)U[j + n];          /* u_{j+n} */
        u64 u1 = (u64)U[j + n - 1];      /* u_{j+n-1} */
        u64 u0 = (j + n >= 2) ? (u64)U[j + n - 2] : 0;

        /* qhat = ((u2*b + u1) / v1), rhat = ((u2*b + u1) % v1) */
        const u64 B = 1ull << 32;
        u64 num = (u2 << 32) | u1;
        u64 qhat = num / v1;
        u64 rhat = num - qhat * v1;
        if (qhat >= B) { qhat = B - 1; rhat = num - qhat * v1; }

        /* correction: while qhat*v0 > rhat*B + u0 -> qhat--, rhat += v1 */
        while (qhat * (u64)v0 > ((rhat << 32) | u0)) {
            --qhat;
            rhat += v1;
            if (rhat >= B) break;
        }

        /* subtract qhat * V from U segment starting at j */
        u32 q = (u32)qhat;
        u32 borrow = submul_1(&U[j], V, n, q);

        if (borrow) {
            /* add back V (i.e., U[j..j+n] += V[0..n-1]) and decrement q */
            u32 c = add_n(&U[j], V, n);
            U[j + n] += c;
            /* q-- ; not needed for remainder */
        }
        /* U[j + n] holds the next word; loop continues */
    }

    /* remainder sits in U[0..n-1], denormalize */
    if (s) rshift_bits(U, n, s);
    memcpy(R, U, sizeof(u32) * n);
}

/* ============================================================
 *  Public: in-place fast mul_mod
 *  a = (a * b) % m
 * ============================================================ */
static inline void cigint_mul_mod_ref_fast(Cigint *a, const Cigint *b, const Cigint *m) {
    const size_t N = CIGINT_N;

    /* Early outs */
    /* if m == 0 -> undefined; assert */
    int mz = 1; for (size_t i = 0; i < N; ++i) if (m->data[i]) { mz = 0; break; }
    assert(!mz && "mod by zero");

    /* LSW scratch */
    u32 *A = (u32*)alloca(sizeof(u32) * N);
    u32 *B = (u32*)alloca(sizeof(u32) * N);
    u32 *M = (u32*)alloca(sizeof(u32) * N);
    u32 *U = (u32*)alloca(sizeof(u32) * (2 * N));
    u32 *R = (u32*)alloca(sizeof(u32) * N);

    msw_to_lsw(a->data, A, N);
    msw_to_lsw(b->data, B, N);
    msw_to_lsw(m->data, M, N);

    /* full product: U = A * B (2N limbs, LSW-first) */
    mul_full_2n(A, B, U, N);

    /* remainder R = U mod M */
    mod_knuth_rem(U, M, R, N);

    /* write back MSW-first */
    lsw_to_msw(R, a->data, N);
}

Cigint cigint_mul_2_mod(Cigint lhs, Cigint rhs, Cigint m) {
    cigint_mul_mod_ref_fast(&lhs, &rhs, &m);
    return lhs;
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
            // Cigint prod2 = cigint_mul_2_mod(A.c[i], B.c[j], m);
            // if (cigint_cmp(prod, prod2) != 0) {
                // printf("Check point");
            // }
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
    // Cigint n = cigint_from_dec("7669373798138436444425917457853162445620122597663674808454541043947311884377"); //up to 76digit
    Cigint n = cigint_from_dec("771458402252315418489174122874340496560720991135040845033515625723932549727443357205751720461233037696609749135759573352834404574203828395431821751181407"); //up to 154 digit
    // Cigint n = cigint_from_dec("12598680682438752944055149498662865476737841742262510750565369004945449194176119763492007280023309651826310840991041734855800903377926420446403888376677867"); //up to 154 digit

    // Cigint a = cigint_from_dec("229911617100");
    // Cigint b = cigint_from_dec("1199165601554601993");
    // Cigint r = cigint_mul_2_mod(a, b, n);
    // std::string r_str = r.toDecStr();
    // printf("r = %s\n", r_str.c_str());
    printf("n = ");
    cigint_print10(n);
    bool ok = aks_like_cigint(n);
    printf(" -> %s\n", ok ? "maybe prime" : "composite-ish");
    return 0;
}
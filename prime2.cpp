#include <array>
#include <vector>
#include <cstdio>
#include <cstdint>

// Cigint configuration
#define CIGINT_IMPLEMENTATION
#define CIGINT_STRIP_PREFIX
#define CIGINT_N (1024 / 32)
#include "cigint.h"

// ===== Parameters =====
// keep the same small r as you used in Python
static constexpr int R = 14;

// ===== Small helpers on Cigint =====
static inline Cigint make_zero() {
    Cigint z = CIGINT_ZERO();
    return z;
}
static inline Cigint make_one() {
    Cigint o = CIGINT_ZERO();
    o.data[CIGINT_N - 1] = 1u;
    return o;
}
static inline bool is_zero_c(const Cigint &x) {
    for (size_t i = 0; i < CIGINT_N; ++i) if (x.data[i] != 0u) return false;
    return true;
}
static inline void reduce_mod(Cigint *x, const Cigint *m) {
    // x = x % m  (in place)
    cigint_divmod_ref(x, m, NULL, x);
}
static inline void add_mod_inplace(Cigint *dst, const Cigint *src, const Cigint *m) {
    cigint_add_ref(dst, src);
    // ensure 0 <= dst < m
    if (cigint_cmp(*dst, *m) >= 0) {
        cigint_sub_ref(dst, m);
        // still possible to be >= m if src was large; ensure with a full reduce:
        if (cigint_cmp(*dst, *m) >= 0) reduce_mod(dst, m);
    }
}

// safe word-wise (a * b) % m that never truncates the 2N product
static void mul_mod_wordwise(Cigint *out, const Cigint *a, const Cigint *b, const Cigint *m) {
    // A' = a % m, B' = b % m  (keep operands bounded)
    Cigint A = *a, B = *b;
    reduce_mod(&A, m);
    reduce_mod(&B, m);

    Cigint res = make_zero();

    // accumulate each limb of B: res += (A * bj) << (32*j)   (mod m)
    for (size_t j = 0; j < CIGINT_N; ++j) {
        uint32_t bj = B.data[CIGINT_N - 1 - j]; // LSW is at highest index
        if (!bj) continue;

        // T = A * bj (1024-bit by 32-bit), column-wise with 64-bit carry
        Cigint T = make_zero();
        uint64_t carry = 0;
        for (size_t i = 0; i < CIGINT_N; ++i) {
            size_t idx = CIGINT_N - 1 - i;
            uint64_t acc = (uint64_t)A.data[idx] * (uint64_t)bj + (uint64_t)T.data[idx] + carry;
            T.data[idx] = (uint32_t)acc;
            carry = acc >> 32;
        }
        // shift by whole limbs (32*j bits)
        Cigint Tshift = cigint_shl(T, (uint32_t)(32 * j));

        // res = (res + Tshift) % m
        cigint_add_ref(&res, &Tshift);
        if (cigint_cmp(res, *m) >= 0) cigint_sub_ref(&res, m);
        reduce_mod(&res, m);
    }

    *out = res;
}

// ===== Polynomials modulo (x^R - 1, m) =====
struct Poly {
    std::array<Cigint, R> c; // coefficients 0..R-1
};

static inline Poly poly_zero() {
    Poly P;
    for (int i = 0; i < R; ++i) P.c[i] = make_zero();
    return P;
}

static inline Poly poly_one() {
    Poly P = poly_zero();
    P.c[0] = make_one();
    return P;
}

static void poly_mul_mod(Poly *A, const Poly *B, const Cigint *m) {
    Poly C = poly_zero();

    for (int i = 0; i < R; ++i) {
        // small optimization: skip zero rows
        bool ai_zero = is_zero_c(A->c[i]);
        if (ai_zero) continue;

        for (int j = 0; j < R; ++j) {
            if (is_zero_c(B->c[j])) continue;

            const int idx = (i + j) % R;

            // term = (A[i] * B[j]) % m  (do NOT mutate A/B in-place)
            Cigint term;
            mul_mod_wordwise(&term, &A->c[i], &B->c[j], m);

            // C[idx] = (C[idx] + term) % m
            add_mod_inplace(&C.c[idx], &term, m);
        }
    }

    *A = C;
}

// exponentiate P(x) = (1 + x)^n modulo (x^R - 1, m)
static Poly poly_pow_1x(const Cigint &n, const Cigint &m) {
    Poly base = poly_zero();
    base.c[0] = make_one();
    base.c[1] = make_one();

    Poly acc = poly_one();

    // iterate over bits of n: LSB -> MSB
    uint32_t bits = cigint_highest_order(n);
    for (uint32_t k = 0; k < bits; ++k) {
        if (cigint_get_bit_ref(&n, (int)k)) {
            poly_mul_mod(&acc, &base, &m);
        }
        // base = base^2 mod (x^R - 1, m)
        Poly tmp = base;         // copy because poly_mul_mod overwrites lhs
        poly_mul_mod(&tmp, &base, &m);
        base = tmp;
    }

    return acc;
}

// n mod R using small-div (remainder only)
static uint32_t cigint_mod_small(const Cigint &n, uint32_t mod) {
    uint32_t r = 0;
    cigint_sdivmod(n, mod, NULL, &r);
    return r;
}

// check P(x) == 1 + x^{n mod R}  (coefficients reduced mod n)
static bool poly_is_xn_plus_1(const Poly &P, const Cigint &n, const Cigint &m) {
    (void)m; // unused, kept for symmetry
    const uint32_t r = cigint_mod_small(n, R);
    const Cigint ZERO = make_zero();
    const Cigint ONE  = make_one();

    for (int i = 0; i < R; ++i) {
        if (i == 0) {
            if (cigint_cmp(P.c[i], ONE) != 0) return false;
        } else if (i == (int)r) {
            if (cigint_cmp(P.c[i], ONE) != 0) return false;
        } else {
            if (!is_zero_c(P.c[i])) return false;
        }
    }
    return true;
}

// top-level: "AKS-like" check with small fixed R
static bool aks_like_cigint(const Cigint &n) {
    // 1) Compute P(x) = (1+x)^n mod (n, x^R-1)
    // 2) Check P(x) == 1 + x^{n mod R}
    Poly p = poly_pow_1x(n, n);
    return poly_is_xn_plus_1(p, n, n);
}

int main() {
	// Cigint n = cigint_from_dec("7669373798138436444425917457853162445620122597663674808454541043947311884377"); //up to 76digit
    Cigint n = cigint_from_dec("592050202721842132167314177912814106098300752505128394576383610271382547200796999123188872375399778064159260533195923283990719135880903");
    std::printf("n = ");
    cigint_print10(n);
    bool ok = aks_like_cigint(n);
    std::printf(" -> %s\n", ok ? "maybe prime" : "composite-ish");
    return 0;
}

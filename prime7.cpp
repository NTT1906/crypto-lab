#include <array>
#include <vector>
#include <cstdio>

#define CIGINT_IMPLEMENTATION
#define CIGINT_STRIP_PREFIX
#define CIGINT_N (512 / 32)
#include "cigint.h"

// we use the same small r as your Python
static constexpr int R = 14;

struct Poly{
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
	}
	if (cigint_cmp(res, m) >= 0) {
		cigint_sub_ref(&res, &m);
	}
}

// static Cigint mul_mod_safe(Cigint a, Cigint b, const Cigint &m) {
// 	// cigint_divmod_ref(&a, &m, NULL, &a);
// 	// cigint_divmod_ref(&b, &m, NULL, &b);
// 	Cigint res = z0();
// 	u32 bits = cigint_highest_order(b);
// 	for (u32 i = 0; i < bits; ++i) {
// 		if (cigint_get_bit(b, i)) {
// 			add_mod(res, a, m);
// 		}
// 		add_mod(a, a, m);
// 	}
// 	return res;
// }

static inline u32 cigint_add_wrap_carry(Cigint *lhs, const Cigint *rhs) {
    u64 carry = 0;
    for (size_t i = CIGINT_N; i-- > 0;) {
        u64 t = (u64)lhs->data[i] + (u64)rhs->data[i] + carry;
        lhs->data[i] = (u32)t;
        carry = t >> 32;
    }
    return (u32)carry;
}

/* res = (res + b) mod m
 * PRE: 0 <= res < m, 0 <= b < m
 * COST: one limb add + at most one subtract
 */
static inline void add_mod_fast(Cigint &res, const Cigint &b, const Cigint &m) {
    u32 c = cigint_add_wrap_carry(&res, &b);
    if (c || cigint_cmp(res, m) >= 0) cigint_sub_ref(&res, &m);
}

/* x = 2x mod m
 * PRE: 0 <= x < m
 */
static inline void dbl_mod_fast(Cigint &x, const Cigint &m) {
    u32 carry = 0;
    for (size_t i = CIGINT_N; i-- > 0;) {
        u32 v = x.data[i];
        u32 nv = (v << 1) | carry;
        carry = v >> 31;
        x.data[i] = nv;
    }
    if (carry || cigint_cmp(x, m) >= 0) cigint_sub_ref(&x, &m);
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
	u32 k = cigint_mod_small(n, R); // exact n % R
	for (int i = 0; i < R; ++i) {
		if (i == 0 || i == (int)k) {
			if (cigint_cmp(P.c[i], z1()) != 0) return false;
		} else {
			if (!cigint_is0(P.c[i])) return false;
		}
	}
	return true;
}


/* a,b < m required in the trusted version (use safe version if unsure) */
static inline Cigint mul_mod_win4_trusted(Cigint a, const Cigint &b, const Cigint &m) {
	if (cigint_is0(a) || cigint_is0(b)) return CIGINT_ZERO();

	Cigint T[16]; for (int i = 0; i < 16; ++i) T[i] = CIGINT_ZERO();
	T[1] = a;
	T[2] = a; dbl_mod_fast(T[2], m);
	T[4] = T[2]; dbl_mod_fast(T[4], m);
	T[8] = T[4]; dbl_mod_fast(T[8], m);
	for (int d = 3; d < 16; ++d) if (d != 4 && d != 8) {
		Cigint acc = CIGINT_ZERO();
		if (d & 1) add_mod_fast(acc, T[1], m);
		if (d & 2) add_mod_fast(acc, T[2], m);
		if (d & 4) add_mod_fast(acc, T[4], m);
		if (d & 8) add_mod_fast(acc, T[8], m);
		T[d] = acc;
	}

	Cigint res = CIGINT_ZERO();
	u32 bits = cigint_highest_order(b);
	if (!bits) return res;
	u32 total_nibbles = (bits + 3u) / 4u;
	u32 skip = 8u * (u32)CIGINT_N - total_nibbles;

	u32 seen = 0; bool first = true;
	for (size_t i = 0; i < CIGINT_N; ++i) {
		u32 v = b.data[i];
		for (int k = 0; k < 8; ++k) {
			if (seen++ < skip) continue;
			u32 d = (v >> (28 - 4*k)) & 0xFu;
			if (first) { if (d) res = T[d]; first = false; }
			else {
				dbl_mod_fast(res, m); dbl_mod_fast(res, m);
				dbl_mod_fast(res, m); dbl_mod_fast(res, m);
				if (d) add_mod_fast(res, T[d], m);
			}
		}
	}
	return first ? CIGINT_ZERO() : res;
}

/* C[k] = sum_{i+j≡k (mod R)} A[i]*B[j] (mod m), with every coeff kept < m */
static inline int add_modR(int t) { return (t >= R) ? (t - R) : t; }

static Poly poly_mul_fast(const Poly &A, const Poly &B, const Cigint &m) {
	Poly C; for (int i = 0; i < R; ++i) C.c[i] = CIGINT_ZERO();

	for (int i = 0; i < R; ++i) {
		if (cigint_is0(A.c[i])) continue;
		for (int j = 0; j < R; ++j) {
			if (cigint_is0(B.c[j])) continue;
			const int idx = add_modR(i + j);
			Cigint prod = mul_mod_win4_trusted(A.c[i], B.c[j], m);
			add_mod_fast(C.c[idx], prod, m);
		}
	}
	return C;
}

#include <cmath>
using namespace std;

// === NTT + CRT option for poly convolution modulo (x^R - 1, n) ===
// Requirements: C++17, your cigint.h already included/compiled.
// Correct IFF product of chosen primes M >= R*(n-1)^2. Otherwise we fall back.

namespace nttcrt {

// ---------- small mod helpers ----------
static inline uint32_t umod(uint64_t x, uint32_t m) { return (uint32_t)(x % m); }
static inline uint32_t addm(uint32_t a, uint32_t b, uint32_t m){ uint32_t s=a+b; return s>=m? s-m : s; }
static inline uint32_t subm(uint32_t a, uint32_t b, uint32_t m){ return a>=b? a-b : a+ (m-b); }

static uint32_t powmod(uint32_t a, uint32_t e, uint32_t m){
    uint64_t r=1, x=a;
    while(e){ if(e&1) r=(r*x)%m; x=(x*x)%m; e>>=1; }
    return (uint32_t)r;
}
static uint32_t invmod(uint32_t a, uint32_t m){ return powmod(a, m-2, m); } // m prime

// ---------- bit-reverse iterative NTT (length power-of-two) ----------
static void ntt(std::vector<uint32_t>& a, bool invert, uint32_t mod, uint32_t primitive_root) {
    const int n = (int)a.size();
    // bit-reverse
    for (int i=1,j=0;i<n;i++){
        int bit = n>>1;
        for(; j&bit; bit>>=1) j^=bit;
        j^=bit;
        if(i<j) std::swap(a[i], a[j]);
    }
    // layers
    for (int len=2; len<=n; len<<=1) {
        uint32_t wlen = powmod(primitive_root, (mod-1)/len, mod);
        if(invert) wlen = invmod(wlen, mod);
        for (int i=0;i<n;i+=len) {
            uint32_t w = 1;
            for (int j=0;j<len/2;j++){
                uint32_t u=a[i+j];
                uint32_t v=(uint32_t)((uint64_t)a[i+j+len/2]*w % mod);
                uint32_t x = u+v; if(x>=mod) x-=mod;
                uint32_t y = u>=v? u-v : u+mod-v;
                a[i+j] = x;
                a[i+j+len/2] = y;
                w = (uint32_t)((uint64_t)w*wlen % mod);
            }
        }
    }
    if (invert) {
        uint32_t ninv = invmod((uint32_t)n, mod);
        for (uint32_t& x: a) x = (uint32_t)((uint64_t)x * ninv % mod);
    }
}

// linear convolution via NTT
static std::vector<uint32_t> conv_ntt(std::vector<uint32_t> A, std::vector<uint32_t> B,
                                      uint32_t mod, uint32_t g) {
    int need = (int)A.size() + (int)B.size() - 1;
    int n=1; while(n<need) n<<=1;
    A.resize(n,0); B.resize(n,0);
    ntt(A,false,mod,g); ntt(B,false,mod,g);
    for (int i=0;i<n;i++) A[i] = (uint32_t)((uint64_t)A[i]*B[i] % mod);
    ntt(A,true,mod,g);
    A.resize(need);
    return A;
}

// ---------- bridge to Cigint ----------

// x % p  (p <= 2^31 fits u32)
static inline uint32_t mod_u32(const Cigint& x, uint32_t p){
    Cigint q; uint32_t r=0; cigint_sdivmod(x,p,&q,&r); return r;
}

// dst = (dst + add) mod modn   (all Cigint)
static inline void add_mod_big(Cigint& dst, const Cigint& add, const Cigint& modn) {
    // dst += add
    u64 carry = 0;
    for (size_t i = CIGINT_N; i-- > 0; ) {
        u64 t = (u64)dst.data[i] + (u64)add.data[i] + carry;
        dst.data[i] = (u32)t;
        carry = t >> 32;
    }
    // reduce (at most one subtract)
    if (carry || cigint_cmp(dst, modn) >= 0) cigint_sub_ref(&dst, &modn);
}

// out = (a * c) mod modn, where c is u32
static inline Cigint mul_u32_mod(const Cigint& a, uint32_t c, const Cigint& modn) {
    Cigint tmp = CIGINT_ZERO();
    u64 carry = 0;
    for (size_t i = CIGINT_N; i-- > 0; ) {
        u64 t = (u64)a.data[i] * c + carry;
        tmp.data[i] = (u32)t;
        carry = t >> 32;
    }
    cigint_divmod_ref(&tmp, &modn, nullptr, &tmp); // tmp %= modn
    return tmp;
}

// out = (a * c)  (no mod), used to update Mn before mod
static inline void mul_u32_ip(Cigint& a, uint32_t c) {
    u64 carry = 0;
    for (size_t i = CIGINT_N; i-- > 0; ) {
        u64 t = (u64)a.data[i] * c + carry;
        a.data[i] = (u32)t;
        carry = t >> 32;
    }
}

// ---------- CRT-to-mod-n via Garner, without big intermediates ----------
// Given residues r[i] (mod p[i]) for a single coefficient:
// return x ≡ CRT(r[i] mod p[i]) reduced modulo modn.
static Cigint crt_to_modn(const std::vector<uint32_t>& r,
                          const std::vector<uint32_t>& p,
                          const Cigint& modn)
{
    // Xn = current result modulo n
    // Mn = current modulus modulo n
    Cigint Xn = CIGINT_ZERO();
    Cigint Mn = CIGINT_ZERO();
    Mn.data[CIGINT_N - 1] = 1; // 1

    // Garner digits k_j; also need r0
    uint32_t r0 = r[0];
    // Initialize Xn with r0
    if (r0) {
        Cigint t = CIGINT_ZERO(); t.data[CIGINT_N - 1] = r0;
        cigint_divmod_ref(&t, &modn, nullptr, &t);
        add_mod_big(Xn, t, modn);
    }
    // Mn = p0 mod n
    if (p[0] != 1) {
        Cigint t = CIGINT_ZERO(); t.data[CIGINT_N - 1] = p[0];
        cigint_divmod_ref(&t, &modn, nullptr, &t);
        Mn = t;
    }

    std::vector<uint32_t> kj; kj.reserve(p.size());

    for (size_t i=1; i<p.size(); ++i) {
        uint32_t pi = p[i];

        // compute X (current) modulo pi using Garner eval:
        uint64_t xmod = r0 % pi;
        uint64_t prod = p[0] % pi;
        for (size_t j=0;j<kj.size();++j){
            xmod = (xmod + prod * kj[j]) % pi;
            prod = (prod * p[j+1]) % pi;
        }

        uint32_t Mi_mod_pi = (uint32_t)prod;
        uint32_t invMi = invmod(Mi_mod_pi, pi);
        uint32_t diff = (r[i] >= xmod) ? (uint32_t)(r[i] - xmod) : (uint32_t)(r[i] + pi - xmod);
        uint32_t ki = (uint32_t)((uint64_t)diff * invMi % pi);
        kj.push_back(ki);

        // Xn = Xn + Mn * ki  (mod n)
        if (ki) {
            Cigint add = mul_u32_mod(Mn, ki, modn);
            add_mod_big(Xn, add, modn);
        }
        // Mn = Mn * pi  (mod n)
        if (pi != 1) {
            // multiply Mn by pi, then reduce mod n
            MulAgain:
            {
                Cigint tmp = Mn; mul_u32_ip(tmp, pi);
                cigint_divmod_ref(&tmp, &modn, nullptr, &tmp);
                Mn = tmp;
            }
        }
    }
    return Xn; // in [0, n)
}

// Return Σ log2(p_i)
static double bits_of(const std::vector<uint32_t>& ps){
    double s=0; for (uint32_t q: ps) s += std::log2((double)q); return s;
}

// ---------- main: cyclic conv modulo (x^R-1, modn) via multi-prime NTT ----------
static Poly poly_mul_ntt_crt(const Poly& A, const Poly& B, const Cigint& modn) {
    // Candidate NTT primes with primitive root 3 and huge 2^k factor.
    // Extend this list if you really want guaranteed correctness for 512-bit n.
    struct P { uint32_t mod, g; };
    static const P pool[] = {
        {998244353u,    3u},
        {1004535809u,   3u},
        {469762049u,    3u},
        {167772161u,    3u},
        {1224736769u,   3u},
        {2013265921u,   31u},  // primitive root for this one is 31
        // add many more 2^k * c + 1 primes as needed ...
    };

    // How many primes do we need?
    uint32_t nbits = cigint_highest_order(modn);           // bitlen(n)
    double needBits = 2.0 * nbits + std::ceil(std::log2((double)R));
    std::vector<P> chosen;
    double acc=0;
    for (const auto &pp: pool) {
        chosen.push_back(pp);
        acc += std::log2((double)pp.mod);
        if (acc >= needBits + 2) break; // +2 safety
    }
    // If not enough headroom, fall back to your existing schoolbook path.
    if (acc < needBits) {
        // Fall back: exact and simple
        return poly_mul_fast(A, B, modn);
    }

    // Prepare linear-conv length N >= 2R-1
    int lin_len = 2*R - 1;
    int N = 1; while (N < lin_len) N <<= 1;  // e.g., for R=14 => N=32

    // For each prime: do linear conv and fold to cyclic (x^R-1)
    const int K = (int)chosen.size();
    std::vector<std::array<uint32_t, R>> residues(K);

    // Prepare A,B reduced modulo each prime
    for (int t=0;t<K;++t) {
        uint32_t mod = chosen[t].mod, g = chosen[t].g;

        std::vector<uint32_t> a(R), b(R);
        for (int i=0;i<R;++i){ a[i] = mod_u32(A.c[i], mod); b[i] = mod_u32(B.c[i], mod); }

        // linear conv (size >= 2R-1) via NTT
        auto aa = a; auto bb = b;
        auto lin = conv_ntt(aa, bb, mod, g); // size lin_len

        // fold to cyclic of length R: C[k] = lin[k] + lin[k+R]
        for (int k=0;k<R;++k){
            uint32_t v = (k < (int)lin.size()? lin[k] : 0);
            if (k+R < (int)lin.size()) v = addm(v, lin[k+R], mod);
            residues[t][k] = v;
        }
    }

    // CRT → modulo n, coefficient-wise
    Poly C; for (int i=0;i<R;++i) C.c[i] = CIGINT_ZERO();
    std::vector<uint32_t> mods; mods.reserve(K);
    for (int t=0;t<K;++t) mods.push_back(chosen[t].mod);

    for (int k=0;k<R;++k) {
        std::vector<uint32_t> rs; rs.reserve(K);
        for (int t=0;t<K;++t) rs.push_back(residues[t][k]);
        C.c[k] = crt_to_modn(rs, mods, modn);
    }
    return C;
}

} // namespace nttcrt

// ---- Swap-in hook inside your code ----
// Replace poly_mul_fast(...) calls with POLY_MUL(...)
// to enable NTT+CRT when safe (enough primes), else fall back.
static inline Poly POLY_MUL(const Poly& A, const Poly& B, const Cigint& m) {
    return nttcrt::poly_mul_ntt_crt(A, B, m); // auto-falls back if not enough bits
}


/* (1+x)^n mod (x^R-1, n) with 4-bit window on polynomials */
static Poly poly_pow_1x_win4(const Cigint &n, const Cigint &modn) {
	// base = 1 + x
	Poly base; for (int i = 0; i < R; ++i) base.c[i] = CIGINT_ZERO();
	base.c[0] = z1(); base.c[1] = z1();

	// precompute P[1..15] = base^d
	Poly P[16];
	P[0] = Poly{}; for (int i = 0; i < R; ++i) P[0].c[i] = CIGINT_ZERO();
	P[1] = base;
	for (int d = 2; d <= 15; ++d) P[d] = POLY_MUL(P[d - 1], base, modn);

	// res = 1
	Poly res; for (int i = 0; i < R; ++i) res.c[i] = CIGINT_ZERO();
	res.c[0] = z1();

	// scan exponent n by nibbles (MSW → LSW)
	u32 bits = cigint_highest_order(n);
	if (bits == 0) return res;
	u32 n_nibbles = (bits + 3u) / 4u;
	u32 skip = 8u * (u32)CIGINT_N - n_nibbles;

	bool first = true;
	u32 seen = 0;

	for (size_t i = 0; i < CIGINT_N; ++i) {
		u32 v = n.data[i];
		for (int k = 0; k < 8; ++k) {
			if (seen++ < skip) continue;
			u32 d = (v >> (28 - 4*k)) & 0xFu;

			if (first) {
				if (d) res = P[d];
				first = false;
			} else {
				// 4 squarings (res = res^(16))
				res = POLY_MUL(res, res, modn);
				res = POLY_MUL(res, res, modn);
				res = POLY_MUL(res, res, modn);
				res = POLY_MUL(res, res, modn);
				if (d) res = POLY_MUL(res, P[d], modn);
			}
		}
	}
	return first ? res : res;
}

static bool aks_like_cigint(const Cigint &n) {
	// reject 0,1
	if (cigint_is0(n)) return false;
	if (cigint_cmp(n, z1()) == 0) return false;

	Poly p = poly_pow_1x_win4(n, n);
	return poly_is_xn_plus_1(p, n, n);
}

#include <random>

void cigint_fill_more_random(Cigint *cig) {
	std::random_device rd;             // Random device
	std::mt19937 gen(rd());            // Mersenne Twister engine
	std::uniform_int_distribution<u32> dist(0, UINT32_MAX);  // Distribution for random u32
	size_t randomLimb = 1 + gen() % CIGINT_N;
	for (u32 &i : cig->data) i = 0;
	// Generate random values for cig->data
	for (size_t i = randomLimb; i-- > 0;) {
		cig->data[i] = dist(gen);  // Using the Mersenne Twister to get random values
	}
}

Cigint get_random_odd() {
	Cigint cig;
	cigint_fill_more_random(&cig);
	cigint_set_bit_ref(&cig, 0, 1);
	return cig;
}

static bool has_small_factor(const Cigint &n) {
	static const int SMALL_PRIMES[] = {
		2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97
		// extend to 1000 if you like
	};
	for (int p : SMALL_PRIMES) {
		u32 r = 0;
		cigint_sdivmod(n, (u32)p, NULL, &r);
		if (r == 0) return cigint_cmp(n, cigint_from_u32((u32)p)) != 0;
	}
	return false;
}

Cigint gen_prime() {
	Cigint cig;
	do {
		cig = get_random_odd();
	} while (has_small_factor(cig) || !aks_like_cigint(cig));
	return cig;
}

int main() {
	// Cigint n = cigint_from_dec("7669373798138436444425917457853162445620122597663674808454541043947311884377"); //up to 76digit
	// Cigint n = cigint_from_dec("771458402252315418489174122874340496560720991135040845033515625723932549727443357205751720461233037696609749135759573352834404574203828395431821751181407"); //up to 154 digit
	Cigint n = cigint_from_dec("12598680682438752944055149498662865476737841742262510750565369004945449194176119763492007280023309651826310840991041734855800903377926420446403888376677867"); //up to 154 digit

	// Cigint n = gen_prime();
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
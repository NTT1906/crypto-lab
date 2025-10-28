#define CIGINT_IMPLEMENTATION
#define CIGINT_STRIP_PREFIX
#define CIGINT_N (512 / 32)
// bench_shifts.cpp
#define CIGINT_IMPLEMENTATION
#include <chrono>
#include <cmath>
#include <random>
#include <cstring>
#include <cassert>
#include <iostream>
#include "cigint.h"

using u32 = uint32_t;
static inline uint64_t now_ns() {
  using namespace std::chrono;
  return duration_cast<nanoseconds>(steady_clock::now().time_since_epoch()).count();
}

// -------- Two-phase versions (shl2/shr2) --------

Cigint shl2(Cigint lhs, uint amnt) {
  Cigint res = (Cigint){0};
  size_t limb_off = amnt / SIZEOF_UINT;
  size_t bits = amnt % SIZEOF_UINT;
  if (limb_off >= CIGINT_N) return res;

  if (bits == 0) {
    for (size_t i = 0; i + limb_off < CIGINT_N; ++i)
      res.data[i] = lhs.data[i + limb_off];
    return res;
  }
  size_t rbits = SIZEOF_UINT - bits;
  size_t i = 0;
  for (; i + limb_off + 1 < CIGINT_N; ++i)
    res.data[i] = (lhs.data[i + limb_off] << bits) |
                  (lhs.data[i + limb_off + 1] >> rbits);
  if (i + limb_off < CIGINT_N)
    res.data[i] = (lhs.data[i + limb_off] << bits);
  return res;
}

Cigint shr2(Cigint lhs, uint amnt) {
  Cigint res = (Cigint){0};
  size_t limb_off = amnt / SIZEOF_UINT;
  size_t bits = amnt % SIZEOF_UINT;
  if (limb_off >= CIGINT_N) return res;

  if (bits == 0) {
    for (size_t i = 0; i + limb_off < CIGINT_N; ++i)
      res.data[CIGINT_N - 1 - i] = lhs.data[CIGINT_N - 1 - (i + limb_off)];
    return res;
  }
  size_t lbits = SIZEOF_UINT - bits;
  size_t i = 0;
  for (; i + limb_off + 1 < CIGINT_N; ++i) {
    uint hi = lhs.data[CIGINT_N - 1 - (i + limb_off)];
    uint lo = lhs.data[CIGINT_N - 1 - (i + limb_off + 1)];
    res.data[CIGINT_N - 1 - i] = (hi >> bits) | (lo << lbits);
  }
  if (i + limb_off < CIGINT_N) {
    uint hi = lhs.data[CIGINT_N - 1 - (i + limb_off)];
    res.data[CIGINT_N - 1 - i] = (hi >> bits);
  }
  return res;
}

// -------- Your old buggy versions (optional) --------
// // Keep only if you want to see them fail correctness on word-multiple shifts.
// static inline Cigint shl_old(Cigint lhs, uint amnt) {
//   Cigint res = lhs;
//   size_t i = 0;
//   size_t offset = (amnt + SIZEOF_UINT - 1) / SIZEOF_UINT;
//   size_t rshift_amnt = SIZEOF_UINT - (amnt % SIZEOF_UINT);
//   while (i < CIGINT_N - offset) {
//     res.data[i] = u1_shl(lhs.data[i], amnt);
//     res.data[i] |= u1_shr(lhs.data[i + offset], rshift_amnt);
//     i++;
//   }
//   while (i < CIGINT_N) {
//     res.data[i] = u1_shl(lhs.data[i], amnt);
//     i++;
//   }
//   return res;
// }
// static inline Cigint shr_old(Cigint lhs, uint amnt) {
//   Cigint res = lhs;
//   size_t i = 0;
//   size_t offset = (amnt + SIZEOF_UINT - 1) / SIZEOF_UINT;
//   size_t bits_amnt = amnt % SIZEOF_UINT;
//   size_t lshift_amnt = SIZEOF_UINT - bits_amnt;
//   while (i < CIGINT_N - offset) {
//     res.data[CIGINT_N - 1 - i] = u1_shr(lhs.data[CIGINT_N - 1 - i], amnt);
//     uint last_bits = u1_get_last_nbits(lhs.data[CIGINT_N - 1 - i - offset], bits_amnt);
//     res.data[CIGINT_N - 1 - i] |= u1_shl(last_bits, lshift_amnt);
//     i++;
//   }
//   while (i < CIGINT_N) {
//     res.data[CIGINT_N - 1 - i] = u1_shr(lhs.data[CIGINT_N - 1 - i], amnt);
//     i++;
//   }
//   return res;
// }

// ----- helpers -----
static inline bool eq(const Cigint& a, const Cigint& b) {
  for (size_t i=0;i<CIGINT_N;++i) if (a.data[i]!=b.data[i]) return false;
  return true;
}

static volatile Cigint sink;

int main(){
  std::mt19937_64 rng(12345);
  auto rand_cig = [&]{
    Cigint x = {0};
    for (size_t i=0;i<CIGINT_N;++i) x.data[i] = (u32)rng();
    // ensure non-zero LSW sometimes
    x.data[CIGINT_N-1] |= 1u;
    return x;
  };

  const int ITERS = 20000;
  uint64_t t_shl2=0, t_shl1=0, t_shl_old=0;
  uint64_t t_shr2=0, t_shr1=0, t_shr_old=0;

  // correctness + warmup
  for (int i=0;i<2000;++i){
    auto x = rand_cig();
    uint sh = (uint)(rng() % (CIGINT_N * SIZEOF_UINT));
    auto a = shl2(x,sh);
    // assert(eq(a,b)); // must match
  }

  // measure
  for (int i=0;i<ITERS;++i){
    auto x = rand_cig();
    uint sh = (uint)(rng() % (CIGINT_N * SIZEOF_UINT));

    auto t0=now_ns(); sink=shl2(x,sh); auto t1=now_ns(); t_shl2 += (t1-t0);
    t0=now_ns(); sink=shr2(x,sh); t1=now_ns(); t_shr2 += (t1-t0);

    // Optional: old buggy timings (don’t assert equality, they can be wrong)
    sh += 1;
    t0=now_ns(); sink=shl_old(x,sh); t1=now_ns(); t_shl_old += (t1-t0);
    t0=now_ns(); sink=shr_old(x,sh); t1=now_ns(); t_shr_old += (t1-t0);
  }

  std::cout << "avg ns/op (lower is better)\n";
  std::cout << "SHL two-phase: " << (t_shl2/ITERS) << "\n";
  std::cout << "SHL old (bug): " << (t_shl_old/ITERS) << " (not reliable)\n";
  std::cout << "SHR two-phase: " << (t_shr2/ITERS) << "\n";
  std::cout << "SHR old (bug): " << (t_shr_old/ITERS) << " (not reliable)\n";
}

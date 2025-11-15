#include <iostream>
#include <chrono>
#include <array>
#include <random>

#define BI_N 16
#include "bigint.h"

static void bui_randomize(bui &x) {
    std::random_device rd; std::mt19937 gen(rd());
    std::uniform_int_distribution<u32> dist(0, UINT32_MAX);
    size_t limbs = 1 + gen() % BI_N;
    for (u32 &i : x) i = 0;
    for (size_t i = limbs; i-- > 0;) x[i] = dist(gen);
}

template<typename Func>
long long benchmark(Func func, const bui& a, const bui& b, int iterations = 1000) {
    using namespace std::chrono;
    auto start = high_resolution_clock::now();

    for (int i = 0; i < iterations; ++i) {
        func(a, b);  // Call the function to be benchmarked
    }

    auto end = high_resolution_clock::now();
    return duration_cast<nanoseconds>(end - start).count() / iterations;
}

int main() {
    // Generate random test inputs
    bui a, b;
    bui_randomize(a);
    bui_randomize(b);

    // Benchmark `mul_low`
    auto time_low = benchmark(mul_low, a, b);
    std::cout << "mul_low average time: " << time_low << " ns\n";

    // Benchmark `mul_low_fast`
    auto time_fast = benchmark(mul_low_fast, a, b);
    std::cout << "mul_low_fast average time: " << time_fast << " ns\n";

    return 0;
}

#ifndef BENCHMARK_H
#define BENCHMARK_H
#include <chrono>
#include <atomic>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace std::chrono;

constexpr int ITER = 10000;
// constexpr int ITER = 1'000;

// your sink
template<class T>
inline void sink_volatile(const T& v) {
    static volatile T sink;
    sink = v;
}

// ------------------------------------------------------------
// generic bench
//  - name: label
//  - body: lambda that does EXACTLY ONE iteration
// Example:
// bench("add_fast_1", [&]{
//     Cigint r = cigint_add_fast_1(a, b);
//     consume(r);
// });
//
// bench("add_fast_1 (body)", [&]{
//     Cigint r = cigint_add_fast_1(a, b);
//     consume(r);
// });
//
// bench("add_ref (const rhs)", [&]{
//     Cigint x = a;
//     cigint_add_ref(&x, &b);
//     consume(x);
// });
//
// bench("add_ref_old (mut)", [&]{
//     Cigint x = a;
//     Cigint y = b;
//     cigint_add_ref_old(&x, &y);
//     consume(x);
// });
// ------------------------------------------------------------
template<typename Body>
void bench(const string& name, Body&& body, int iters = ITER) {
    atomic_signal_fence(memory_order_seq_cst);
    auto t0 = steady_clock::now();
    body();
    auto t1 = steady_clock::now();
    for (int i = 1; i < iters; ++i)
        body();
    auto t2 = steady_clock::now();
    atomic_signal_fence(memory_order_seq_cst);

    auto first = duration_cast<nanoseconds>(t1 - t0).count();
    auto total_loop = duration_cast<nanoseconds>(t2 - t0).count();
    double mean_loop = double(total_loop) / double(iters);

    long long total_body = 0;
    atomic_signal_fence(memory_order_seq_cst);
    for (int i = 0; i < iters; ++i) {
        auto t3 = steady_clock::now();
        body();
        auto t4 = steady_clock::now();
        total_body += duration_cast<nanoseconds>(t4 - t3).count();
    }
    atomic_signal_fence(memory_order_seq_cst);
    double mean_body = double(total_body) / double(iters);

    cout << left << setw(25) << name
         << " | first: " << right << setw(10) << first << " ns"
         << " | total_loop: " << setw(12) << total_loop << " ns"
         << " | mean_loop: "  << setw(12) << fixed << setprecision(2) << mean_loop << " ns"
         << " | total_body: " << right << setw(12) << total_body << " ns"
         << " | mean_body: "  << setw(12) << fixed << setprecision(2) << mean_body << " ns\n";
}
#endif // BENCHMARK_H

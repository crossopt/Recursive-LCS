#include <iostream>
#include <algorithm>
#include <chrono>
#include <string>

#include "lcs_kernel.h"

using Clock = std::chrono::steady_clock;
using std::chrono::time_point;
using std::chrono::duration_cast;
using std::chrono::milliseconds;

std::string generate_random_string(unsigned length) {
    auto randchar = []() -> char {
        const size_t max_index = 26;
        return 'a' + rand() % max_index;
    };
    std::string result(length, 0);
    std::generate_n(result.begin(), length, randchar);
    return result;
}

double time_dp(const std::string &a, const std::string &b) {
    time_point<Clock> start = Clock::now();
    LCS::kernel::dp_lcs(a, b);
    time_point<Clock> end = Clock::now();
    milliseconds diff = duration_cast<milliseconds>(end - start);
    return diff.count();
}

double time_recursive(const std::string &a, const std::string &b) {
    time_point<Clock> start = Clock::now();
    auto kernel = LCS::kernel::RecursiveLCS(a, b);
    kernel.lcs_whole_a(0, b.size());
    time_point<Clock> end = Clock::now();
    milliseconds diff = duration_cast<milliseconds>(end - start);
    return diff.count();
}

double time_recursive_base(const std::string &a, const std::string &b, unsigned recursion_base) {
    time_point<Clock> start = Clock::now();
    auto kernel = LCS::kernel::RecursiveLCS(a, b, recursion_base);
    kernel.lcs_whole_a(0, b.size());
    time_point<Clock> end = Clock::now();
    milliseconds diff = duration_cast<milliseconds>(end - start);
    return diff.count();
}

double time_iterative(const std::string &a, const std::string &b) {
    time_point<Clock> start = Clock::now();
    auto kernel = LCS::kernel::IterativeLCS(a, b);
    kernel.lcs_whole_a(0, b.size());
    time_point<Clock> end = Clock::now();
    milliseconds diff = duration_cast<milliseconds>(end - start);
    return diff.count();
}

void test_strings(const std::string &a, const std::string &b) {
    std::cout << "Times for LCS calculation for string lengths " << a.size() << " and " << b.size() << ":" << std::endl;
    std::cout << "Time for dynamic programming is " << time_dp(a, b) << "ms" << std::endl;
    std::cout << "Time for iterative kernel is " << time_iterative(a, b) << "ms" << std::endl;
    std::cout << "Time for recursive kernel is " << time_recursive(a, b) << "ms" << std::endl;
}

void test_recursive_bases(const std::string &a, const std::string &b) {
    std::cout << "Times for Recursive LCS calculation for string lengths " << a.size() << " and " << b.size() << ":" << std::endl;
    std::cout << "Time for base 4 is " << time_recursive_base(a, b, 4) << "ms" << std::endl;
    std::cout << "Time for base 8 is " << time_recursive_base(a, b, 8) << "ms" << std::endl;
    std::cout << "Time for base 16 is " << time_recursive_base(a, b, 16) << "ms" << std::endl;
    std::cout << "Time for base 32 is " << time_recursive_base(a, b, 32) << "ms" << std::endl;
}

int main() {
    test_recursive_bases(generate_random_string(100), generate_random_string(200));
    test_recursive_bases(generate_random_string(500), generate_random_string(500));

    test_strings(generate_random_string(100), generate_random_string(200));
    test_strings(generate_random_string(500), generate_random_string(500));
    test_strings(generate_random_string(1000), generate_random_string(1000));
    test_strings(generate_random_string(2000), generate_random_string(2000));
    return 0;
}

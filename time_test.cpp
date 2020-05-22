#include <iostream>
#include <algorithm>
#include <chrono>
#include <string>

#include "lcs_kernel.h"
#include "fibonacci.h"

using Clock = std::chrono::steady_clock;
using std::chrono::time_point;
using std::chrono::duration_cast;
using std::chrono::milliseconds;

std::string generate_random_string(unsigned length) {
    auto randchar = []() -> char {
        const size_t max_index = 26;
        return 'A' + rand() % max_index;
    };
    std::string result(length, 0);
    std::generate_n(result.begin(), length, randchar);
    return result;
}

std::string generate_fib_string(unsigned length) {
    if (length == 0) {
        return "A";
    }
    std::string prev = "A", current = "AB";
    for (unsigned n = 2; n <= length; ++n) {
        std::string next = current + prev;
        prev = current;
        current = next;
    }
    return current;
}

double time_dp(const std::string &a, const std::string &b) {
    time_point<Clock> start = Clock::now();
    std::cout << "dp returned " << LCS::kernel::dp_lcs(a, b) << std::endl;
    time_point<Clock> end = Clock::now();
    milliseconds diff = duration_cast<milliseconds>(end - start);
    return diff.count();
}

double time_recursive(const std::string &a, unsigned b_number) {
    time_point<Clock> start = Clock::now();
    std::cout << "recursive calculation returned " << LCS::fkernel::FibonacciKernel(a, b_number).lcs << std::endl;
    time_point<Clock> end = Clock::now();
    milliseconds diff = duration_cast<milliseconds>(end - start);
    return diff.count();
}

void test_strings(const std::string &a, unsigned b_number) {
    std::string b = generate_fib_string(b_number);
    std::cout << "Times for LCS calculation for string lengths " << a.size() << " and " << b.size() << ":" << std::endl;
    auto dp_time = time_dp(a, b);
    std::cout << "Time for dynamic programming is " << dp_time << "ms" << std::endl;
    auto recursive_time = time_recursive(a, b_number);
    std::cout << "Time for recursive kernel is " << recursive_time << "ms" << std::endl;
}


int main() {
    test_strings(generate_random_string(3), 10);
    test_strings(generate_random_string(5), 10);
    test_strings(generate_random_string(3), 30);
    test_strings(generate_random_string(5), 30);
    return 0;
}

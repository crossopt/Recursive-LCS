#include <iostream>
#include <algorithm>
#include <chrono>
#include <string>

#include "lcs_kernel.h"
#include "grammar_compressed.h"

using Clock = std::chrono::steady_clock;
using std::chrono::time_point;
using std::chrono::duration_cast;
using std::chrono::milliseconds;

std::string generate_random_string(unsigned length) {
    auto randchar = []() -> char {
        const size_t max_index = 5;
        return 'A' + rand() % max_index;
    };
    std::string result(length, 0);
    std::generate_n(result.begin(), length, randchar);
    return result;
}

LCS::gc::GrammarCompressed generate_fib_string(unsigned length) {
    if (length == 0) {
        return LCS::gc::GrammarCompressed(length + 1, 'A');
    } else if (length == 1) {
        return LCS::gc::GrammarCompressed(length + 1, generate_fib_string(0), LCS::gc::GrammarCompressed(0, 'B'));
    } else {
        LCS::gc::GrammarCompressed prev = generate_fib_string(length - 1);
        return LCS::gc::GrammarCompressed(length + 1, prev, *prev.first_symbol);
    }
}

double time_dp(const std::string &a, const std::string &b) {
    time_point<Clock> start = Clock::now();
    std::cout << "dp returned " << LCS::kernel::dp_lcs(a, b) << std::endl;
    time_point<Clock> end = Clock::now();
    milliseconds diff = duration_cast<milliseconds>(end - start);
    return diff.count();
}

double time_recursive(const std::string &a, const LCS::gc::GrammarCompressed &b_number) {
    time_point<Clock> start = Clock::now();
    std::cout << "recursive calculation returned " << LCS::gc::GCKernel(a, b_number).lcs << std::endl;
    time_point<Clock> end = Clock::now();
    milliseconds diff = duration_cast<milliseconds>(end - start);
    return diff.count();
}

void test_strings(const std::string &a, unsigned b_number) {
    LCS::gc::GrammarCompressed b = generate_fib_string(b_number);
    std::string b_string = b.decompress();
    std::cout << "Times for LCS calculation for string lengths " << a.size() << " and " << b_string.size() << ":" << std::endl;
    auto dp_time = time_dp(a, b_string);
    auto recursive_time = time_recursive(a, b);
    std::cout << "Time for dynamic programming is " << dp_time << "ms" << std::endl;
    std::cout << "Time for recursive kernel is " << recursive_time << "ms" << std::endl;
}


void test_recursive(const std::string &a, unsigned b_number) {
    LCS::gc::GrammarCompressed b = generate_fib_string(b_number);
    std::cout << "LCS calculation for string length " << a.size() << " and " << b_number << " Fibonacci string" << std::endl;
    auto recursive_time = time_recursive(a, b);
    std::cout << "Time for recursive kernel is " << recursive_time << "ms" << std::endl;
}


int main() {
    test_strings(generate_random_string(4), 10);
    test_strings(generate_random_string(16), 10);
    test_strings(generate_random_string(64), 10);
    test_strings(generate_random_string(256), 10);

    test_strings(generate_random_string(4), 20);
    test_strings(generate_random_string(16), 20);
    test_strings(generate_random_string(64), 20);
    test_strings(generate_random_string(256), 20);

    test_strings(generate_random_string(4), 30);
    test_strings(generate_random_string(16), 30);
    test_strings(generate_random_string(64), 30);
    test_strings(generate_random_string(256), 30);

    test_recursive(generate_random_string(4), 40);
    test_recursive(generate_random_string(16), 40);
    test_recursive(generate_random_string(64), 40);
    test_recursive(generate_random_string(256), 40);

    test_recursive(generate_random_string(4), 50);
    test_recursive(generate_random_string(16), 50);
    test_recursive(generate_random_string(64), 50);
    test_recursive(generate_random_string(256), 50);

    test_recursive(generate_random_string(4), 60);
    test_recursive(generate_random_string(16), 60);
    test_recursive(generate_random_string(64), 60);
    test_recursive(generate_random_string(256), 60);
    return 0;
}

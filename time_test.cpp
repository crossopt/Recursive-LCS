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

std::string generate_random_a_string(unsigned length) {
    auto randchar = []() -> char {
        return 'A';
    };
    std::string result(length, 0);
    std::generate_n(result.begin(), length, randchar);
    return result;
}


unsigned fib_string_storage_index(unsigned length, LCS::gc::GrammarCompressedStorage &gc_storage) {
    if (length == 0) {
        gc_storage.add_rule(LCS::gc::GrammarCompressed(gc_storage, length + 1, 'A'));
        return gc_storage.rules.size() - 1;
    } else if (length == 1) {
        unsigned index = fib_string_storage_index(0, gc_storage);
        gc_storage.add_rule(LCS::gc::GrammarCompressed(gc_storage, 0, 'B'));
        unsigned b_index = gc_storage.rules.size() - 1;
        gc_storage.add_rule(LCS::gc::GrammarCompressed(gc_storage, length + 1, index, b_index));
        return gc_storage.rules.size() - 1;
    } else {
        unsigned index = fib_string_storage_index(length - 1, gc_storage);
        gc_storage.add_rule(LCS::gc::GrammarCompressed(gc_storage, length + 1, index, (gc_storage.rules[index]).first_symbol));
        return gc_storage.rules.size() - 1;
    }
}

LCS::gc::GrammarCompressedStorage generate_fib_string(unsigned length) {
    LCS::gc::GrammarCompressedStorage storage = LCS::gc::GrammarCompressedStorage();
    unsigned index = fib_string_storage_index(length, storage);
    storage.final_rule = index;
    return storage;
}

double time_dp(const std::string &a, const std::string &b) {
    time_point<Clock> start = Clock::now();
    std::cout << "dp returned " << LCS::kernel::dp_lcs(a, b) << std::endl;
    time_point<Clock> end = Clock::now();
    milliseconds diff = duration_cast<milliseconds>(end - start);
    return diff.count();
}

double time_recursive(const std::string &a, const LCS::gc::GrammarCompressedStorage &b) {
    time_point<Clock> start = Clock::now();
    std::cout << "recursive calculation returned " << LCS::gc::GCKernel(a, b).lcs << std::endl;
    time_point<Clock> end = Clock::now();
    milliseconds diff = duration_cast<milliseconds>(end - start);
    return diff.count();
}

void test_strings(const std::string &a, unsigned b_number) {
    LCS::gc::GrammarCompressedStorage b = generate_fib_string(b_number);
    std::string b_string = b.rules[b.final_rule].decompress();
    std::cout << "Times for LCS calculation for string lengths " << a.size() << " and " << b_string.size() << ":" << std::endl;
    auto dp_time = time_dp(a, b_string);
    auto recursive_time = time_recursive(a, b);
    std::cout << "Time for dynamic programming is " << dp_time << "ms" << std::endl;
    std::cout << "Time for recursive kernel is " << recursive_time << "ms" << std::endl;
}


void test_recursive(const std::string &a, unsigned b_number) {
    LCS::gc::GrammarCompressedStorage b = generate_fib_string(b_number);
    std::cout << "LCS calculation for string length " << a.size() << " and " << b_number << " Fibonacci string" << std::endl;
    auto recursive_time = time_recursive(a, b);
    std::cout << "Time for recursive kernel is " << recursive_time << "ms" << std::endl;
}

void test_lz(const std::string &p, const std::string &t) {
    LCS::gc::GrammarCompressedStorage compress_78 = LCS::gc::LZ78(t);
    LCS::gc::GrammarCompressedStorage compress_w = LCS::gc::LZW(t);
    auto lz78_time = time_recursive(p, compress_78);
    auto lzw_time = time_recursive(p, compress_w);
    auto dp_time = time_dp(p, t);
    std::cout << "Time for dynamic programming is " << dp_time << "ms" << std::endl;
    std::cout << "Time for lz78 recursive kernel is " << lz78_time << "ms" << std::endl;
    std::cout << "Time for lzw recursive kernel is " << lzw_time << "ms" << std::endl;
}


int main() {
    test_lz(generate_random_string(4), generate_random_string(4096));
    test_lz(generate_random_string(16), generate_random_string(4096));
    test_lz(generate_random_string(64), generate_random_string(4096));
    test_lz(generate_random_string(256), generate_random_string(4096));

    test_lz(generate_random_string(4), generate_random_string(16384));
    test_lz(generate_random_string(16), generate_random_string(16384));
    test_lz(generate_random_string(64), generate_random_string(16384));
    test_lz(generate_random_string(256), generate_random_string(16384));

    test_lz(generate_random_string(4), generate_random_a_string(16384));
    test_lz(generate_random_string(16), generate_random_a_string(16384));
    test_lz(generate_random_string(64), generate_random_a_string(16384));
    test_lz(generate_random_string(256), generate_random_a_string(16384));

    test_lz(generate_random_string(4), generate_random_a_string(16384 * 16));
    test_lz(generate_random_string(16), generate_random_a_string(16384 * 16));
    test_lz(generate_random_string(64), generate_random_a_string(16384 * 16));
    test_lz(generate_random_string(256), generate_random_a_string(16384 * 16));

    test_lz(generate_random_string(4), generate_random_a_string(16384 * 16 * 16));
    test_lz(generate_random_string(16), generate_random_a_string(16384 * 16 * 16));
    test_lz(generate_random_string(64), generate_random_a_string(16384 * 16 * 16));
    test_lz(generate_random_string(256), generate_random_a_string(16384 * 16 * 16));
    return 0;
}

#include <iostream>
#include <algorithm>
#include <chrono>
#include <string>
#include <fstream>
#include <sstream>

#include "lcs_kernel.h"
#include "grammar_compressed.h"

using Clock = std::chrono::steady_clock;
using std::chrono::time_point;
using std::chrono::duration_cast;
using std::chrono::milliseconds;

std::string generate_random_abc_string(unsigned length) {
    auto randchar = []() -> char {
        const size_t max_index = 3;
        return 'A' + rand() % max_index;
    };
    std::string result(length, 0);
    std::generate_n(result.begin(), length, randchar);
    return result;
}

std::string generate_random_alpha_string(unsigned length) {
    auto randchar = []() -> char {
        const size_t max_index = 26;
        return 'A' + rand() % max_index;
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

double time_dp(const std::string &a, const std::string &b, bool dbg) {
    time_point<Clock> start = Clock::now();
    auto res = LCS::kernel::dp_lcs(a, b);
    if (dbg) {
        std::cout << "dp returned " << res << std::endl;
    }
    time_point<Clock> end = Clock::now();
    milliseconds diff = duration_cast<milliseconds>(end - start);
    return diff.count();
}

double time_recursive(const std::string &a, const LCS::gc::GrammarCompressedStorage &b, bool dbg) {
    time_point<Clock> start = Clock::now();
    auto res = LCS::gc::GCKernel(a, b).lcs;
    if (dbg) {
        std::cout << "recursive calculation returned " << res << std::endl;
    }
    time_point<Clock> end = Clock::now();
    milliseconds diff = duration_cast<milliseconds>(end - start);
    return diff.count();
}

void test_fibonacci(const std::string &a, unsigned b_number, bool dbg) {
    LCS::gc::GrammarCompressedStorage b = generate_fib_string(b_number);
    std::string b_string = b.rules[b.final_rule].decompress(b);
    if (dbg) {
        std::cout << "Times for LCS calculation for string lengths " << a.size() << " and " << b_string.size() << ":" << std::endl;
    }
    auto dp_time = time_dp(a, b_string, dbg);
    auto recursive_time = time_recursive(a, b, dbg);
    if (dbg) {
        std::cout << "Time for dynamic programming is " << dp_time << "ms" << std::endl;
        std::cout << "Time for recursive kernel is " << recursive_time << "ms" << std::endl;
    }
    // to-latex-format: pattern length, grammar length = fib number, uncompressed string length, compressed time, dp time
    if (!dbg) {
        std::cout << a.size() << '&' << b.final_rule + 1 << '&' << b_string.size() << '&' << recursive_time << '&' << dp_time << "\\\\" << std::endl;
    }
}

void test_aa(const std::string &a, unsigned long long b_number, bool dbg) {
    LCS::gc::GrammarCompressedStorage b = LCS::gc::get_aaaa(b_number);
    auto recursive_time = time_recursive(a, b, dbg);
    if (dbg) {
        std::cout << "Time for recursive kernel is " << recursive_time << "ms" << std::endl;
    }
    // to-latex-format: pattern length, grammar length = fib number, uncompressed string length, compressed time, dp time
    if (!dbg) {
        std::cout << a.size() << '&' << b.final_rule + 1 << '&' << b_number << '&' << recursive_time << "\\\\" << std::endl;
    }
}

void test_recursive(const std::string &a, unsigned b_number) {
    LCS::gc::GrammarCompressedStorage b = generate_fib_string(b_number);
    std::cout << "LCS calculation for string length " << a.size() << " and " << b_number << " Fibonacci string" << std::endl;
    auto recursive_time = time_recursive(a, b, 1);
    std::cout << "Time for recursive kernel is " << recursive_time << "ms" << std::endl;
}

void test_lz(const std::string &p, unsigned lz_number, bool dbg) {
    std::string t = LCS::gc::get_lz78_grammar_string(lz_number);
    LCS::gc::GrammarCompressedStorage compress_78 = LCS::gc::LZ78(t);
    auto lz78_time = time_recursive(p, compress_78, dbg);
    auto dp_time = time_dp(p, t, dbg);
    if (dbg) {
        std::cout << "Time for dynamic programming is " << dp_time << "ms" << std::endl;
        std::cout << "Time for lz78 recursive kernel is " << lz78_time << "ms" << std::endl;
    }
    // to-latex-format: pattern length, grammar length for lz78, uncompressed string length, lz78 time, dp time
    if (!dbg) {
        std::cout << p.size() << '&' << compress_78.final_rule + 3 << '&' <<
        t.size() << '&' << lz78_time << '&' << dp_time << "\\\\" << std::endl;
    }
}

void test_lzw(const std::string &p, unsigned lz_number, bool dbg) {
    std::string t = LCS::gc::get_lzw_grammar_string(lz_number);
    LCS::gc::GrammarCompressedStorage compress_w = LCS::gc::LZW(t);
    auto lzw_time = time_recursive(p, compress_w, dbg);
    auto dp_time = time_dp(p, t, dbg);
    if (dbg) {
        std::cout << "Time for dynamic programming is " << dp_time << "ms" << std::endl;
        std::cout << "Time for lzw recursive kernel is " << lzw_time << "ms" << std::endl;
    }
    // to-latex-format: pattern length, grammar length for lzw, uncompressed string length, lzw time, dp time
    if (!dbg) {
        std::cout << p.size() << '&' << compress_w.final_rule + 3 << '&' <<
        t.size() << '&' << lzw_time << '&' << dp_time << "\\\\" << std::endl;
    }
}

void test_unix_compress(const std::string &p, const std::string &file_name, bool dbg) {
    std::string t = LCS::gc::get_uncompress_string(file_name);
    LCS::gc::GrammarCompressedStorage compress_w = LCS::gc::get_compress_string(file_name);
    auto lzw_time = time_recursive(p, compress_w, dbg);
    auto dp_time = time_dp(p, t, dbg);
    // std::cout << "STRING " << t << std::endl;
    if (dbg) {
        std::cout << "Time for dynamic programming is " << dp_time << "ms" << std::endl;
        std::cout << "Time for lzw recursive kernel is " << lzw_time << "ms" << std::endl;
    }
    // to-latex-format: pattern length, compressed string length, uncompressed string length, lz time, dp time
    if (!dbg) {
        std::cout << p.size() << '&' << compress_w.final_rule + 3 << '&' <<
        t.size() << '&' << lzw_time << '&' << dp_time << "\\\\" << std::endl;
    }
}

void test_graph(unsigned int pattern_size, unsigned lz_number, unsigned int repeats, bool dbg) {
    double timedp = 0, timelzw = 0; int grammar_length = 0;
    for (unsigned int i = 0; i < repeats; ++i) {
        std::string p = generate_random_abc_string(pattern_size);
        std::string t = LCS::gc::get_lzw_grammar_string(lz_number, i);
        LCS::gc::GrammarCompressedStorage compress_w = LCS::gc::LZW(t);
        grammar_length = compress_w.final_rule + 3;
        auto lzw_time = time_recursive(p, compress_w, dbg);
        auto dp_time = time_dp(p, t, dbg);
        if (dbg) {
            std::cout << "Time for dynamic programming is " << dp_time << "ms" << std::endl;
            std::cout << "Time for lzw recursive kernel is " << lzw_time << "ms" << std::endl;
        }
        timedp += dp_time; timelzw += lzw_time;
    }
    // to-python-format: pattern length (const), grammar length for lzw, lzw time, dp time
    if (!dbg) {
        std::cout << pattern_size << '&' << grammar_length << '&' << 
        timelzw / repeats << '&' << timedp / repeats << std::endl;
    }
}

void test_graph_lz78(unsigned int pattern_size, unsigned lz_number, unsigned int repeats, bool dbg) {
    double timedp = 0, timelz78 = 0; int grammar_length = 0;
    for (unsigned int i = 0; i < repeats; ++i) {
        std::string p = generate_random_abc_string(pattern_size);
        std::string t = LCS::gc::get_lz78_grammar_string(lz_number, i);
        LCS::gc::GrammarCompressedStorage compress_w = LCS::gc::LZ78(t);
        grammar_length = compress_w.final_rule + 3;
        auto lz78_time = time_recursive(p, compress_w, dbg);
        auto dp_time = time_dp(p, t, dbg);
        if (dbg) {
            std::cout << "Time for dynamic programming is " << dp_time << "ms" << std::endl;
            std::cout << "Time for lz78 recursive kernel is " << lz78_time << "ms" << std::endl;
        }
        timedp += dp_time; timelz78 += lz78_time;
    }
    // to-python-format: pattern length (const), grammar length for lz78, lz78 time, dp time
    if (!dbg) {
        std::cout << pattern_size << '&' << grammar_length << '&' << 
        timelz78 / repeats << '&' << timedp / repeats << std::endl;
    }
}

void test_bad_graph(unsigned int pattern_size, const std::string &tlist, bool dbg) {
    double timedp = 0, timelzw = 0;
    int grammar_length = 0;
    unsigned int repeats = 20;
    unsigned int tlen = tlist.size() / 20;
    for (unsigned int i = 0; i < repeats; ++i) {
        std::string t = tlist.substr(i * tlen, tlen);
        std::string p = generate_random_abc_string(pattern_size);
        LCS::gc::GrammarCompressedStorage compress_w = LCS::gc::LZW(t); // TODO LZW2
        grammar_length = compress_w.final_rule + 3;
        auto lzw_time = time_recursive(p, compress_w, dbg);
        auto dp_time = time_dp(p, t, dbg);
        if (dbg) {
            std::cout << "Time for dynamic programming is " << dp_time << "ms" << std::endl;
            std::cout << "Time for lzw recursive kernel is " << lzw_time << "ms" << std::endl;
        }        
        timedp += dp_time; timelzw += lzw_time;
    }
    // to-python-format: pattern length (const), grammar length for lzw, lzw time, dp time
    if (!dbg) {
        std::cout << pattern_size << '&' << grammar_length << '&' << 
        timelzw / repeats << '&' << timedp / repeats<< std::endl;
    }
}

void test_bad_graph_lz78(unsigned int pattern_size, const std::string &tlist, bool dbg) {
    double timedp = 0, timelz78 = 0;
    int grammar_length = 0;
    unsigned int repeats = 20;
    unsigned int tlen = tlist.size() / 20;
    for (unsigned int i = 0; i < repeats; ++i) {
        std::string t = tlist.substr(i * tlen, tlen);
        std::string p = generate_random_abc_string(pattern_size);
        LCS::gc::GrammarCompressedStorage compress_w = LCS::gc::LZ78(t);
        grammar_length = compress_w.final_rule + 3;
        auto lz78_time = time_recursive(p, compress_w, dbg);
        auto dp_time = time_dp(p, t, dbg);
        if (dbg) {
            std::cout << "Time for dynamic programming is " << dp_time << "ms" << std::endl;
            std::cout << "Time for lz78 recursive kernel is " << lz78_time << "ms" << std::endl;
        }        
        timedp += dp_time; timelz78 += lz78_time;
    }
    // to-python-format: pattern length (const), grammar length for lz78, lz78 time, dp time
    if (!dbg) {
        std::cout << pattern_size << '&' << grammar_length << '&' << 
        timelz78 / repeats << '&' << timedp / repeats << std::endl;
    }
}



int main() {
    // All time tests that have been run for this code.
    // This is not meant to be run simultaneously!

    // test_fibonacci(generate_random_abc_string(4), 14, 0);
    // test_fibonacci(generate_random_abc_string(16), 14, 0);
    // test_fibonacci(generate_random_abc_string(64), 14, 0);
    // test_fibonacci(generate_random_abc_string(256), 14, 0);

    // test_fibonacci(generate_random_abc_string(4), 14 + 8, 0);
    // test_fibonacci(generate_random_abc_string(16), 14 + 8, 0);
    // test_fibonacci(generate_random_abc_string(64), 14 + 8, 0);
    // test_fibonacci(generate_random_abc_string(256), 14 + 8, 0);

    // test_fibonacci(generate_random_abc_string(4), 30, 0);
    // test_fibonacci(generate_random_abc_string(16), 30, 0);
    // test_fibonacci(generate_random_abc_string(64), 30, 0);
    // test_fibonacci(generate_random_abc_string(256), 30, 0);
    
    // test_lz(generate_random_alpha_string(4), 64, 0);
    // test_lz(generate_random_alpha_string(16), 64, 0);
    // test_lz(generate_random_alpha_string(64), 64, 0);
    // test_lz(generate_random_alpha_string(256), 64, 0);

    // test_lz(generate_random_alpha_string(4), 64 * 8, 0);
    // test_lz(generate_random_alpha_string(16), 64 * 8, 0);
    // test_lz(generate_random_alpha_string(64), 64 * 8, 0);
    // test_lz(generate_random_alpha_string(256), 64 * 8, 0);

    // test_lz(generate_random_alpha_string(4), 64 * 64, 0);
    // test_lz(generate_random_alpha_string(16), 64 * 64, 0);
    // test_lz(generate_random_alpha_string(64), 64 * 64, 0);
    // test_lz(generate_random_alpha_string(256), 64 * 64, 0);
    
    // test_lzw(generate_random_alpha_string(4), 64, 0);
    // test_lzw(generate_random_alpha_string(16), 64, 0);
    // test_lzw(generate_random_alpha_string(64), 64, 0);
    // test_lzw(generate_random_alpha_string(256), 64, 0);

    // // test_lzw(generate_random_alpha_string(4), 64 * 8, 0);
    // test_lzw(generate_random_alpha_string(16), 64 * 8, 0);
    // test_lzw(generate_random_alpha_string(64), 64 * 8, 0);
    // test_lzw(generate_random_alpha_string(256), 64 * 8, 0);

    // // test_lzw(generate_random_alpha_string(4), 64 * 64, 0);
    // test_lzw(generate_random_alpha_string(16), 64 * 64, 0);
    // test_lzw(generate_random_alpha_string(64), 64 * 64, 0);
    // test_lzw(generate_random_alpha_string(256), 64 * 64, 0);

    // test_unix_compress(generate_random_alpha_string(16), "../test_files/t10.Z", 1);
    // test_unix_compress(generate_random_alpha_string(64), "../test_files/t10.Z", 1);
    // test_unix_compress(generate_random_alpha_string(256), "../test_files/t10.Z", 1);

    // test_unix_compress(generate_random_alpha_string(16), "../test_files/t9.Z", 1);
    // test_unix_compress(generate_random_alpha_string(64), "../test_files/t9.Z", 1);
    // test_unix_compress(generate_random_alpha_string(256), "../test_files/t9.Z", 1);

    // test_unix_compress(generate_random_alpha_string(16), "../test_files/t8.Z", 1);
    // test_unix_compress(generate_random_alpha_string(64), "../test_files/t8.Z", 1);
    // test_unix_compress(generate_random_alpha_string(256), "../test_files/t8.Z", 1);

    // test_aa(generate_random_abc_string(30), 8388608, 0);
    // test_aa(generate_random_abc_string(30), 1ll * 536870000, 0);
    // test_aa(generate_random_abc_string(30), 1ll * 536800000 * 10, 0);

    srand(time(0));

    // LZW & LZ78 generated runs
    // for (unsigned int i = 1; i < 100000; i += 100) {
    //     test_graph(16, i, 20, 0);
    // }

    // std::ifstream t("../test_files/rand_large");
    // std::stringstream buffer;
    // buffer << t.rdbuf();
    // std::string all = buffer.str();

    // // std::cerr << all.size() << '\n';
    // for (unsigned int i = 100000; 20 * i < all.size(); i += 300000) {
    //     std::string ith = all.substr(0, 20 * i);
    //     test_bad_graph_lz78(16, ith, 0);
    // }

    return 0;
}

#include <string>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "gtest/gtest.h"
#include "monge_matrix.h"
#include "lcs_kernel.h"
#include "grammar_compressed.h"

namespace LCS {
namespace gc {
namespace {

std::string fib_string(int n) {
	if (n == 0) {
		return "A";
	} else if (n == 1) {
		return "AB";
	} else {
		return fib_string(n - 1) + fib_string(n - 2);
	}
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

/*
LCS::gc::GrammarCompressedStorage gc_fib_string(unsigned length) {
    LCS::gc::GrammarCompressedStorage storage = LCS::gc::GrammarCompressedStorage();
    unsigned index = fib_string_storage_index(length, storage);
    storage.final_rule = index;
    return storage;
}

TEST(GrammarCompressedTest, FibonacciLengthOneTest) {
    std::string s = "APATTERNA";
    auto gc_string = gc_fib_string(0);
    auto fkernel = GCKernel(s, gc_string);
    ASSERT_EQ(gc_string.rules[gc_string.final_rule].decompress(gc_string), fib_string(0));
    ASSERT_EQ(kernel::dp_lcs(s, fib_string(0)), fkernel.lcs);
}

TEST(GrammarCompressedTest, FibonacciLengthTwoTest) {
    std::string s = "APATBTERNAB";
    auto gc_string = gc_fib_string(1);
    auto fkernel = GCKernel(s, gc_string);
    ASSERT_EQ(gc_string.rules[gc_string.final_rule].decompress(gc_string), fib_string(1));
    ASSERT_EQ(kernel::dp_lcs(s, fib_string(1)), fkernel.lcs);
}

TEST(GrammarCompressedTest, FibonacciLengthThreeTest) {
    std::string s = "APATBTERNAB";
    auto gc_string = gc_fib_string(2);
    auto fkernel = GCKernel(s, gc_string);
    ASSERT_EQ(gc_string.rules[gc_string.final_rule].decompress(gc_string), fib_string(2));
    ASSERT_EQ(kernel::dp_lcs(s, fib_string(2)), fkernel.lcs);
}



TEST(GrammarCompressedTest, FibonacciSameFibStringTest) {
    std::string s = "ABAABABAABAAB";
    auto gc_string = gc_fib_string(5);
    auto fkernel = GCKernel(s, gc_string);
    ASSERT_EQ(gc_string.rules[gc_string.final_rule].decompress(gc_string), fib_string(5));
    ASSERT_EQ(kernel::dp_lcs(s, fib_string(5)), fkernel.lcs);
}


TEST(GrammarCompressedTest, FibonacciDifferentFibStringTest) {
    std::string s = "ABACABABDAABAAAB";
    auto gc_string = gc_fib_string(8);
    auto fkernel = GCKernel(s, gc_string);
    ASSERT_EQ(gc_string.rules[gc_string.final_rule].decompress(gc_string), fib_string(8));
    ASSERT_EQ(kernel::dp_lcs(s, fib_string(8)), fkernel.lcs);
}

TEST(GrammarCompressedTest, FibonacciManyABsTest) {
    std::string s = "ABABABABABA";
    auto gc_string = gc_fib_string(5);
    auto fkernel = GCKernel(s, gc_string);
    ASSERT_EQ(gc_string.rules[gc_string.final_rule].decompress(gc_string), fib_string(5));
    ASSERT_EQ(kernel::dp_lcs(s, fib_string(5)), fkernel.lcs);
}

TEST(GrammarCompressedTest, FibonacciNoMatchesTest) {
    std::string s = "QWERTY";
    auto gc_string = gc_fib_string(2);
    auto fkernel = GCKernel(s, gc_string);
    ASSERT_EQ(gc_string.rules[gc_string.final_rule].decompress(gc_string), fib_string(2));
    ASSERT_EQ(kernel::dp_lcs(s, fib_string(2)), fkernel.lcs);
}

TEST(GrammarCompressedTest, FibonacciOneMatchATest) {
    std::string s = "QWAERTY";
    auto gc_string = gc_fib_string(7);
    auto fkernel = GCKernel(s, gc_string);
    ASSERT_EQ(gc_string.rules[gc_string.final_rule].decompress(gc_string), fib_string(7));
    ASSERT_EQ(kernel::dp_lcs(s, fib_string(7)), fkernel.lcs);
}


TEST(GrammarCompressedTest, FibonacciOneMatchBTest) {
    std::string s = "QWERTBY";
    auto gc_string = gc_fib_string(7);
    auto fkernel = GCKernel(s, gc_string);
    ASSERT_EQ(gc_string.rules[gc_string.final_rule].decompress(gc_string), fib_string(7));
    ASSERT_EQ(kernel::dp_lcs(s, fib_string(7)), fkernel.lcs);
}

TEST(GrammarCompressedTest, LZ78ComputesCorrectlyTest) {
    std::string s = "ABACABA";
    auto gc_lz78_string = LZ78(s);
    // Test the general structure of the compressed string.
    ASSERT_EQ(gc_lz78_string.rules[gc_lz78_string.final_rule].decompress(gc_lz78_string), s);
    auto left = gc_lz78_string.rules[gc_lz78_string.final_rule].first_symbol;
    auto right = gc_lz78_string.rules[gc_lz78_string.final_rule].second_symbol;
    ASSERT_EQ(gc_lz78_string.rules[left].decompress(gc_lz78_string), "ABACAB");
    ASSERT_EQ(gc_lz78_string.rules[right].decompress(gc_lz78_string), "A");
    auto lower_left = gc_lz78_string.rules[left].first_symbol;
    auto lower_right = gc_lz78_string.rules[left].second_symbol;
    ASSERT_EQ(gc_lz78_string.rules[lower_left].decompress(gc_lz78_string), "ABAC");
    ASSERT_EQ(gc_lz78_string.rules[lower_right].decompress(gc_lz78_string), "AB");
    auto final_left = gc_lz78_string.rules[lower_left].first_symbol;
    auto final_right = gc_lz78_string.rules[lower_left].second_symbol;
    ASSERT_EQ(gc_lz78_string.rules[final_left].decompress(gc_lz78_string), "AB");
    ASSERT_EQ(gc_lz78_string.rules[final_right].decompress(gc_lz78_string), "AC");
    // Test that the decompression works correctly.
    ASSERT_EQ(gc_lz78_string.rules[gc_lz78_string.final_rule].number, 11);
    ASSERT_EQ(gc_lz78_string.rules[left].number, 9);
    ASSERT_EQ(gc_lz78_string.rules[right].number, 10);
    ASSERT_EQ(gc_lz78_string.rules[lower_left].number, 6);
    ASSERT_EQ(gc_lz78_string.rules[lower_right].number, 8);
    ASSERT_EQ(gc_lz78_string.rules[final_left].number, 3);
    ASSERT_EQ(gc_lz78_string.rules[final_right].number, 5);
} 

TEST(GrammarCompressedTest, LZWComputesCorrectlyTest) {
    std::string s = "ABACABA";
    auto gc_lzw_string = LZW(s);
    // Test the general structure of the compressed string.
    ASSERT_EQ(gc_lzw_string.rules[gc_lzw_string.final_rule].decompress(gc_lzw_string), s);
    auto left = gc_lzw_string.rules[gc_lzw_string.final_rule].first_symbol;
    auto right = gc_lzw_string.rules[gc_lzw_string.final_rule].second_symbol;
    ASSERT_EQ(gc_lzw_string.rules[left].decompress(gc_lzw_string), "ABAC");
    ASSERT_EQ(gc_lzw_string.rules[right].decompress(gc_lzw_string), "ABA");
    auto lower_left = gc_lzw_string.rules[left].first_symbol;
    auto lower_right = gc_lzw_string.rules[left].second_symbol;
    ASSERT_EQ(gc_lzw_string.rules[lower_left].decompress(gc_lzw_string), "AB");
    ASSERT_EQ(gc_lzw_string.rules[lower_right].decompress(gc_lzw_string), "AC");
    // Test that the decompression works correctly.
    ASSERT_EQ(gc_lzw_string.rules[gc_lzw_string.final_rule].number, 27 + 5);
    ASSERT_EQ(gc_lzw_string.rules[left].number, 27 + 3);
    ASSERT_EQ(gc_lzw_string.rules[right].number, 27 + 4);
    ASSERT_EQ(gc_lzw_string.rules[lower_left].number, 27 + 1);
    ASSERT_EQ(gc_lzw_string.rules[lower_right].number, 27 + 2);
}

TEST(GrammarCompressedTest, LZ78GrammarStringIsCalculatedCorrectlyTest) {
    ASSERT_EQ(get_lz78_grammar_string(0), "A");
    ASSERT_EQ(get_lz78_grammar_string(1), "AAB");
    ASSERT_EQ(get_lz78_grammar_string(2), "AABABC");
    ASSERT_EQ(get_lz78_grammar_string(3), "AABABCABCD");
    ASSERT_EQ(get_lz78_grammar_string(4), "AABABCABCDABCDE");
    ASSERT_EQ(get_lz78_grammar_string(5), "AABABCABCDABCDEABCDEF");
    ASSERT_EQ(get_lz78_grammar_string(6), "AABABCABCDABCDEABCDEFABCDEFG");
}

TEST(GrammarCompressedTest, LZWGrammarStringIsCalculatedCorrectlyTest) {
    ASSERT_EQ(get_lzw_grammar_string(0), "AA");
    ASSERT_EQ(get_lzw_grammar_string(1), "AAAAB");
    ASSERT_EQ(get_lzw_grammar_string(2), "AAAABAABC");
    ASSERT_EQ(get_lzw_grammar_string(3), "AAAABAABCAABCD");
    ASSERT_EQ(get_lzw_grammar_string(4), "AAAABAABCAABCDAABCDE");
    ASSERT_EQ(get_lzw_grammar_string(5), "AAAABAABCAABCDAABCDEAABCDEF");
    ASSERT_EQ(get_lzw_grammar_string(6), "AAAABAABCAABCDAABCDEAABCDEFAABCDEFG");
}

TEST(GrammarCompressedTest, LZGrammarStringIsCalculatedCorrectlyTest) {
    ASSERT_EQ(get_lz_grammar_string(0), "ABAA");
    ASSERT_EQ(get_lz_grammar_string(1), "ABAAAAB");
    ASSERT_EQ(get_lz_grammar_string(2), "ABAAAABAABC");
    ASSERT_EQ(get_lz_grammar_string(3), "ABAAAABAABCAABCD");
}

TEST(GrammarCompressedTest, LZ78GrammarStringIsCompressedQuadraticallyByLZ78Compressions) {
    std::vector <unsigned int> number_list = {5, 179, 300};
    for (auto number: number_list) {
        std::string s = get_lz78_grammar_string(number);
        ASSERT_EQ(s.size(), (number + 1) * (number + 2) / 2);
        auto gc_lz78_string = LZ78(s);
        // Exactly n * 3: one rule per symbol, one per substring, one for concatenation.
        ASSERT_EQ(gc_lz78_string.final_rule, number * 3);   
    }
}

TEST(GrammarCompressedTest, LZWGrammarStringIsCompressedQuadraticallyByLZWCompressions) {
    std::vector <unsigned int> number_list = {5, 179, 300};
    for (auto number: number_list) {
        std::string s = get_lzw_grammar_string(number);
        ASSERT_EQ(s.size(), (number + 1) * (number + 2) / 2 + number + 1);
        auto gc_lzw_string = LZW(s);
        // Exactly n * 2 + 27: the alphabet rules, one rule per substring, one for concatenation.
        ASSERT_EQ(gc_lzw_string.final_rule, number * 2 + 27);   
    }
}

TEST(GrammarCompressedTest, LZGrammarStringIsCompressedQuadraticallyByLZ78Compressions) {
    std::vector <unsigned int> number_list = {5, 179, 300};
    for (auto number: number_list) {
        std::string s = get_lz_grammar_string(number);
        ASSERT_EQ(s.size(), (number + 1) * (number + 2) / 2  + (number + 1) + 2);
        auto gc_lz78_string = LZ78(s);
        // 5 original rules for 'ABAA': 'A', 'B', 'AB', 'A', AA'.
        // Exactly n * 3: one rule per symbol, one per substring, one for concatenation.
        ASSERT_EQ(gc_lz78_string.final_rule, 5 + number * 3);   
    }
}

TEST(GrammarCompressedTest, LZGrammarStringIsCompressedQuadraticallyByLZWCompressions) {
    std::vector <unsigned int> number_list = {5, 179, 300};
    for (auto number: number_list) {
        std::string s = get_lz_grammar_string(number);
        ASSERT_EQ(s.size(), (number + 1) * (number + 2) / 2 + (number + 1) + 2);
        auto gc_lzw_string = LZW(s);
        // 2 original rules for 'ABAA': 'AB', 'AA'.
        // Exactly n * 2 + 27: the alphabet rules, one rule per substring, one for concatenation.
        ASSERT_EQ(gc_lzw_string.final_rule, 2 + number * 2 + 27);   
    }
}

TEST(GrammarCompressedTest, StringDecompressReturnsCorrectStringTest) {
    ASSERT_EQ(get_uncompress_string("../test_files/f1.Z"), "aaaaaaaa\n");
    ASSERT_EQ(get_uncompress_string("../test_files/f2.Z"), "This is a test file!\n");
    ASSERT_EQ(get_uncompress_string("../test_files/f3.Z"), "aaaabaabcaabcd");
    ASSERT_EQ(get_uncompress_string("../test_files/f4.Z"), "aaaaaaaaaaaaaaaaaaaa");
}

TEST(GrammarCompressedTest, StringAaaaReturnsCorrectGrammarTest) {
    std::string aa = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
    GrammarCompressedStorage gcs = get_aaaa(40);
    ASSERT_EQ(gcs.rules[gcs.final_rule].decompress(gcs), aa);

    std::string aa2 = "aaaa";
    GrammarCompressedStorage gcs2 = get_aaaa(4);
    ASSERT_EQ(gcs2.rules[gcs2.final_rule].decompress(gcs2), aa2);


    GrammarCompressedStorage gcs3 = get_aaaa(536800000);
    std::cout << gcs3.final_rule << '\n';
} */

TEST(GrammarCompressedTest, GrammarDecompressReturnsCorrectStringTest) {
    // auto gcs1 = get_compress_string("../test_files/f1.Z");
    // ASSERT_EQ(gcs1.rules[gcs1.final_rule].decompress(gcs1), "aaaaaaaa\n");

    // auto gcs2 = get_compress_string("../test_files/f2.Z");
    // ASSERT_EQ(gcs2.rules[gcs2.final_rule].decompress(gcs2), "This is a test file!\n");

    auto str3 = get_uncompress_string("../test_files/f3.Z");
    std::cout << "break\n\n\n\n\n\n";
    auto gcs3 = get_compress_string("../test_files/f3.Z");
    ASSERT_EQ(gcs3.rules[gcs3.final_rule].decompress(gcs3), "aaaabaabcaabcd");
    std::cout << "break\n\n\n\n\n\n";

    // auto gcs4 = get_compress_string("../test_files/f4.Z");
    // auto str = get_uncompress_string("../test_files/f4.Z");
    // ASSERT_EQ(gcs4.rules[gcs4.final_rule].decompress(gcs4), "aaaaaaaaaaaaaaaaaaaa"); // 20
} 

}  // namespace
}  // namespace gc
}  // namespace LCS

#include <iostream>
#include <algorithm>
#include <string>

#include "lcs_kernel.h"
#include "grammar_compressed.h"

int main() {
    std::string pattern_string;
    unsigned long long n;
    std::cin >> n >> pattern_string;
    auto compress_w = LCS::gc::get_compress_string("../file.Z");
    auto res = LCS::gc::GCKernel(pattern_string, compress_w).lcs;
    std::cerr << res << '\n';
    return 0;
}

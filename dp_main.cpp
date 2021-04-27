#include <iostream>
#include <algorithm>
#include <string>

#include "lcs_kernel.h"
#include "grammar_compressed.h"

int main() {
    std::string pattern_string, file_name;
    cin >> pattern_string >> file_name;
    std::string t = LCS::gc::get_uncompress_string("../file.Z");
    auto res = LCS::kernel::dp_lcs(pattern_string, t);
    std::cerr << res << '\n';
    return 0;
}

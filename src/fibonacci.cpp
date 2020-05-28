#include "fibonacci.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>

namespace LCS {
namespace gc {

GCKernel::GCKernel(const std::string &p, const GrammarCompressed &t): lcs(calculate_lcs(p, t)) {}


matrix::Permutation GCKernel::calculate_char_kernel(const std::string &p, char c) {
    unsigned last_row = p.size();
    std::vector <unsigned> last_col(p.size());
    for (unsigned i = 0; i < p.size(); ++i) {
        last_col[i] = p.size() - i - 1;
        if (p[i] == c || last_col[i] > last_row) {
            std::swap(last_col[i], last_row);
        }
    }
    std::vector <unsigned> result(p.size() + 1);
    last_col.push_back(last_row);
    for (unsigned i = 0; i <= p.size(); ++i) {
        result[last_col[i]] = p.size() + 1 - i;
    }
    if (last_row == p.size()) {  // from top to bottom
        result.pop_back();
    }
    return matrix::Permutation(result);
}


// return permutation split into stuff touching left and not
std::pair <matrix::Permutation, matrix::Permutation> get_left(const matrix::Permutation &p, unsigned m) {
    // from left => row index is from 1 to m
    return p.split_row(m);
}

// return permutation split into stuff touching right and not
std::pair <matrix::Permutation, matrix::Permutation> get_right(const matrix::Permutation &p, unsigned min_m) {
    // to right => col index is from n + 1 to n + m
    return p.split_col(min_m);
}

matrix::Permutation combine(matrix::Permutation &left_side,
                            matrix::Permutation &both_sides,
                            matrix::Permutation &right_side,
                            unsigned add) {
    // left side needs no fixing
    for (auto &i : both_sides.rows) {
        i.second += add;
    }
    for (auto &i : both_sides.cols) {
        i.first += add;
    }
    for (auto &i : right_side.rows) {
        i.first += add;
        i.second += add;
    }
    for (auto &i : right_side.cols) {
        i.first += add;
        i.second += add;
    }
    std::vector <std::pair <unsigned, unsigned> > part_combined_rows, all_combined_rows;
    std::vector <std::pair <unsigned, unsigned> > part_combined_cols, all_combined_cols;
    std::merge(left_side.cols.begin(), left_side.cols.end(),
               both_sides.cols.begin(), both_sides.cols.end(),
               std::back_inserter(part_combined_cols));
    std::merge(right_side.cols.begin(), right_side.cols.end(),
               part_combined_cols.begin(), part_combined_cols.end(),
               std::back_inserter(all_combined_cols));
    std::merge(left_side.rows.begin(), left_side.rows.end(),
               both_sides.rows.begin(), both_sides.rows.end(),
               std::back_inserter(part_combined_rows), std::greater<std::pair <unsigned, unsigned>>());
    std::merge(right_side.rows.begin(), right_side.rows.end(),
               part_combined_rows.begin(), part_combined_rows.end(),
               std::back_inserter(all_combined_rows), std::greater<std::pair <unsigned, unsigned>>());
    return matrix::Permutation{all_combined_rows, all_combined_cols};
}


matrix::Permutation GCKernel::calculate_count_kernel(const std::string &p, const GrammarCompressed &t, int &fm) {
    if (t.is_base) { // t is a single symbol
        fm = 1;
        return calculate_char_kernel(p, t.value);
    } else { // t = uv
        int x1 = 0, x2 = 0;
        auto first_half = calculate_count_kernel(p, *t.first_symbol, x1);
        auto second_half = calculate_count_kernel(p, *t.second_symbol, x2);
        fm = x1 + x2;

        auto to_right = get_right(first_half, x1);
        auto from_left = get_left(second_half, p.size());
        auto intersection = to_right.second * from_left.first;
        return combine(to_right.first, intersection, from_left.second, x1);
    }
}

unsigned GCKernel::calculate_lcs(const std::string &p, const GrammarCompressed &t) {
    int size = 0;
    auto kernel = calculate_count_kernel(p, t, size);
    unsigned count_dom = 0;
    for (auto i: kernel.rows) {
        count_dom += (i.first <= p.size() && i.second > (unsigned)size);
    }
    return p.size() - count_dom;
}

}  // namespace fkernel
}  // namespace LCS

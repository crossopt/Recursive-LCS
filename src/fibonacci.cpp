#include "fibonacci.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>

namespace LCS {
namespace fkernel {

FibonacciKernel::FibonacciKernel(const std::string &p, int fn): p(p), fn(fn), lca(calculate_lca(p, fn)) {}


matrix::Permutation FibonacciKernel::calculate_char_kernel(const std::string &p, char c) {
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
        result.erase(result.begin());
    }
    return matrix::Permutation(result);
}


// return permutation split into stuff touching left and not
std::pair <matrix::Permutation, matrix::Permutation> get_left(const matrix::Permutation &p, unsigned m) {
    // from left => row index is from 1 to m
    unsigned split_value = m;
    std::vector <std::pair <unsigned, unsigned> > row_first, row_second;
    std::vector <std::pair <unsigned, unsigned> > col_first, col_second;
    for (const auto &matched_pair: p.rows) {
        if (matched_pair.first <= split_value) {
            row_first.push_back(matched_pair);
        }
        if (matched_pair.first > split_value) {
            row_second.push_back(matched_pair);
        }
    }
    for (const auto &matched_pair: p.cols) {
        if (matched_pair.second <= split_value) {
            col_first.push_back(matched_pair);
        }
        if (matched_pair.second > split_value) {
            col_second.push_back(matched_pair);
        }
    }
    return {matrix::Permutation{row_first, col_first}, matrix::Permutation{row_second, col_second}};
}

// return permutation split into stuff touching right and not
std::pair <matrix::Permutation, matrix::Permutation> get_right(const matrix::Permutation &p, unsigned min_m) {
    // to right => col index is from n + 1 to n + m
    unsigned split_value = min_m;
    std::vector <std::pair <unsigned, unsigned> > row_first, row_second;
    std::vector <std::pair <unsigned, unsigned> > col_first, col_second;
    for (const auto &matched_pair: p.rows) {
        if (matched_pair.second <= split_value) {
            row_first.push_back(matched_pair);
        }
        if (matched_pair.second > split_value) {
            row_second.push_back(matched_pair);
        }
    }
    for (const auto &matched_pair: p.cols) {
        if (matched_pair.first <= split_value) {
            col_first.push_back(matched_pair);
        }
        if (matched_pair.first > split_value) {
            col_second.push_back(matched_pair);
        }
    }
    return {matrix::Permutation{row_first, col_first}, matrix::Permutation{row_second, col_second}};
}

matrix::Permutation combine(matrix::Permutation &left_side,
                            matrix::Permutation &both_sides,
                            matrix::Permutation &right_side,
                            unsigned add) {
    // left side needs no fixing
    for (auto &i : both_sides.rows) {
        i.second += add; // only if something?
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


matrix::Permutation FibonacciKernel::calculate_count_kernel(const std::string &p, int fn, int &fm) {
    if (fn == 0) { // string A
        fm = 1;
        return calculate_char_kernel(p, 'A');
    } else if (fn == 1) { // string AB
        int x1 = 0;
        std::cout << "first\n";
        auto first_half = calculate_count_kernel(p, fn - 1, x1);
        std::cout << "second\n";
        auto second_half = calculate_char_kernel(p, 'B');
        std::cout << "got\n";
        fm = x1 + 1;

        std::cout << "right\n";
        auto to_right = get_right(first_half, x1);
        std::cout << "left\n";
        auto from_left = get_left(second_half, p.size());
        std::cout << "int\n";
        auto intersection = to_right.second * from_left.first;
        std::cout << "comb\n";
        return combine(to_right.first, intersection, from_left.second, x1);
    } else {  // F_{fn} = F_{fn - 1} + F_{fn - 2}
        int x1 = 0, x2 = 0;
        auto first_half = calculate_count_kernel(p, fn - 1, x1);
        auto second_half = calculate_count_kernel(p, fn - 2, x2);
        fm = x1 + x2;

        auto to_right = get_right(first_half, x1);
        auto from_left = get_left(second_half, p.size());
        auto intersection = to_right.second * from_left.first;
        return combine(to_right.first, intersection, from_left.second, x1);
    }
}

unsigned FibonacciKernel::calculate_lca(const std::string &p, int fn) {
    int size = 0;
    auto kernel = calculate_count_kernel(p, fn, size);
    for (auto i : kernel.rows) {
        std::cout << i.first << ',' << i.second << ' ';
    } std::cout << '\n';
    // Count lca from the kernel
    //  return n - 0 - kernel_sum(a.size(), n);
    unsigned count_dom = 0;
    for (auto element: kernel.rows) {
        count_dom += (element.first < p.size());// && element.second < (unsigned)size);
    }
    std::cout << count_dom << ' ' << p.size() - count_dom << '\n';

    return 0;
}

}  // namespace fkernel
}  // namespace LCS

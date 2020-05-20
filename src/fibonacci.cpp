#include "fibonacci.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>

namespace LCS {
namespace fkernel {

FibonacciKernel::FibonacciKernel(const std::string &p, int fn): p(p), fn(fn), kernel(calculate_kernel(p, fn)) {}



// std::pair <matrix::Permutation, matrix::Permutation> FibonacciKernel::calculate_char_kernel(const std::string &p, char c) {
//     std::cout << "start " << p << ' ' << c << '\n';
//     unsigned last_row = p.size();
//     std::vector <unsigned> last_col(p.size());
//     for (unsigned i = 0; i < p.size(); ++i) {
//         last_col[i] = p.size() - i - 1;
//         if (p[i] == c || last_col[i] > last_row) {
//             std::swap(last_col[i], last_row);
//         }
//     }
//     std::vector <unsigned> result(p.size() + 1);
//     last_col.push_back(last_row);
//         std::cout << "lc ";
//     for (auto i : last_col) std::cout << i << ' ';
//         std::cout << '\n';
//     for (unsigned i = 0; i <= p.size(); ++i) {
//         result[last_col[i]] = p.size() + 1 - i;
//     }
//     result.pop_back();
//     return {matrix::Permutation({{last_col.back() + 1, 1}},
//                                 {{1, last_col.back() + 1}}),
//                                 matrix::Permutation(result)};
// }


std::pair <matrix::Permutation, matrix::Permutation> FibonacciKernel::calculate_char_kernel(const std::string &p, char c) {
    std::cout << "start " << p << ' ' << c << '\n';
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
        std::cout << "lc ";
    for (auto i : last_col) std::cout << i << ' ';
        std::cout << '\n';
    for (unsigned i = 0; i <= p.size(); ++i) {
        result[last_col[i]] = i + 1;
    }
    result.pop_back();
    return {matrix::Permutation({{last_col[0] + 1, 1}},
                                {{1, last_col[0] + 1}}),
                                matrix::Permutation(result)};
}

// Returns pair of sub.s(0) and sub.s(n) kernels
std::pair <matrix::Permutation, matrix::Permutation> FibonacciKernel::calculate_kernel(const std::string &p, int fn) {
    if (fn == 0) { // string A
        return calculate_char_kernel(p, 'A');
    } else if (fn == 1) { // string AB
        auto first_half = calculate_kernel(p, fn - 1);
        std::cout << "SIZE " << p.size() << '\n';
        std::cout << "subs0 1: "; for (auto i : first_half.first.rows) { std::cout << i.first << ',' << i.second << ' ';} std::cout << '\n';
        std::cout << "subsn 1: "; for (auto i : first_half.second.rows) { std::cout << i.first << ',' << i.second << ' ';} std::cout << '\n';
        auto second_half = calculate_char_kernel(p, 'B');
        std::cout << "subs0 1: "; for (auto i : second_half.first.rows) { std::cout << i.first << ',' << i.second << ' ';} std::cout << '\n';
        std::cout << "subsn 1: "; for (auto i : second_half.second.rows) { std::cout << i.first << ',' << i.second << ' ';} std::cout << '\n';
        
        // matrix::Permutation test1 = first_half.first * second_half.second;
        // std::cout << "test 1: "; for (auto i : test1.rows) { std::cout << i.first << ',' << i.second << ' ';} std::cout << '\n';
        // matrix::Permutation test2 = second_half.first * first_half.second;
        // std::cout << "test 2: "; for (auto i : test2.rows) { std::cout << i.first << ',' << i.second << ' ';} std::cout << '\n';
        // matrix::Permutation test3 = first_half.second * second_half.first;
        // std::cout << "test 3: "; for (auto i : test3.rows) { std::cout << i.first << ',' << i.second << ' ';} std::cout << '\n';
        // matrix::Permutation test4 = second_half.second * first_half.first;
        // std::cout << "test 4: "; for (auto i : test4.rows) { std::cout << i.first << ',' << i.second << ' ';} std::cout << '\n';
        
        // first_half.first.grow_front(p.size() + 1);
        first_half.second.grow_back(p.size() + 1);
        // second_half.first.grow_back(p.size() + 1);
        second_half.second.grow_front(p.size() + 1);
        // test1.grow_front(p.size() + 1); test2.grow_front(p.size() + 1);
        // test3.grow_front(p.size() + 1); test4.grow_front(p.size() + 1);
        // std::cout << "test 1: "; for (auto i : test1.rows) { std::cout << i.first << ',' << i.second << ' ';} std::cout << '\n';
        // std::cout << "test 2: "; for (auto i : test2.rows) { std::cout << i.first << ',' << i.second << ' ';} std::cout << '\n';
        // std::cout << "test 3: "; for (auto i : test3.rows) { std::cout << i.first << ',' << i.second << ' ';} std::cout << '\n';
        // std::cout << "test 4: "; for (auto i : test4.rows) { std::cout << i.first << ',' << i.second << ' ';} std::cout << '\n';
        // matrix::Permutation ttet1 = test1 * test3;
        // std::cout << "ttet 1: "; for (auto i : ttet1.rows) { std::cout << i.first << ',' << i.second << ' ';} std::cout << '\n';
        // matrix::Permutation ttet2 = test2 * test4;
        // std::cout << "ttet 2: "; for (auto i : ttet2.rows) { std::cout << i.first << ',' << i.second << ' ';} std::cout << '\n';
        // matrix::Permutation ttet3 = test3 * test1;
        // std::cout << "ttet 3: "; for (auto i : ttet3.rows) { std::cout << i.first << ',' << i.second << ' ';} std::cout << '\n';
        // matrix::Permutation ttet4 = test4 * test2;
        // std::cout << "ttet 4: "; for (auto i : ttet4.rows) { std::cout << i.first << ',' << i.second << ' ';} std::cout << '\n';
        // now join, n = 2
        matrix::Permutation subs0 = first_half.first * second_half.first; // ideally: (10, 12), (12, 13)
        matrix::Permutation subsn = first_half.second * second_half.second;

        return {subs0, subsn};
    } else {  // F_{fn} = F_{fn - 2} + F_{fn - 1}
        auto first_half = calculate_kernel(p, fn - 2);
        auto second_half = calculate_kernel(p, fn - 1);
        return first_half;
        // first_half.grow_back(sum_length);
        // second_half.grow_front(sum_length);
        // return first_half * second_half;
    }
}


}  // namespace fkernel
}  // namespace LCS

#include "grammar_compressed.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <numeric>

namespace LCS {
namespace gc {

const int ALPHABET_SIZE = 26;

// Compress the string s with alphabet characters 'A'-'Z' using LZ78 compression.
GrammarCompressed LZ78(const std::string &s) {
    std::vector <GrammarCompressed> gcs(1);
    int current_entry = 0;  // The entry corresponding to the current buffer.
    std::vector <std::vector <int>> next_entry(1, std::vector <int> (ALPHABET_SIZE, 0));
    int last_string_entry = 0;  // The last entry corresponding to the piece of a string.
    for (unsigned int i = 0; i < s.size(); ++i) {
        int c = s[i] - 'A';
        if (next_entry[current_entry][c] != 0 && i + 1 != s.size()) {
            // If current prefix + c is in the dictionary, and the string has not ended.
            current_entry = next_entry[current_entry][c];
        } else {
            // Add current character as string to grammar.
            int dict_char = gcs.size();
            int dict_entry = dict_char;
            gcs.push_back(GrammarCompressed(dict_char, s[i]));
            next_entry.push_back(std::vector <int> (ALPHABET_SIZE, 0));

            // Add the new string (current_entry + c) to the dictionary.
            if (!current_entry) {
                next_entry[current_entry][c] = dict_char;
            } else {  // A previous non-empty entry existed.
                dict_entry = gcs.size();
                gcs.push_back(GrammarCompressed(dict_entry, gcs[current_entry], gcs[dict_char]));
                next_entry.push_back(std::vector <int> (ALPHABET_SIZE, 0));
                next_entry[current_entry][c] = dict_entry;
                current_entry = 0;
            }

            // Concatenate two dictionary strings, if necessary.
            if (!last_string_entry) {
                last_string_entry = dict_entry;
            } else {
                int string_entry = gcs.size();
                gcs.push_back(GrammarCompressed(string_entry, gcs[last_string_entry], gcs[dict_entry]));
                next_entry.resize(next_entry.size() + 1);
                last_string_entry = string_entry;
            }
        }
    }
    return gcs[last_string_entry];
}

// Compress the string s with alphabet characters 'A'-'Z' using LZW compression.
GrammarCompressed LZW(const std::string &s) {
    std::vector <GrammarCompressed> gcs(1);
    int current_entry = 0;  // The entry corresponding to the current buffer.
    std::vector <std::vector <int>> next_entry(ALPHABET_SIZE + 1, std::vector <int> (ALPHABET_SIZE, 0));
    // Initialize LZW alphabet.
    for (unsigned int i = 0; i < ALPHABET_SIZE; ++i) {
        next_entry[0][i] = i + 1;
        gcs.push_back(GrammarCompressed(i + 1, (char)('A' + i)));
    }
    int last_string_entry = 0;  // The last entry corresponding to the piece of a string.
    for (unsigned int i = 0; i < s.size(); ++i) {
        int c = s[i] - 'A';
        if (next_entry[current_entry][c] != 0 && i + 1 != s.size()) {
            // If current prefix + c is in the dictionary, and the string has not ended.
            current_entry = next_entry[current_entry][c];
        } else {
            // Add the new string (current_entry + c) to the dictionary.
            int dict_char = c + 1, dict_entry = gcs.size();
            gcs.push_back(GrammarCompressed(dict_entry, gcs[current_entry], gcs[dict_char]));
            next_entry.push_back(std::vector <int> (ALPHABET_SIZE, 0));
            next_entry[current_entry][c] = dict_entry;
            current_entry = 0;

            // Concatenate two dictionary strings, if necessary.
            if (!last_string_entry) {
                last_string_entry = dict_entry;
            } else {
                int string_entry = gcs.size();
                gcs.push_back(GrammarCompressed(string_entry, gcs[last_string_entry], gcs[dict_entry]));
                next_entry.resize(next_entry.size() + 1);
                last_string_entry = string_entry;
            }
        }
    }
    return gcs[last_string_entry];
}


GCKernel::GCKernel(const std::string &p, const GrammarCompressed &t): lcs(calculate_lcs(p, t)) {}

// Returns permutation split into strings touching the left side and not.
std::pair <matrix::Permutation, matrix::Permutation> get_left(const matrix::Permutation &p,
                                                              unsigned left, unsigned right) {
    // from left => row index is from 1 to m
    auto result = p.split_row(left);
    return {result.first, result.second.split_col(right).second};
}

// Returns permutation split into strings touching the right side and not.
std::pair <matrix::Permutation, matrix::Permutation> get_right(const matrix::Permutation &p,
                                                               unsigned left, unsigned right) {
    // to right => col index is from n + 1 to n + m
    auto result = p.split_col(right);
    return {result.first.split_row(left).first, result.second};
}

// Returns permutation with compressed coordinates (values from 1 to x).
matrix::Permutation compress(const matrix::Permutation &uncompressed) {
    std::vector <int> row_values(uncompressed.rows.size());
    std::vector <int> col_values(uncompressed.rows.size());
    for (unsigned i = 0; i < row_values.size(); ++i) {
        row_values[i] = uncompressed.rows[i].first;
        col_values[i] = uncompressed.rows[i].second;
    }
    std::sort(col_values.begin(), col_values.end());
    std::reverse(row_values.begin(), row_values.end());
    std::unordered_map <unsigned, unsigned> row_compression;
    std::unordered_map <unsigned, unsigned> col_compression;
    for (unsigned i = 0; i < row_values.size(); ++i) {
        row_compression[row_values[i]] = i + 1;
    }
    for (unsigned i = 0; i < col_values.size(); ++i) {
        col_compression[col_values[i]] = i + 1;
    }
    std::vector <std::pair <unsigned, unsigned> > compressed_rows(uncompressed.rows.size());
    std::vector <std::pair <unsigned, unsigned> > compressed_cols(uncompressed.cols.size());
    for (unsigned i = 0; i < uncompressed.rows.size(); ++i) {
        compressed_rows[i].first = row_compression[uncompressed.rows[i].first];
        compressed_rows[i].second = col_compression[uncompressed.rows[i].second];
    }
    for (unsigned i = 0; i < uncompressed.cols.size(); ++i) {
        compressed_cols[i].first = col_compression[uncompressed.cols[i].first];
        compressed_cols[i].second = row_compression[uncompressed.cols[i].second];
    }
    return matrix::Permutation{compressed_rows, compressed_cols};
}

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
    return compress(matrix::Permutation(result));
}


// Collects a permutation from the strings adjacent to the left side, right side and both.
matrix::Permutation combine(matrix::Permutation &left_side,
                            matrix::Permutation &both_sides,
                            matrix::Permutation &right_side,
                            unsigned row_add, unsigned col_add) {
    for (auto &i : both_sides.rows) {
        i.second += col_add;
    }
    for (auto &i : both_sides.cols) {
        i.first += col_add;
    }
    for (auto &i : right_side.rows) {
        i.first += row_add;
        i.second += col_add;
    }
    for (auto &i : right_side.cols) {
        i.first += col_add;
        i.second += row_add;
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
    return compress(matrix::Permutation(all_combined_rows, all_combined_cols));
}


void GCKernel::calculate_gc_kernel(std::vector<matrix::Permutation> &calculated, 
                                                  const std::string &p, const GrammarCompressed &t) {
    if (calculated[t.number].get_nonzero_amount()) {
    } else if (t.is_base) { // t is a single symbol
        calculated[t.number] = calculate_char_kernel(p, t.value);
    } else { // t = uv
        auto first_half = *t.first_symbol;
        auto second_half = *t.second_symbol;
        if (!calculated[first_half.number].get_nonzero_amount()) {
            calculate_gc_kernel(calculated, p, first_half);
        }
        if (!calculated[second_half.number].get_nonzero_amount()) {
            calculate_gc_kernel(calculated, p, second_half);
        }

        auto to_right = get_right(calculated[first_half.number], p.size(),
                                    calculated[first_half.number].cols.size() - p.size());
        auto from_left = get_left(calculated[second_half.number], p.size(),
                                    calculated[second_half.number].cols.size() - p.size());
        auto intersection = to_right.second * from_left.first;

        calculated[t.number] = combine(to_right.first, intersection, from_left.second,
                                       calculated[first_half.number].rows.size() - p.size(),
                                       calculated[first_half.number].cols.size() - p.size());
    }
}

unsigned GCKernel::calculate_lcs(const std::string &p, const GrammarCompressed &t) {
    std::vector <matrix::Permutation> kernels(t.number + 1);
    calculate_gc_kernel(kernels, p, t);
    auto kernel = kernels.back();
    int size = kernel.cols.size() - p.size();
    unsigned count_dom = 0;
    for (auto i: kernel.rows) {
        count_dom += (i.first <= p.size() && i.second > (unsigned)size);
    }
    return p.size() - count_dom;
}

}  // namespace gc
}  // namespace LCS

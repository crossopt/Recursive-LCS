#include "grammar_compressed.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <numeric>

namespace LCS {
namespace gc {

const int ALPHABET_SIZE = 26;
const int ASCII_SIZE = 256;


// Returns the nth LZ78 grammar string. A LZ78 grammar string is a specially constructed string
// that is compressed quadratically by the LZ78 compression.
std::string get_lz78_grammar_string(unsigned int number) {
    std::vector <std::string> trie_entries(number + 1);
    trie_entries[0] = "A";
    for (unsigned int i = 1; i <= number; ++i) {
        trie_entries[i] = trie_entries[i - 1] + (char)('A' + i % ALPHABET_SIZE);
    }
    std::string result;
    for (unsigned int i = 0; i <= number; ++i) {
        result += trie_entries[i];
    }
    return result;
}

// Returns the nth LZW grammar string. A LZW grammar string is a specially constructed string
// that is compressed quadratically by the LZW compression.
std::string get_lzw_grammar_string(unsigned int number) {
    std::vector <std::string> trie_entries(number + 1);
    trie_entries[0] = "AA";
    for (unsigned int i = 1; i <= number; ++i) {
        trie_entries[i] = trie_entries[i - 1] + (char)('A' + i % ALPHABET_SIZE);
    }
    std::string result;
    for (unsigned int i = 0; i <= number; ++i) {
        result += trie_entries[i];
    }
    return result;
}

// Returns the nth LZ grammar string. A LZW grammar string is a specially constructed string
// that is compressed quadratically by the LZW compression.
std::string get_lz_grammar_string(unsigned int number) {
    std::vector <std::string> trie_entries(number + 1);
    std::string add = "ABAA";
    trie_entries[0] = "AAB";
    for (unsigned int i = 0; i + 1 < number; ++i) {
        trie_entries[i + 1] = trie_entries[i] + (char)('A' + (i + 2) % ALPHABET_SIZE);
    }
    std::string result = add;
    for (unsigned int i = 0; i < number; ++i) {
        result += trie_entries[i];
    }
    return result;
}

// Compress the string s with alphabet characters 'A'-'Z' using LZ78 compression.
GrammarCompressedStorage LZ78(const std::string &s) {
    GrammarCompressedStorage gcs = GrammarCompressedStorage();
    std::vector <unsigned int> gcs_index(1);
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
            int dict_char = gcs_index.size();
            int dict_entry = dict_char;
            gcs_index.push_back(gcs.rules.size());
            gcs.add_rule(GrammarCompressed(gcs, dict_char, s[i]));
            next_entry.push_back(std::vector <int> (ALPHABET_SIZE, 0));

            // Add the new string (current_entry + c) to the dictionary.
            if (!current_entry) {
                next_entry[current_entry][c] = dict_char;
            } else {  // A previous non-empty entry existed.
                dict_entry = gcs_index.size();
                gcs_index.push_back(gcs.rules.size());
                gcs.add_rule(GrammarCompressed(gcs, dict_entry, gcs_index[current_entry], gcs_index[dict_char]));
                next_entry.push_back(std::vector <int> (ALPHABET_SIZE, 0));
                next_entry[current_entry][c] = dict_entry;
                current_entry = 0;
            }

            // Concatenate two dictionary strings, if necessary.
            if (!last_string_entry) {
                last_string_entry = dict_entry;
            } else {
                int string_entry = gcs_index.size();
                gcs_index.push_back(gcs.rules.size());
                gcs.add_rule(GrammarCompressed(gcs, string_entry, gcs_index[last_string_entry], gcs_index[dict_entry]));
                next_entry.resize(next_entry.size() + 1);
                last_string_entry = string_entry;
            }
        }
    }
    gcs.final_rule = gcs_index[last_string_entry];
    return gcs;
}

// Compress the string s with alphabet characters 'A'-'Z' using LZW compression.
GrammarCompressedStorage LZW(const std::string &s) {
    GrammarCompressedStorage gcs = GrammarCompressedStorage();
    std::vector <unsigned int> gcs_index(1);
    int current_entry = 0;  // The entry corresponding to the current buffer.
    std::vector <std::vector <int>> next_entry(ALPHABET_SIZE + 1, std::vector <int> (ALPHABET_SIZE, 0));
    // Initialize LZW alphabet.
    for (unsigned int i = 0; i < ALPHABET_SIZE; ++i) {
        next_entry[0][i] = i + 1;
        gcs_index.push_back(gcs.rules.size());
        gcs.add_rule(GrammarCompressed(gcs, i + 1, (char)('A' + i)));
    }
    int last_string_entry = 0;  // The last entry corresponding to the piece of a string.
    for (unsigned int i = 0; i < s.size(); ++i) {
        int c = s[i] - 'A';
        if (next_entry[current_entry][c] != 0 && i + 1 != s.size()) {
            // If current prefix + c is in the dictionary, and the string has not ended.
            current_entry = next_entry[current_entry][c];
        } else {
            // Add the new string (current_entry + c) to the dictionary.
            int dict_char = c + 1, dict_entry = gcs_index.size();
            gcs_index.push_back(gcs.rules.size());
            gcs.add_rule(GrammarCompressed(gcs, dict_entry, gcs_index[current_entry], gcs_index[dict_char]));
            next_entry.push_back(std::vector <int> (ALPHABET_SIZE, 0));
            next_entry[current_entry][c] = dict_entry;
            current_entry = 0;

            // Concatenate two dictionary strings, if necessary.
            if (!last_string_entry) {
                last_string_entry = dict_entry;
            } else {
                int string_entry = gcs_index.size();
                gcs_index.push_back(gcs.rules.size());
                gcs.add_rule(GrammarCompressed(gcs, string_entry, gcs_index[last_string_entry], gcs_index[dict_entry]));
                next_entry.resize(next_entry.size() + 1);
                last_string_entry = string_entry;
            }
        }
    }
    gcs.final_rule = gcs_index[last_string_entry];
    return gcs;
}

GrammarCompressedStorage get_aaaa(unsigned long long number) {
    GrammarCompressedStorage gcs = GrammarCompressedStorage();
    std::vector <unsigned long long> gcs_index(1);
    unsigned long long current_entry = 0;  // The entry corresponding to the current buffer.
    std::vector <unsigned long long> next_entry(ASCII_SIZE + 1, 0);
    // Initialize LZW alphabet.
    next_entry[0] = 'a' + 1;
    for (unsigned int i = 0; i < ASCII_SIZE; ++i) {
        gcs_index.push_back(gcs.rules.size());
        gcs.add_rule(GrammarCompressed(gcs, i + 1, i));
    }
    unsigned long long last_string_entry = 0;  // The last entry corresponding to the piece of a string.
    for (unsigned long long i = 0; i < number; ++i) {
        // int c = 'a';
        if (next_entry[current_entry] != 0 && i + 1 != number) {
            // If current prefix + c is in the dictionary, and the string has not ended.
            current_entry = next_entry[current_entry];
        } else {
            // Add the new string (current_entry + c) to the dictionary.
            // int dict_char = c + 1;
            int dict_entry = gcs_index.size();
            gcs_index.push_back(gcs.rules.size());
            gcs.add_rule(GrammarCompressed(gcs, dict_entry, gcs_index[current_entry], 'a'));
            next_entry.push_back(0);
            next_entry[current_entry] = dict_entry;
            current_entry = 0;

            // Concatenate two dictionary strings, if necessary.
            if (!last_string_entry) {
                last_string_entry = dict_entry;
            } else {
                int string_entry = gcs_index.size();
                gcs_index.push_back(gcs.rules.size());
                gcs.add_rule(GrammarCompressed(gcs, string_entry, gcs_index[last_string_entry], gcs_index[dict_entry]));
                next_entry.resize(next_entry.size() + 1);
                last_string_entry = string_entry;
            }
        }
    }
    gcs.final_rule = gcs_index[last_string_entry];
    return gcs;
}

int intify(char c) {
    return (int)c >= 0 ? (int)c : 256 + (int)c;
}

std::string read_file_contents(const std::string &file_name) {
    std::ifstream t(file_name);
    std::stringstream buffer;
    buffer << t.rdbuf();
    return buffer.str();
}

  
/*
  Uncompress file compressed with UNIX Z-compress.
  This function and the one below were adapted from Mark Adler's C implementation
  on Stackoverflow. Original copyright notice below.

  unlzw version 1.4, 22 August 2015
  Copyright (C) 2014, 2015 Mark Adler
  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.
  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:
  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.
  Mark Adler
  madler@alumni.caltech.edu
*/
std::string get_uncompress_string(const std::string &file_name) {
    std::string input = read_file_contents(file_name);
    unsigned int inlen = input.size();
    std::vector <unsigned int> prefix(65536, 0), suffix(65536, 0);

    unsigned int flags = intify(input[2]);
    unsigned int mx = flags & 0x1f;
    if (mx == 9) mx = 10;

    flags = flags & 0x80;
    unsigned int bits = 9;
    unsigned int mask = 0x1ff;
    unsigned int fend = flags ? 256 : 255;

    if (inlen == 3) {
        return "";
    }

    unsigned int buf = intify(input[3]);
    buf += intify(input[4]) << 8;
    unsigned int prev = buf & mask;

    unsigned int fin = prev;
    buf >>= bits;

    // std::cout << "buf" << ' ' << buf <<  " fin " << fin << '\n';
    unsigned int left = 16 - bits;
    std::string put(1, fin);

    unsigned int mark = 3, nxt = 5;
    while (nxt < inlen) {
         if ((fend >= mask) && (bits < mx)) {
            unsigned int rem = (nxt - mark) % bits;

            if (rem) {
                rem = bits - rem;
                if (rem >= inlen - nxt) {
                    break;
                }
                nxt += rem;
            }

            buf = 0;
            left = 0;

            mark = nxt;

            bits += 1;
            mask <<= 1;
            mask += 1;
         }
        buf += (intify(input[nxt]) << left);
        nxt += 1;
        left += 8;
        if (left < bits) {
            buf += (intify(input[nxt]) << left);
            nxt += 1;
            left += 8;
        }
        unsigned int code = buf & mask;
        // std::cout << "bmc " << buf << ' ' << mask << ' ' << code << '\n';
        buf >>= bits;
        left -= bits;

        if ((code == 256) && flags) {
            unsigned int rem = (nxt - mark) % bits;
            if (rem) {
                rem = bits - rem;
                if (rem > inlen - nxt) {
                    break;
                }
                nxt += rem;
            }
            buf = 0;
            left = 0;

            mark = nxt;

            bits = 9;
            mask = 0x1ff;
            fend = 255;
            continue;
        }

        unsigned int temp = code;
        // std::cout << "code " << code << '\n';
        std::string stck;

        if (code > fend) {
            stck.push_back(fin);
            code = prev;
        }

        while (code >= 256) {
            stck.push_back(suffix[code]);
            code = prefix[code];
        }

        stck.push_back(code);
        fin = code;

        if (fend < mask) {
            fend += 1;
            prefix[fend] = prev;
            suffix[fend] = fin;
            // std::cout << "ADD " << fend << ' ' << prev << ' ' << fin << '\n';
        }
        prev = temp;
        for (int i = (int)stck.size() - 1; i >= 0; --i) {
            put.push_back(stck[i]);
        }
        // std::cout << "stack " << stck << std::endl;
    }
    return put;
}

// Convert file compressed with UNIX Z-compress without decompressing.
GrammarCompressedStorage get_compress_string(const std::string &file_name) {
    std::string input = read_file_contents(file_name);
    unsigned int inlen = input.size();
    std::vector <unsigned int> prefix(65536, 0), suffix(65536, 0), rule_num(65536, 0);

    unsigned int flags = intify(input[2]);
    unsigned int mx = flags & 0x1f;
    if (mx == 9) mx = 10;

    flags = flags & 0x80;
    unsigned int bits = 9;
    unsigned int mask = 0x1ff;
    unsigned int fend = flags ? 256 : 255;

    if (inlen == 3) {
        return GrammarCompressedStorage();
    }
    GrammarCompressedStorage gcs = GrammarCompressedStorage();
    std::vector <unsigned int> gcs_index(1);

    for (unsigned int i = 0; i < ASCII_SIZE; ++i) {
        gcs_index.push_back(gcs.rules.size());
        gcs.add_rule(GrammarCompressed(gcs, i + 1, i));
        // std::cout << "RULE" << i + 1 << "  ( actual " << gcs.rules.size() - 1  << ") \n";
        
        //<< gcs.rules.back().decompress() << '\n';
    }
    // unsigned int last_rule = ASCII_SIZE;

    unsigned int buf = intify(input[3]);
    buf += intify(input[4]) << 8;
    unsigned int prev = buf & mask;

    unsigned int fin = prev;
    buf >>= bits;

    // std::cout << "buf" << ' ' << buf <<  " fin " << fin << '\n';
    unsigned int left = 16 - bits;
    // std::string put(1, fin);
    std::vector <unsigned int> fin_list;
    fin_list.push_back(fin + 1);

    unsigned int mark = 3, nxt = 5;
    while (nxt < inlen) {
         if ((fend >= mask) && (bits < mx)) {
            unsigned int rem = (nxt - mark) % bits;

            if (rem) {
                rem = bits - rem;
                if (rem >= inlen - nxt) {
                    break;
                }
                nxt += rem;
            }

            buf = 0;
            left = 0;

            mark = nxt;

            bits += 1;
            mask <<= 1;
            mask += 1;
         }
        buf += (intify(input[nxt]) << left);
        nxt += 1;
        left += 8;
        if (left < bits) {
            buf += (intify(input[nxt]) << left);
            nxt += 1;
            left += 8;
        }
        unsigned int code = buf & mask;
        // std::cout << "bmc " << buf << ' ' << mask << ' ' << code << '\n';
        buf >>= bits;
        left -= bits;

        if ((code == 256) && flags) {
            unsigned int rem = (nxt - mark) % bits;
            if (rem) {
                rem = bits - rem;
                if (rem > inlen - nxt) {
                    break;
                }
                nxt += rem;
            }
            buf = 0;
            left = 0;

            mark = nxt;

            bits = 9;
            mask = 0x1ff;
            fend = 255;
            continue;
        }

        unsigned int temp = code;
        // std::cout << "code " << code << '\n';
        // std::string stck;

        if (code > fend) {
            // stck.push_back(fin);
            code = prev;
        }

        while (code >= 256) {
            // stck.push_back(suffix[code]);
            code = prefix[code];
        }

        // stck.push_back(code);
        fin = code;

        // перечитать дяц и все переисать на интедкс
        if (fend < mask) {
            fend += 1;
            prefix[fend] = prev;
            suffix[fend] = fin;
            gcs_index.push_back(gcs.rules.size());
            rule_num[fend] = gcs_index.back();

            if (rule_num[prev] != 0 && rule_num[prev] != 1) {
                // std::cout << "SANITY CHECK " << suffix[fend] << ' ' << code << '\n';
                gcs.add_rule(GrammarCompressed(gcs, rule_num[fend] + 1, rule_num[prev], suffix[fend]));
                // std::cout << "RULE ASS GCS   " << rule_num[fend] + 1 << ' ' << rule_num[prev] << ' ' << suffix[fend] << "   "  << gcs.rules.size() << "    " << gcs.rules.back().decompress() << '\n';
            } else {
                gcs.add_rule(GrammarCompressed(gcs, rule_num[fend] + 1, suffix[fend]));
                // std::cout << "RULE   " << rule_num[fend] << ' ' << suffix[fend] << '\n';
            }
            fin_list.push_back(gcs_index.back());
            // std::cout << "add fin: " << gcs_index.back() << '\n';
            // std::cout << "ADD " << fend << ' ' << prev << ' ' << fin << '\n';
        }
        // else {
        //    std::cout << "ELSE FEND " << fend << ' ' << rule_num[fend] << '\n';
        // }
        prev = temp;
        // for (int i = (int)stck.size() - 1; i >= 0; --i) {
        //     put.push_back(stck[i]);
        // }
        // std::cout << "stack " << stck << std::endl;
    }
    // for (auto i: fin_list) {
    //     std::cout << "qqqqq " << gcs.rules[i].decompress() << '\n';
    // }
    unsigned int last_fin_rule = gcs_index[fin_list[0]];
    for (unsigned int i = 1; i < fin_list.size(); ++i) {
        gcs_index.push_back(gcs.rules.size());
        gcs.add_rule(GrammarCompressed(gcs, gcs_index.back() + 1, last_fin_rule, fin_list[i]));
        // std::cout << "RULE CONC GCS   " << gcs_index.back() + 1 << "   " <<  last_fin_rule << ' ' << fin_list[i] << "   "  << '\n';
        last_fin_rule = gcs_index.back();
    }

    gcs.final_rule = last_fin_rule;
    // std::cout << "FINAL RULE IS " << gcs.final_rule << ' ' << gcs.rules.size() << '\n';
    auto res = gcs.rules[gcs.final_rule].decompress();
    // std::cout << "RES:  " <<  res << '\n';
    for (unsigned int i = 0; i < gcs.rules.size(); ++i) {
        // std::cout << 'r' << i << ' ' << gcs.rules[i].is_base << ' ' << gcs.rules[i].value << ' ' << gcs.rules[i].first_symbol << ' ' << gcs.rules[i].second_symbol << '\n';
    }
    return gcs;
}


GCKernel::GCKernel(const std::string &p, const GrammarCompressedStorage &t): lcs(calculate_lcs(p, t)) {}

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
                                                  const std::string &p, const GrammarCompressedStorage &t, unsigned int index) {
    if (calculated[t.rules[index].number].get_nonzero_amount()) {
    } else if (t.rules[index].is_base) { // t is a single symbol
        calculated[t.rules[index].number] = calculate_char_kernel(p, t.rules[index].value);
    } else { // t = uv
        unsigned int first_half = t.rules[index].first_symbol;
        unsigned int second_half = t.rules[index].second_symbol;
        unsigned int first_number = t.rules[first_half].number;
        unsigned int second_number = t.rules[second_half].number;
        if (!calculated[first_number].get_nonzero_amount()) {
            calculate_gc_kernel(calculated, p, t, first_half);
        }
        if (!calculated[second_number].get_nonzero_amount()) {
            calculate_gc_kernel(calculated, p, t, second_half);
        }

        auto to_right = get_right(calculated[first_number], p.size(), calculated[first_number].cols.size() - p.size());
        auto from_left = get_left(calculated[second_number], p.size(), calculated[second_number].cols.size() - p.size());
        auto intersection = to_right.second * from_left.first;

        calculated[t.rules[index].number] = combine(to_right.first, intersection, from_left.second,
                                                    calculated[first_number].rows.size() - p.size(),
                                                    calculated[first_number].cols.size() - p.size());
    }
}

unsigned GCKernel::calculate_lcs(const std::string &p, const GrammarCompressedStorage &t) {
    std::vector <matrix::Permutation> kernels(t.rules[t.final_rule].number + 1);
    calculate_gc_kernel(kernels, p, t, t.final_rule);
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

#ifndef INC_FIBONACCI_H_
#define INC_FIBONACCI_H_

#include <string>
#include <iostream>
#include <memory>

#include "lcs_kernel.h"
#include "monge_matrix.h"

namespace LCS {
namespace gc {

class GrammarCompressed;

// Class that stores a context-free grammar for GC-string generation.
class GrammarCompressedStorage {
public:
	unsigned final_rule; // The outer rule index for the rule that sets the whole string.
	std::vector <GrammarCompressed> rules; // All rules in the grammar.

	void add_rule(const GrammarCompressed &rule) {
		rules.push_back(rule);
	}
};

// Class that stores a context-free grammar rule for GC-string generation.
class GrammarCompressed {
public:
	GrammarCompressedStorage &gc_storage;
	const int number;  // the number of the current rule in the grammar
	const bool is_base;  // is the current rule a single alphabet value
	const char value;  // the alphabet value for base rules
	unsigned int first_symbol;  // the index of the left rule for non-bases
	unsigned int second_symbol;  // the index of the right rule for non-bases

    GrammarCompressed(GrammarCompressedStorage &gc_storage): gc_storage(gc_storage), number(0), is_base(0), value(0), first_symbol(0), second_symbol(0) {}

	explicit GrammarCompressed(GrammarCompressedStorage &gc_storage,
		                       int number, char c): gc_storage(gc_storage),
													number(number),
													is_base(true),
													value(c),
													first_symbol(0),
													second_symbol(0) {}
	GrammarCompressed(GrammarCompressedStorage &gc_storage, int number, unsigned int first_symbol, unsigned int second_symbol):
	    gc_storage(gc_storage),
		number(number),
		is_base(false),
		value(0), 
		first_symbol(first_symbol),
		second_symbol(second_symbol) {}

	// Returns the decompressed string for the grammar.
	std::string decompress() const {
		// if (is_base) {
		// 	std::cout << "decompressing char rule " << number << " " << value << " " << std::string(1, value) << "\n";
		// } else {
		// 	std::cout << "decompressing rule " << number << " " <<  first_symbol << " " <<  second_symbol << "\n";
		// }
		return is_base ? std::string(1, value) : gc_storage.rules[first_symbol].decompress() 
		                                       + gc_storage.rules[second_symbol].decompress();
	}
};

GrammarCompressedStorage LZ78(const std::string &s);
GrammarCompressedStorage LZW(const std::string &s);
GrammarCompressedStorage LZW2(const std::string &s);
GrammarCompressedStorage LZWASCII(const std::string &s);
std::string get_lz78_grammar_string(unsigned int number, unsigned int repeat_number = 0);
std::string get_lzw_grammar_string(unsigned int number, unsigned int repeat_number = 0);
std::string get_lz_grammar_string(unsigned int number);

std::string get_uncompress_string(const std::string &file_name);
GrammarCompressedStorage get_compress_string(const std::string &file_name);

// Time agrep, delete later.
GrammarCompressedStorage get_aaaa(unsigned long long number);

// Class that calculates the LCS kernel to solve the semi-local LCS problem
// for a plain pattern and a grammar-compressed text.
class GCKernel {
public:
    // Initialize the LCS kernel for pattern p and text t.
    GCKernel(const std::string &p, const GrammarCompressedStorage &t);
    const unsigned lcs;
private:
	// Returns the lcs for pattern p and text t.
    unsigned calculate_lcs(const std::string &p,  const GrammarCompressedStorage &t);
    // Returns the compressed kernel for pattern p and character c.
    matrix::Permutation calculate_char_kernel(const std::string &p, char c);
	// Recursively calculates the compressed kernel for pattern p and text t.
    void calculate_gc_kernel(std::vector<matrix::Permutation> &calculated,
    						 				const std::string &p, const GrammarCompressedStorage &t, unsigned int index);
};


}  // namespace gc
}  // namespace LCS

#endif

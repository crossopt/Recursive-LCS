#ifndef INC_FIBONACCI_H_
#define INC_FIBONACCI_H_

#include <string>
#include <iostream>
#include <memory>

#include "lcs_kernel.h"
#include "monge_matrix.h"

namespace LCS {
namespace gc {

// Class that stores a context-free grammar for GC-string generation.
class GrammarCompressed {
public:
	const int number;  // the number of the current rule in the grammar
	const bool is_base;  // is the current rule a single alphabet value
	const char value;  // the alphabet value for base rules
	const std::shared_ptr<GrammarCompressed> first_symbol;  // the left rule for non-bases
	const std::shared_ptr<GrammarCompressed> second_symbol;  // the right rule for non-bases

	explicit GrammarCompressed(int number, char c): number(number),
													is_base(true),
													value(c),
													first_symbol(0),
													second_symbol(0) {}
	GrammarCompressed(int number, const GrammarCompressed &first_symbol,
								  const GrammarCompressed &second_symbol):
		number(number),
		is_base(false),
		value(0), 
		first_symbol(std::make_shared<GrammarCompressed>(first_symbol)),
		second_symbol(std::make_shared<GrammarCompressed>(second_symbol)) {}

	// Returns the decompressed string for the grammar.
	std::string decompress() const {
		return is_base ? std::string(1, value) : first_symbol->decompress() + second_symbol->decompress();
	}
};

// Class that calculates the LCS kernel to solve the semi-local LCS problem
// for a plain pattern and a grammar-compressed text.
class GCKernel {
public:
    // Initialize the LCS kernel for pattern p and text t.
    GCKernel(const std::string &p, const GrammarCompressed &t);
    const unsigned lcs;
private:
	// Returns the lcs for pattern p and text t.
    unsigned calculate_lcs(const std::string &p,  const GrammarCompressed &t);
    // Returns the compressed kernel for pattern p and character c.
    matrix::Permutation calculate_char_kernel(const std::string &p, char c);
	// Recursively calculates the compressed kernel for pattern p and text t.
    void calculate_gc_kernel(std::vector<matrix::Permutation> &calculated,
    						 				const std::string &p, const GrammarCompressed &t);
};


}  // namespace gc
}  // namespace LCS

#endif

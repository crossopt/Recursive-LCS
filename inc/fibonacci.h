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
	const bool is_base;
	const char value;
	const std::shared_ptr<GrammarCompressed> first_symbol;
	const std::shared_ptr<GrammarCompressed> second_symbol;

	explicit GrammarCompressed(char c): is_base(true), value(c), first_symbol(0), second_symbol(0) {}
	GrammarCompressed(const GrammarCompressed &first_symbol, const GrammarCompressed &second_symbol):
		is_base(false), value(0), 
		first_symbol(std::make_shared<GrammarCompressed>(first_symbol)),
		second_symbol(std::make_shared<GrammarCompressed>(second_symbol)) {}

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

    unsigned calculate_lcs(const std::string &p,  const GrammarCompressed &t);
    matrix::Permutation calculate_count_kernel(const std::string &p, const GrammarCompressed &t, int &fm);
    matrix::Permutation calculate_char_kernel(const std::string &p, char c);
};


}  // namespace fkernel
}  // namespace LCS

#endif

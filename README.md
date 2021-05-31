# Recursive-LCS
An implementation of an algorithm for computing semi-local LCS using recursive combing, also its generalization for grammar-compressed texts.


### Build project
Build using cmake, from root:

    mkdir build
    cd build
    cmake ..
    make

### Build targets:
* ./lcs_test: tests for everything, including semi-local LCS and grammar-compressed LCS
* ./time_test: various util functions used to time grammar-compressed LCS
* ./main: run recursive LCS, not used lately

### Files & Bachelor's relevant code:
* src/monge_matrix: utils for working with permutations & monge matrices, steady ant algorithm
* src/lcs_kernel: recursive lcs for uncompressed strings
* src/grammar_compressed: all code relevant for GC-compressed strings: LCS calculation & compressed formats
   * generation of strongly compressed strings: get_lz78_grammar_string, get_lzw_grammar_string, get_aaaa
   * compression: LZW, LZ78 (as compression is written decompression is not necessary here)
   * decompression (for UNIX-compress): get_uncompress_string, get_compress_string

### Graph & results generation:
* Most results were run using time_test, some UNIX-compress specific things with agrep required run_all.sh (check versions here?)
* Graphs were generated via utils in graph/main.py
* Used test files are in test_files directory, so is tokenization

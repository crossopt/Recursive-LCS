#ifndef INC_LCS_H_
#define INC_LCS_H_

#include <iostream>
#include <string>
#include <algorithm>

#include "monge_matrix.h"

namespace LCS {
namespace kernel {

class LCSKernel {
public:
    const std::string &a;
    const std::string &b;
    matrix::SubpermutationMatrix kernel;
    matrix::MongeMatrix kernel_sum;

    LCSKernel(const std::string &a, const std::string &b);
private:
    matrix::SubpermutationMatrix calculate_kernel(unsigned a_l, unsigned a_r,
                                                  unsigned b_l, unsigned b_r);
};

}  // namespace kernel
}  // namespace LCS


#endif
#include <iostream>
#include <string>

#include "lcs_kernel.h"

int main() {
	std::cout << "Enter two strings in separate lines:" << std::endl;
	std::string a, b;
	std::getline(std::cin, a);
	std::getline(std::cin, b);
	auto kernel = LCS::kernel::RecursiveLCS(a, b);
	unsigned lcs_length = kernel.lcs_whole_a(0, b.size());
	std::cout << "The LCS length for the strings is " << lcs_length << std::endl;
	return 0;
}

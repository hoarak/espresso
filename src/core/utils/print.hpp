#ifndef CORE_UTILS_PRINT_HPP
#define CORE_UTILS_PRINT_HPP

#include <iostream>

namespace Utils {

void print() { std::cout << '\n'; }

/**
 * @brief Python style print function.
 */
template <typename Arg, typename... Args> void print(Arg v, Args... args) {
  std::cout << v << " ";
  print(args...);
}

} /* namespace Utils */

#endif

#ifndef CORE_UTILS_PRINT_HPP
#define CORE_UTILS_PRINT_HPP

#include <iostream>

namespace Utils {

/**
 * @brief Python style print function.
 */
template <typename Arg, typename... Args> void print(Arg v, Args... args) {
  std::cout << v << " ";
  print(args...);
}

template <typename T> void print(T v) { std::cout << v << '\n'; }

} /* namespace Utils */

#endif

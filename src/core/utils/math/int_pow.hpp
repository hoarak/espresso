#ifndef CORE_UTILS_INT_POW_HPP
#define CORE_UTILS_INT_POW_HPP

namespace Utils {
/**
 * @brief Calculate integer powers.
 *
 * This functions calculates x^n, where
 * n is a positive integer that is known
 * at compile time. It uses exponentiation by
 * squaring to construct an efficient function.
 */
template <unsigned n, typename T> inline T int_pow(T x) {
  switch (n) {
  case 0:
    return T(1);
  case 1:
    return x;
  default:
    /** Even branch */
    if (n % 2 == 0) {
      return int_pow<n / 2, T>(x * x);
    } else {
      return x * int_pow<(n - 1) / 2, T>(x * x);
    }
  }
}
}

#endif

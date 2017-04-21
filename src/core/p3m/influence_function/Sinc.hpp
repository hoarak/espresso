#ifndef CORE_P3M_INFLUENCE_FUNCTION_SINC_HPP
#define CORE_P3M_INFLUENCE_FUNCTION_SINC_HPP

#include <boost/math/constants/constants.hpp>

namespace P3M {
namespace InfluenceFunction {

/**
 * @brief sin(pi*x)/(pi*x) functor.
 *
 * Close to zero the Taylor expansion of sinc is used
 * to avoid numerical problems with very small divisors.
 */
template <typename T> class Sinc {
  static constexpr T eps = 0.1;
  static constexpr T c2 = -0.1666666666667e-0;
  static constexpr T c4 = 0.8333333333333e-2;
  static constexpr T c6 = -0.1984126984127e-3;
  static constexpr T c8 = 0.2755731922399e-5;

public:
  T operator()(T x) const {
    auto const pi_x = boost::math::constants::pi<T>() * x;

    if (std::abs(x) > eps) {
      return std::sin(pi_x) / pi_x;
    } else {
      auto const pi_x_2 = pi_x * pi_x;
      return 1.0 + pi_x_2 * (c2 + pi_x_2 * (c4 + pi_x_2 * (c6 + pi_x_2 * c8)));
    }
  }
};
}
}
#endif

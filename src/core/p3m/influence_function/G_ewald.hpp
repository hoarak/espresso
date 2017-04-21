#ifndef CORE_P3M_INFLUENCE_FUNCTION_G_EWALD_HPP
#define CORE_P3M_INFLUENCE_FUNCTION_G_EWALD_HPP

#include <boost/math/constants/constants.hpp>

namespace P3M {
namespace InfluenceFunction {

template <typename T> class G_ewald {
  const T f;
  static constexpr T limit{30};

  T sqr(T x) const { return x * x; }

public:
  explicit G_ewald(T alpha) : f(sqr(boost::math::constants::pi<T>() / alpha)) {}
  T operator()(T kx, T ky, T kz) const {
    auto const k2 = sqr(kx) + sqr(ky) + sqr(kz);
    auto const exponent = f * k2;

    return (exponent < limit) ? exp(-exponent) / k2 : 0.0;
  }
};
}
}

#endif

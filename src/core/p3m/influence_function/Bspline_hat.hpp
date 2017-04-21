#ifndef CORE_P3M_INFLUENCE_FUNCTION_BSPLINE_HAT_HPP
#define CORE_P3M_INFLUENCE_FUNCTION_BSPLINE_HAT_HPP

#include <cmath>

#include "utils/math/int_pow.hpp"
#include "Sinc.hpp"

namespace P3M {
namespace InfluenceFunction {
/**
 * @brief sinc(x)^p
 */
template <typename T, std::size_t order> class BSpline_hat {
  Sinc<T> sinc;

public:
  T operator()(T x) const { return Utils::int_pow<order>(sinc(x)); }
};
}
}

#endif

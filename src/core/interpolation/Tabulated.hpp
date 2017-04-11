#ifndef __INERPOLATION_INTERPOLATION_HPP
#define __INERPOLATION_INTERPOLATION_HPP

#include <array>
#include <cassert>

namespace Interpolation {

template <typename T> inline static int int_floor(const T x) {
  const int i = static_cast<int>(x); /* truncate */
  const int n = (x != static_cast<T>(i));
  const int g = (x < 0);
  return i - (n & g); /* i-1 if x<0 and x!=i */
}

template <unsigned n_interpolation, unsigned order, typename T = double>
class Tabulated {
public:
  template <typename W> explicit Tabulated(W w) {
    for (int i = -n_interpolation; i <= n_interpolation; i++) {
      const T x = i / ((2.0 * n_interpolation) + 1.0);
      for (unsigned j = 0; j < order; j++) {
        data[n_interpolation + i][j] = w(j, x);
      }
    }
  }

  /** Get interpolation weight for position x \in [-0.5, 0.5] for interpolation
   * point ip. */
  T operator()(const int ip, const T x) {
    const unsigned int ind = int_floor((x + 0.5) * 2.0 * n_interpolation);
    assert((ind < data.size()) && (ip < order));

    return data[ind][ip];
  }

private:
  std::array<std::array<T, order>, 2 * n_interpolation + 1> data;
};
}

#endif

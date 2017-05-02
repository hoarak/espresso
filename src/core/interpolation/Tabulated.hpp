#ifndef __INERPOLATION_INTERPOLATION_HPP
#define __INERPOLATION_INTERPOLATION_HPP

#include <array>
#include <cassert>
#include <vector>

namespace Interpolation {

template <typename T> inline static int int_floor(T x) {
  const int i = static_cast<int>(x); /* truncate */
  const int n = (x != (T)i);
  const int g = (x < 0);
  return i - (n & g); /* i-1 if x<0 and x!=i */
}

template <typename T, int n_interpolation, int cao> class Tabulated {
public:
  template <typename W> explicit Tabulated(W w) {

    for (int i = -n_interpolation; i <= n_interpolation; i++) {
      const double x = i / ((T(2) * n_interpolation) + T{1});
      for (int j = 0; j < cao; j++) {
        data[n_interpolation + i][j] = w(j, x);
      }
    }
  }

  /** Get interpolation weight for position x \in [-0.5, 0.5] for interpolation
   * point ip. */
  T operator()(const int ip, const T x) const {
    const unsigned int ind = int_floor((x + T(0.5)) * T(2) * n_interpolation);
    assert((ind < data.size()) && (ip < cao));

    return data[ind][ip];
  }

private:
  std::array<std::array<T, cao>, 2 * n_interpolation + 1> data;
};
}

#endif

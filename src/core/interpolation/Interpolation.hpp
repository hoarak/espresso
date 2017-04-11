#ifndef CORE_INTERPOLATION_INTERPOLATION_HPP
#define CORE_INTERPOLATION_INTERPOLATION_HPP

#include <array>
#include <cmath>

namespace Interpolation {

template <typename T, unsigned order, typename Weights> class Interpolation {
  Weights weights;
  /* Should be constexpr, but constexpr if requires c++14 */
  static const T pos_shift() {
    if (order % 2 == 0)
      return 0.5;
    else
      return 0.0;
  }

public:
  using index_t = std::array<int, 3>;

  explicit Interpolation(Weights weights = Weights())
      : weights(std::forward<Weights>(weights)) {}

  template <typename Particles, typename Offset, typename Kernel>
  void operator()(Particles const &parts, Offset const &offset, T h,
                  Kernel kernel) const {
    using array_t = std::array<T, 3>;

    const T hi = 1. / h;
    for (auto const &p : parts) {
      /* Distance to the nearest mesh point in units of h \in [-0.5, 0.5) */
      array_t dist;
      /* Index of the lower left corner of the assignment cube */
      index_t ll;

      auto const &pos = p.position();

      for (int dim = 0; dim < 3; dim++) {
        const T nmp_pos = (pos[dim] - offset[dim]) * hi + pos_shift();
        const T nmp_ind = static_cast<int>(std::floor(nmp_pos + 0.5));
        dist[dim] = nmp_pos - nmp_ind;
        ll[dim] = nmp_ind - order / 2;
      }

      /* Precompute the weights that are reused */
      std::array<T, order> inter_y;
      std::array<T, order> inter_z;

      for (int i = 0; i < order; i++) {
        inter_y[i] = weights(i, dist[1]);
        inter_z[i] = weights(i, dist[2]);
      }

      /* Loop the cube */
      index_t ind;
      for (int i = 0; i < order; i++) {
        ind[0] = ll[0] + i;
        const T wx = weights(i, dist[0]);
        for (int j = 0; j < order; j++) {
          ind[1] = ll[1] + j;
          const T wxy = wx * inter_y[j];
          for (int k = 0; k < order; k++) {
            ind[2] = ll[2] + k;
            kernel(p, ind, wxy * inter_z[k]);
          }
        }
      }
    }
  }
};

} /* namespace Interpolation */

#endif

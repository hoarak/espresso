#ifndef CORE_INTERPOLATION_CACHED_HPP
#define CORE_INTERPOLATION_CACHED_HPP

#include <array>
#include <cmath>
#include <vector>

namespace Interpolation {

template <typename T, int order> class Cached {
  /* Should be constexpr, but constexpr if requires c++14 */
  static const T pos_shift() {
    if (order % 2 == 0)
      return 0.5;
    else
      return 0.0;
  }

  std::vector<std::array<std::array<T, order>, 3>> m_weights;
  std::vector<std::array<int, 3>> m_ll;

public:
  using index_t = std::array<int, 3>;

  template <typename Particles, typename Kernel>
  void interpolate_from_cache(Particles const &parts,
                              Kernel const &kernel) const {
    auto weight_it = m_weights.cbegin();
    auto ll_it = m_ll.cbegin();

    for (auto const &p : parts) {
      auto const &w = *weight_it;
      auto const &ll = *ll_it;

      /* Pin the z weights, because they are reused */
      const std::array<T, order> wy = w[1];
      const std::array<T, order> wz = w[2];

      index_t ind;
      for (int i = 0; i < order; i++) {
        ind[0] = ll[0] + i;
        const T wx = w[0][i];
        for (int j = 0; j < order; j++) {
          ind[1] = ll[1] + j;
          const T wxy = wx * wy[j];
          for (int k = 0; k < order; k++) {
            ind[2] = ll[2] + k;
            kernel(p, ind, wxy * wz[k]);
          }
        }
      }

      ++weight_it;
      ++ll_it;
    }
  }

  template <typename Particles, typename Offset, typename Weights,
            typename Kernel>
  void operator()(Particles const &parts, Weights const &weights,
                  Offset const &offset, T h, Kernel const &kernel) {
    using array_t = std::array<T, 3>;

    m_weights.resize(parts.size());
    m_ll.resize(parts.size());
    auto weight_it = m_weights.begin();
    auto ll_it = m_ll.begin();

    const T hi = 1. / h;
    for (auto const &p : parts) {
      /* Distance to the nearest mesh point in units of h \in [-0.5, 0.5) */
      array_t dist;
      /* Index of the lower left corner of the assinment cube */
      index_t &ll = *ll_it;

      auto const &pos = p.position();

      for (int dim = 0; dim < 3; dim++) {
        const T nmp_pos = (pos[dim] - offset[dim]) * hi + pos_shift();
        const T nmp_ind = static_cast<int>(std::floor(nmp_pos + 0.5));
        dist[dim] = nmp_pos - nmp_ind;
        ll[dim] = nmp_ind - order / 2;
      }

      auto &w = *weight_it;

      for (int dim = 0; dim < 3; dim++) {
        for (int i = 0; i < order; i++) {
          w[dim][i] = weights(i, dist[dim]);
        }
      }

      index_t ind;
      for (int i = 0; i < order; i++) {
        ind[0] = ll[0] + i;
        const T wx = w[0][i];
        for (int j = 0; j < order; j++) {
          ind[1] = ll[1] + j;
          const T wxy = wx * w[1][j];
          for (int k = 0; k < order; k++) {
            ind[2] = ll[2] + k;
            kernel(p, ind, wxy * w[2][k]);
          }
        }
      }

      ++weight_it;
      ++ll_it;
    }
  };
};

} /* namespace Interpolation */

#endif

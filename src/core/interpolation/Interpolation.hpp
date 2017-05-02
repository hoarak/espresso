#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

#include <array>
#include <cmath>
#include <tuple>
#include <type_traits>

namespace Interpolation {

template <typename Scalar, int order> class Plain {
  /* Should be constexpr, but constexpr if requires c++14 */
  static const Scalar pos_shift() {
    if (order % 2 == 0)
      return 0.5;
    else
      return 0.0;
  }

public:
  using index_t = std::array<int, 3>;

  template <typename Particles, typename Offset, typename Weights,
            typename Kernel>
  void operator()(Particles const &parts, Weights const &weights,
                  Offset const &offset, Scalar h, Kernel const &kernel) const {
    using array_t = std::array<Scalar, 3>;

    const Scalar hi = 1. / h;
    for (auto const &p : parts) {
      /* Distance to the nearest mesh point in units of h \in [-0.5, 0.5) */
      array_t dist;
      /* Index of the lower left corner of the assinment cube */
      index_t ll;

      auto const &pos = p.position();

      for (int dim = 0; dim < 3; dim++) {
        const Scalar nmp_pos = (pos[dim] - offset[dim]) * hi + pos_shift();
        const Scalar nmp_ind = static_cast<int>(std::floor(nmp_pos + 0.5));

        dist[dim] = nmp_pos - nmp_ind;
        ll[dim] = nmp_ind - order / 2;
      }

      index_t ind;
      for (int i = 0; i < order; i++) {
        ind[0] = ll[0] + i;
        const Scalar wx = weights(i, dist[0]);
        for (int j = 0; j < order; j++) {
          ind[1] = ll[1] + j;
          const Scalar wxy = wx * weights(j, dist[1]);
          for (int k = 0; k < order; k++) {
            ind[2] = ll[2] + k;
            kernel(p, ind, wxy * weights(k, dist[2]));
          }
        }
      }
    }
  };
};

template <typename Scalar, int order> class CacheYZ {
  /* Should be constexpr, but constexpr if requires c++14 */
  static const Scalar pos_shift() {
    if (order % 2 == 0)
      return 0.5;
    else
      return 0.0;
  }

public:
  template <typename Particles, typename Offset, typename Weights,
            typename Kernel>
  void operator()(Particles const &parts, Weights const &weights,
                  Offset const &offset, Scalar h, Kernel const &kernel) const {
    using array_t = std::array<Scalar, 3>;
    using index_t = std::array<int, 3>;

    const Scalar hi = 1. / h;
    for (auto const &p : parts) {
      /* Distance to the nearest mesh point in units of h \in [-0.5, 0.5) */
      array_t dist;
      /* Index of the lower left corner of the assinment cube */
      index_t ll;

      auto const &pos = p.position();

      for (int dim = 0; dim < 3; dim++) {
        const Scalar nmp_pos = (pos[dim] - offset[dim]) * hi + pos_shift();
        const Scalar nmp_ind = static_cast<int>(std::floor(nmp_pos + 0.5));
        dist[dim] = nmp_pos - nmp_ind;
        ll[dim] = nmp_ind - order / 2;
      }

      std::array<Scalar, order> inter_y;
      std::array<Scalar, order> inter_z;

      for (int i = 0; i < order; i++) {
        inter_y[i] = weights(i, dist[1]);
        inter_z[i] = weights(i, dist[2]);
      }

      index_t ind;
      for (int i = 0; i < order; i++) {
        ind[0] = ll[0] + i;
        const Scalar wx = weights(i, dist[0]);
        for (int j = 0; j < order; j++) {
          ind[1] = ll[1] + j;
          const Scalar wxy = wx * inter_y[j];
          for (int k = 0; k < order; k++) {
            ind[2] = ll[2] + k;
            kernel(p, ind, wxy * inter_z[k]);
          }
        }
      }
    }
  };
};

template <typename Scalar, int order> class CacheXYZ {
  /* Should be constexpr, but constexpr if requires c++14 */
  static const Scalar pos_shift() {
    if (order % 2 == 0)
      return 0.5;
    else
      return 0.0;
  }

public:
  template <typename Particles, typename Offset, typename Weights,
            typename Kernel>
  void operator()(Particles const &parts, Weights const &weights,
                  Offset const &offset, Scalar h, Kernel const &kernel) const {
    using array_t = std::array<Scalar, 3>;
    using index_t = std::array<int, 3>;

    const Scalar hi = 1. / h;
    for (auto const &p : parts) {
      /* Distance to the nearest mesh point in units of h \in [-0.5, 0.5) */
      array_t dist;
      /* Index of the lower left corner of the assinment cube */
      index_t ll;

      auto const &pos = p.position();

      for (int dim = 0; dim < 3; dim++) {
        const Scalar nmp_pos = (pos[dim] - offset[dim]) * hi + pos_shift();
        const Scalar nmp_ind = static_cast<int>(std::floor(nmp_pos + 0.5));
        dist[dim] = nmp_pos - nmp_ind;
        ll[dim] = nmp_ind - order / 2;
      }

      alignas(16) std::array<Scalar, order> inter_x;
      alignas(16) std::array<Scalar, order> inter_y;
      alignas(16) std::array<Scalar, order> inter_z;

      for (int i = 0; i < order; i++) {
        inter_x[i] = weights(i, dist[0]);
        inter_y[i] = weights(i, dist[1]);
        inter_z[i] = weights(i, dist[2]);
      }

      index_t ind;
      for (int i = 0; i < order; i++) {
        ind[0] = ll[0] + i;
        const Scalar wx = inter_x[i];
        for (int j = 0; j < order; j++) {
          ind[1] = ll[1] + j;
          const Scalar wxy = wx * inter_y[j];
          for (int k = 0; k < order; k++) {
            ind[2] = ll[2] + k;
            kernel(p, ind, wxy * inter_z[k]);
          }
        }
      }
    }
  };
};

template <typename Scalar, int order> class Cached {
  /* Should be constexpr, but constexpr if requires c++14 */
  static const Scalar pos_shift() {
    if (order % 2 == 0)
      return 0.5;
    else
      return 0.0;
  }

  std::vector<std::array<std::array<Scalar, order>, 3>> m_weights;
  std::vector<std::array<int, 3>> m_ll;

public:
  template <typename Particles, typename Kernel>
  void interpolate_from_cache(Particles const &parts,
                              Kernel const &kernel) const {
    using array_t = std::array<Scalar, 3>;
    using index_t = std::array<int, 3>;

    auto weight_it = m_weights.cbegin();
    auto ll_it = m_ll.cbegin();

    for (auto const &p : parts) {
      auto const &w = *weight_it;
      auto const &ll = *ll_it;

      /* Pin the z weights, because they are reused */
      const std::array<Scalar, order> wy = w[1];
      const std::array<Scalar, order> wz = w[2];

      index_t ind;
      for (int i = 0; i < order; i++) {
        ind[0] = ll[0] + i;
        const Scalar wx = w[0][i];
        for (int j = 0; j < order; j++) {
          ind[1] = ll[1] + j;
          const Scalar wxy = wx * wy[j];
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
                  Offset const &offset, Scalar h, Kernel const &kernel) {
    using array_t = std::array<Scalar, 3>;
    using index_t = std::array<int, 3>;

    m_weights.resize(parts.size());
    m_ll.resize(parts.size());
    auto weight_it = m_weights.begin();
    auto ll_it = m_ll.begin();

    const Scalar hi = 1. / h;
    for (auto const &p : parts) {
      /* Distance to the nearest mesh point in units of h \in [-0.5, 0.5) */
      array_t dist;
      /* Index of the lower left corner of the assinment cube */
      index_t &ll = *ll_it;

      auto const &pos = p.position();

      for (int dim = 0; dim < 3; dim++) {
        const Scalar nmp_pos = (pos[dim] - offset[dim]) * hi + pos_shift();
        const Scalar nmp_ind = static_cast<int>(std::floor(nmp_pos + 0.5));
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
        const Scalar wx = w[0][i];
        for (int j = 0; j < order; j++) {
          ind[1] = ll[1] + j;
          const Scalar wxy = wx * w[1][j];
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

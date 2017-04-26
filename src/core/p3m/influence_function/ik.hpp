#ifndef CORE_P3M_INFLUENCE_FUNCTION_IK_HPP
#define CORE_P3M_INFLUENCE_FUNCTION_IK_HPP

#include <array>
#include <cmath>
#include <vector>

#include <boost/math/constants/constants.hpp>

#include "utils/math/sqr.hpp"

#include "AliasingSum.hpp"

namespace P3M {
namespace InfluenceFunction {

struct Transposed {
  constexpr static std::size_t RX = 0;
  constexpr static std::size_t RY = 1;
  constexpr static std::size_t RZ = 2;

  constexpr static std::size_t KX = 2;
  constexpr static std::size_t KY = 0;
  constexpr static std::size_t KZ = 1;
};

struct NotTransposed {
  constexpr static std::size_t RX = 0;
  constexpr static std::size_t RY = 1;
  constexpr static std::size_t RZ = 2;

  constexpr static std::size_t KX = 0;
  constexpr static std::size_t KY = 1;
  constexpr static std::size_t KZ = 2;
};

/**
 * @brief Implementation of doi: 10.1063/1.477414, Eq. (31).
 */
template <typename T, typename G_hat, typename W_hat, int m_max = 0,
          typename index_t = std::array<unsigned, 3>>
class IK {
  constexpr static T pi = boost::math::constants::pi<T>();
  constexpr static std::size_t RX = 0;
  constexpr static std::size_t RY = 1;
  constexpr static std::size_t RZ = 2;

  constexpr static std::size_t KX = 2;
  constexpr static std::size_t KY = 0;
  constexpr static std::size_t KZ = 1;

  std::vector<T> m_data;
  index_t m_mesh;
  std::array<std::vector<int>, 3> m_dop;
  std::array<T, 3> m_box;
  G_hat g_hat;

  AliasingSum<T, W_hat, m_max, index_t> m_aliasing_sum;

  /** Calculates the Fourier transformed differential operator.
   *  Remark: This is done on the level of n-vectors and not k-vectors,
   *           i.e. the prefactor i*2*PI/L is missing! */
  void calc_dop() {
    for (int i = 0; i < 3; i++) {
      auto &dop = m_dop[i];
      dop.resize(m_mesh[i]);

      dop[0] = 0;
      dop[m_mesh[i] / 2] = 0;
      for (int j = 1; j < m_mesh[i] / 2; j++) {
        dop[j] = j;
        dop[m_mesh[i] - j] = -j;
      }
    }
  }

  std::pair<std::array<T, 3>, T> aliasing_sums_force(index_t const &n) const {
    std::pair<std::array<T, 3>, T> ret{};
    using Utils::sqr;

    m_aliasing_sum(n, [&ret, this](int nmx, int nmy, int nmz, T w) {
      auto const f =
          w * g_hat(nmx / m_box[RX], nmy / m_box[RY], nmz / m_box[RZ]);

      ret.first[RX] += f * nmx / m_box[RX];
      ret.first[RY] += f * nmy / m_box[RY];
      ret.first[RZ] += f * nmz / m_box[RZ];

      ret.second += w;
    });

    return ret;
  }

  T aliasing_sums_energy(index_t const &n) const {
    T num{0}, denum{0};

    m_aliasing_sum(n, [&num, &denum, this](int nmx, int nmy, int nmz, T w) {
      num += w * g_hat(nmx / m_box[RX], nmy / m_box[RY], nmz / m_box[RZ]);
      denum += w;
    });

    return num / (denum * denum);
  }

public:
  IK(index_t mesh, std::array<T, 3> box, G_hat g_hat, W_hat w_hat)
      : g_hat(g_hat), m_mesh(mesh), m_box(box), m_aliasing_sum(mesh, w_hat) {
    calc_dop();
  }

  T energy(index_t const &n) const {
    if ((n[KX] % (m_mesh[RX] / 2) == 0) && (n[KY] % (m_mesh[RY] / 2) == 0) &&
        (n[KZ] % (m_mesh[RZ] / 2) == 0)) {
      return 0.0;
    } else {
      return aliasing_sums_energy(n) / pi;
    }
  }

  T force(index_t const &n) const {
    using Utils::sqr;

    if ((n[KX] % (m_mesh[RX] / 2) == 0) && (n[KY] % (m_mesh[RY] / 2) == 0) &&
        (n[KZ] % (m_mesh[RZ] / 2) == 0)) {
      return 0.0;
    } else {
      auto const as = aliasing_sums_force(n);

      auto const f1 = m_dop[RX][n[KX]] * as.first[RX] / m_box[RX] +
                      m_dop[RY][n[KY]] * as.first[RY] / m_box[RY] +
                      m_dop[RZ][n[KZ]] * as.first[RZ] / m_box[RZ];

      auto const f2 = sqr(m_dop[RX][n[KX]] / m_box[RX]) +
                      sqr(m_dop[RY][n[KY]] / m_box[RY]) +
                      sqr(m_dop[RZ][n[KZ]] / m_box[RZ]);

      auto const f3 = f1 / (f2 * sqr(as.second));

      return 2 * f3 / pi;
    }
  }
};
} /* namespace InfluenceFunction */
} /* namespace P3M */

#endif

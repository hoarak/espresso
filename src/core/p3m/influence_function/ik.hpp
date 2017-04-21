#ifndef CORE_P3M_INFLUENCE_FUNCTION_IK_HPP
#define CORE_P3M_INFLUENCE_FUNCTION_IK_HPP

#include <array>
#include <cmath>
#include <vector>

#include <boost/math/constants/constants.hpp>

#include "utils/math/sqr.hpp"

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
  std::array<std::vector<int>, 3> m_shift;
  std::array<std::vector<int>, 3> m_dop;
  std::array<T, 3> m_box;
  G_hat g_hat;
  W_hat w_hat;

  void calc_shift() {
    for (int i = 0; i < 3; i++) {
      auto &shift = m_shift[i];
      shift.resize(m_mesh[i]);

      shift[0] = 0;
      for (int j = 1; j <= m_mesh[i] / 2; j++) {
        shift[j] = j;
        shift[m_mesh[i] - j] = -j;
      }
    }
  }

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

  std::pair<std::array<T, 3>, T> aliasing_sums(index_t const &n) const {
    std::pair<std::array<T, 3>, T> ret{};
    using Utils::sqr;

    for (int mx = -m_max; mx <= m_max; ++mx) {
      const int nmx = m_shift[RX][n[KX]] + m_mesh[RX] * mx;
      auto const sx = sqr(w_hat(nmx / static_cast<T>(m_mesh[RX])));
      for (int my = -m_max; my <= m_max; ++my) {
        int const nmy = m_shift[RY][n[KY]] + m_mesh[RY] * my;
        auto const sxy = sx * sqr(w_hat(nmy / static_cast<T>(m_mesh[RY])));
        for (int mz = -m_max; mz <= m_max; ++mz) {
          int const nmz = m_shift[RZ][n[KZ]] + m_mesh[RZ] * mz;
          auto const sxyz = sxy * sqr(w_hat(nmz / static_cast<T>(m_mesh[RZ])));

          auto const f = sxyz * g_hat(nmx / m_box[RX], nmy / m_box[RY],
                                              nmz / m_box[RZ]);

          ret.first[RX] += f * nmx / m_box[RX];
          ret.first[RY] += f * nmy / m_box[RY];
          ret.first[RZ] += f * nmz / m_box[RZ];

          ret.second += sxyz;
        }
      }
    }

    return ret;
  }

public:
  IK(index_t mesh, std::array<T, 3> box, G_hat g_hat, W_hat w_hat)
      : g_hat(g_hat), w_hat(w_hat), m_mesh(mesh), m_box(box) {
    calc_shift();
    calc_dop();
  }

  T operator()(index_t const &n) const {
    using Utils::sqr;

    if ((n[KX] % (m_mesh[RX] / 2) == 0) && (n[KY] % (m_mesh[RY] / 2) == 0) &&
        (n[KZ] % (m_mesh[RZ] / 2) == 0)) {
      return 0.0;
    } else {
      auto const as = aliasing_sums(n);

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

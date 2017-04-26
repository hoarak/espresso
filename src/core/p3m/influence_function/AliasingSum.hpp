#ifndef CORE_P3M_INFLUENCE_ALIASING_SUM_HPP
#define CORE_P3M_INFLUENCE_ALIASING_SUM_HPP

#include <array>

#include <boost/math/constants/constants.hpp>

#include "utils/math/sqr.hpp"

namespace P3M {
namespace InfluenceFunction {

template <typename T, typename W_hat, int m_max = 0,
          typename index_t = std::array<unsigned, 3>>
class AliasingSum {
  constexpr static T pi = boost::math::constants::pi<T>();
  constexpr static unsigned RX = 0;
  constexpr static unsigned RY = 1;
  constexpr static unsigned RZ = 2;

  constexpr static unsigned KX = 2;
  constexpr static unsigned KY = 0;
  constexpr static unsigned KZ = 1;

  index_t m_mesh;
  std::array<std::vector<int>, 3> m_shift;

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

public:
  explicit AliasingSum(index_t mesh, W_hat w_hat) : m_mesh(mesh), w_hat(w_hat) {
    calc_shift();
  }

  template <typename Kernel>
  void operator()(index_t const &n, Kernel kernel) const {
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

          kernel(nmx, nmy, nmz, sxyz);
        }
      }
    }
  }
};
}
}

#endif

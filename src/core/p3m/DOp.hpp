#ifndef CORE_P3M_DOP_HPP
#define CORE_P3M_DOP_HPP

#include <array>
#include <vector>

namespace P3M {

class DOp {
  std::array<std::vector<int>, 3> m_dop;

public:
  DOp() = default;
  template <typename Mesh> explicit DOp(Mesh mesh) {
    for (int i = 0; i < 3; i++) {
      auto &dop = m_dop[i];
      dop.resize(mesh[i]);

      dop[0] = 0;
      dop[mesh[i] / 2] = 0;
      for (int j = 1; j < mesh[i] / 2; j++) {
        dop[j] = j;
        dop[mesh[i] - j] = -j;
      }
    }
  }

  std::vector<int> const &operator[](int dim) const { return m_dop[dim]; }
  int operator()(int dim, int n) const { return m_dop[dim][n]; }
};
}

#endif

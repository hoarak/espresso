#ifndef CORE_UNIT_TESTS_MOCK_PARTICLE_HPP
#define CORE_UNIT_TESTS_MOCK_PARTICLE_HPP

#include <array>
#include <limits>

namespace Testing {
class Particle {
  unsigned m_id;
  std::array<double, 3> m_position;

public:
  Particle()
      : m_id(-1), m_position({{std::numeric_limits<double>::infinity(),
                               std::numeric_limits<double>::infinity(),
                               std::numeric_limits<double>::infinity()}}) {}
  explicit Particle(int id) : Particle() {}
  Particle(int id, std::array<double, 3> position)
      : m_id(id), m_position(position) {}

  int identity() const { return m_id; }
  std::array<double, 3> const &position() const { return m_position; }
};
}

#endif

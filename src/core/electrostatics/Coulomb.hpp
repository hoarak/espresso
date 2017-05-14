#ifndef CORE_ELECTROSTATICS_COULOMB_HPP
#define CORE_ELECTROSTATICS_COULOMB_HPP

#include <array>
#include <stdexcept>

template <typename ParticleRange, typename T = double> class Coulomb {
public:
  virtual T energy(ParticleRange const &) {
    throw std::runtime_error("Energy calculation not supported.");
  }

  virtual void add_forces(ParticleRange &) {
    throw std::runtime_error("Force calculation not supported.");
  }

  virtual std::array<T, 9> stress(ParticleRange const &) {
    std::runtime_error("Force calculation not supported.");
  }
};

#endif

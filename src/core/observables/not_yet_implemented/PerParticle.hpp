#ifndef OBSERABLES_PER_PARTICLE_HPP
#define OBSERABLES_PER_PARTICLE_HPP

#include "Observable.hpp"

namespace Observables {

template<class Kernel>
class PerParticle : public Observable {
  PerParticle(Kernel const& kernel) : m_kernel(kernel) {}

  int actual_calculate() override {

  }

  private:
    Kernel const& m_kernel;
};

} /* namespace Observables */

#endif

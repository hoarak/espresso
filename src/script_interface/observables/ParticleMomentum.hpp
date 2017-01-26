#ifndef SCRIPT_INTERFACE_OBSERVABLES_PARTICLE_MOMENTUM_HPP
#define SCRIPT_INTERFACE_OBSERVABLES_PARTICLE_MOMENTUM_HPP

#include "PidObservable.hpp"

namespace ScriptInterface {
namespace Observables {
class ParticleMomenta
    : public PidObservable<::Observables::ParticleProperties::Momentum> {
public:
  using PidObservable<::Observables::ParticleProperties::Momentum>::observable;
  const std::string name() const override { return "Observables::Momentum"; }
};
}
}
#endif

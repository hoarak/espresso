#ifndef SCRIPT_INTERFACE_OBSERVABLES_PARTICLE_POSITIONS_HPP
#define SCRIPT_INTERFACE_OBSERVABLES_PARTICLE_POSITIONS_HPP

#include "PidObservable.hpp"

namespace ScriptInterface {
namespace Observables {
class ParticlePositions
    : public PidObservable<::Observables::ParticleProperties::Position> {
public:
  using PidObservable<::Observables::ParticleProperties::Position>::observable;
  const std::string name() const override {
    return "Observables::ParticlePositions";
  }
};
}
}
#endif

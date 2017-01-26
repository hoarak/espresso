#ifndef OBSERVABLES_PARTICLE_PROPERTIES_HPP
#define OBSERVABLES_PARTICLE_PROPERTIES_HPP

#include "core/integrate.hpp"
#include "core/particle_data.hpp"
#include "core/rotation.hpp"

#include "operators.hpp"
#include "types.hpp"

namespace Observables {
namespace ParticleProperties {

struct Constant : public Scalar {
  Constant(double val) : m_value(val) {}

  double operator()(size_t i) const { return m_value; }

private:
  const double m_value;
};

template <size_t numerator, size_t denominator>
struct RationalConstant : Scalar {
  template <typename... Args> RationalConstant(Args...) {}

  double operator()(size_t) const {
    return static_cast<double>(numerator) / denominator;
  }
};

struct Property {
  Property(Particle const &p) : m_p(p) {}

protected:
  Particle const &m_p;
};

struct Mass : Scalar, Property {
  using Property::Property;
  double operator()(size_t) const { return m_p.p.mass; }
};

struct Charge : Scalar, Property {
  using Property::Property;
  double operator()(size_t) const { return m_p.p.q; }
};

struct Position : Vector<3>, Property {
  using Property::Property;
  double operator()(size_t i) const { return m_p.r.p[i]; }
};

struct Velocity : Vector<3>, Property {
  using Property::Property;
  double operator()(size_t i) const { return m_p.m.v[i] / time_step; }
};

struct Omega : Vector<3>, Property {
  using Property::Property;
  double operator()(size_t i) const { return m_p.m.omega[i]; }
};

struct Force : Vector<3>, Property {
  using Property::Property;
  double operator()(size_t i) const {
    return m_p.f.f[i] / time_step / time_step * 2.;
  }
};

struct RotationMatrix : Matrix<3, 3> {
  RotationMatrix(Particle const &p) {
    define_rotation_matrix(&p, m_matrix.data());
  }

  double operator()(size_t i) const { return m_matrix[i]; }

private:
  std::array<double, 9> m_matrix;
};

template <typename T> using Flux = OuterProduct<T, Velocity>;

using Current = Flux<Charge>;
using Momentum = Flux<Mass>;
using StressTensor = Flux<Momentum>;
using Power = ScalarProduct<Velocity, Force>;
using AngularMomentum = MatrixProduct<RotationMatrix, Omega>;
using KineticEnergy = ComponentwiseProduct<RationalConstant<1, 2>,
                                           ScalarProduct<Velocity, Velocity>>;
}
}

#endif

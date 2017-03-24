#ifndef OBSERVABLES_PIDOBSERVABLE_HPP
#define OBSERVABLES_PIDOBSERVABLE_HPP

#include "Observable.hpp"
#include "core/particle_data.hpp"

namespace Observables {
template <class Property> class ParticleProperty : public Observable {
public:
  Tensor calculate() const override {
    std::array<size_t, Property::rank() + 1> extent;

    std::copy(Property::extent().begin(), Property::extent().end(),
              extent.begin() + 1);
    extent[0] = 0;

    Tensor ret(extent);

    updatePartCfg(0);

    auto out = ret.begin();
    for (int i = 0; i < n_part; i++) {
      auto extent = ret.extents();
      extent.front()++;
      ret.resize(extent);

      auto const property = Property(partCfg[i]);
      for (size_t j = 0; j < Property::size(); ++j) {
        *out++ = property(j);
      }
    }

  return ret;
}
};
}

#endif

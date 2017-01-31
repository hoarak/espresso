#ifndef OBSERVABLES_PIDOBSERVABLE_HPP
#define OBSERVABLES_PIDOBSERVABLE_HPP

#include "Observable.hpp"
#include "core/particle_data.hpp"
#include "particle_properties.hpp"

namespace Observables {

// Observable which acts on a given list of particle ids
template <class ParticleProperty> class PidObservable : public Observable {
public:
  PidObservable() {
    auto const property_extent = ParticleProperty::extent();
    m_extent.resize(property_extent.size() + 1, 0);
    std::copy(property_extent.begin(), property_extent.end(),
              m_extent.begin() + 1);
  }

  template <typename Container> void set_ids(Container const &ids) {
    m_ids.clear();
    std::copy(std::begin(ids), std::end(ids), std::back_inserter(m_ids));
    m_extent[0] = m_ids.size();
  }

  std::vector<int> const &ids() const { return m_ids; }

  Tensor calculate() const override {
    Tensor ret(m_extent);

    if (!sortPartCfg()) {
      throw std::runtime_error("could not sort partCfg");
    }

    auto out = ret.begin();
    for (auto const &i : m_ids) {
      auto const property = ParticleProperty(partCfg[i]);
      for (size_t j = 0; j < ParticleProperty::size(); ++j) {
        *out++ = property(j);
      }
    }

    return ret;
  }

  std::vector<size_type> extent() const override { return m_extent; }

private:
  std::vector<int> m_ids;
  std::vector<size_t> m_extent;
};
}

#endif

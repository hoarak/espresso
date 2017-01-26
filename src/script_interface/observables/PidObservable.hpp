/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
  Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef SCRIPT_INTERFACE_OBSERVABLES_PIDOBSERVABLE_HPP
#define SCRIPT_INTERFACE_OBSERVABLES_PIDOBSERVABLE_HPP

#include "ScriptInterface.hpp"

#include "Observable.hpp"
#include "core/observables/PidObservable.hpp"

namespace ScriptInterface {
namespace Observables {

template <typename ParticleProperty> class PidObservable : public Observable {
  using obs_t = ::Observables::PidObservable<ParticleProperty>;

public:
  PidObservable() : m_observable(new obs_t) {}

  VariantMap get_parameters() const override {
    auto params = Observable::get_parameters();
    params["ids"] = m_observable->ids();

    return params;
  }

  ParameterMap valid_parameters() const override {
    auto params = Observable::valid_parameters();
    params["ids"] = {ParameterType::INT_VECTOR, true};
    return params;
  }

  void set_parameter(std::string const &name, Variant const &value) override {
    if (name == "ids") {
      m_observable->set_ids(get_value<std::vector<int>>(value));
      return;
    }
    Observable::set_parameter(name, value);
  }

  virtual std::shared_ptr<::Observables::Observable>
  observable() const override {
    return std::static_pointer_cast<::Observables::Observable>(m_observable);
  }

private:
  std::shared_ptr<obs_t> m_observable;
};

} /* namespace Observables */
} /* namespace ScriptInterface */

#endif

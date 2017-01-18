/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project

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
#ifndef OBSERVABLES_OBSERVABLE_HPP
#define OBSERVABLES_OBSERVABLE_HPP

#include "config.hpp"
#include <stdexcept>
#include <vector>

namespace Observables {

class Observable {
public:
  Observable();
  int calculate();

  virtual int actual_update(){};

  virtual int n_values() const { return 0; };
  std::vector<double> last_value;

protected:
  double last_update;

private:
  virtual int actual_calculate() {
    throw std::runtime_error(
        "Observable did not override actual_caucluate()\n");
  };
};

} // Namespace Observables
#endif

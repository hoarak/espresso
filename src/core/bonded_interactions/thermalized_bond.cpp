/*
  Copyright (C) 2010-2018 The ESPResSo project
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
/** \file
 *
 *  Implementation of \ref thermalized_bond.hpp
 */

#include "thermalized_bond.hpp"
#include "bonded_interaction_data.hpp"
#include "communication.hpp"
#include "global.hpp"

#include <utils/constants.hpp>

int n_thermalized_bonds = 0;

int thermalized_bond_set_params(int bond_type, double temp, double gamma, double r_cut) {
  if (bond_type < 0)
    return ES_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].p.thermalized_bond.temp = temp;
  bonded_ia_params[bond_type].p.thermalized_bond.gamma =gamma;
  bonded_ia_params[bond_type].p.thermalized_bond.r_cut = r_cut;
  bonded_ia_params[bond_type].p.thermalized_bond.pref1 = gamma;
  bonded_ia_params[bond_type].p.thermalized_bond.pref2 =
      sqrt(24.0 * gamma / time_step * temp);
  bonded_ia_params[bond_type].type = BONDED_IA_THERMALIZED_DIST;
  bonded_ia_params[bond_type].num = 1;

  n_thermalized_bonds += 1;
  mpi_bcast_ia_params(bond_type, -1);
  mpi_bcast_parameter(FIELD_THERMALIZEDBONDS);

  return ES_OK;
}

void thermalized_bond_heat_up() {
  double pref_scale = sqrt(3);
  thermalized_bond_update_params(pref_scale);
}

void thermalized_bond_cool_down() {
  double pref_scale = 1.0 / sqrt(3);
  thermalized_bond_update_params(pref_scale);
}

void thermalized_bond_init() {

  for (auto &bonded_ia_param : bonded_ia_params) {
    if (bonded_ia_param.type == BONDED_IA_THERMALIZED_DIST) {
      Thermalized_bond_parameters &t = bonded_ia_param.p.thermalized_bond;
      t.pref1_dist = t.gamma;
      t.pref2 =
          sqrt(24.0 * t.gamma / time_step * t.temp);
    }
  }
}

void thermalized_bond_update_params(double pref_scale) {
  for (auto &bonded_ia_param : bonded_ia_params) {
    if (bonded_ia_param.type == BONDED_IA_THERMALIZED_DIST) {
      Thermalized_bond_parameters &t = bonded_ia_param.p.thermalized_bond;
      t.pref2 *= pref_scale;
    }
  }
}

/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
#ifndef _THERMOSTAT_H
#define _THERMOSTAT_H
/** \file thermostat.hpp 

*/

#include "config.hpp"

#include <cmath>
#include "debug.hpp"
#include "particle_data.hpp"
#include "random.hpp"
#include "global.hpp"
#include "integrate.hpp"
#include "cells.hpp"
#include "lb.hpp"
#include "dpd.hpp"
#include "virtual_sites.hpp"

/** \name Thermostat switches*/
/************************************************************/
/*@{*/

#define THERMO_OFF        0
#define THERMO_LANGEVIN   1
#define THERMO_DPD        2
#define THERMO_NPT_ISO    4
#define THERMO_LB         8
#define THERMO_INTER_DPD  16
#define THERMO_GHMC       32
#define THERMO_CPU        64
/*@}*/

namespace Thermostat {
  auto noise = []() { return (d_random() - 0.5); };
}

/************************************************
 * exported variables
 ************************************************/

/** Switch determining which thermostat to use. This is a or'd value
    of the different possible thermostats (defines: \ref THERMO_OFF,
    \ref THERMO_LANGEVIN, \ref THERMO_DPD \ref THERMO_NPT_ISO). If it
    is zero all thermostats are switched off and the temperature is
    set to zero.  */
extern int thermo_switch;

/** temperature. */
extern double temperature;

#ifndef PARTICLE_ANISOTROPY
/** Langevin friction coefficient gamma. */
extern double langevin_gamma;
#else
extern double langevin_gamma[3];
#endif

#ifndef ROTATIONAL_INERTIA
/** Langevin friction coefficient gamma. */
extern double langevin_gamma_rotation;
#else
extern double langevin_gamma_rotation[3];
#endif

/** Langevin for translations */
extern bool langevin_trans;

/** Langevin for rotations */
extern bool langevin_rotate;

/** Friction coefficient for nptiso-thermostat's inline-function friction_therm0_nptiso */
extern double nptiso_gamma0;
/** Friction coefficient for nptiso-thermostat's inline-function friction_thermV_nptiso */
extern double nptiso_gammav;

/** Number of NVE-MD steps in GHMC Cycle*/
extern int ghmc_nmd;
/** Phi parameter for GHMC partial momenum update step */
extern double ghmc_phi;

/************************************************
 * functions
 ************************************************/


/** initialize constants of the thermostat on
    start of integration */
void thermo_init();

/** very nasty: if we recalculate force when leaving/reentering the integrator,
    a(t) and a((t-dt)+dt) are NOT equal in the vv algorithm. The random
    numbers are drawn twice, resulting in a different variance of the random force.
    This is corrected by additional heat when restarting the integrator here.
    Currently only works for the Langevin thermostat, although probably also others
    are affected.
*/
void thermo_heat_up();

/** pendant to \ref thermo_heat_up */
void thermo_cool_down();
/** Get current temperature for CPU thermostat */
int get_cpu_temp();

/** Start the CPU thermostat */
void set_cpu_temp(int temp);

#endif

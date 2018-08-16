//
// Created by florian on 17.08.18.
//

#include "coulomb.hpp"

#include "debye_hueckel.hpp"
#include "reaction_field.hpp"

Coulomb_parameters coulomb = {
        0.0, COULOMB_NONE,
        0.0, DIPOLAR_NONE,
};

Debye_hueckel_params dh_params = {0.0, 0.0};
Reaction_field_params rf_params = {0.0, 0.0};

/** Induced field (for const. potential feature) **/
double field_induced;
/** Applied field (for const. potential feature) **/
double field_applied;
//
// Created by florian on 17.08.18.
//

#include "coulomb.hpp"

Coulomb_parameters coulomb = {
#ifdef ELECTROSTATICS
        0.0, COULOMB_NONE,
#endif
#ifdef DIPOLES
        0.0, DIPOLAR_NONE,
#endif
};

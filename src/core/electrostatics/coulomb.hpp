//
// Created by florian on 17.08.18.
//

#ifndef ESPRESSO_COULOMB_HPP
#define ESPRESSO_COULOMB_HPP


#include <core/nonbonded_interactions/cos2.hpp>
#include "../../../cmake-build-debug/myconfig.hpp"

enum CoulombMethod {
  COULOMB_NONE,      //< Coulomb interaction switched off (NONE)
  COULOMB_DH,        //< Coulomb method is Debye-Hueckel
  COULOMB_P3M,       //< Coulomb method is P3M
  COULOMB_MMM1D,     //< Coulomb method is one-dimensional MMM
  COULOMB_MMM2D,     //< Coulomb method is two-dimensional MMM
  COULOMB_MAGGS,     //< Coulomb method is "Maggs"
  COULOMB_ELC_P3M,   //< Coulomb method is P3M plus ELC
  COULOMB_RF,        //< Coulomb method is Reaction-Field
  COULOMB_P3M_GPU,   //< Coulomb method is P3M with GPU based long range part
                     // calculation
  COULOMB_MMM1D_GPU, //< Coulomb method is one-dimensional MMM running on GPU
  COULOMB_EK,        //< Coulomb method is electrokinetics
  COULOMB_SCAFACOS,  //< Coulomb method is scafacos
};
/** field containing the interaction parameters for
 *  the coulomb  interaction.  */
struct Coulomb_parameters {

#ifdef ELECTROSTATICS
  /** bjerrum length times temperature. */
  double prefactor;

  /** Method to treat coulomb interaction. */
  CoulombMethod method;
#endif

#ifdef DIPOLES
  double Dprefactor;
  DipolarInteraction Dmethod;
#endif

};

/** Structure containing the coulomb parameters. */
extern Coulomb_parameters coulomb;

#endif //ESPRESSO_COULOMB_HPP

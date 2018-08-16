//
// Created by florian on 17.08.18.
//

#ifndef ESPRESSO_COULOMB_HPP
#define ESPRESSO_COULOMB_HPP

#include "config.hpp"

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

/** \name Type codes for the type of dipolar interaction
  Enumeration of implemented methods for the magnetostatic
  interaction.
 */
/************************************************************/
/*@{*/

enum DipolarInteraction {
    /** dipolar interation switched off (NONE). */
            DIPOLAR_NONE = 0,
    /** dipolar method is P3M. */
            DIPOLAR_P3M,
    /** Dipolar method is P3M plus DLC. */
            DIPOLAR_MDLC_P3M,
    /** Dipolar method is all with all and no replicas */
            DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA,
    /** Dipolar method is magnetic dipolar direct sum */
            DIPOLAR_DS,
    /** Dipolar method is direct sum plus DLC. */
            DIPOLAR_MDLC_DS,
    /** Direct summation on gpu */
            DIPOLAR_DS_GPU,
#ifdef DIPOLAR_BARNES_HUT
    /** Direct summation on gpu by Barnes-Hut algorithm */
  DIPOLAR_BH_GPU,
#endif
    /** Scafacos library */
            DIPOLAR_SCAFACOS
};

/** field containing the interaction parameters for
 *  the coulomb  interaction.  */
struct Coulomb_parameters {
  /** bjerrum length times temperature. */
  double prefactor;

  /** Method to treat coulomb interaction. */
  CoulombMethod method;

  double Dprefactor;
  DipolarInteraction Dmethod;
};

/** Structure containing the coulomb parameters. */
extern Coulomb_parameters coulomb;

/** Induced field (for const. potential feature). **/
extern double field_induced;
/** Applied field (for const. potential feature) **/
extern double field_applied;

#endif //ESPRESSO_COULOMB_HPP

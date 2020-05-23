/*

 $Id: eval.h,v 1.18 2007/04/27 06:01:48 garrett Exp $

 AutoDock

 Copyright (C) 1989-2007,  Scott Halliday, Rik Belew, Garrett M. Morris, David S. Goodsell, Ruth Huey, Arthur J. Olson,
 All Rights Reserved.

 AutoDock is a Trade Mark of The Scripps Research Institute.

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

 */

/********************************************************************
    The header file for the eval class

                                rsh 09/06/95
********************************************************************/
#ifndef _EVAL_H
#define _EVAL_H

#include <stdio.h>
#include "structs.h"
#include "rep.h"
#include "qmultiply.h"
#include "cnv_state_to_coords.h"
#include "trilinterp.h"
#include "eintcal.h"
#include "energy.h"

#ifdef CUDA_READY
typedef float (*evdWHbtype)[ATOM_MAPS][ATOM_MAPS];
typedef unsigned short (*UStorProfilestype)[MAX_TORS][NTORDIVS];
typedef float (*torstype)[MAX_TORS];
typedef Boole (*istorconstrainedtype)[MAX_TORS];
typedef float (*qspabschargestype)[MAX_ATOMS];
typedef int (*ignore_interstype)[MAX_ATOMS];
typedef float (*mapstype)[MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS];
typedef int (*typestype)[MAX_ATOMS];
typedef float (*ABSchargestype)[MAX_ATOMS];
typedef float (*chargestype)[MAX_ATOMS];
typedef float (*crdstype)[MAX_ATOMS][SPACE];
typedef int (*tlisttype)[MAX_ATOMS];
typedef Real (*crdtype)[SPACE];
typedef Real (*vttype)[SPACE];
typedef Real (*crdpdbtype)[SPACE];
typedef Real (*maptype)[MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS];
typedef unsigned short (*US_torProfiletype)[NTORDIVS];
#endif

#ifdef DEBUG
extern FILE *logFile;
#endif

#if defined(USING_COLINY)
void make_state_from_rep(double *x, int n, State *now);
#endif

// Variables required for energy calculations
//added extern to variables. Everton Mendonça 03/04/2016
/*extern int *B_outsidesgpu;
extern float *energiesgpu;
extern unsigned int *evalflagsgpu;
extern int cpunatoms;
extern int nBlocks;
extern int blocksize;
extern unsigned int num_individuals;
extern float *nb_group_energycpu;

extern float *float_arraygpu;
extern int *int_arraygpu;

//tril params
extern float *crdsgpu;
extern float *chargesgpu;
extern float *ABSchargesgpu;
extern int   *typesgpu;
extern float *gridinfosgpu;
extern int   *ignore_intersgpu;
extern float *cudamap;

//eint cal params
extern float *nonbondlistsgpu;
extern Boole *incelecgpu;
extern Boole inc14interactgpu;
extern float *scale14sgpu;
extern Boole *usenonbondcutsgpu;
extern Boole haveflexresiduesgpu;
extern float *unboundinternalFEsgpu;
extern float *evdWHbgpu;
extern float *solfngpu;
extern float *epsilonfngpu;
extern float *repsilonfngpu;
extern int *nnb_arraygpu;
extern float *nb_group_energygpu;

#define BLOCK_SIZE 128*/

// Global variables required for energy calculations
extern int ElecMap;
extern int DesolvMap;
extern Real nb_group_energy[3];
extern int Nnb_array[3];

extern void CHECK_ERROR(int num);

void make_state_from_rep(Representation **rep, State *stateNow);

class Eval
{
   private:
      UnsignedFourByteLong num_evals;
      int natom, Nnb;
      GridMapSetInfo *info;
      Real eval_elec[MAX_ATOMS]; // gmm added 21-Jan-1998, for writePDBQState
      Real eval_emap[MAX_ATOMS]; // gmm added 21-Jan-1998, for writePDBQState
      Boole B_calcIntElec, B_isGaussTorCon, B_ShowTorE;
      State stateNow;
      unsigned short *US_TorE, (*US_torProfile)[NTORDIVS];
      int *type, (*tlist)[MAX_ATOMS];
      NonbondParam *nonbondlist;
      Real *charge, *abs_charge, *qsp_abs_charge;
      Real (*crd)[SPACE], (*vt)[SPACE], (*crdpdb)[SPACE], (*crdreo)[SPACE];
      EnergyTables *ptr_ad_energy_tables;
      Real (*map)[MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS];
      Boole *B_isTorConstrained;
      Molecule mol;
      int ignore_inter[MAX_ATOMS]; // gmm 2002-05-21, for CA, CB in flexible sidechains
      Boole         B_include_1_4_interactions; // gmm 2005-01-8, for scaling 1-4 nonbonds
      Real scale_1_4;                  // gmm 2005-01-8, for scaling 1-4 nonbonds
      ParameterEntry *parameterArray;
      Real  unbound_internal_FE;
      Boole B_compute_intermol_energy; // use for computing unbound state
      Boole B_use_non_bond_cutoff;  // set this to FALSE if we are computing unbound extended conformations
      Boole B_have_flexible_residues;

   public:
      Eval(void);
      void setup( Real init_crd[MAX_ATOMS][SPACE],
                  Real  init_charge[MAX_ATOMS],
                  Real  init_abs_charge[MAX_ATOMS],
                  Real  init_qsp_abs_charge[MAX_ATOMS],
                  int            init_type[MAX_ATOMS], int init_natom,
                  Real  init_map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],

                  Real  init_elec[MAX_ATOMS], // gmm added 21-Jan-1998, for writePDBQState
                  Real  init_emap[MAX_ATOMS], // gmm added 21-Jan-1998, for writePDBQState

                  NonbondParam *init_nonbondlist,
                  EnergyTables   *init_ptr_ad_energy_tables,
                  int init_Nnb,
                  Boole          init_B_calcIntElec,
                  Boole          init_B_isGaussTorCon, Boole init_B_isTorConstrained[MAX_TORS],
                  Boole          init_B_ShowTorE, unsigned short init_US_TorE[MAX_TORS],
                  unsigned short init_US_torProfile[MAX_TORS][NTORDIVS],
                  Real  init_vt[MAX_TORS][SPACE], int init_tlist[MAX_TORS][MAX_ATOMS],
                  Real  init_crdpdb[MAX_ATOMS][SPACE],
                  Real  init_crdreo[MAX_ATOMS][SPACE],
                  State stateInit, Molecule molInit,
                  int            init_ignore_inter[MAX_ATOMS],
                  Boole          init_B_include_1_4_interactions, // gmm 2005-01-8, for scaling 1-4 nonbonds
                  Real  init_scale_1_4,                   // gmm 2005-01-8, for scaling 1-4 nonbonds
                  ParameterEntry init_parameterArray[MAX_MAPS],
                  Real  init_unbound_internal_FE,
                  GridMapSetInfo *init_info,
                  Boole  init_B_use_non_bond_cutoff,  // set this to FALSE if we are computing unbound extended conformations
                  Boole  init_B_have_flexible_residues
                  );
      void update_crds( Real init_crdreo[MAX_ATOMS][SPACE],
                        Real init_vt[MAX_TORS][SPACE] );

      double operator()(Representation **);
      double operator()(Representation **, int); // GMM - allows calculation of a particular term of the total energy
#if defined(USING_COLINY)
      double operator()(double*, int);
#endif
      double eval();    // WEH - a basic change that facilitates the use of Coliny
      double eval(int); // GMM - allows calculation of a particular term of the total energy
      UnsignedFourByteLong evals(void);
      void reset(void);
      int write(FILE *out_file, Representation **rep);
      void compute_intermol_energy(Boole init_B_compute_intermol_energy); // for computing unbound state
#ifdef CUDA_READY
      void evalPrep(
        unsigned int, //pop_size
        double *,    //penergies
        int,         //indv
//        int *,       //pB_outsides
//      Boole *,     //pb_comp_intermol
        /* int *,*/  //pnatoms
        float *,     //pcrds
        float *,
        int *);
//        float *,     //pcharges
//        float *,     //pABScharges
//        int *,       //ptypes
        /*maptype ,*///map
//        float *,     //pgridinfos
//        int *,       //pignore_inters
//        int *,       //pnnbs
//        Boole *,     //pincelec
//        Boole *,     //pinc14interact
//        float *,     //pscale14s
        /*float *,*/ //pqspabscharges
        /*ParameterEntry *,*/ // pparameterarrays
//        Boole *,     //pusenonbondcuts
//        Boole *,     //phaveflexresidueses
//        float *);    // punboundinternalFEs
/*, Boole *, int *, istorconstrainedtype, torstype, Boole *, UStorProfilestype);*/
      Real getMapCrd(int k, int l, int m, int n);
      State *getPTRState();
      State getState();
      vttype getVT();
      tlisttype getTList();
//      int (*getTList())[MAX_ATOMS];
      crdpdbtype getCRDPDB();
      crdtype getCRD();
      int getNAtom();
      int getStateNTOR();
      int getCRDCoord(int x, int y);
      GridMapSetInfo *getGridInfo();
      Boole getBCompIntEn();
      Real *getCharge();
      Real *getABSCharge();
      int *getType();
      maptype getMap();
      int *getIgnoreInter();
      NonbondParam *getNonBondList();
      EnergyTables *getPtrAdEnTbl();
      int getNNB();
      Boole getBCalcIntElec();
      Boole getBInc14Interact();
      Real getScale14();
      Real *getQSPABSCharge();
      ParameterEntry *getParamArray();
      Boole getBUseNonBondCut();
      Boole getBHaveFlexResid();
      Real getUnboundInternalFE();
      Boole getBIsGaussTorCon();
      Boole getBIsTorConstrainedBoole(int i);
      Boole *getBIsTorConstrained();
      float getStateTOR(int i);
      Boole getBShowTorE();
      unsigned short getUSTorEIdx(int i);
      unsigned short getUSTorProfileIdx(int i, int j);
      unsigned short *getUSTorE();
      US_torProfiletype getUSTorProfile();
      Real getChargeIdx(int i);
      void incNumEval(UnsignedFourByteLong i);
      void setUSTorEIdx(int indx, unsigned short info);
      NonbondParam getNonBondListIdx(int);
      Real getETblsolfn(int);
      Real getETblepsilonfn(int);
      Real getETblrepsilonfn(int);
#endif
      Real getETblevdWHb(int , int, int);
};
inline Real Eval::getETblevdWHb(int i, int j, int k)
{
    return ptr_ad_energy_tables->e_vdW_Hb[i][j][k];
}


#ifdef CUDA_READY

inline Real Eval::getETblsolfn(int i)
{
    return ptr_ad_energy_tables->sol_fn[i];
}

inline Real Eval::getETblepsilonfn(int i)
{
    return ptr_ad_energy_tables->epsilon_fn[i];
}
inline Real Eval::getETblrepsilonfn(int i)
{
    return ptr_ad_energy_tables->r_epsilon_fn[i];
}

inline Real Eval::getMapCrd(int k, int l, int m, int n)
{
    return map[k][l][m][n];
}

inline void Eval::setUSTorEIdx(int indx, unsigned short info)
{
    US_TorE[indx] = info;
}

inline void Eval::incNumEval(UnsignedFourByteLong i)
{
    num_evals+=i;
}
inline Real Eval::getChargeIdx(int i)
{
    return charge[i];
}

inline US_torProfiletype Eval::getUSTorProfile(void)
{
    return US_torProfile;
}

inline unsigned short *Eval::getUSTorE(void)
{
    return US_TorE;
}

inline unsigned short Eval::getUSTorProfileIdx(int i, int j)
{
    return US_torProfile[i][j];
}

inline unsigned short Eval::getUSTorEIdx(int i)
{
    return US_TorE[i];
}

inline Boole Eval::getBShowTorE(void)
{
    return B_ShowTorE;
}

inline float Eval::getStateTOR(int i)
{
    return stateNow.tor[i];
}

inline Boole *Eval::getBIsTorConstrained(void)
{
    return B_isTorConstrained;
}

inline Boole Eval::getBIsTorConstrainedBoole(int i)
{
    return B_isTorConstrained[i];
}
inline Boole Eval::getBIsGaussTorCon(void)
{
    return B_isGaussTorCon;
}

inline Real Eval::getUnboundInternalFE(void)
{
    return unbound_internal_FE;
}

inline Boole Eval::getBHaveFlexResid(void)
{
    return B_have_flexible_residues;
}

inline Boole Eval::getBUseNonBondCut(void)
{
    return B_use_non_bond_cutoff;
}

inline ParameterEntry *Eval::getParamArray(void)
{
    return parameterArray;
}

inline Real *Eval::getQSPABSCharge(void)
{
    return qsp_abs_charge;
}

inline Real Eval::getScale14(void)
{
    return scale_1_4;
}

inline Boole Eval::getBInc14Interact(void)
{
    return B_include_1_4_interactions;
}

inline Boole Eval::getBCalcIntElec(void)
{
    return B_calcIntElec;
}

inline int Eval::getNNB(void)
{
    return Nnb;
}

inline EnergyTables *Eval::getPtrAdEnTbl(void)
{
    return ptr_ad_energy_tables;
}

inline NonbondParam Eval::getNonBondListIdx(int idx)
{
    return nonbondlist[idx];
}

inline NonbondParam *Eval::getNonBondList(void)
{
    return nonbondlist;
}

inline int *Eval::getIgnoreInter(void)
{
    return ignore_inter;
}

inline maptype Eval::getMap(void)
{
    return map;
}

inline int *Eval::getType(void)
{
    return type;
}

inline Real *Eval::getABSCharge(void)
{
    return abs_charge;
}
inline Real *Eval::getCharge(void)
{
    return charge;
}

inline Boole Eval::getBCompIntEn(void)
{
    return B_compute_intermol_energy;
}
inline GridMapSetInfo * Eval::getGridInfo(void)
{
    return info;
}

inline int Eval::getCRDCoord(int x, int y)
{
    return crd[x][y];
}

inline int Eval::getNAtom(void)
{
    return natom;
}
inline crdtype Eval::getCRD(void)
{
    return crd;
}
inline crdpdbtype Eval::getCRDPDB(void)
{
    return crdpdb;
}

//inline int (*Eval::getTList(void))[MAX_ATOMS]
inline tlisttype Eval::getTList(void)
{
    return tlist;
}

inline vttype Eval::getVT(void)
{
    return vt;
}

inline State Eval::getState(void)
{
    return stateNow;
}

inline State *Eval::getPTRState(void)
{
    return &stateNow;
}

inline int Eval::getStateNTOR(void)
{
    return stateNow.ntor;
}
#endif

inline Eval::Eval(void)
: num_evals(0)
{
}

inline void Eval::setup(Real init_crd[MAX_ATOMS][SPACE],
                        Real init_charge[MAX_ATOMS],
                        Real init_abs_charge[MAX_ATOMS],
                        Real init_qsp_abs_charge[MAX_ATOMS],
                        int init_type[MAX_ATOMS],
                        int init_natom,
                        Real init_map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],

                        Real init_elec[MAX_ATOMS], // gmm added 21-Jan-1998, for writePDBQState
                        Real init_emap[MAX_ATOMS], // gmm added 21-Jan-1998, for writePDBQState
                        NonbondParam *init_nonbondlist,
                        EnergyTables   *init_ptr_ad_energy_tables,
                        int init_Nnb,
                        Boole init_B_calcIntElec,
                        Boole init_B_isGaussTorCon,
                        Boole init_B_isTorConstrained[MAX_TORS],
                        Boole init_B_ShowTorE,
                        unsigned short init_US_TorE[MAX_TORS],
                        unsigned short init_US_torProfile[MAX_TORS][NTORDIVS],
                        Real init_vt[MAX_TORS][SPACE],
                        int init_tlist[MAX_TORS][MAX_ATOMS],
                        Real init_crdpdb[MAX_ATOMS][SPACE],
                        Real init_crdreo[MAX_ATOMS][SPACE],
                        State stateInit,
                        Molecule molInit,

                        int init_ignore_inter[MAX_ATOMS],

                        Boole init_B_include_1_4_interactions,
                        Real init_scale_1_4,

                        ParameterEntry init_parameterArray[MAX_MAPS],

                        Real init_unbound_internal_FE,
                        GridMapSetInfo *init_info,
                        Boole init_B_use_non_bond_cutoff,  // set this to FALSE if we are computing unbound extended conformations
                        Boole init_B_have_flexible_residues
                       )

{
    register int i;

    crd = init_crd;
    charge = init_charge;
    abs_charge = init_abs_charge;
    qsp_abs_charge = init_qsp_abs_charge;
    type = init_type;
    natom = init_natom;
    map = init_map;

    nonbondlist = init_nonbondlist;
    ptr_ad_energy_tables = init_ptr_ad_energy_tables;
    Nnb = init_Nnb;
    B_calcIntElec = init_B_calcIntElec;
    B_isGaussTorCon = init_B_isGaussTorCon;
    B_isTorConstrained = init_B_isTorConstrained;
    B_ShowTorE = init_B_ShowTorE;
    US_TorE = init_US_TorE;
    US_torProfile = init_US_torProfile;
    vt = init_vt;
    tlist = init_tlist;
    crdpdb = init_crdpdb;
    crdreo = init_crdreo;
    // set all of the components of the State, one at a time...
    copyState( &stateNow, stateInit );
#ifdef DEBUG
    pr(logFile, "\n\nstateNow:\n");
    printState( logFile, stateNow, 2 );
#endif
    num_evals = 0;
    for (i=0; i<MAX_ATOMS; i++) {
       init_elec[i] = init_emap[i] = 0.0;
       ignore_inter[i] = init_ignore_inter[i];
    }
    mol = molInit;

    B_include_1_4_interactions = init_B_include_1_4_interactions;
    scale_1_4 = init_scale_1_4;

    parameterArray = init_parameterArray;

    unbound_internal_FE = init_unbound_internal_FE;

    info = init_info;
    B_compute_intermol_energy = TRUE; // default is "Yes, calculate the intermolecular energy".

    B_use_non_bond_cutoff = init_B_use_non_bond_cutoff;  // set this to FALSE if we are computing unbound extended conformations
    B_have_flexible_residues = init_B_have_flexible_residues;
}

inline void Eval::update_crds( Real init_crdreo[MAX_ATOMS][SPACE],
                               Real init_vt[MAX_TORS][SPACE] )
{
    crdreo = init_crdreo;
    vt = init_vt;
}

inline void Eval::compute_intermol_energy(Boole init_B_compute_intermol_energy)
    // For computing the conformation and the internal energy of the unbound state.
{
    B_compute_intermol_energy = init_B_compute_intermol_energy;
}


inline UnsignedFourByteLong Eval::evals(void)
{
   return(num_evals);
}

inline void Eval::reset(void)
{
   num_evals = 0;
}

#ifdef CUDA_READY
void cpu_alloc(int, int);
void cpu_free(void);
#endif

#endif

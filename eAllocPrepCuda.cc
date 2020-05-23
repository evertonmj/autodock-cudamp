/*
 * ePrepCuda - Prepares population for energy calculations on CUDA card
 * @author TEAM E51
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "autocomm.h"
#include "typedefs.h"
#include "support.h"
#include "assert.h"
#include "gs.h"
#include "eval.h"
#include "structs.h"
//#include "constants.h"
#include "trilinterp.h"
#include "eintcal.h"


extern FILE *logFile;
extern class Eval evaluate;

#define is_out_grid_infocuda(x,y,z) (((x)<=(::evaluate.getGridInfo()->lo[X])) || ((x)>=(::evaluate.getGridInfo()->hi[X])) || ((y)<=(::evaluate.getGridInfo()->lo[Y])) || ((y)>=(::evaluate.getGridInfo()->hi[Y])) || ((z)<=(::evaluate.getGridInfo()->lo[Z])) || ((z)>=(::evaluate.getGridInfo()->hi[Z])))

void Genetic_Algorithm::eAllocPrepCuda(Population &original_population, double invdiffwa)
{
    int i, j, B_outside = 0, I_tor = 0, indx = 0;
    float energy = 0.0L;
/*
// Arrays for trilinterp params
// MAP might be able to be done outside of function if each individual doesn't prepare it.
    bool *b_comp_intermol = (bool *)malloc(sizeof(bool) * original_population.num_individuals());
    int *natoms = (int *)malloc(sizeof(int) * original_population.num_individuals());
    crdtype *crds = (crdtype *)malloc(sizeof(crdtype) * original_population.num_individuals());
    Real **charges = (Real **)malloc(sizeof(Real *) * original_population.num_individuals());
    Real **ABScharges = (Real **)malloc(sizeof(Real *) * original_population.num_individuals());
    int **types = (int **)malloc(sizeof(int *) * original_population.num_individuals());
    maptype *maps = (maptype *)malloc(sizeof(maptype) * original_population.num_individuals());
    GridMapSetInfo **gridinfos = (GridMapSetInfo **)malloc(sizeof(GridMapSetInfo *) * original_population.num_individuals());
    int *atomsInsideOrOut = (int *)malloc(sizeof(int) * original_population.num_individuals());

    int **ignore_inters = (int **)malloc(sizeof(int *) * original_population.num_individuals());

    //eintcal params
    NonbondParam **nonbondlists = (NonbondParam **)malloc(sizeof(NonbondParam *) * original_population.num_individuals());

    EnergyTables **ptrsToEnergyTables = (EnergyTables**)malloc(sizeof(EnergyTables *) * original_population.num_individuals());

    int *nnbs = (int *)malloc(sizeof(int) * original_population.num_individuals());

    bool *incelec = (bool *)malloc(sizeof(bool) * original_population.num_individuals());

    bool *inc14interact = (bool *)malloc(sizeof(bool) * original_population.num_individuals());

    Real *scale14s = (Real *)malloc(sizeof(Real) * original_population.num_individuals());
    Real **qspabscharges = (Real **)malloc(sizeof(Real *) * original_population.num_individuals());

    ParameterEntry **parameterarrays = (ParameterEntry **)malloc(sizeof(ParameterEntry *) * original_population.num_individuals());



    bool *usenonbondcuts = (bool *)malloc(sizeof(bool) * original_population.num_individuals());

    bool *haveflexresidueses = (bool *)malloc(sizeof(bool) * original_population.num_individuals());

    Real *unboundinternalFEs = (Real *)malloc(sizeof(Real) * original_population.num_individuals());
*/

    for (i = 0; i < original_population.num_individuals(); i++)
    {
        switch(e_mode)
        {
            case Normal_Eval:
//                original_population[i].phenotyp.setEvalFlag(2);
                if (original_population[i].phenotyp.getEvalFlag() != 1 )
                {
//                fprintf(stdout, "Looking at normal_eval flag = %u\n", original_population[i].phenotyp.getEvalFlag());
                }
                if (!original_population[i].phenotyp.getEvalFlag())
                {
                  //  fprintf(stderr," Flag is false\n");
                    //prepare for ownage
                    make_state_from_rep(original_population[i].phenotyp.getRep(), ::evaluate.getPTRState());
                    cnv_state_to_coords(::evaluate.getState(),
                                        ::evaluate.getVT(),
                                        ::evaluate.getTList(),
                                        ::evaluate.getStateNTOR(),
                                        ::evaluate.getCRDPDB(),
                                        ::evaluate.getCRD(),
                                        ::evaluate.getNAtom());

                    for (j = 0; (j< ::evaluate.getNAtom() )&&(!B_outside); j++)
                    {
                          B_outside = is_out_grid_infocuda(
                                        ::evaluate.getCRDCoord(j,0),
                                        ::evaluate.getCRDCoord(j,1),
                                        ::evaluate.getCRDCoord(j,2));
                    }

/*
                    // Trilinterp params
                    // store flags for every individual
                    b_comp_intermol[i] = ::evaluate.getBCompIntEn();

                    //1st param of trilinterp
                    natoms[i] = ::evaluate.getNAtom();

                    //2nd param of trilinterp
                    memcpy(crds[i], ::evaluate.getCRD(), sizeof(crdtype));

                    //3rd param of trilinterp
                    memcpy(charges[i], ::evaluate.getCharge(), sizeof(Real *));

                    //4th param of tril
                    memcpy(ABScharges[i], ::evaluate.getABSCharge(), sizeof(Real *));

                    //5th param of tril
                    memcpy(types[i], ::evaluate.getType(), sizeof(int *));

                    //6th param of tril
                    memcpy(maps[i], ::evaluate.getMap(), sizeof(maptype));

                    //7th param of tril
                    memcpy(gridinfos[i], ::evaluate.getGridInfo(), sizeof(GridMapSetInfo *));

                    //8th param of tril
                    atomsInsideOrOut[i] = B_outside?SOME_ATOMS_OUTSIDE_GRID:ALL_ATOMS_INSIDE_GRID;

                    //9th param of tril
                    ignore_inters[i] = ::evaluate.getIgnoreInter();
                    //10th, 11th, 12th, and 13th params are constant

*/

                    if (::evaluate.getBCompIntEn())
                    {
                      //#pragma omp critical
                      //{
                        energy = trilinterp(0, ::evaluate.getNAtom(),
                        ::evaluate.getCRD(),
                        ::evaluate.getCharge(),
                        ::evaluate.getABSCharge(),
                        ::evaluate.getType(),
                        ::evaluate.getMap(),
                        ::evaluate.getGridInfo(),
                        B_outside?SOME_ATOMS_OUTSIDE_GRID:ALL_ATOMS_INSIDE_GRID,
                        ::evaluate.getIgnoreInter(),
                        NULL_ELEC,
                        NULL_EVDW,
                        NULL_ELEC_TOTAL,
                        NULL_EVDW_TOTAL);
                      //}

                    }
/*
                    //1st param of eintcal
                    memcpy(nonbondlists[i], ::evaluate.getNonBondList(), sizeof(NonbondParam *));

                    //2nd param of eintcal
                    memcpy(ptrsToEnergyTables[i], ::evaluate.getPtrAdEnTbl(), sizeof(EnergyTables *));

                    //3rd param of eintcal
                    //already stored from above

                    //4th param of eintcal
                    nnbs[i] = ::evaluate.getNNB();

                    //5th param of eintcal
                    incelec[i] = ::evaluate.getBCalcIntElec();

                    //6th param of eintcal
                    inc14interact[i] = ::evaluate.getBInc14Interact();

                    //7th param of eintcal
                    scale14s[i] = ::evaluate.getScale14();

                    //8th param of eintcal
                    memcpy(qspabscharges[i], ::evaluate.getQSPABSCharge(), sizeof(Real *));

                    //9th param of eintcal
                    memcpy(parameterarrays[i], ::evaluate.getParamArray(), sizeof(ParameterEntry *));

                    //10th param of eintcal
                    usenonbondcuts[i] = ::evaluate.getBUseNonBondCut();

                    //11th param of eintcal
                    haveflexresidueses[i] = ::evaluate.getBHaveFlexResid();

                    //12th param of eintcal
                    unboundinternalFEs[i] = ::evaluate.getUnboundInternalFE();
*/
                      energy += eintcal(
                                    ::evaluate.getNonBondList(),
                                    ::evaluate.getPtrAdEnTbl(),
                                    ::evaluate.getCRD(),
                                    ::evaluate.getNNB(),
                                    ::evaluate.getBCalcIntElec(),
                                    ::evaluate.getBInc14Interact(),
                                    ::evaluate.getScale14(),
                                    ::evaluate.getQSPABSCharge(),
                                    ::evaluate.getParamArray(),
                                    ::evaluate.getBUseNonBondCut(),
                                    ::evaluate.getBHaveFlexResid()) -
                                    ::evaluate.getUnboundInternalFE();

                        if (::evaluate.getBIsGaussTorCon())
                        {
                            for(I_tor = 0; I_tor <= ::evaluate.getStateNTOR(); I_tor++)
                            {
                                if (::evaluate.getBIsTorConstrainedBoole(I_tor)
                                    == 1)
                                {
                                    indx = RadiansToDivs( WrpModRad(
                                               ::evaluate.getStateTOR(I_tor) ));
                                   if (::evaluate.getBShowTorE())
                                   {
                                      ::evaluate.setUSTorEIdx(I_tor, ::evaluate.getUSTorProfileIdx(I_tor,indx));
                                      energy += (float)(::evaluate.getUSTorEIdx(I_tor));
                                }
                                else
                                {
                                    energy += (float)::evaluate.getUSTorProfileIdx(I_tor, indx);
                                }
                            }
                        }
                    }

                    ::evaluate.incNumEval(1);
                    if (!finite(energy))
                    {
                        (void)fprintf(logFile, "eval.cc: ERROR! energ is infinite!\n\n");
                        for (j = 0; j < ::evaluate.getNAtom(); j++)
                        {
                            (void)fprintf(logFile, FORMAT_PDBQ_ATOM_RESSTR, "",
                                    j+1, "C   INF     1", ::evaluate.getCRDCoord(j,X), ::evaluate.getCRDCoord(i,Y), ::evaluate.getCRDCoord(i,Z), 0.0, 0.0, ::evaluate.getChargeIdx(j));
                            (void)fprintf(logFile, "\n");
                        }

                    }

                    if (ISNAN(energy))
                    {
                        (void)fprintf( logFile, "eval.cc:  ERROR!  energy is not a number!\n\n");
                        for (j=0; j< ::evaluate.getNAtom(); j++)
                        {
                            (void)fprintf(logFile, FORMAT_PDBQ_ATOM_RESSTR, "",
                                    j+1, "C   NaN     1", ::evaluate.getCRDCoord(j,X), ::evaluate.getCRDCoord(i,Y), ::evaluate.getCRDCoord(i,Z), 0.0, 0.0, ::evaluate.getChargeIdx(j));
                            (void)fprintf(logFile, "\n");
                        }
                    }

                    original_population[i].phenotyp.setValue(energy);
                    original_population[i].phenotyp.setEvalFlag(1);
                    break;
                }
            break;

            default:
            fprintf(logFile, "Unknown Evaluation Mode!\n");
            break;
        }

//                    fprintf(stderr, "energy = %E\n",energy);
//        fprintf(stderr, "worst = %E, energy = %E, invdiffwa = %E\n", worst, energy, invdiffwa);
//get value
        alloc[i] = (worst - original_population[i].phenotyp.getValue()) * invdiffwa;
        assert(finite(energy));
        assert(finite(alloc[i]));
        assert(!ISNAN(energy));
        assert(!ISNAN(alloc[i]));

    }
}

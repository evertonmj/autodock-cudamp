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
#include "cuda_wrapper.h"
#include "eval_wrapper.h"
#include "eMapPrepCuda.h"

/*extern "C" void eval_wrapper(
    Boole  b_comp_intermol,
    float *crds,
    float *cpuenergies,
    float *float_arraycpu,
    int *int_arraycpu);*/

extern FILE *logFile;
extern class Eval evaluate;
extern int Nnb_array[3];

float *float_arraycpu;
int *int_arraycpu;
Boole b_comp_intermol;
int natoms;
float *crds;
double *energies;
float *cpuenergies;

/**
 * Separate cuda preparation function, such that it checks if energies are equal at end
 * @param original_population, original population
 * @param e_mode energy mode (only supports Always_Eval and Normal_Eval
 * @param firstEnergy first energy (needed for comparrison)
 * @return whether the evals are equal or not
 */
int eEqualPrepCuda(Population &original_population, int e_mode, double firstEnergy)
{
    int i;
    for (i = 0; i < original_population.num_individuals(); i++)
    {
        e_prep(e_mode, i, original_population);
    }

    eval_wrapper(
      b_comp_intermol,
      crds,
      cpuenergies,
      float_arraycpu,
      int_arraycpu);


    return e_complete_evalEqual(original_population, ::evaluate.getCharge(), firstEnergy);

}

/**
 * Mapping evaluation cuda preparation (does necessary functionality that mapping does)
 * @param original_population original_population
 * @param e_mode energy mode (only support Normal_Eval and Always_Eval)
 */
void eMapPrepCuda(Population &original_population, int e_mode)
{
    int k,l,m,n;
    int j;
    int i;
    int indv;
    int     I_tor = 0;
    int     indx = 0;
    double energy = 0.0L;

     #pragma omp critical
     {
      //int individual = ignlgi_t(omp_get_thread_num())%original_population.num_individuals();
      //int gene_number = ignlgi_t(omp_get_thread_num())%original_population[individual].genotyp.num_genes();

      for (i = 0; i < original_population.num_individuals(); i++)
      {
          original_population[i].phenotyp.write(*(original_population[i].genotyp.vread(0)), 0);
          original_population[i].phenotyp.write(*(original_population[i].genotyp.vread(1)), 1);
          original_population[i].phenotyp.write(*(original_population[i].genotyp.vread(2)), 2);
          original_population[i].phenotyp.write(*(original_population[i].genotyp.vread(3)), 3);
          original_population[i].phenotyp.write(*(original_population[i].genotyp.vread(4)), 4);

          e_prep(e_mode, i, original_population);
          /*original_population[individual].phenotyp.write(*(original_population[individual].genotyp.vread(0)), 0);
          original_population[individual].phenotyp.write(*(original_population[individual].genotyp.vread(1)), 1);
          original_population[individual].phenotyp.write(*(original_population[individual].genotyp.vread(2)), 2);
          original_population[individual].phenotyp.write(*(original_population[individual].genotyp.vread(3)), 3);
          original_population[individual].phenotyp.write(*(original_population[individual].genotyp.vread(4)), 4);

          e_prep(e_mode, individual, original_population);*/
    }

    eval_wrapper(
      b_comp_intermol,
      crds,
      cpuenergies,
      float_arraycpu,
      int_arraycpu);
   }


    e_complete_eval(original_population, ::evaluate.getCharge());
}

/*
 * Prepares data for gpu calculations
 * @param pop_size population size
 * @param penergies energies array
 * @param indv which individual we are looking at
 * @param pcrds crds
 * @param pfloat_arraycpu float array of variables
 * @param pint_arraycpu int array of variables
 */
void Eval::evalPrep(unsigned int pop_size,
                    double * penergies,
                        int indv,
                        float *pcrds,
                        float *pfloat_arraycpu,
                        int *pint_arraycpu)
{
    register int i;
    int     B_outside= 0;
    int     I_tor = 0;
    int     indx = 0;
    double energy = 0.0L;

      cnv_state_to_coords(stateNow, vt, tlist, stateNow.ntor, crdpdb, crd, natom);

      for (i=0; (i<natom)&&(!B_outside); i++) {
        B_outside = is_out_grid_info(crd[i][0], crd[i][1], crd[i][2]);
      }

      int j,k,l,m;

      //trilinterp setup
      pint_arraycpu[INTBOUTS * pop_size + indv] = (int)B_outside;
      for (j = 0; j < natom; j++)
      {
        for(k = 0; k < SPACE; k++)
        {
          pcrds[indv * natom * SPACE + j * SPACE + k] = (float)crd[j][k];
        }
      }

      pfloat_arraycpu[FLOATINFO * pop_size + 0] = (float)info->center[0];
      pfloat_arraycpu[FLOATINFO * pop_size + 1] = (float)info->center[1];
      pfloat_arraycpu[FLOATINFO * pop_size + 2] = (float)info->center[2];
      pfloat_arraycpu[FLOATINFO * pop_size + 3] = (float)info->lo[0];
      pfloat_arraycpu[FLOATINFO * pop_size + 4] = (float)info->lo[1];
      pfloat_arraycpu[FLOATINFO * pop_size + 5] = (float)info->lo[2];
      pfloat_arraycpu[FLOATINFO * pop_size + 6] = (float)info->inv_spacing;
      pfloat_arraycpu[FLOATINFO * pop_size + 7] = (float)info->hi[0];
      pfloat_arraycpu[FLOATINFO * pop_size + 8] = (float)info->hi[1];
      pfloat_arraycpu[FLOATINFO * pop_size + 9] = (float)info->hi[2];


      //eintcal setup
      pint_arraycpu[INTINCELEC * pop_size + indv] = (int)B_calcIntElec;
      pfloat_arraycpu[FLOATSCALE14 * pop_size + indv] = (float)scale_1_4;
      pint_arraycpu[INTNONBONDCUT * pop_size + indv] = (int)B_use_non_bond_cutoff;

      pfloat_arraycpu[FLOATUNBOUNDINTERNAL * pop_size + indv] = (float)unbound_internal_FE;

      penergies[indv] = 0.0L;

      if (B_isGaussTorCon) {
        for (I_tor = 0; I_tor <= stateNow.ntor; I_tor++)
        {
          if (B_isTorConstrained[I_tor] == 1) {
            indx = RadiansToDivs( WrpModRad(stateNow.tor[I_tor]) );
            if (B_ShowTorE) {
              penergies[indv] += (double)(US_TorE[I_tor] = US_torProfile[I_tor][indx]);
            } else {
              penergies[indv] += (double)US_torProfile[I_tor][indx];
            }
          }
        } // I_tor
      }//if

      num_evals++;
}

/**
 * Allocates memory for the CPU
 * @param num_individuals number of individuals
 * @param natoms number of atoms
 */
void cpu_alloc(int num_individuals, int natoms)
{
    float_arraycpu = (float *)malloc(sizeof(float) * FLOATSIZE * num_individuals);
    int_arraycpu = (int *)malloc(sizeof(int) * INTSIZE * num_individuals);
    b_comp_intermol = ::evaluate.getBCompIntEn();
    natoms = ::evaluate.getNAtom();
    crds = (float *)malloc(sizeof(float) * natoms * SPACE * num_individuals);
    energies = (double *)malloc(sizeof(double) * num_individuals);
    cpuenergies = (float *)malloc(sizeof(float) * num_individuals);

}

/**
 * Frees memory for the CPU
 */
void cpu_free(void)
{
    free(float_arraycpu);
    free(int_arraycpu);
    free(crds);
    free(energies);
    free(cpuenergies);
}

/**
 * Helper function for preparing evaluations on GPU
 * @param e_mode energy mode
 * @param i individual being looked at
 * @param original_population original population
 */
void e_prep(int e_mode, int i, Population &original_population)
{
    //int individual = ignlgi_t(omp_get_thread_num())%original_population.num_individuals();
    //gene_number = ignlgi_t(omp_get_thread_num())%pure[individual].genotyp.num_genes();
    int index = INTEVALFLAG * original_population.num_individuals() + i;


    switch(e_mode)
    {
        case Always_Eval:

            //int_arraycpu[INTEVALFLAG * original_population.num_individuals() + i] = 0;
            int_arraycpu[index] = 0;
            //int_arraycpu[INTEVALFLAG * individual + i] = 0;
              make_state_from_rep(original_population[i].phenotyp.getRep(), ::evaluate.getPTRState());
            ::evaluate.evalPrep(
              original_population.num_individuals(),
              energies,
              i,
              crds,
              float_arraycpu,
              int_arraycpu);

              original_population[i].phenotyp.setEvalFlag(1);


            break;

        case Normal_Eval:
            //int_arraycpu[INTEVALFLAG * original_population.num_individuals() + i] = (int)original_population[i].phenotyp.getEvalFlag();

            int_arraycpu[index] = (int)original_population[i].phenotyp.getEvalFlag();

            //int_arraycpu[INTEVALFLAG * individual + i] = (int)original_population[i].phenotyp.getEvalFlag();
            if (!original_population[i].phenotyp.getEvalFlag())
            {
                make_state_from_rep(original_population[i].phenotyp.getRep(), ::evaluate.getPTRState());
                ::evaluate.evalPrep(
                    original_population.num_individuals(),
                    energies,
                    i,
                    crds,
                    float_arraycpu,
                    int_arraycpu);

                  original_population[i].phenotyp.setEvalFlag(1);
            }
        break;
        default:
        fprintf(logFile, "Unknown Evaluation Mode!\n");
        break;
    }
}

/**
 * Completes the evaluation (checks for finite and if answers are not numbers
 * Also stores remaining values into class
 * @param original_population original population
 * @param charges charges
 */
void e_complete_eval(Population &original_population, Real *charges)
{
    int indv;
    int i;
    for (indv = 0; indv < original_population.num_individuals(); indv++)
    {
        if(!int_arraycpu[INTEVALFLAG * original_population.num_individuals() + indv])
        {

            energies[indv] += (double)cpuenergies[indv];

            if (!finite(energies[indv]))
            {
                (void)fprintf( logFile, "eval.cc:  ERROR!  energy is infinite!\n\n");
                for (i=0; i<natoms; i++)
                {
                    (void)fprintf(logFile, FORMAT_PDBQ_ATOM_RESSTR, "", i+1, "C   INF     1", crds[indv * natoms * SPACE + i * SPACE + X], crds[indv * natoms * SPACE + i * SPACE + Y], crds[indv * natoms * SPACE + i * SPACE + Z], 0.0, 0.0, charges[i]);
                    (void)fprintf(logFile, "\n");
                } // i
            }
            if (ISNAN(energies[indv]))
            {
                (void)fprintf( logFile, "eval.cc:  ERROR!  energy is not a number!\n\n");
                for (i=0; i<natoms; i++)
                {
                    (void)fprintf(logFile, FORMAT_PDBQ_ATOM_RESSTR, "", i+1, "C   NaN     1", crds[indv * natoms * SPACE + i * SPACE + X], crds[indv * natoms * SPACE + i * SPACE + Y], crds[indv * natoms * SPACE + i * SPACE + Z], 0.0, 0.0, charges[i]);
                    (void)fprintf(logFile, "\n");
               } // i
            }

            original_population[indv].phenotyp.setValue(energies[indv]);
        }
    }
}


/**
 * Completes the evaluation (checks for finite and if answers are not numbers
 * Also stores remaining values into class and does check for energies equal
 * @param original_population original population
 * @param charges charges
 * @return if energies are equal
 */
int e_complete_evalEqual(Population &original_population, Real *charges, double firstEnergy)
{
    int indv;
    int i;
    int energiesEqual = 1;
    for (indv = 1; indv < original_population.num_individuals(); indv++)
    {
        if(!int_arraycpu[INTEVALFLAG * original_population.num_individuals() + indv])
        {

//            fprintf(stderr, "energies = %E\n", energies[indv]);

            energies[indv] += (double)cpuenergies[indv];
            if (!finite(energies[indv]))
            {
                (void)fprintf( logFile, "eval.cc:  ERROR!  energy is infinite!\n\n");
                for (i=0; i<natoms; i++)
                {
                    (void)fprintf(logFile, FORMAT_PDBQ_ATOM_RESSTR, "", i+1, "C   INF     1", crds[indv * natoms * SPACE + i * SPACE + X], crds[indv * natoms * SPACE + i * SPACE + Y], crds[indv * natoms * SPACE + i * SPACE + Z], 0.0, 0.0, charges[i]);
                    (void)fprintf(logFile, "\n");
                } // i
            }
            if (ISNAN(energies[indv]))
            {
                (void)fprintf( logFile, "eval.cc:  ERROR!  energy is not a number!\n\n");
                for (i=0; i<natoms; i++)
                {
                    (void)fprintf(logFile, FORMAT_PDBQ_ATOM_RESSTR, "", i+1, "C   NaN     1", crds[indv * natoms * SPACE + i * SPACE + X], crds[indv * natoms * SPACE + i * SPACE + Y], crds[indv * natoms * SPACE + i * SPACE + Z], 0.0, 0.0, charges[i]);
                    (void)fprintf(logFile, "\n");
               } // i
            }

            original_population[indv].phenotyp.setValue(energies[indv]);
        }
        energiesEqual = energiesEqual && (original_population[indv].phenotyp.getValue() == firstEnergy);
    }

    return energiesEqual;
}

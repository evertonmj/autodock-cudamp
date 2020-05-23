/*
 * Wrapper for selection allocation
 * Compiled with Cuda compiler.
 */

// includes, system
#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include "typedefs.h"

//C++ #defines
#include "autocomm.h"
#include "grid.h"
#include "eval.h"
#include "constants.h"
#include "trilinterp.h"
#include "eintcal.h"
#include "distdepdiel.h"
#include "cuda_wrapper.h"
#include "eval_wrapper.h"

texture<float, 1> tex;

// Variables required for energy calculations
//added extern to variables. Everton Mendonça 03/04/2016
int *B_outsidesgpu;
float *energiesgpu;
unsigned int *evalflagsgpu;
int cpunatoms;
int nBlocks;
int blocksize;
unsigned int num_individuals;
float *nb_group_energycpu;

float *float_arraygpu;
int *int_arraygpu;

//tril params
float *crdsgpu;
float *chargesgpu;
float *ABSchargesgpu;
int   *typesgpu;
float *gridinfosgpu;
int   *ignore_intersgpu;
float *cudamap;

//eint cal params
float *nonbondlistsgpu;
Boole *incelecgpu;
Boole inc14interactgpu;
float *scale14sgpu;
Boole *usenonbondcutsgpu;
Boole haveflexresiduesgpu;
float *unboundinternalFEsgpu;
float *evdWHbgpu;
float *solfngpu;
float *epsilonfngpu;
float *repsilonfngpu;
int *nnb_arraygpu;
float *nb_group_energygpu;

#define BLOCK_SIZE 128

// Global variables required for energy calculations
/*int ElecMap = 0;
int DesolvMap = 0;
Real nb_group_energy[3];
int Nnb_array[3];

void CHECK_ERROR(int num);

/**
 * Trilinterp GPU kernel, does trilinterp energy calculations for each
 * individual in the population.
 * @param num_individualsgpu number of individuals in population
 * @param penergiesgpu array of energies used to store individual's energy
 * @param b_comp_intermolgpu flag (used in cpu trilinterp)
 * @param natomsgpu number of atoms (used in cpu trilinterp)
 * @param crdsgpu (used in cpu trilinterp)
 * @param chargesgpu (used in cpu trilinterp)
 * @param ABSchargesgpu (used in cpu trilinterp)
 * @param typesgpu (used in cpu trilinterp)
 * @param ignore_intersgpu (used in cpu trilinterp)
 * @param p_elec_total (used in cpu trilinterp)
 * @param p_emap_total (used in cpu trilinterp)
 * @param elecMap (used in cpu trilinterp)
 * @param desolvMap (used in cpu trilinterp)
 * @param SomeAtomsOutside (used in cpu trilinterp)
 * @param AllAtomsInside (used in cpu trilinterp)
 * @param pfloat_arraygpu array of float variables used in cpu trilinterp
 * @param pint_arraygpu array of integer varibales used in cpu trilinterp
 */

__global__ void eval_tril_kernel(unsigned int num_individualsgpu,
                                float *penergiesgpu,
                                Boole b_comp_intermolgpu,
                                int natomsgpu,
                                float *crdsgpu,
                   float *chargesgpu,
                   float *ABSchargesgpu,
                   int *typesgpu,
                   int * ignore_intersgpu,
                   float *p_elec_total,
                   float *p_emap_total,
                   int elecMap,
                   int desolvMap,
                   int SomeAtomsOutside,
                   int AllAtomsInside,
                   float *pfloat_arraygpu,
                   int *pint_arraygpu)
{
    int idx = blockIdx.x  * blockDim.x + threadIdx.x;

    if (idx < num_individualsgpu)
    {
        int some_atoms_outside_grid;

        if ((int)pint_arraygpu[INTBOUTS * num_individualsgpu + idx]) //B_outsidesgpu[idx])
        {
            some_atoms_outside_grid = SomeAtomsOutside;
        }
        else
        {
            some_atoms_outside_grid = AllAtomsInside;
        }

        if (!(unsigned int)pint_arraygpu[INTEVALFLAG * num_individualsgpu + idx])//!evalflagsgpu[idx])
        {
            if (b_comp_intermolgpu)
            {
//            fprintf(stderr, "woot\n");
                float elec_total=0.0f, emap_total=0.0f;
                int i;

                for (i=0; i<natomsgpu;i++)
                {
                    float e, m, d;
                    float u, v, w;
                    float p0u, p0v, p0w;
                    float p1u, p1v, p1w;
                    int AtomType;
                    int u0, v0, w0;
                    int u1, v1, w1;
                    if (ignore_intersgpu[i])
                    {
                        //if (elec != NULL) elec[i] = 0;
                        //if (emap != NULL) emap[i] = 0;
                        continue;
                    }
                    if (some_atoms_outside_grid)
                    {
                        float x,y,z;
                        x = crdsgpu[idx * natomsgpu * SPACE + i * SPACE + X];
                        y = crdsgpu[idx * natomsgpu * SPACE + i * SPACE + Y];
                        z = crdsgpu[idx * natomsgpu * SPACE + i * SPACE + Z];

                        if (((x)<=((float)pfloat_arraygpu[FLOATINFO * num_individualsgpu + 3])) || ((x)>=((float)pfloat_arraygpu[FLOATINFO * num_individualsgpu + 7])) || ((y)<=((float)pfloat_arraygpu[FLOATINFO * num_individualsgpu + 4])) ||
                                ((y)>=((float)pfloat_arraygpu[FLOATINFO * num_individualsgpu + 8])) || ((z)<=((float)pfloat_arraygpu[FLOATINFO * num_individualsgpu + 5])) || ((z)>=((float)pfloat_arraygpu[FLOATINFO * num_individualsgpu + 9])))
                        {
                            float epenalty;
                            x -= (float)pfloat_arraygpu[FLOATINFO * num_individualsgpu + 0];
                            y -= (float)pfloat_arraygpu[FLOATINFO * num_individualsgpu + 1];
                            z -= (float)pfloat_arraygpu[FLOATINFO * num_individualsgpu + 2];

                            epenalty = sqhypotenuse(x,y,z) * ENERGYPENALTY;
                            //if (elec != NULL) elec[i] = epenalty;
                            //if (emap != NULL) emap[i] = epenalty;
                            elec_total += epenalty;
                            emap_total += epenalty;
                            continue;
                        }
                    }


                    AtomType = typesgpu[i];

                    u1  = (u0 = (int) (u = ((float)crdsgpu[idx * natomsgpu * SPACE + i * SPACE + X]-(float)pfloat_arraygpu[FLOATINFO * num_individualsgpu + 3]) * (float)pfloat_arraygpu[FLOATINFO * num_individualsgpu + 6])) + 1;
                    p1u = 1.0f - (p0u = u - (float) u0);

                    v1  = (v0 = (int) (v = ((float)crdsgpu[idx * natomsgpu * SPACE + i * SPACE + Y]-(float)pfloat_arraygpu[FLOATINFO * num_individualsgpu + 4]) * (float)pfloat_arraygpu[FLOATINFO * num_individualsgpu + 6])) + 1;
                    p1v = 1.0f - (p0v = v - (float) v0);

                    w1  = (w0 = (int) (w = ((float)crdsgpu[idx * natomsgpu * SPACE + i * SPACE + Z]-(float)pfloat_arraygpu[FLOATINFO * num_individualsgpu + 5]) * (float)pfloat_arraygpu[FLOATINFO * num_individualsgpu + 6])) + 1;
                    p1w = 1.0f - (p0w = w - (float) w0);


            #ifdef MINPOINT
                    int ix,iy,iz;                      //MINPOINT
                    ix = (p0u < p1u)? u0 : u1;                  //MINPOINT
                    iy = (p0v < p1v)? v0 : v1;                  //MINPOINT
                    iz = (p0w < p1w)? w0 : w1;                  //MINPOINT


                    e = tex1Dfetch(tex, iz * MAX_GRID_PTS * MAX_GRID_PTS * MAX_MAPS + iy * MAX_GRID_PTS * MAX_MAPS + ix * MAX_MAPS + elecMap);               //MINPOINT
                    m = tex1Dfetch(tex, iz * MAX_GRID_PTS * MAX_GRID_PTS * MAX_MAPS + iy * MAX_GRID_PTS * MAX_MAPS + ix * MAX_MAPS + AtomType);              //MINPOINT
                    d = tex1Dfetch(tex, iz * MAX_GRID_PTS * MAX_GRID_PTS * MAX_MAPS + iy * MAX_GRID_PTS * MAX_MAPS + ix * MAX_MAPS + desolvMap);             //MINPOINT
            #else

                    e = m = d = 0.0f;

                    e += p1u * p1v * p1w * tex1Dfetch(tex, w0 * MAX_GRID_PTS * MAX_GRID_PTS * MAX_MAPS + v0 * MAX_GRID_PTS * MAX_MAPS + u0 * MAX_MAPS + elecMap);
                    m += p1u * p1v * p1w * tex1Dfetch(tex, w0 * MAX_GRID_PTS * MAX_GRID_PTS * MAX_MAPS + v0 * MAX_GRID_PTS * MAX_MAPS + u0 * MAX_MAPS + AtomType);
                    d += p1u * p1v * p1w * tex1Dfetch(tex, w0 * MAX_GRID_PTS * MAX_GRID_PTS * MAX_MAPS + v0 * MAX_GRID_PTS * MAX_MAPS + u0 * MAX_MAPS + desolvMap);

                    d += p0u * p1v * p1w * tex1Dfetch(tex, w0 * MAX_GRID_PTS * MAX_GRID_PTS * MAX_MAPS + v0 * MAX_GRID_PTS * MAX_MAPS + u1 * MAX_MAPS + desolvMap);
                    m += p0u * p1v * p1w * tex1Dfetch(tex, w0 * MAX_GRID_PTS * MAX_GRID_PTS * MAX_MAPS + v0 * MAX_GRID_PTS * MAX_MAPS + u1 * MAX_MAPS + AtomType);
                    e += p0u * p1v * p1w * tex1Dfetch(tex, w0 * MAX_GRID_PTS * MAX_GRID_PTS * MAX_MAPS + v0 * MAX_GRID_PTS * MAX_MAPS + u1 * MAX_MAPS + elecMap);

                    e += p1u * p0v * p1w * tex1Dfetch(tex, w0 * MAX_GRID_PTS * MAX_GRID_PTS * MAX_MAPS + v1 * MAX_GRID_PTS * MAX_MAPS + u0 * MAX_MAPS + elecMap);
                    m += p1u * p0v * p1w * tex1Dfetch(tex, w0 * MAX_GRID_PTS * MAX_GRID_PTS * MAX_MAPS + v1 * MAX_GRID_PTS * MAX_MAPS + u0 * MAX_MAPS + AtomType);
                    d += p1u * p0v * p1w * tex1Dfetch(tex, w0 * MAX_GRID_PTS * MAX_GRID_PTS * MAX_MAPS + v1 * MAX_GRID_PTS * MAX_MAPS + u0 * MAX_MAPS + desolvMap);

                    d += p0u * p0v * p1w * tex1Dfetch(tex, w0 * MAX_GRID_PTS * MAX_GRID_PTS * MAX_MAPS + v1 * MAX_GRID_PTS * MAX_MAPS + u1 * MAX_MAPS + desolvMap);
                    m += p0u * p0v * p1w * tex1Dfetch(tex, w0 * MAX_GRID_PTS * MAX_GRID_PTS * MAX_MAPS + v1 * MAX_GRID_PTS * MAX_MAPS + u1 * MAX_MAPS + AtomType);
                    e += p0u * p0v * p1w * tex1Dfetch(tex, w0 * MAX_GRID_PTS * MAX_GRID_PTS * MAX_MAPS + v1 * MAX_GRID_PTS * MAX_MAPS + u1 * MAX_MAPS + elecMap);

                    e += p1u * p1v * p0w * tex1Dfetch(tex, w1 * MAX_GRID_PTS * MAX_GRID_PTS * MAX_MAPS + v0 * MAX_GRID_PTS * MAX_MAPS + u0 * MAX_MAPS + elecMap);
                    m += p1u * p1v * p0w * tex1Dfetch(tex, w1 * MAX_GRID_PTS * MAX_GRID_PTS * MAX_MAPS + v0 * MAX_GRID_PTS * MAX_MAPS + u0 * MAX_MAPS + AtomType);
                    d += p1u * p1v * p0w * tex1Dfetch(tex, w1 * MAX_GRID_PTS * MAX_GRID_PTS * MAX_MAPS + v0 * MAX_GRID_PTS * MAX_MAPS + u0 * MAX_MAPS + desolvMap);

                    d += p0u * p1v * p0w * tex1Dfetch(tex, w1 * MAX_GRID_PTS * MAX_GRID_PTS * MAX_MAPS + v0 * MAX_GRID_PTS * MAX_MAPS + u1 * MAX_MAPS + desolvMap);
                    m += p0u * p1v * p0w * tex1Dfetch(tex, w1 * MAX_GRID_PTS * MAX_GRID_PTS * MAX_MAPS + v0 * MAX_GRID_PTS * MAX_MAPS + u1 * MAX_MAPS + AtomType);
                    e += p0u * p1v * p0w * tex1Dfetch(tex, w1 * MAX_GRID_PTS * MAX_GRID_PTS * MAX_MAPS + v0 * MAX_GRID_PTS * MAX_MAPS + u1 * MAX_MAPS + elecMap);

                    e += p1u * p0v * p0w * tex1Dfetch(tex, w1 * MAX_GRID_PTS * MAX_GRID_PTS * MAX_MAPS + v1 * MAX_GRID_PTS * MAX_MAPS + u0 * MAX_MAPS + elecMap);
                    m += p1u * p0v * p0w * tex1Dfetch(tex, w1 * MAX_GRID_PTS * MAX_GRID_PTS * MAX_MAPS + v1 * MAX_GRID_PTS * MAX_MAPS + u0 * MAX_MAPS + AtomType);
                    d += p1u * p0v * p0w * tex1Dfetch(tex, w1 * MAX_GRID_PTS * MAX_GRID_PTS * MAX_MAPS + v1 * MAX_GRID_PTS * MAX_MAPS + u0 * MAX_MAPS + desolvMap);

                    d += p0u * p0v * p0w * tex1Dfetch(tex, w1 * MAX_GRID_PTS * MAX_GRID_PTS * MAX_MAPS + v1 * MAX_GRID_PTS * MAX_MAPS + u1 * MAX_MAPS + desolvMap);
                    m += p0u * p0v * p0w * tex1Dfetch(tex, w1 * MAX_GRID_PTS * MAX_GRID_PTS * MAX_MAPS + v1 * MAX_GRID_PTS * MAX_MAPS + u1 * MAX_MAPS + AtomType);
                    e += p0u * p0v * p0w * tex1Dfetch(tex, w1 * MAX_GRID_PTS * MAX_GRID_PTS * MAX_MAPS + v1 * MAX_GRID_PTS * MAX_MAPS + u1 * MAX_MAPS + elecMap);

            #endif // not MINPOINT


                    elec_total += e * chargesgpu[i];
                    emap_total += m + d * ABSchargesgpu[i];

                    //if (elec != NULL)
                    //{
                        //elec[i] = chargesgpu[idx * natomsgpu + i]; // e
                    //}
    //                if (emap != NULL) emap[i] =ABSchargesgpu[idx * natomsgpu + i];// m + d * ABSchargesgpu[idx * natomsgpu + i];



                }

    //            if (p_elec_total != NULL) *p_elec_total = elec_total;
    //            if (p_emap_total != NULL) *p_emap_total = emap_total;
                //energiesgpu[idx] = (float)(elec_total+emap_total);
                penergiesgpu[idx] = (float)(elec_total+emap_total);
//                fprintf(stderr," energy = %f, energy gpu = %f\n", (float)(elec_total +emap_total), pfloat_arraygpu[FLOATENERGIES * num_individualsgpu + idx]);
            }

        }
    }
}

/**
 * eintcal GPU kernel, does eintcal energy calculations for each
 * individual in the population.
 * @param num_individualsgpu number of individuals in population
 * @param natomsgpu number of atoms
 * @param penergiesgpu array of energies used to store individual's energy
 * @param nonbondlist (used in cpu eintcal)
 * @param tcoord (used in cpu eintcal)
 * @param B_include_1_4_interactions (used in cpu eintcal)
 * @param B_have_flexible_residues (used in cpu eintcal)
 * @param nnb_array (used in cpu eintcal)
 * @param Nb_group_energy (used in cpu eintcal)
 * @param stre_vdW_Hb (used in cpu eintcal)
 * @param strsol_fn (used in cpu eintcal)
 * @param strepsilon_fn (used in cpu eintcal)
 * @param strr_epsilon_fn (used in cpu eintcal)
 * @param b_comp_intermolgpu (used in cpu eintcal)
 * @param pfloat_arraygpu array of float variables used in cpu trilinterp
 * @param pint_arraygpu array of integer varibales used in cpu trilinterp
 */
__global__ void eintcal_kernel(
                        unsigned int num_individualsgpu,
                        int natomsgpu,
                        float *penergiesgpu,
                        float *nonbondlist,
                        float *tcoord,
                        Boole B_include_1_4_interactions,
                        Boole B_have_flexible_residues,
                        int *nnb_array,
                        float *Nb_group_energy,
                        float *stre_vdW_Hb,
                        float *strsol_fn,
                        float *strepsilon_fn,
                        float *strr_epsilon_fn,
                        Boole b_comp_intermolgpu,
                        float *pfloat_arraygpu,
                        int *pint_arraygpu)
{

    int idx = blockIdx.x  * blockDim.x + threadIdx.x;

    if (idx < num_individualsgpu)
    {

        if (!pint_arraygpu[INTEVALFLAG * num_individualsgpu + idx])//!evalflagsgpu[idx])
        {

    #ifndef EINTCALPRINT
    #   ifndef NOSQRT
            float r = 0.0f;
//            float nbc = B_use_non_bond_cutoff[idx] ? NBC : 999;
            float nbc = (Boole)pint_arraygpu[INTNONBONDCUT * num_individualsgpu + idx] ? NBC : 999;
    #   else
//            float nbc2 = B_use_non_bond_cutoff[idx] ? NBC2 : 999 * 999;
//            float nbc2 = (Boole)pint_arraygpu[INTNONBONDCUT * num_individualsgpu + idx] ? NBC2 : 999 * 999;
            float nbc = (Boole)pint_arraygpu[INTNONBONDCUT * num_individualsgpu + idx] ? NBC2 : 999 * 999;
    #   endif

    #else
    #   ifndef NOSQRT
            float d = 0.0f;
//            float nbc = B_use_non_bond_cutoff[idx] ? NBC : 999;
            float nbc = (Boole)pint_arraygpu[INTNONBONDCUT * num_individualsgpu + idx] ? NBC : 999;
    #   else
//            float nbc2 = B_use_non_bond_cutoff[idx] ? NBC2 : 999 * 999;
            float nbc = (Boole)pint_arraygpu[INTNONBONDCUT * num_individualsgpu + idx] ? NBC2 : 999 * 999;
    #   endif
    #endif

            float dx = 0.0f, dy = 0.0f, dz = 0.0f;
            float r2 = 0.0f;

            float total_e_internal = 0.0f;

            float e_elec = 0.0f;

    #ifdef EINTCALPRINT
            float total_e_elec = 0.0f;
            float total_e_vdW_Hb = 0.0f;
            float e_vdW_Hb = 0.0f;
            float total_e_desolv = 0.0f;
    #endif

            int inb = 0;
            int a1 = 0, a2 = 0;
            int t1 = 0, t2 = 0;
            int nonbond_type = 0;

            int index_1t_NEINT = 0;
            int index_1t_NDIEL = 0;
            int nb_group = 0;
            int inb_from = 0;
            int inb_to = 0;
            int nb_group_max = 1;

            if (B_have_flexible_residues)
            {
                nb_group_max = 3;
            }

            for (nb_group = 0; nb_group < nb_group_max; nb_group++)
            {
    #ifdef EINTCALPRINT
                if (nb_group ==0)
                {
                    //prints stuff
                }
                if (nb_group == 1)
                {
                    //prints stuff
                }
                if (nb_group == 2)
                {
                    //prints stuff
                }
                if ((Boole)pint_arraygpu[INTINCELEC * num_individualsgpu + idx])//B_calcIntElec[idx])
                {
                    //prints stuff
                } else {
                    //prints stuff
                }
    #endif


                if (nb_group == 0)
                {
                    inb_from = 0;
                } else {
                    inb_from = nnb_array[nb_group-1];
                }
                inb_to = nnb_array[nb_group];

                for (inb = inb_from; inb < inb_to; inb++)
                {

                    float e_internal = 0.0f;
                    float e_desolv = 0.0f;

                    a1 = (int)nonbondlist[inb * 7 + 0];
                    a2 = (int)nonbondlist[inb * 7 + 1];
                    t1 = (int)nonbondlist[inb * 7 + 2];
                    t2 = (int)nonbondlist[inb * 7 + 3];

                    nonbond_type = (int)nonbondlist[inb * 7 + 4];
                    float nb_desolv = nonbondlist[inb * 7  + 5];
                    float q1q2 = nonbondlist[inb * 7 + 6];


                    dx = tcoord[idx * natomsgpu * SPACE + a1 * SPACE + X] - tcoord[idx * natomsgpu * SPACE + a2 * SPACE + X];
                    dy = tcoord[idx * natomsgpu * SPACE + a1 * SPACE + Y] - tcoord[idx * natomsgpu * SPACE + a2 * SPACE + Y];
                    dz = tcoord[idx * natomsgpu * SPACE + a1 * SPACE + Z] - tcoord[idx * natomsgpu * SPACE + a2 * SPACE + Z];

    #ifndef NOSQRT
                    r = clamp(hypotenuse(dx,dy,dz), RMIN_ELEC);
                    r2 = r*r;
                    int index = Ang_to_index(r);

    #else
                    r2 = sqhypotenuse(dx,dy,dz);
                    r2 = clamp(r2, RMIN_ELEC2);
                    int index = SqAng_to_index(r2);
    #endif

                    index_1t_NEINT = BoundedNeint(index);
                    index_1t_NDIEL = BoundedNdiel(index);

                    if ((Boole)pint_arraygpu[INTINCELEC * num_individualsgpu + idx])//B_calcIntElec[idx])
                    {
                        float r_dielectric = strr_epsilon_fn[index_1t_NDIEL];
                        e_elec = q1q2 * r_dielectric;
                        e_internal = e_elec;

                    }

                    //if (r2 < nbc2)
                    if (r2 < nbc)
                    {
                        e_desolv = strsol_fn[index_1t_NEINT] * nb_desolv;
                        int myidx;
                        if (B_include_1_4_interactions != 0 && nonbond_type == 4)
                        {
                            myidx = index_1t_NEINT * ATOM_MAPS * ATOM_MAPS + t2 * ATOM_MAPS + t1;
                            if (myidx == NEINT * ATOM_MAPS * ATOM_MAPS)
                            {
//                                e_internal += scale_1_4[idx] * (stre_vdW_Hb[myidx-1] + e_desolv);
                                e_internal += pfloat_arraygpu[FLOATSCALE14 * num_individualsgpu + idx] * (stre_vdW_Hb[myidx-1] + e_desolv);
                            }
                            else
                            {
//                                e_internal += scale_1_4[idx] * (stre_vdW_Hb[myidx] + e_desolv);
                                e_internal += pfloat_arraygpu[FLOATSCALE14 * num_individualsgpu + idx] * (stre_vdW_Hb[myidx] + e_desolv);
                            }
                        } else {
//                            fprintf(stderr," stre_vdW_Hb[%d][%d][%d] = %f\n", index_1t_NEINT, t2, t1, stre_vdW_Hb[index_1t_NEINT * ATOM_MAPS * ATOM_MAPS + t2 * ATOM_MAPS + t1]);
//i                            e_internal += stre_vdW_Hb[index_1t_NEINT * ATOM_MAPS * ATOM_MAPS + t2 * ATOM_MAPS + t1] + e_desolv;
                            myidx = index_1t_NEINT * ATOM_MAPS * ATOM_MAPS + t2 * ATOM_MAPS + t1;
                            if (myidx == NEINT * ATOM_MAPS * ATOM_MAPS)
                            {
                                e_internal += stre_vdW_Hb[myidx-1] + e_desolv;
//                                fprintf(stderr,"NEINT = %d, index = %d, t2 = %d, t1 = %d\n", NEINT, index_1t_NEINT, t2, t1);

                            }
                            else
                            {
                                e_internal += stre_vdW_Hb[myidx] + e_desolv;
                            }

                        }



                    }
                    total_e_internal += e_internal;

    #ifdef EINTCALPRINT
            total_e_desolv  += e_desolv;
            total_e_elec    += e_elec;
            float dielectric = strepsilon_fn[index_1t_NDIEL];

            if ((Boole)pint_arraygpu[INTINCELEC * num_individualsgpu + idx])//B_calcIntElec[idx])
            {
                e_vdW_Hb = e_internal - e_desolv - e_elec;
                // print stuff
            } else {
                e_vdW_Hb = e_internal - e_desolv;
                // print stuff
            }

            total_e_vdW_Hb += e_vdW_Hb;

    #endif

                }

                if (nb_group == INTRA_LIGAND)
                {
                    Nb_group_energy[INTRA_LIGAND] = total_e_internal;
                } else if (nb_group == INTER) {
                    Nb_group_energy[INTER] = total_e_internal - Nb_group_energy[INTRA_LIGAND];
                } else if (nb_group == INTRA_RECEPTOR) {
                    Nb_group_energy[INTRA_RECEPTOR] = total_e_internal - Nb_group_energy[INTRA_LIGAND] - Nb_group_energy[INTER];
                }

            }

    #ifdef EINTCALPRINT
            if((Boole)pint_arraygpu[INTINCELEC * num_individualsgpu + idx])//B_calcIntElec[idx])
            {
                //print stuff
            } else {
                //print stuff
            }
            //print stuff
    #endif
            if(b_comp_intermolgpu)
            {
                //energiesgpu[idx] += ((float)total_e_internal - (float)unboundinternalFEs[idx]);
                penergiesgpu[idx] += ((float)total_e_internal - pfloat_arraygpu[FLOATUNBOUNDINTERNAL * num_individualsgpu + idx]);//(float)unboundinternalFEs[idx]);
            }
            else
            {
                //energiesgpu[idx] = ((float)total_e_internal - (float)unboundinternalFEs[idx]);
                penergiesgpu[idx] = ((float)total_e_internal - pfloat_arraygpu[FLOATUNBOUNDINTERNAL * num_individualsgpu + idx]);//(float)unboundinternalFEs[idx]);
            }

        }

    }

}


////////////////////////////////////////////////////////////////////////////////
//! Entry point for Cuda function
//! @param b_comp_intermol
//! @param crds
//! @param energiescpu array of energies for returning to cpu
//! @param float_arraycpu
//! @param int_arraycpu
////////////////////////////////////////////////////////////////////////////////
extern "C" void eval_wrapper(
                       Boole b_comp_intermol,
                       float *crds,
                       float *energiescpu,
                       float *float_arraycpu,
                       int *int_arraycpu
                        )

{
    cudaMemcpy(crdsgpu, crds, sizeof(float) * cpunatoms * SPACE * num_individuals, cudaMemcpyHostToDevice);
    CHECK_ERROR(19);

    cudaMemcpy(int_arraygpu, int_arraycpu, sizeof(int) * INTSIZE * num_individuals, cudaMemcpyHostToDevice);
    CHECK_ERROR(18);
    cudaMemcpy(float_arraygpu, float_arraycpu, sizeof(float) * FLOATSIZE * num_individuals, cudaMemcpyHostToDevice);
    CHECK_ERROR(17);

    // execute trilinterp kernel
    eval_tril_kernel<<< nBlocks, blocksize >>>(num_individuals,
                energiesgpu,
                b_comp_intermol, cpunatoms,
                crdsgpu, chargesgpu, ABSchargesgpu, typesgpu,
                ignore_intersgpu,
                (float *)NULL_ELEC_TOTAL,
                (float *)NULL_EVDW_TOTAL,
                ElecMap,
                DesolvMap, SOME_ATOMS_OUTSIDE_GRID,
                ALL_ATOMS_INSIDE_GRID,
                float_arraygpu,
                int_arraygpu);

    CHECK_ERROR(666);
    cudaThreadSynchronize();
    CHECK_ERROR(888);

    // execute eintcal kernel
    eintcal_kernel<<< nBlocks, blocksize >>>(
        num_individuals,
        cpunatoms,
        energiesgpu,
        nonbondlistsgpu,
        crdsgpu,
        inc14interactgpu,
        haveflexresiduesgpu,
        nnb_arraygpu,
        nb_group_energygpu,
        evdWHbgpu,
        solfngpu,
        epsilonfngpu,
        repsilonfngpu,
        b_comp_intermol,
        float_arraygpu,
        int_arraygpu);

    CHECK_ERROR(888);
    cudaThreadSynchronize();
    CHECK_ERROR(23);

    cudaMemcpy(energiescpu, energiesgpu, sizeof(float) * num_individuals, cudaMemcpyDeviceToHost);
    CHECK_ERROR(33);

}


/**
 * Helper function to check if cuda had an error, and outputs the error
 * @param num error number (user defined)
 */
void CHECK_ERROR(int num)
{
    cudaError_t kerr;

    kerr = cudaGetLastError();

    if (kerr != cudaSuccess)
    {
        fprintf(stderr, "################################\nCUDA ERROR %d = %s\n################################\n", num, cudaGetErrorString(kerr));
    }

}

/**
 * Allocates memory on the graphics card and copies consistent variables for later use in the gen alg, also
 * store the map into texture memory for fast access.
 * @param natom number of atoms
 * @param nnum_individuals number of individuals in population
 * @param map the map representation
 * @param nonbondlist nonbonded list
 * @param etbl energie table
 * @param charge charges array
 * @param ABScharge ABScharge array
 * @param type type array
 * @param ignore_inter ingnore_inter array
 * @param inc14intereact inc14interact flag
 * @param haveflexresidues haveflexresidues flag
 */
extern "C" void cuda_alloc_wrapper(int natom, int nnum_individuals, maptype map, NonbondParam * nonbondlist, EnergyTables *etbl, Real *charge, Real *ABScharge, int *type, int *ignore_inter, Boole inc14interact, Boole haveflexresidues, Real *entable_solfn, Real *entable_epsilon_fn, Real *entable_r_epsilon_fn,Real e_vdW_Hb[NEINT][ATOM_MAPS][ATOM_MAPS])
{
    int memsize = 0;
    int i,j,k,l;

    cpunatoms = natom;
    num_individuals = nnum_individuals;

    nb_group_energycpu = (float *)malloc(sizeof(float) * 3);
    float *cudamapcpu = (float *)malloc(sizeof(float) *MAX_GRID_PTS * MAX_GRID_PTS * MAX_GRID_PTS * MAX_MAPS);
    //printf("estou em cuda_alloc_wrapper 2\n");
    for (i = 0; i < MAX_GRID_PTS; i++)
    {
        //printf("estou no for 1\n");
        for (j = 0; j < MAX_GRID_PTS; j++)
        {
            //printf("estou no for 2\n");
            for (k = 0; k < MAX_GRID_PTS; k++)
            {
                //printf("estou no for 3\n");
                for (l = 0; l < MAX_MAPS; l++)
                {
                    //printf("estou no for 4\n");
                    cudamapcpu[i * MAX_GRID_PTS * MAX_GRID_PTS * MAX_MAPS + j * MAX_GRID_PTS * MAX_MAPS + k * MAX_MAPS + l] = map[i][j][k][l];
                }
            }
        }
    }

    cudaMalloc((void**)&cudamap,  sizeof(float) * MAX_GRID_PTS * MAX_GRID_PTS * MAX_GRID_PTS * MAX_MAPS);
    CHECK_ERROR(75);

    memsize += sizeof(float) * MAX_GRID_PTS * MAX_GRID_PTS * MAX_GRID_PTS * MAX_MAPS;

    cudaMemcpy(cudamap, cudamapcpu, sizeof(float) * MAX_GRID_PTS * MAX_GRID_PTS * MAX_GRID_PTS * MAX_MAPS, cudaMemcpyHostToDevice);
    CHECK_ERROR(76);

    tex.addressMode[0] = cudaAddressModeWrap;
    tex.filterMode = cudaFilterModePoint;
    tex.normalized = false;
    cudaBindTexture(0, tex, cudamap);
    free(cudamapcpu);
    CHECK_ERROR(77);

    cudaMalloc((void**)&crdsgpu, sizeof(float) * cpunatoms * SPACE * num_individuals);
    CHECK_ERROR(5);

    memsize+= sizeof(float) * cpunatoms * SPACE * num_individuals;

    //get gpu memory for first kernel
    cudaMalloc((void**)&evalflagsgpu,sizeof(unsigned int) * num_individuals);
    CHECK_ERROR(0);
    cudaMalloc((void**)&energiesgpu, sizeof(float) * num_individuals);
    CHECK_ERROR(1);

    float *energiescpus = (float *)malloc(num_individuals * sizeof(float));
    bzero(energiescpus, sizeof(float) * num_individuals);

    cudaMemcpy(energiesgpu, energiescpus, sizeof(float) * num_individuals, cudaMemcpyHostToDevice);

    free(energiescpus);

    memsize += sizeof(unsigned int) * num_individuals + sizeof(float) * num_individuals;

    int nb_group_max = 3;
    int nb_group = 0;

    int largest = 0;
    for (nb_group = 0; nb_group < nb_group_max; nb_group++)
    {
        if (largest < Nnb_array[nb_group])
        {
            largest = Nnb_array[nb_group];
        }
    }

    float *nonbondlists = (float *)malloc(sizeof(float) * 7 * largest);

    for (k = 0; k < largest; k++)
    {
        nonbondlists[k * 7 + 0] = (float)nonbondlist[k].a1;
        nonbondlists[k * 7 + 1] = (float)nonbondlist[k].a2;
        nonbondlists[k * 7 + 2] = (float)nonbondlist[k].t1;
        nonbondlists[k * 7 + 3] = (float)nonbondlist[k].t2;
        nonbondlists[k * 7 + 4] = (float)nonbondlist[k].nonbond_type;
        nonbondlists[k * 7 + 5] = (float)nonbondlist[k].desolv;
        nonbondlists[k * 7 + 6] = (float)nonbondlist[k].q1q2;
    }

    cudaMalloc((void**)&nonbondlistsgpu, sizeof(float) * 7 * largest);
    CHECK_ERROR(34);

    memsize+= sizeof(float) * 7 * largest;

    cudaMemcpy(nonbondlistsgpu, nonbondlists, sizeof(float) * 7 * largest, cudaMemcpyHostToDevice);
    CHECK_ERROR(55);

    free(nonbondlists);

    blocksize = BLOCK_SIZE;
    nBlocks = num_individuals/BLOCK_SIZE + ((num_individuals%BLOCK_SIZE==0)?0:1);
    // tril params
    cudaMalloc((void**)&B_outsidesgpu, sizeof(int) * num_individuals);
    CHECK_ERROR(2);
    CHECK_ERROR(3);
    CHECK_ERROR(4);
    cudaMalloc((void**)&chargesgpu, sizeof(float) * cpunatoms);
    CHECK_ERROR(6);
    cudaMalloc((void**)&ABSchargesgpu, sizeof(float) * cpunatoms);
    CHECK_ERROR(7);
    cudaMalloc((void**)&typesgpu,sizeof(int) * cpunatoms  * num_individuals);
    CHECK_ERROR(8);
    CHECK_ERROR(9);
    cudaMalloc((void**)&gridinfosgpu, sizeof(float) * 10);
    CHECK_ERROR(10);
    cudaMalloc((void**)&ignore_intersgpu, sizeof(int) * cpunatoms);
    CHECK_ERROR(11);

    memsize += sizeof(float) * cpunatoms * num_individuals;
    memsize += sizeof(float) * cpunatoms * num_individuals;
    memsize += sizeof(int) * cpunatoms * num_individuals + sizeof(float) * 10 + sizeof(int) * cpunatoms * num_individuals;

    //eintcal params
    cudaMalloc((void**)&incelecgpu, sizeof(Boole) * num_individuals);
    CHECK_ERROR(36);
    CHECK_ERROR(37);
    cudaMalloc((void**)&scale14sgpu, sizeof(float) * num_individuals);
    CHECK_ERROR(38);
    CHECK_ERROR(39);
    cudaMalloc((void**)&usenonbondcutsgpu,sizeof(Boole) * num_individuals);
    CHECK_ERROR(40);
    CHECK_ERROR(41);
    cudaMalloc((void**)&unboundinternalFEsgpu,sizeof(float) * num_individuals);
    CHECK_ERROR(42);
    cudaMalloc((void**)&nnb_arraygpu, sizeof(int) * 3);
    cudaMalloc((void **)&nb_group_energygpu, sizeof(float) * 3);

    memsize += sizeof(Boole) * num_individuals + sizeof(float) * num_individuals + sizeof(Boole) * num_individuals + sizeof(float) * num_individuals + sizeof(int) * 3 + sizeof(float) * 3;

    cudaMemcpy(nnb_arraygpu, Nnb_array, sizeof(int) * 3, cudaMemcpyHostToDevice);

    float *chargescpu = (float *)malloc(sizeof(float) * cpunatoms);
    float *ABSchargescpu = (float *)malloc(sizeof(float) * cpunatoms);
    for (k = 0; k < cpunatoms; k++)
    {
        chargescpu[k] = (float)charge[k];
        ABSchargescpu[k] = (float)ABScharge[k];
    }

    cudaMemcpy(ignore_intersgpu, ignore_inter, sizeof(int) * cpunatoms, cudaMemcpyHostToDevice);
    CHECK_ERROR(12);

    cudaMemcpy(typesgpu, type, sizeof(int) * cpunatoms, cudaMemcpyHostToDevice);
    CHECK_ERROR(20);
    cudaMemcpy(chargesgpu, chargescpu, sizeof(float) * cpunatoms, cudaMemcpyHostToDevice);
    CHECK_ERROR(16);
    cudaMemcpy(ABSchargesgpu, ABSchargescpu, sizeof(float) * cpunatoms, cudaMemcpyHostToDevice);
    CHECK_ERROR(15);
    free(chargescpu);
    haveflexresiduesgpu = haveflexresidues;
    inc14interactgpu = inc14interact;
    int m;

    float *evdWHb;
    float *solfn;
    float *epsilonfn;
    float *repsilonfn;

    evdWHb = (float *) malloc(sizeof(float) * NEINT * ATOM_MAPS * ATOM_MAPS);
    solfn = (float *) malloc(sizeof(float) * NEINT);
    epsilonfn = (float *)malloc(sizeof(float) * NDIEL);
    repsilonfn = (float *)malloc(sizeof(float) * NDIEL);

    for (j = 0; j < NEINT; j++)
    {
        //solfn[j] = (float) etbl->sol_fn[j];
        solfn[j] = (float) entable_solfn[j];
        for (k = 0; k < ATOM_MAPS; k++)
        {
            for(m = 0; m < ATOM_MAPS; m++)
            {
                //printf("evdw %d - %d - %d: %f\n\n\n", j, k, m, e_vdW_Hb[j][k][m]);
                //evdWHb[j * ATOM_MAPS * ATOM_MAPS + k * ATOM_MAPS + m] = (float)etbl->e_vdW_Hb[j][k][m];
                evdWHb[j * ATOM_MAPS * ATOM_MAPS + k * ATOM_MAPS + m] = (float) e_vdW_Hb[j][k][m];
//                fprintf(stderr,"evdW_Hb[%d][%d][%d] = %f\n", j,k,m,::evaluate.getETblevdWHb(j,k,m));
            }
        }
    }

    for (k = 0; k < NDIEL; k++)
    {
        // epsilonfn[k] = (float)etbl->epsilon_fn[k];
        // repsilonfn[k] = (float)etbl->r_epsilon_fn[k];
        epsilonfn[k] = (float) entable_epsilon_fn[k];
        repsilonfn[k] = (float) entable_r_epsilon_fn[k];
    }

    cudaMalloc((void**)&evdWHbgpu, sizeof(float) * NEINT * ATOM_MAPS * ATOM_MAPS);
    CHECK_ERROR(43);
    cudaMalloc((void**)&solfngpu,sizeof(float) * NEINT);
    CHECK_ERROR(44);
    cudaMalloc((void**)&epsilonfngpu,sizeof(float) * NDIEL);
    CHECK_ERROR(45);
    cudaMalloc((void**)&repsilonfngpu,sizeof(float) *NDIEL);
    CHECK_ERROR(46);

    memsize+= sizeof(float) * NEINT * ATOM_MAPS * ATOM_MAPS + sizeof(float) * NEINT + sizeof(float) * NDIEL + sizeof(float) * NDIEL;

    cudaMemcpy(evdWHbgpu, evdWHb, sizeof(float) * NEINT * ATOM_MAPS * ATOM_MAPS, cudaMemcpyHostToDevice);
    cudaMemcpy(solfngpu, solfn, sizeof(float) * NEINT, cudaMemcpyHostToDevice);
    CHECK_ERROR(56);
    cudaMemcpy(epsilonfngpu, epsilonfn, sizeof(float) * NDIEL, cudaMemcpyHostToDevice);
    CHECK_ERROR(57);
    cudaMemcpy(repsilonfngpu, repsilonfn, sizeof(float) * NDIEL, cudaMemcpyHostToDevice);
    CHECK_ERROR(58);

    nb_group_energycpu[0] = (float)nb_group_energy[0];
    nb_group_energycpu[1] = (float)nb_group_energy[1];
    nb_group_energycpu[2] = (float)nb_group_energy[2];

    cudaMemcpy(nb_group_energygpu, nb_group_energycpu, sizeof(float) * 3, cudaMemcpyHostToDevice);

    CHECK_ERROR(54);

    free(evdWHb);
    free(solfn);
    free(epsilonfn);
    free(repsilonfn);

    cudaMalloc((void**)&int_arraygpu, sizeof(int) * INTSIZE * num_individuals);
    cudaMalloc((void**)&float_arraygpu, sizeof(float) * FLOATSIZE * num_individuals);


}

/**
 * Free's graphics memory
 */
extern "C" void cuda_free_wrapper(void)
{
    cudaFree(int_arraygpu);
    cudaFree(float_arraygpu);
    //#pragma omp barrier
    //{
      //free(nb_group_energycpu);
    //}
    cudaFree(cudamap);
    CHECK_ERROR(78);
    cudaFree(crdsgpu);
    CHECK_ERROR(79);

    cudaFree(nonbondlistsgpu);
    CHECK_ERROR(59);
    cudaFree(evdWHbgpu);
    CHECK_ERROR(68);
    cudaFree(solfngpu);
    CHECK_ERROR(69);
    cudaFree(epsilonfngpu);
    CHECK_ERROR(70);
    cudaFree(repsilonfngpu);
    CHECK_ERROR(71);

    //Free memory for trilinterp params
    cudaFree(B_outsidesgpu);
    CHECK_ERROR(24);
    cudaFree(chargesgpu);
    CHECK_ERROR(25);
    cudaFree(ABSchargesgpu);
    CHECK_ERROR(26);
    cudaFree(typesgpu);
    CHECK_ERROR(27);
    CHECK_ERROR(28);
    cudaFree(gridinfosgpu);
    CHECK_ERROR(29);
    cudaFree(ignore_intersgpu);
    CHECK_ERROR(30);
    CHECK_ERROR(32);

    //free eintcal
    cudaFree(incelecgpu);
    CHECK_ERROR(61);
    CHECK_ERROR(62);
    cudaFree(scale14sgpu);
    CHECK_ERROR(63);
    CHECK_ERROR(64);
    cudaFree(usenonbondcutsgpu);
    CHECK_ERROR(65);
    CHECK_ERROR(66);
    cudaFree(unboundinternalFEsgpu);
    CHECK_ERROR(67);
    cudaFree(nnb_arraygpu);
    cudaFree(nb_group_energygpu);

    //Free others
    cudaFree(evalflagsgpu);
    CHECK_ERROR(80);
    cudaFree(energiesgpu);
    CHECK_ERROR(81);
}

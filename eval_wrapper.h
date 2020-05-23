//eval_wrapper.h
extern "C" void eval_wrapper(Boole b_comp_intermol, float *crds, float *energiescpu, float *float_arraycpu, int *int_arraycpu);

extern "C" void cuda_alloc_wrapper(int natom, int nnum_individuals, maptype map, NonbondParam * nonbondlist, EnergyTables *etbl, Real *charge, Real *ABScharge, int *type, int *ignore_inter, Boole inc14interact, Boole haveflexresidues, Real *entable_solfn, Real *entable_epsilon_fn, Real *entable_r_epsilon_fn,Real e_vdW_Hb[NEINT][ATOM_MAPS][ATOM_MAPS]);

void CHECK_ERROR(int num);

extern "C" void cuda_free_wrapper(void);

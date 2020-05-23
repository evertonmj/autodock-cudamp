//eMapPrepCuda.h

int eEqualPrepCuda(Population &original_population, int e_mode, double firstEnergy);

void eMapPrepCuda(Population &original_population, int e_mode);
/*
void Eval::evalPrep(unsigned int pop_size,
                    double * penergies,
                        int indv,
                        float *pcrds,
                        float *pfloat_arraycpu,
                        int *pint_arraycpu);*/

void cpu_alloc(int num_individuals, int natoms);

void cpu_free(void);

void e_prep(int e_mode, int i, Population &original_population);

void e_complete_eval(Population &original_population, Real *charges);

int e_complete_evalEqual(Population &original_population, Real *charges, double firstEnergy);

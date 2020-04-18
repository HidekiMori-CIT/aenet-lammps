
#ifdef PAIR_CLASS

PairStyle(aenet,PairAENET)

#else

#ifndef LMP_PAIR_AENET_H
#define LMP_PAIR_AENET_H

#include "pair.h"

extern "C"{
void aenet_init(int ntypes, char* atom_types[], int* stat);
void aenet_final(int* stat);
void aenet_load_potential(int type_id, char* filename, int* stat);

void aenet_atomic_energy_and_forces(double coo_i[3], int type_i, int index_i,
                                    int n_j, double coo_j[], int type_j[],
                                    int index_j[], int natoms, double* E_i,
                                    double F[], int* stat);
}

extern int aenet_nnb_max;
extern double aenet_Rc_max;

extern int aenet_sfb_ver;

namespace LAMMPS_NS {

class PairAENET : public Pair {
 public:
  PairAENET(class LAMMPS *);
  ~PairAENET();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);

 protected:
  virtual void allocate();
  double cutmax;                // max cutoff for all elements
  int nelements;                // # of unique elements
  char **elements;              // names of unique elements

  int *map;                     // mapping from atom types to elements

};

}

#endif
#endif

/* ERROR/WARNING messages:


*/

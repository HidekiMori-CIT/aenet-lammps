
#ifdef PAIR_CLASS

PairStyle(aenet/zbl,PairAENET_ZBL)

#else

#ifndef LMP_PAIR_AENET_ZBL_H
#define LMP_PAIR_AENET_ZBL_H

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

class PairAENET_ZBL : public Pair {
 public:
  PairAENET_ZBL(class LAMMPS *);
  ~PairAENET_ZBL();
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
  
  double cut_global, cut_inner,cut_range;
  double cut_globalsq, cut_innersq;
  double ialpha;
  double E0;
  int dw_off;
  double *z;
  double **d1a,**d2a,**d3a,**d4a,**zze;
  double **sw1,**sw2,**sw3,**sw4,**sw5;

  double e_zbl(double, int, int);
  double dzbldr(double, int, int);
  double d2zbldr2(double, int, int);
  void set_coeff(int, int, double, double);
  
  double weight(double);
  double dweight(double);

};

namespace PairAENET_ZBLConstants {

  // ZBL constants

  static const double pzbl = 0.23;
  static const double a0 = 0.46850;
  static const double c1 = 0.02817;
  static const double c2 = 0.28022;
  static const double c3 = 0.50986;
  static const double c4 = 0.18175;
  static const double d1 = 0.20162;
  static const double d2 = 0.40290;
  static const double d3 = 0.94229;
  static const double d4 = 3.19980;
}

}

#endif
#endif

/* ERROR/WARNING messages:


*/

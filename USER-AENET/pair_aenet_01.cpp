
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_aenet.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

PairAENET::PairAENET(LAMMPS *lmp) : Pair(lmp)
{

  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;

  nelements = 0;
  elements = NULL;

}

/* ----------------------------------------------------------------------
   free all arrays
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairAENET::~PairAENET()
{
  //if(copymode) return;
  int i, stat; 

  aenet_final(&stat);
  if (stat != 0)error->all(FLERR,"Error: aenet finalization failed");
  
  for (i = 0; i < nelements; i++) delete [] elements[i];
  delete [] elements;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(map);
  }
}

/* ---------------------------------------------------------------------- */

void PairAENET::compute(int eflag, int vflag)
{
  int i,j,k,ii,jj,kk,inum,jnum,jtmp;
  int *ilist,*jlist,*numneigh,**firstneigh;
  
  ev_init(eflag,vflag);
  //if (eflag || vflag) ev_setup(eflag,vflag);
  //else evflag = vflag_fdotr = vflag_atom = 0;
  
  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  
  double E_i;
  int itype;
  int index_i = 1;
  double x_i[3];
  
  
  int njmax = 1;
  for (ii = 0; ii < inum; ii++){
     if(numneigh[ilist[ii]] > njmax)njmax = numneigh[ilist[ii]];
  }
  
  int *jtype = new int[njmax];
  int *index_j = new int[njmax];
  for (jj = 0; jj < njmax; jj++)index_j[jj] = jj + 2;
  double *x_j = new double[njmax*3];
  
  int nf_num = njmax+1;
  int nfc_num = nf_num*3;
  double *f_ann = new double[nfc_num];
  
  
  int stat;
  for (ii = 0; ii < inum; ii++) {
    
    E_i = 0.0;
    for (kk = 0; kk < nfc_num; kk++)f_ann[kk] = 0.0;
    
    i = ilist[ii];    
    itype = map[type[i]]+1;
    for (kk = 0; kk < 3; kk++)x_i[kk] = x[i][kk];
    
    jlist = firstneigh[i];
    jnum = numneigh[i];
    for (jj = 0; jj < jnum; jj++) {
      
      j = jlist[jj];
      j &= NEIGHMASK;
      
      jtype[jj] = map[type[j]]+1;
      
      jtmp = 3*jj;
      for (kk = 0; kk < 3; kk++)x_j[jtmp + kk] = x[j][kk];
      
    }
    
    aenet_nnb_max = jnum;
    aenet_atomic_energy_and_forces(&x_i[0], itype, index_i, jnum, &x_j[0], &jtype[0], &index_j[0], nf_num, &E_i, &f_ann[0], &stat);
    
    if (eflag) {
      if (eflag_global) eng_vdwl += E_i;
      if (eflag_atom) eatom[i] += E_i;
    }
    
    for (kk = 0; kk < 3; kk++)f[i][kk] += f_ann[kk];
    
    for (jj = 0; jj < jnum; jj++) {
      
      j = jlist[jj];
      j &= NEIGHMASK;
      
      jtmp = 3*(jj+1);
      for (kk = 0; kk < 3; kk++)f[j][kk] += f_ann[jtmp + kk];
      
    }
    
  }
  
  if (vflag_fdotr) virial_fdotr_compute();
  
  delete [] jtype;
  delete [] index_j;
  delete [] x_j;
  delete [] f_ann;

}

/* ---------------------------------------------------------------------- */

void PairAENET::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(map,n+1,"pair:map");
  
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairAENET::settings(int narg, char **/*arg*/)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairAENET::coeff(int narg, char **arg)
{
  int i,j,m,n;
  int stat;
  int itype;
  char netFile[128];
  
  std::string wc = "**";
  int len = wc.length();

  if (narg < 6) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // read ANN element names between 2 filenames
  // nelements = # of ANN elements
  // elements = list of unique element names

  if (nelements) {
    for (i = 0; i < nelements; i++) delete [] elements[i];
    delete [] elements;
  }
  nelements = narg - 4 - atom->ntypes;
  if (nelements < 1) error->all(FLERR,"Incorrect args for pair coefficients");
  elements = new char*[nelements];

  for (i = 0; i < nelements; i++) {
    n = strlen(arg[i+3]) + 1;
    elements[i] = new char[n];
    strcpy(elements[i],arg[i+3]);
  }
  
  aenet_init(nelements, &elements[0], &stat);
  if(stat != 0)error->all(FLERR,"Error: aenet initialization failed");
  
  aenet_sfb_ver = atoi((std::string(arg[2])).substr(1).c_str());
  std::string file_mask = std::string(arg[3+nelements]);
  for (i = 0; i < nelements; i++){
    itype = i + 1;
    std::string ele_str = std::string(elements[i]);
    std::string ith_file_name = file_mask;
    if (ith_file_name.find(wc) != std::string::npos) {
      ith_file_name.replace(ith_file_name.find(wc),len,ele_str);
    } else 
      ith_file_name = ele_str+"."+ith_file_name;
    snprintf(netFile, 128, "%-s", ith_file_name.c_str());
    aenet_load_potential(itype, netFile, &stat);
    if (stat != 0)error->all(FLERR,"Error: could not load ANN potentials");
  }

  // read args that map atom types to ANN elements
  // map[i] = which element the Ith atom type is, -1 if not mapped

  for (i = 4 + nelements; i < narg; i++) {
    m = i - (4+nelements) + 1;
    for (j = 0; j < nelements; j++)
      if (strcmp(arg[i],elements[j]) == 0) break;
    if (j < nelements) map[m] = j;
    else if (strcmp(arg[i],"NULL") == 0) map[m] = -1;
    else error->all(FLERR,"Incorrect args for pair coefficients");
  }

  // clear setflag since coeff() called once with I,J = * *

  n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements

  int count = 0;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        count++;
      }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");

}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairAENET::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style AENET requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style AENET requires newton pair on");

  // need a full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  
}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairAENET::init_one(int /*i*/, int /*j*/)
{
  return aenet_Rc_max;
}

/* ---------------------------------------------------------------------- */




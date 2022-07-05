
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include "pair_aenet_zbl.h"
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
using namespace PairAENET_ZBLConstants;

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

PairAENET_ZBL::PairAENET_ZBL(LAMMPS *lmp) : Pair(lmp)
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

PairAENET_ZBL::~PairAENET_ZBL()
{
  //if(copymode) return;
  int i, stat;
  char Error_message[256];

  aenet_final(&stat);
  if (stat != 0){
    snprintf(Error_message, 256, "Aenet finalization failed, error_ID from aenet (stat) is %2d",stat);
    error->all(FLERR,Error_message);
  }
  
  for (i = 0; i < nelements; i++) delete [] elements[i];
  delete [] elements;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(map);
    
    memory->destroy(z);
    memory->destroy(d1a);
    memory->destroy(d2a);
    memory->destroy(d3a);
    memory->destroy(d4a);
    memory->destroy(zze);
    memory->destroy(sw1);
    memory->destroy(sw2);
    memory->destroy(sw3);
    memory->destroy(sw4);
    memory->destroy(sw5);
  }
}

/* ---------------------------------------------------------------------- */

void PairAENET_ZBL::compute(int eflag, int vflag)
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
  
  double delx,dely,delz,evdwl,fpair;
  double rsq,r,t,fswitch,eswitch;
  int itz,jtz;
  
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
  
  double *f_zbl = new double[3];
  double *f_part = new double[nfc_num];
  double *f_fact = new double[nfc_num];
  
  double sum_part,sum_fact;
  double part,fact,fpart,ffact;
  double sigma,dwsi_f,dwsi_p;
  double inv_r,wi,dwi;
  double e0_zbl,de0_zbl,cf_zbl, dcf_zbl;
  
  double wi_zbl,wi_ann;
  
  int stat;
  for (ii = 0; ii < inum; ii++) {
    
    E_i = 0.0;
    for (kk = 0; kk < nfc_num; kk++)f_ann[kk] = 0.0;
    
    evdwl = 0.0;
    sum_part = 0.0;
    sum_fact = 0.0;
    for (kk = 0; kk < 3; kk++){f_zbl[kk]=0.0;}
    for (kk = 0; kk < nfc_num; kk++){f_part[kk]=0.0; f_fact[kk]=0.0;}
    
    i = ilist[ii];
    itype = map[type[i]]+1;
    itz = map[type[i]];
    for (kk = 0; kk < 3; kk++)x_i[kk] = x[i][kk];
    
    jlist = firstneigh[i];
    jnum = numneigh[i];
    
    if (jnum == 0)continue;
    
    for (jj = 0; jj < jnum; jj++) {
      
      j = jlist[jj];
      j &= NEIGHMASK;
      
      jtype[jj] = map[type[j]]+1;
      jtz = map[type[j]];
      
      jtmp = 3*jj;
      for (kk = 0; kk < 3; kk++)x_j[jtmp + kk] = x[j][kk];
      
      delx = x_i[0] - x[j][0];
      dely = x_i[1] - x[j][1];
      delz = x_i[2] - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      r = sqrt(rsq);
      inv_r = (-1.0/r);
      
      if (dw_off == 0){
        e0_zbl  = e_zbl(r, itz, jtz);
        de0_zbl = dzbldr(r, itz, jtz);
        cf_zbl  = weight(r);    
        dcf_zbl = dweight(r);
        evdwl += 0.5*e0_zbl*cf_zbl;
        fpair = (cf_zbl*de0_zbl+dcf_zbl*e0_zbl)*inv_r;
      } else{
        evdwl += 0.5*e_zbl(r, itz, jtz);
        fpair = dzbldr(r, itz, jtz)*inv_r;
      }
      
      f_zbl[0] += delx*fpair;
      f_zbl[1] += dely*fpair;
      f_zbl[2] += delz*fpair;
      
      part = exp(-r*ialpha);
      fact = r*part;
      fpart = -ialpha*part*inv_r;
      ffact = part*inv_r + r*fpart;
      
      f_part[0] += delx*fpart;
      f_part[1] += dely*fpart;
      f_part[2] += delz*fpart;
      
      f_fact[0] += delx*ffact;
      f_fact[1] += dely*ffact;
      f_fact[2] += delz*ffact;
      
      f_part[jtmp+3+0] = -delx*fpart;
      f_part[jtmp+3+1] = -dely*fpart;
      f_part[jtmp+3+2] = -delz*fpart;
      
      f_fact[jtmp+3+0] = -delx*ffact;
      f_fact[jtmp+3+1] = -dely*ffact;
      f_fact[jtmp+3+2] = -delz*ffact;
      
      sum_part += part;
      sum_fact += fact;
      
    }
    
    //std::cout <<"dwi= "<<dwi<<" by mori\n";
    
    aenet_nnb_max = jnum;
    aenet_atomic_energy_and_forces(&x_i[0], itype, index_i, jnum, &x_j[0], &jtype[0], &index_j[0], nf_num, &E_i, &f_ann[0], &stat);
    
    E_i -= E0;
    
    sigma = sum_fact/sum_part;
    wi    = weight(sigma);
    dwi   = dweight(sigma);
    
    if(dw_off == 0){
      dwsi_f = -dwi*(1.0/sum_part)*E_i;
      dwsi_p = -dwi*(-sum_fact/(sum_part*sum_part))*E_i;
      wi_zbl = 1.0;
      wi_ann = 1.0 - wi;
    } else {
      dwsi_f = dwi*(1.0/sum_part)*(evdwl-E_i);
      dwsi_p = dwi*(-sum_fact/(sum_part*sum_part))*(evdwl-E_i);
      wi_zbl = wi;
      wi_ann = 1.0 - wi;
    }
    
    for (kk = 0; kk < 3; kk++){
      f[i][kk] += wi_zbl*f_zbl[kk] + wi_ann*f_ann[kk] + dwsi_f*f_fact[kk] + dwsi_p*f_part[kk];
    }
            
    for (jj = 0; jj < jnum; jj++) {
      
      j = jlist[jj];
      j &= NEIGHMASK;
      
      jtmp = 3*(jj+1);
      for (kk = 0; kk < 3; kk++){
        int ktmp = jtmp + kk;
        f[j][kk] += wi_ann*f_ann[ktmp] + dwsi_f*f_fact[ktmp] + dwsi_p*f_part[ktmp];
      }
      
    }
    
    if (eflag) {
      if (eflag_global) eng_vdwl += wi_zbl*evdwl + wi_ann*E_i;
      if (eflag_atom)   eatom[i] += wi_zbl*evdwl + wi_ann*E_i;
    }
    
  }
  
  if (vflag_fdotr) virial_fdotr_compute();
  
  delete [] jtype;
  delete [] index_j;
  delete [] x_j;
  delete [] f_ann;
  delete [] f_zbl;
  delete [] f_part;
  delete [] f_fact;

}

/* ---------------------------------------------------------------------- */

void PairAENET_ZBL::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(map,n+1,"pair:map");

  memory->create(z,n+1,"pair:z");
  memory->create(d1a,n+1,n+1,"pair:d1a");
  memory->create(d2a,n+1,n+1,"pair:d2a");
  memory->create(d3a,n+1,n+1,"pair:d3a");
  memory->create(d4a,n+1,n+1,"pair:d4a");
  memory->create(zze,n+1,n+1,"pair:zze");
  memory->create(sw1,n+1,n+1,"pair:sw1");
  memory->create(sw2,n+1,n+1,"pair:sw2");
  memory->create(sw3,n+1,n+1,"pair:sw3");
  memory->create(sw4,n+1,n+1,"pair:sw4");
  memory->create(sw5,n+1,n+1,"pair:sw5");
  
  
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairAENET_ZBL::settings(int narg, char **arg)
{
  if (narg == 4 ) {
    cut_inner  = utils::numeric(FLERR,arg[0],false,lmp);
    cut_global = utils::numeric(FLERR,arg[1],false,lmp);
    ialpha     = utils::numeric(FLERR,arg[2],false,lmp);
    //E0 = -4475.52435072301;
    E0         = utils::numeric(FLERR,arg[3],false,lmp);
    dw_off = 0;
  } else if (narg == 5 ) {
    cut_inner  = utils::numeric(FLERR,arg[0],false,lmp);
    cut_global = utils::numeric(FLERR,arg[1],false,lmp);
    ialpha     = utils::numeric(FLERR,arg[2],false,lmp);
    E0         = utils::numeric(FLERR,arg[3],false,lmp);
    dw_off     = utils::numeric(FLERR,arg[4],false,lmp);
  } else error->all(FLERR,"Illegal pair_style command");

  if (cut_inner <= 0.0 )
    error->all(FLERR,"Illegal pair_style command");
  if (cut_inner > cut_global)
    error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairAENET_ZBL::coeff(int narg, char **arg)
{
  int i,j,m,n;
  int stat;
  int itype;
  char netFile[128];
  char Error_message[256];
  
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
  if (nelements%2 != 0) error->all(FLERR,"Incorrect args for pair coefficients");
  if (nelements < 1) error->all(FLERR,"Incorrect args for pair coefficients");
  
  nelements = nelements/2;
  elements = new char*[nelements];

  for (i = 0; i < nelements; i++) {
    n = strlen(arg[2*i+3]) + 1;
    elements[i] = new char[n];
    strcpy(elements[i],arg[2*i+3]);
    z[i] = utils::numeric(FLERR,arg[(2*i+1)+3],false,lmp);
  }
  
  aenet_init(nelements, &elements[0], &stat);
  if(stat != 0){
      snprintf(Error_message, 256, "Aenet initialization failed, error_ID from aenet (stat) is %2d",stat);
	  error->all(FLERR,Error_message);
  }
  
  aenet_sfb_ver = atoi((std::string(arg[2])).substr(1).c_str());
  std::string file_mask = std::string(arg[3+2*nelements]);
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
    if (stat != 0){
        snprintf(Error_message, 256, "Could not load %-s, error_ID from aenet (stat) is%2d", ith_file_name.c_str(),stat);
	    error->all(FLERR,Error_message);
    }
  }

  // read args that map atom types to ANN elements
  // map[i] = which element the Ith atom type is, -1 if not mapped

  for (i = 4 + nelements*2; i < narg; i++) {
    m = i - (4+nelements*2) + 1;
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
        set_coeff(map[i], map[j], z[map[i]], z[map[j]]);
        count++;
      }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");

}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairAENET_ZBL::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style AENET requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style AENET requires newton pair on");

  // need a full neighbor list

  //int irequest = neighbor->request(this,instance_me);
  //neighbor->requests[irequest]->half = 0;
  //neighbor->requests[irequest]->full = 1;
  // request a full neighbor list
  neighbor->add_request(this, NeighConst::REQ_FULL);

  
  cut_innersq = cut_inner * cut_inner;
  cut_globalsq = cut_global * cut_global;
  cut_range = cut_global - cut_inner;
  
}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairAENET_ZBL::init_one(int /*i*/, int /*j*/)
{
  return std::max(aenet_Rc_max,cut_global);
}

/* ----------------------------------------------------------------------
   compute ZBL pair energy
------------------------------------------------------------------------- */

double PairAENET_ZBL::e_zbl(double r, int i, int j) {

  double d1aij = d1a[i][j];
  double d2aij = d2a[i][j];
  double d3aij = d3a[i][j];
  double d4aij = d4a[i][j];
  double zzeij = zze[i][j];
  double rinv = 1.0/r;

  double sum = c1*exp(-d1aij*r);
  sum += c2*exp(-d2aij*r);
  sum += c3*exp(-d3aij*r);
  sum += c4*exp(-d4aij*r);

  double result = zzeij*sum*rinv;

  return result;
}


/* ----------------------------------------------------------------------
   compute ZBL first derivative
------------------------------------------------------------------------- */

double PairAENET_ZBL::dzbldr(double r, int i, int j) {

  double d1aij = d1a[i][j];
  double d2aij = d2a[i][j];
  double d3aij = d3a[i][j];
  double d4aij = d4a[i][j];
  double zzeij = zze[i][j];
  double rinv = 1.0/r;

  double e1 = exp(-d1aij*r);
  double e2 = exp(-d2aij*r);
  double e3 = exp(-d3aij*r);
  double e4 = exp(-d4aij*r);

  double sum = c1*e1;
  sum += c2*e2;
  sum += c3*e3;
  sum += c4*e4;

  double sum_p = -c1*d1aij*e1;
  sum_p -= c2*d2aij*e2;
  sum_p -= c3*d3aij*e3;
  sum_p -= c4*d4aij*e4;

  double result = zzeij*(sum_p - sum*rinv)*rinv;

  return result;
}

/* ----------------------------------------------------------------------
   compute ZBL second derivative
------------------------------------------------------------------------- */

double PairAENET_ZBL::d2zbldr2(double r, int i, int j) {

  double d1aij = d1a[i][j];
  double d2aij = d2a[i][j];
  double d3aij = d3a[i][j];
  double d4aij = d4a[i][j];
  double zzeij = zze[i][j];
  double rinv = 1.0/r;

  double e1 = exp(-d1aij*r);
  double e2 = exp(-d2aij*r);
  double e3 = exp(-d3aij*r);
  double e4 = exp(-d4aij*r);

  double sum = c1*e1;
  sum += c2*e2;
  sum += c3*e3;
  sum += c4*e4;

  double sum_p = c1*e1*d1aij;
  sum_p += c2*e2*d2aij;
  sum_p += c3*e3*d3aij;
  sum_p += c4*e4*d4aij;

  double sum_pp = c1*e1*d1aij*d1aij;
  sum_pp += c2*e2*d2aij*d2aij;
  sum_pp += c3*e3*d3aij*d3aij;
  sum_pp += c4*e4*d4aij*d4aij;

  double result = zzeij*(sum_pp + 2.0*sum_p*rinv +
                         2.0*sum*rinv*rinv)*rinv;

  return result;
}

/* ----------------------------------------------------------------------
   calculate the i,j entries in the various coeff arrays
------------------------------------------------------------------------- */

void PairAENET_ZBL::set_coeff(int i, int j, double zi, double zj)
{
  double ainv = (pow(zi,pzbl) + pow(zj,pzbl))/(a0*force->angstrom);
  d1a[i][j] = d1*ainv;
  d2a[i][j] = d2*ainv;
  d3a[i][j] = d3*ainv;
  d4a[i][j] = d4*ainv;
  zze[i][j] = zi*zj*force->qqr2e*force->qelectron*force->qelectron;
  //std::cout << "zze= " << force->qelectron << std::endl;

  if (i != j){
    d1a[j][i] = d1a[i][j];
    d2a[j][i] = d2a[i][j];
    d3a[j][i] = d3a[i][j];
    d4a[j][i] = d4a[i][j];
    zze[j][i] = zze[i][j];
  };

  // e =  t^3 (sw3 + sw4*t) + sw5
  //   = A/3*t^3 + B/4*t^4 + C
  // sw3 = A/3
  // sw4 = B/4
  // sw5 = C

  // dedr = t^2 (sw1 + sw2*t)
  //      = A*t^2 + B*t^3
  // sw1 = A
  // sw2 = B

  // de2dr2 = 2*A*t + 3*B*t^2

  // Require that at t = tc:
  // e = -Fc
  // dedr = -Fc'
  // d2edr2 = -Fc''

  // Hence:
  // A = (-3Fc' + tc*Fc'')/tc^2
  // B = ( 2Fc' - tc*Fc'')/tc^3
  // C = -Fc + tc/2*Fc' - tc^2/12*Fc''

  double tc = cut_global - cut_inner;
  double fc = e_zbl(cut_global, i, j);
  double fcp = dzbldr(cut_global, i, j);
  double fcpp = d2zbldr2(cut_global, i, j);

  double swa = (-3.0*fcp + tc*fcpp)/(tc*tc);
  double swb = ( 2.0*fcp - tc*fcpp)/(tc*tc*tc);
  double swc = -fc + (tc/2.0)*fcp - (tc*tc/12.0)*fcpp;

  sw1[i][j] = swa;
  sw2[i][j] = swb;
  sw3[i][j] = swa/3.0;
  sw4[i][j] = swb/4.0;
  sw5[i][j] = swc;

  if (i != j){
    sw1[j][i] = sw1[i][j];
    sw2[j][i] = sw2[i][j];
    sw3[j][i] = sw3[i][j];
    sw4[j][i] = sw4[i][j];
    sw5[j][i] = sw5[i][j];
  }
}

double PairAENET_ZBL::weight(double sigma) {

  double wi;  

  if (sigma > cut_global){
    wi = 0.0;
  }else if (sigma < cut_inner){
    wi = 1.0;
  }else {
    double ui  = (sigma - cut_inner)/cut_range;
    double ui3 = ui*ui*ui;
    wi = 1.0 - ui3*(6.0*ui*ui - 15.0*ui + 10.0);
  }
  return wi;
}

double PairAENET_ZBL::dweight(double sigma) {

  double dwi;

  if (sigma > cut_global){
    dwi = 0.0;
  }else if (sigma < cut_inner){
    dwi = 0.0;
  }else {
    double ui  = (sigma - cut_inner)/cut_range;
    double ui2 = ui*ui;
    dwi = -ui2*(30.0*ui*ui - 60.0*ui + 30.0)/cut_range;
  }
  return dwi;
}


/* ---------------------------------------------------------------------- */




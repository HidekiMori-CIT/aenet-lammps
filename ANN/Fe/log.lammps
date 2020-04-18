LAMMPS (3 Mar 2020)
atom_style      atomic
units           metal
boundary        p p p

lattice         bcc 2.9
Lattice spacing in x,y,z = 2.9 2.9 2.9
region          region1 block 0 1 0 1 0 1 units lattice
create_box      1 region1
Created orthogonal box = (0 0 0) to (2.9 2.9 2.9)
  1 by 1 by 1 MPI processor grid

create_atoms    1 box
Created 2 atoms
  create_atoms CPU = 0.00150323 secs

pair_style      aenet
pair_coeff      * * v01 Fe 10tw-10tw.ann Fe
mass            1 55.845

thermo_style    custom step etotal temp lx vol press pxx pyy pxy
thermo          10

fix             f1 all box/relax iso 0.0

minimize        1.0e-20 0.0 1000 100000
WARNING: Using 'neigh_modify every 1 delay 0 check yes' setting during minimization (../min.cpp:190)
Neighbor list info ...
  update every 1 steps, delay 0 steps, check yes
  max neighbors/atom: 2000, page size: 100000
  master list distance cutoff = 8.5
  ghost atom cutoff = 8.5
  binsize = 4.25, bins = 1 1 1
  1 neighbor lists, perpetual/occasional/extra = 1 0 0
  (1) pair aenet, perpetual
      attributes: full, newton on
      pair build: full/bin/atomonly
      stencil: full/bin/3d
      bin: standard
Per MPI rank memory allocation (min/avg/max) = 4.115 | 4.115 | 4.115 Mbytes
Step TotEng Temp Lx Volume Press Pxx Pyy Pxy 
       0   -8959.7099            0          2.9       24.389    -105938.5    -105938.5    -105938.5 -1.1170779e-10 
      10   -8959.7147            0       2.8971    24.315906   -103326.24   -103326.24   -103326.24 2.8675807e-10 
      20   -8959.7193            0       2.8942    24.242958   -100412.17   -100412.17   -100412.17 -2.4967524e-11 
      30   -8959.7238            0       2.8913    24.170157   -97221.184   -97221.184   -97221.184 -1.8323163e-10 
      40   -8959.7281            0       2.8884    24.097501    -93775.55    -93775.55    -93775.55 9.1970435e-11 
      50   -8959.7323            0       2.8855    24.024991   -90094.834   -90094.834   -90094.834 2.3903431e-11 
      60   -8959.7363            0       2.8826    23.952627   -86195.977   -86195.977   -86195.977 5.0953689e-11 
      70   -8959.7401            0       2.8797    23.880408   -82093.435   -82093.435   -82093.435 2.2342379e-10 
      80   -8959.7437            0       2.8768    23.808334    -77799.38    -77799.38    -77799.38 -9.1216019e-11 
      90   -8959.7471            0       2.8739    23.736406   -73323.964   -73323.964   -73323.964 -1.4515696e-10 
     100   -8959.7502            0        2.871    23.664622   -68675.612   -68675.612   -68675.612 -4.8865192e-11 
     110   -8959.7532            0       2.8681    23.592984   -63861.338   -63861.338   -63861.338 -4.1103071e-10 
     120    -8959.756            0       2.8652     23.52149   -58887.053   -58887.053   -58887.053 9.4558695e-11 
     130   -8959.7585            0       2.8623    23.450141    -53757.85    -53757.85    -53757.85 -3.8670426e-10 
     140   -8959.7607            0       2.8594    23.378936   -48478.252   -48478.252   -48478.252 -2.7204645e-10 
     150   -8959.7628            0       2.8565    23.307875   -43052.408   -43052.408   -43052.408 -6.1947457e-11 
     160   -8959.7645            0       2.8536    23.236959   -37484.224   -37484.224   -37484.224 4.6204687e-10 
     170   -8959.7661            0       2.8507    23.166186   -31777.453   -31777.453   -31777.453 -3.791169e-11 
     180   -8959.7673            0       2.8478    23.095558    -25935.73    -25935.73    -25935.73 -3.7660589e-10 
     190   -8959.7684            0       2.8449    23.025073   -19962.606   -19962.606   -19962.606 -4.7963715e-10 
     200   -8959.7691            0        2.842    22.954732   -13861.598   -13861.598   -13861.598 4.1832753e-11 
     210   -8959.7696            0       2.8391    22.884534   -7636.3009   -7636.3009   -7636.3009 1.5919096e-10 
     220   -8959.7698            0       2.8362    22.814479   -1290.6199   -1290.6199   -1290.6199 7.4312161e-11 
     224   -8959.7698            0    2.8356167    22.800405 -3.5037541e-10 1.8053536e-10 2.5647427e-10 1.9037485e-10 
Loop time of 1.58565 on 1 procs for 224 steps with 2 atoms

99.9% CPU use with 1 MPI tasks x no OpenMP threads

Minimization stats:
  Stopping criterion = energy tolerance
  Energy initial, next-to-last, final = 
        -8959.70990256     -8959.76977477     -8959.76977477
  Force two-norm initial, final = 4.83792 2.17542e-14
  Force max component initial, final = 4.83792 1.52981e-14
  Final line search alpha, max atom move = 1 1.52981e-14
  Iterations, force evaluations = 224 226

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Pair    | 1.5808     | 1.5808     | 1.5808     |   0.0 | 99.70
Neigh   | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.0010184  | 0.0010184  | 0.0010184  |   0.0 |  0.06
Output  | 0.00064756 | 0.00064756 | 0.00064756 |   0.0 |  0.04
Modify  | 0          | 0          | 0          |   0.0 |  0.00
Other   |            | 0.003142   |            |       |  0.20

Nlocal:    2 ave 2 max 2 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Nghost:    557 ave 557 max 557 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Neighs:    0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
FullNghs:  360 ave 360 max 360 min
Histogram: 1 0 0 0 0 0 0 0 0 0

Total # of neighbors = 360
Ave neighs/atom = 180
Neighbor list builds = 0
Dangerous builds = 0

Total wall time: 0:00:01

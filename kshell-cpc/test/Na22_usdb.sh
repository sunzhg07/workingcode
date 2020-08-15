#!/bin/sh 
# export OMP_STACKSIZE=1g
export GFORTRAN_UNBUFFERED_PRECONNECTED=y
# ulimit -s unlimited

# ---------- Na22_usdb --------------
echo "start running log_Na22_usdb_m0p.txt ..."
cat > Na22_usdb_0.input <<EOF
&input
  beta_cm = 0.0
  eff_charge = 1.5, 0.5, 
  fn_int = "usdb.snt"
  fn_ptn = "Na22_usdb_p.ptn"
  fn_save_wave = "Na22_usdb_m0p.wav"
  gl = 1.0, 0.0, 
  gs = 3.91, -2.678, 
  hw_type = 2
  is_double_j = .false.
  max_lanc_vec = 200
  maxiter = 300
  mode_lv_hdd = 0
  mtot = 0
  n_block = 0
  n_eigen = 10
  n_restart_vec = 15
&end
EOF
nice ./kshell.exe Na22_usdb_0.input > log_Na22_usdb_m0p.txt 2>&1 

rm -f tmp_snapshot_Na22_usdb_p.ptn_0_* tmp_lv_Na22_usdb_p.ptn_0_* Na22_usdb_0.input 


./collect_logs.py log_*Na22_usdb* > summary_Na22_usdb.txt
echo "Finish computing Na22_usdb.    See summary_Na22_usdb.txt"
echo 


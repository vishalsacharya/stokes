#!/bin/bash

export EXEC=bin/stokes

nodes+=(               1          1           1           1           1 );
cores+=(              16         16          16          16          16 );
mpi_proc+=(            1          1           1           1           1 );
threads+=(            16         16          16          16          16 );
test_case+=(           2          2           2           2           2 );
pt_cnt+=(              1          1           1           1           1 );
pt_rad+=(           0.15       0.15        0.15        0.15        0.15 );
jump_width+=(      1e-10      1e-10       1e-10       1e-10       1e-10 );
rho+=(              1e+7       1e+7        1e+7        1e+7        1e+7 );
ref_tol+=(         1e-10      1e-10       1e-10       1e-10       1e-10 );
min_depth+=(           1          1           1           1           1 );
max_depth+=(           6          6           7           7           8 );
fmm_q+=(              14         14           6           4           2 );
fmm_m+=(              10          6           6           6           6 );
gmres_tol+=(      1.0e-7       1e-7        1e-7        1e-7        1e-7 );
gmres_iter+=(        200        200         200         200         200 );
max_time+=(      3600000    3600000     3600000     3600000     3600000 );

# Export arrays
export      nodes_="$(declare -p      nodes)";
export      cores_="$(declare -p      cores)";
export   mpi_proc_="$(declare -p   mpi_proc)";
export    threads_="$(declare -p    threads)";
export  test_case_="$(declare -p  test_case)";
export     pt_cnt_="$(declare -p     pt_cnt)";
export     pt_rad_="$(declare -p     pt_rad)";
export jump_width_="$(declare -p jump_width)";
export        rho_="$(declare -p        rho)";
export    ref_tol_="$(declare -p    ref_tol)";
export  min_depth_="$(declare -p  min_depth)";
export  max_depth_="$(declare -p  max_depth)";
export      fmm_q_="$(declare -p      fmm_q)";
export      fmm_m_="$(declare -p      fmm_m)";
export  gmres_tol_="$(declare -p  gmres_tol)";
export gmres_iter_="$(declare -p gmres_iter)";
export   max_time_="$(declare -p   max_time)";

export RESULT_SUBDIR=$(basename ${0%.*});
export WORK_DIR=$(dirname ${PWD}/$0)/..
cd ${WORK_DIR}

TERM_WIDTH=$(stty size | cut -d ' ' -f 2)
./scripts/.submit_jobs.sh | cut -b -${TERM_WIDTH}




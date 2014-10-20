#!/bin/bash

export EXEC=bin/stokes

#POROUS MEDIA
nodes+=(             1        9       47      227     1024 );
cores+=(            16       16       16       16       16 );
mpi_proc+=(          1        9       47      227     1024 );
threads+=(          16       16       16       16       16 );
test_case+=(         3        3        3        3        3 );
pt_cnt+=(         1000     1000     1000     1000     1000 );
pt_rad+=(         0.10     0.10     0.10     0.10     0.10 );
jump_width+=(     1e-9     1e-9     1e-9     1e-9     1e-9 );
rho+=(            1e+9     1e+9     1e+9     1e+9     1e+9 );
ref_tol+=(        1e-9     1e-9     1e-9     1e-9     1e-9 );
min_depth+=(         1        1        1        1        1 );
max_depth+=(         5        6        7        8        9 );
fmm_q+=(            14       14       14       14       14 );
fmm_m+=(            10       10       10       10       10 );
gmres_tol+=(      1e-8     1e-8     1e-8     1e-8     1e-8 );
gmres_iter+=(      200      200      200      200      200 );
max_time+=(     360000   360000   360000   360000   360000 );

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


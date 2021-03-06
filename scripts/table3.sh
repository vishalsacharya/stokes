#!/bin/bash

export EXEC=bin/stokes

NODES=1;
THRDS=1;
nodes+=(        $NODES   $NODES   $NODES   $NODES   $NODES );
cores+=(        $THRDS   $THRDS   $THRDS   $THRDS   $THRDS );
mpi_proc+=(     $NODES   $NODES   $NODES   $NODES   $NODES );
threads+=(      $THRDS   $THRDS   $THRDS   $THRDS   $THRDS );
test_case+=(         0        0        0        0        0 );
pt_cnt+=(            0        0        0        0        0 );
pt_rad+=(            0        0        0        0        0 );
jump_width+=(        0        0        0        0        0 );
rho+=(            0e+0     0e+0     0e+0     0e+0     0e+0 ); #
ref_tol+=(        1e-0     1e-1     1e-2     1e-3     1e-4 ); #
min_depth+=(         1        1        1        1        1 );
max_depth+=(        15       15       15       15       15 );
fmm_q+=(             8       10       12       14       15 );
fmm_m+=(             2        4        6        8       10 );
gmres_tol+=(    1.0e-6   1.0e-6   1.0e-6   1.0e-6   1.0e-6 ); #
gmres_iter+=(      400      400      400      400      400 );
max_time+=(    3600000  3600000  3600000  3600000  3600000 );
                                                                                                                       
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


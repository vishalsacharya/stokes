static char help[] = "\n\
-test_case       0 - Constant coefficient\n\
                 1 - Synthetic variable coefficient\n\
                 2 - Spheres\n\
                 3 - Porous media\n\
\n\
-pt_cnt     <Int>    Number of spheres\n\
-pt_rad     <Real>   Sphere radius\n\
-jump_width <Real>   Jump width\n\
-rho        <Real>   Inf norm of \\rho\n\
-ref_tol    <Real>   Tree refinement tolerance\n\
-min_depth  <Int>    Minimum tree depth\n\
-max_depth  <Int>    Maximum tree depth\n\
-fmm_q      <Int>    Chebyshev polynomial degree\n\
-fmm_m      <Int>    Multipole order (+ve even integer)\n\
-gmres_tol  <Real>   GMRES residual tolerance\n\
-gmres_iter <Int>    GMRES maximum iterations\n\
";

#include <petscksp.h>
#include <cassert>
#include <cstring>

// PVFMM INCLUDES
#include <profile.hpp>
#include <fmm_cheb.hpp>
#include <fmm_node.hpp>
#include <fmm_tree.hpp>
#include <cheb_node.hpp>

// TBSLAS INCLUDES
#include <utils/common.h>
#include <utils/fields.h>
#include <tree/utils_tree.h>
#include <tree/semilag_tree.h>
#include <utils.hpp>

typedef pvfmm::FMM_Node<pvfmm::Cheb_Node<double> > FMMNode_t;
typedef pvfmm::FMM_Cheb<FMMNode_t> FMM_Mat_t;
typedef pvfmm::FMM_Tree<FMM_Mat_t> FMM_Tree_t;

double time_ksp;
int    iter_ksp;
size_t num_oct;
double BG_INF_ERR ;
double BG_L2_ERR  ;
double AVG_INF_ERR;
double AVG_L2_ERR ;
double SS_INF_ERR ;
double SS_L2_ERR  ;

// TEST_CASE: 0 - Constant coefficient
// TEST_CASE: 1 - Synthetic variable coefficient
// TEST_CASE: 2 - Sphere
// TEST_CASE: 3 - Porous media
PetscInt   TEST_CASE        = 0;

PetscInt   ERROR_RESOLUTION = 1;
PetscInt   VTK_ORDER        = 8;
PetscInt   INPUT_DOF        = 3;
PetscReal  SCAL_EXP         = 1.0;
PetscBool  PERIODIC         = PETSC_TRUE;
PetscBool  TREE_ONLY        = PETSC_FALSE;

PetscInt   MAXDEPTH         = MAX_DEPTH; // Maximum tree depth
PetscInt   MINDEPTH         = 4;         // Minimum tree depth
PetscReal  TOL              = 1e-3;      // Tolerance
PetscReal  GMRES_TOL        = 1e-6;      // Fine mesh GMRES tolerance

PetscInt   CHEB_DEG         = 14;        // Fine mesh Cheb. order
PetscInt   MUL_ORDER        = 10;        // Fine mesh mult  order

PetscInt   MAX_ITER         = 200;

PetscInt   PT_CNT           = 1;         // Number of spheres
PetscReal  PT_RAD           = 0.15;
PetscReal  JUMP_WIDTH       = 0.001;

PetscReal  f_max            = 0;
PetscReal  rho_             = 1000000;

PetscInt   NTSTEP = 1;

// TBSLAS GLOBAL
void (*fn_con)(const double* , int , double*)=NULL;
double tcurr = 0;


//PetscReal  L                = 500;
std::vector<double>  pt_coord;

template<typename real_t, int sdim>
void
get_gaussian_field_cylinder_atT(const real_t* points_pos,
                                int num_points,
                                real_t* out) {
  real_t xc      = 0.6;
  real_t yc      = 0.5;
  real_t r = sqrt((xc-0.5)*(xc-0.5) + (yc-0.5)*(yc-0.5));
  xc = 0.5+r*cos(tcurr);
  yc = 0.5+r*sin(tcurr);
  const real_t theta   = 0.0;
  const real_t sigma_x = 0.06;
  const real_t sigma_y = 0.06;
  const real_t A       = 1.0;

  tbslas::get_gaussian_field_cylinder<real_t, sdim>(points_pos,
                                                    num_points,
                                                    out,
                                                    xc,
                                                    yc,
                                                    theta,
                                                    sigma_x,
                                                    sigma_y,
                                                    A);
}

template<typename real_t, int sdim>
void
get_slotted_cylinder_atT(const real_t* points_pos,
                         int num_points,
                         real_t* out) {
  real_t xc = 0.5;
  real_t yc = 0.5;
  real_t zc = 0.5;
  real_t R  = 0.3;
  real_t w  = 0.1;
  real_t a  = tcurr;
  tbslas::get_slotted_cylinder<real_t, sdim>(points_pos,
                                             num_points,
                                             out,
                                             xc,
                                             yc,
                                             zc,
                                             R,
                                             w,
                                             a);
}

template<typename real_t>
void get_gaussian_position (real_t xi,  real_t xf,
                            real_t yi,  real_t yf,
                            real_t zi,  real_t zf,
                            real_t &xc, real_t &yc, real_t &zc)
{

  for(size_t i=0;i<PT_CNT;i++){
    if (pt_coord[3*i  ] > xi && pt_coord[3*i  ] < xf  &&
        pt_coord[3*i+1] > yi && pt_coord[3*i+1] < yf  &&
        pt_coord[3*i+2] > zi && pt_coord[3*i+2] < zf){

      xc = pt_coord[3*i  ];
      yc = pt_coord[3*i+1];
      zc = pt_coord[3*i+2];
      break;

    }
  }
}

template<typename real_t, int sdim>
void
get_gaussian_field_3d_wrapper_01(const real_t* points_pos,
                                 const int num_points,
                                 real_t* out) {

  real_t xc;
  real_t yc;
  real_t zc;
  real_t A       = 1.0 ;
  real_t sigma_x = 0.06;
  real_t sigma_y = 0.06;
  real_t sigma_z = 0.06;

  get_gaussian_position <real_t>(0.6, 0.7, 0.4, 0.6, 0.4, 0.6, xc, yc, zc);

  tbslas::get_gaussian_field_3d<real_t,sdim>(points_pos, num_points, out,
                                              xc, yc, zc,
                                              A,
                                              sigma_x, sigma_y, sigma_z);
}


template<typename real_t, int sdim>
void
get_gaussian_field_3d_wrapper_02(const real_t* points_pos,
                                 const int num_points,
                                 real_t* out) {

  real_t xc;
  real_t yc;
  real_t zc;
  real_t A       = 1.0 ;
  real_t sigma_x = 0.06;
  real_t sigma_y = 0.06;
  real_t sigma_z = 0.06;

  get_gaussian_position <real_t>(0.5, 0.7, 0.2, 0.4, 0.6, 0.8, xc, yc, zc);

  tbslas::get_gaussian_field_3d<real_t,sdim>(points_pos, num_points, out,
                                              xc, yc, zc,
                                              A,
                                              sigma_x, sigma_y, sigma_z);
}

template<typename real_t, int sdim>
void
get_gaussian_field_3d_wrapper_03(const real_t* points_pos,
                                 const int num_points,
                                 real_t* out) {

  real_t xc;
  real_t yc;
  real_t zc;
  real_t A       = 1.0 ;
  real_t sigma_x = 0.06;
  real_t sigma_y = 0.06;
  real_t sigma_z = 0.06;

  get_gaussian_position <real_t>(0.5, 0.7, 0.5, 0.85, 0.2, 0.4, xc, yc, zc);

  tbslas::get_gaussian_field_3d<real_t,sdim>(points_pos, num_points, out,
                                             xc, yc, zc,
                                             A,
                                             sigma_x, sigma_y, sigma_z);
}

void rho(const double* coord, int n, double* out){ //Input function
  int dof=1;
  size_t pt_cnt=pt_coord.size()/3;

  // tansfer to the middle of the domain
  // std::vector<double> coord_local(n*COORD_DIM);
  // for (int i=0;i<n*3;i++) {
  //   coord_local[i]=(coord[i]-0.5)*2+0.5;
  // }

  switch (TEST_CASE)
  {
    case 0: // Constant coefficient
      for(int i=0;i<n;i++){
        const double* c=&coord[i*COORD_DIM];
        {
          out[i*dof+0]= 0;
        }
      }
      break;
    case 1: // Synthetic variable coefficient
      for(int i=0;i<n;i++){
        double L=500;
        const double* c=&coord[i*COORD_DIM];
        {
          double r_2=(c[0]-0.5)*(c[0]-0.5)+(c[1]-0.5)*(c[1]-0.5)+(c[2]-0.5)*(c[2]-0.5);
          out[i*dof+0]=rho_*exp(-L*r_2);
        }
      }
      break;
    case 2: // Sphere
      for(int i=0;i<n;i++){
        out[i*dof+0]=1;
        // const double* c=&coord_local[i*COORD_DIM];
        const double* c=&coord[i*COORD_DIM];
        for(size_t j=0;j<pt_cnt;j++){
          double r_2=0;
          r_2+=(c[0]-pt_coord[j*COORD_DIM+0])*(c[0]-pt_coord[j*COORD_DIM+0]);
          r_2+=(c[1]-pt_coord[j*COORD_DIM+1])*(c[1]-pt_coord[j*COORD_DIM+1]);
          r_2+=(c[2]-pt_coord[j*COORD_DIM+2])*(c[2]-pt_coord[j*COORD_DIM+2]);
          double r=sqrt(r_2);
          out[i*dof+0]*=1.0/(1.0+exp(-(4.0/JUMP_WIDTH)*(r-PT_RAD)));
        }
        out[i*dof+0]=rho_*(1-out[i*dof+0]);
        //out[i*dof+0]=rho_*(1-out[i*dof+0])*out[i*dof+0]; // Surface
        for(int i=0;i<n;i++) {
          if(fabs(coord[i*3+0]-0.5)>0.5) out[i*dof+0]=0;
          if(fabs(coord[i*3+1]-0.5)>0.5) out[i*dof+0]=0;
          if(fabs(coord[i*3+2]-0.5)>0.5) out[i*dof+0]=0;
        }
      }
      break;
    case 3: // Porous media
      for(int i=0;i<n;i++){
        out[i*dof+0]=1;
        const double* c=&coord[i*COORD_DIM];
        for(size_t j=0;j<pt_cnt;j++){
          double r_2=0;
          r_2+=(c[0]-pt_coord[j*COORD_DIM+0])*(c[0]-pt_coord[j*COORD_DIM+0]);
          r_2+=(c[1]-pt_coord[j*COORD_DIM+1])*(c[1]-pt_coord[j*COORD_DIM+1]);
          r_2+=(c[2]-pt_coord[j*COORD_DIM+2])*(c[2]-pt_coord[j*COORD_DIM+2]);
          double r=sqrt(r_2);
          out[i*dof+0]*=1.0/(1.0+exp(-(4.0/JUMP_WIDTH)*(r-1.2*PT_RAD)));
        }
        //{
        //  double r=0.49-fabs(c[0]-0.5);
        //  out[i*dof+0]*=1.0/(1.0+exp(-(4.0/JUMP_WIDTH)*r));
        //}
        //{
        //  double r=0.49-fabs(c[1]-0.5);
        //  out[i*dof+0]*=1.0/(1.0+exp(-(4.0/JUMP_WIDTH)*r));
        //}
        //{
        //  double r=0.49-fabs(c[2]-0.5);
        //  out[i*dof+0]*=1.0/(1.0+exp(-(4.0/JUMP_WIDTH)*r));
        //}
        out[i*dof+0]=rho_*out[i*dof+0];
        for(int i=0;i<n;i++) {
          if(fabs(coord[i*3+0]-0.5)>0.5) out[i*dof+0]=0;
          if(fabs(coord[i*3+1]-0.5)>0.5) out[i*dof+0]=0;
          if(fabs(coord[i*3+2]-0.5)>0.5) out[i*dof+0]=0;
        }
      }
      break;
    default:
      break;
  }
  // // tansfer back the domain
  // for (int i=0;i<n*3;i++) {
  //   coord[i] = (coord[i]-0.5)*0.5+0.5;
  // }
}

void fn_input(const double* coord, int n, double* out){ //Input function
  int dof=INPUT_DOF;

  switch (TEST_CASE)
  {
    case 0: // Constant coefficient
      for(int i=0;i<n;i++){
        double L=125;
        const double* c=&coord[i*COORD_DIM];
        {
          double r_2=(c[0]-0.5)*(c[0]-0.5)+(c[1]-0.5)*(c[1]-0.5)+(c[2]-0.5)*(c[2]-0.5);
          if(INPUT_DOF>0) out[i*dof+0]=                                        0;
          if(INPUT_DOF>1) out[i*dof+1]= 4*L*L*(c[2]-0.5)*(5-2*L*r_2)*exp(-L*r_2);
          if(INPUT_DOF>2) out[i*dof+2]=-4*L*L*(c[1]-0.5)*(5-2*L*r_2)*exp(-L*r_2);
        }
      }
      break;
    case 1: // Synthetic variable coefficient
      for(int i=0;i<n;i++){
        double L=500;
        const double* c=&coord[i*COORD_DIM];
        {
          double r_2=(c[0]-0.5)*(c[0]-0.5)+(c[1]-0.5)*(c[1]-0.5)+(c[2]-0.5)*(c[2]-0.5);
          double rho_val=1.0;
          rho(c, 1, &rho_val);
          if(INPUT_DOF>0) out[i*dof+0]=                                        0+                         0*rho_val;
          if(INPUT_DOF>1) out[i*dof+1]= 4*L*L*(c[2]-0.5)*(5-2*L*r_2)*exp(-L*r_2)+2*L*(c[2]-0.5)*exp(-L*r_2)*rho_val;
          if(INPUT_DOF>2) out[i*dof+2]=-4*L*L*(c[1]-0.5)*(5-2*L*r_2)*exp(-L*r_2)-2*L*(c[1]-0.5)*exp(-L*r_2)*rho_val;
        }
      }
      break;
    case 2: // Sphere
      for(int i=0;i<n;i++){
        const double* c=&coord[i*COORD_DIM];
        {
          double r_2=(c[0]-0.1)*(c[0]-0.1)+(c[1]-0.5)*(c[1]-0.5)+(c[2]-0.5)*(c[2]-0.5);
          double rho_val=1.0;
          rho(c, 1, &rho_val);
          if(INPUT_DOF>0) out[i*dof+0]=1*rho_val;
          if(INPUT_DOF>1) out[i*dof+1]=0*rho_val;
          if(INPUT_DOF>2) out[i*dof+2]=0*rho_val;
        }
      }
      break;
    case 3: // Porous media
      for(int i=0;i<n;i++){
        const double* c=&coord[i*COORD_DIM];
        {
          double r_2=(c[0]-0.1)*(c[0]-0.1)+(c[1]-0.5)*(c[1]-0.5)+(c[2]-0.5)*(c[2]-0.5);
          double rho_val=1.0;
          rho(c, 1, &rho_val);
          if(INPUT_DOF>0) out[i*dof+0]=1*rho_val;
          if(INPUT_DOF>1) out[i*dof+1]=0*rho_val;
          if(INPUT_DOF>2) out[i*dof+2]=0*rho_val;
        }
      }
      break;
    default:
      break;
  }
}

void u_ref(const double* coord, int n, double* out){ //Analytical solution
  int dof=INPUT_DOF;

  switch (TEST_CASE)
  {
    case 0: // Constant coefficient
      for(int i=0;i<n;i++){
        double L=125;
        const double* c=&coord[i*COORD_DIM];
        {
          double r_2=(c[0]-0.5)*(c[0]-0.5)+(c[1]-0.5)*(c[1]-0.5)+(c[2]-0.5)*(c[2]-0.5);
          out[i*dof+0]=                          0;
          out[i*dof+1]=-2*L*(c[2]-0.5)*exp(-L*r_2);
          out[i*dof+2]=+2*L*(c[1]-0.5)*exp(-L*r_2);
        }
      }
      break;
    case 1: // Synthetic variable coefficient
      for(int i=0;i<n;i++){
        double L=500;
        const double* c=&coord[i*COORD_DIM];
        {
          double r_2=(c[0]-0.5)*(c[0]-0.5)+(c[1]-0.5)*(c[1]-0.5)+(c[2]-0.5)*(c[2]-0.5);
          out[i*dof+0]=                          0;
          out[i*dof+1]=-2*L*(c[2]-0.5)*exp(-L*r_2);
          out[i*dof+2]=+2*L*(c[1]-0.5)*exp(-L*r_2);
        }
      }
      break;
    case 2: // Sphere
      for(int i=0;i<n;i++){
        const double* c=&coord[i*COORD_DIM];
        double r_2=0;
        r_2+=(c[0]-0.5)*(c[0]-0.5);
        r_2+=(c[1]-0.5)*(c[1]-0.5);
        r_2+=(c[2]-0.5)*(c[2]-0.5);
        double r=sqrt(r_2);
        if(r>PT_RAD){
          double cos_t=(c[0]-0.5)/r;
          double sin_t=sqrt(1-cos_t*cos_t);
          double r0_r=PT_RAD/r;
          double ur= cos_t*(1.0-1.50*r0_r+0.50*r0_r*r0_r*r0_r);
          double ut=-sin_t*(1.0-0.75*r0_r-0.25*r0_r*r0_r*r0_r);
          out[i*dof+0]=cos_t*ur-sin_t*ut-1.0;

          double r_yz=sqrt((c[1]-0.5)*(c[1]-0.5)+(c[2]-0.5)*(c[2]-0.5));
          out[i*dof+1]=(sin_t*ur+cos_t*ut)*(c[1]-0.5)/r_yz;
          out[i*dof+2]=(sin_t*ur+cos_t*ut)*(c[2]-0.5)/r_yz;
        }else{
          out[i*dof+0]=-1.0;
          out[i*dof+1]= 0.0;
          out[i*dof+2]= 0.0;
        }
      }
      break;
    case 3: // Porous media
      for(int i=0;i<n;i++){
        const double* c=&coord[i*COORD_DIM];
        {
          double r_2=(c[0]-0.5)*(c[0]-0.5)+(c[1]-0.5)*(c[1]-0.5)+(c[2]-0.5)*(c[2]-0.5);
          out[i*dof+0]=0;
          out[i*dof+1]=0;
          out[i*dof+2]=0;
        }
      }
      break;
    default:
      break;
  }
}

////////////////////////////////////////////////////////////////////////////////
struct FMMData{
  const pvfmm::Kernel<double>* kernel;
  FMM_Mat_t* fmm_mat;
  FMM_Tree_t* tree;
  PetscInt m,n,M,N;
  std::vector<double> rho;
  std::vector<double> rho1; // for checking error
  std::vector<double> u_ref;
  pvfmm::BoundaryType bndry;
};

int tree2vec(FMMData fmm_data, Vec& Y){
  PetscErrorCode ierr;
  FMM_Tree_t* tree=fmm_data.tree;
  int cheb_deg=fmm_data.fmm_mat->ChebDeg();

  std::vector<FMMNode_t*> nlist;
  { // Get non-ghost, leaf nodes.
    std::vector<FMMNode_t*>& nlist_=tree->GetNodeList();
    for(size_t i=0;i<nlist_.size();i++){
      if(nlist_[i]->IsLeaf() && !nlist_[i]->IsGhost()){
        nlist.push_back(nlist_[i]);
      }
    }
  }
  assert(nlist.size()>0);

  int omp_p=omp_get_max_threads();
  size_t n_coeff3=(cheb_deg+1)*(cheb_deg+2)*(cheb_deg+3)/6;

  {
    PetscInt Y_size;
    ierr = VecGetLocalSize(Y, &Y_size);
    int data_dof=Y_size/(n_coeff3*nlist.size());

    PetscScalar *Y_ptr;
    ierr = VecGetArray(Y, &Y_ptr);

#pragma omp parallel for
    for(size_t tid=0;tid<omp_p;tid++){
      size_t i_start=(nlist.size()* tid   )/omp_p;
      size_t i_end  =(nlist.size()*(tid+1))/omp_p;
      for(size_t i=i_start;i<i_end;i++){
        pvfmm::Vector<double>& coeff_vec=nlist[i]->ChebData();
        double s=std::pow(0.5,COORD_DIM*nlist[i]->Depth()*0.5*SCAL_EXP);

        size_t Y_offset=i*n_coeff3*data_dof;
        for(size_t j=0;j<n_coeff3*data_dof;j++) Y_ptr[j+Y_offset]=coeff_vec[j]*s;
      }
    }
    ierr = VecRestoreArray(Y, &Y_ptr);
  }

  return 0;
}

int vec2tree(Vec& Y, FMMData fmm_data){
  PetscErrorCode ierr;
  FMM_Tree_t* tree=fmm_data.tree;
  const MPI_Comm* comm=tree->Comm();
  int cheb_deg=fmm_data.fmm_mat->ChebDeg();

  std::vector<FMMNode_t*> nlist;
  { // Get non-ghost, leaf nodes.
    std::vector<FMMNode_t*>& nlist_=tree->GetNodeList();
    for(size_t i=0;i<nlist_.size();i++){
      if(nlist_[i]->IsLeaf() && !nlist_[i]->IsGhost()){
        nlist.push_back(nlist_[i]);
      }
    }
  }

  int omp_p=omp_get_max_threads();
  size_t n_coeff3=(cheb_deg+1)*(cheb_deg+2)*(cheb_deg+3)/6;

  {
    PetscInt Y_size;
    ierr = VecGetLocalSize(Y, &Y_size);
    int data_dof=Y_size/(n_coeff3*nlist.size());

    const PetscScalar *Y_ptr;
    ierr = VecGetArrayRead(Y, &Y_ptr);

#pragma omp parallel for
    for(size_t tid=0;tid<omp_p;tid++){
      size_t i_start=(nlist.size()* tid   )/omp_p;
      size_t i_end  =(nlist.size()*(tid+1))/omp_p;
      for(size_t i=i_start;i<i_end;i++){
        pvfmm::Vector<double>& coeff_vec=nlist[i]->ChebData();
        double s=std::pow(2.0,COORD_DIM*nlist[i]->Depth()*0.5*SCAL_EXP);

        size_t Y_offset=i*n_coeff3*data_dof;
        for(size_t j=0;j<n_coeff3*data_dof;j++) coeff_vec[j]=PetscRealPart(Y_ptr[j+Y_offset])*s;
        nlist[i]->DataDOF()=data_dof;
      }
    }
  }

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
int FMM_Init(MPI_Comm& comm, FMMData *fmm_data){
  int myrank, np;
  MPI_Comm_rank(comm, &myrank);
  MPI_Comm_size(comm,&np);

  FMM_Mat_t *fmm_mat=new FMM_Mat_t;
  FMM_Tree_t* tree=new FMM_Tree_t(comm);

  //Kernel function
  const pvfmm::Kernel<double>* kernel;
  // if(INPUT_DOF==1) kernel=&pvfmm::LaplaceKernel<double>::potn_ker();
  //else if(INPUT_DOF==2) kernel=&pvfmm::ker_helmholtz;
  if(INPUT_DOF==1) kernel=&pvfmm::LaplaceKernel<double>::potential();
  else if(INPUT_DOF==3) kernel=&pvfmm::StokesKernel<double>::velocity();

  //Setup FMM data structure.
  int mult_order=MUL_ORDER;
  int cheb_deg=CHEB_DEG;
  bool adap=true;
  double tol=TOL;
  pvfmm::BoundaryType bndry=pvfmm::FreeSpace;
  if(PERIODIC==PETSC_TRUE) bndry=pvfmm::Periodic;

  typename FMMNode_t::NodeData tree_data;
  { // Tree Construction (rho).
    FMM_Tree_t* tree=new FMM_Tree_t(comm);
    //Various parameters.
    tree_data.dim=COORD_DIM;
    tree_data.max_depth=MAXDEPTH;
    tree_data.cheb_deg=cheb_deg;

    //Set input function pointer
    tree_data.input_fn=rho;
    tree_data.data_dof=1;
    tree_data.tol=tol*fabs(rho_);

    std::vector<double> pt_coord;
    { //Set source coordinates.
      size_t NN=ceil(pow((double)np,1.0/3.0));
      NN=std::max<size_t>(NN,pow(2.0,MINDEPTH));
      size_t N_total=NN*NN*NN;
      size_t start= myrank   *N_total/np;
      size_t end  =(myrank+1)*N_total/np;
      for(size_t i=start;i<end;i++){
        pt_coord.push_back(((double)((i/  1    )%NN)+0.5)/NN);
        pt_coord.push_back(((double)((i/ NN    )%NN)+0.5)/NN);
        pt_coord.push_back(((double)((i/(NN*NN))%NN)+0.5)/NN);
      }
      tree_data.pt_coord=pt_coord;
    }
    tree_data.max_pts=1; // Points per octant.

    //Create Tree and initialize with input data.
    tree->Initialize(&tree_data);
    tree->InitFMM_Tree(adap,bndry);

    pt_coord.clear();
    std::vector<FMMNode_t*> nlist=tree->GetNodeList();
    for(size_t i=0;i<nlist.size();i++){
      if(nlist[i]->IsLeaf() && !nlist[i]->IsGhost()){
        double s=pow(0.5,nlist[i]->Depth()+1);
        double* c=nlist[i]->Coord();
        pt_coord.push_back(c[0]+s);
        pt_coord.push_back(c[1]+s);
        pt_coord.push_back(c[2]+s);
      }
    }
    tree_data.pt_coord=pt_coord;

    char out_name_buffer[300];

    snprintf(out_name_buffer, sizeof(out_name_buffer),
             "%s/stokes_rho_%d_", tbslas::get_result_dir().c_str(), 0);

    tree->Write2File(out_name_buffer,VTK_ORDER);
    delete tree;
  }
  { // Tree Construction.
    bool adap=true;

    //Set input function pointer
    tree_data.input_fn=fn_input;
    tree_data.data_dof=kernel->ker_dim[0];
    tree_data.tol=tol*f_max;

    //Create Tree and initialize with input data.
    tree->Initialize(&tree_data);
    tree->InitFMM_Tree(adap,bndry);

    std::vector<double> pt_coord;
    pt_coord.clear();
    std::vector<FMMNode_t*> nlist=tree->GetNodeList();
    for(size_t i=0;i<nlist.size();i++){
      if(nlist[i]->IsLeaf() && !nlist[i]->IsGhost()){
        double s=pow(0.5,nlist[i]->Depth()+1);
        double* c=nlist[i]->Coord();
        pt_coord.push_back(c[0]+s);
        pt_coord.push_back(c[1]+s);
        pt_coord.push_back(c[2]+s);
      }
    }
    tree_data.pt_coord=pt_coord;

    { //Output max tree depth.
      std::vector<size_t> all_nodes(MAXDEPTH+1,0);
      std::vector<size_t> leaf_nodes(MAXDEPTH+1,0);
      std::vector<FMMNode_t*>& nodes=tree->GetNodeList();
      for(size_t i=0;i<nodes.size();i++){
        FMMNode_t* n=nodes[i];
        if(!n->IsGhost()) all_nodes[n->Depth()]++;
        if(!n->IsGhost() && n->IsLeaf()) leaf_nodes[n->Depth()]++;
      }

      if(!myrank) std::cout<<"All  Nodes: ";
      for(int i=0;i<MAXDEPTH;i++){
        int local_size=all_nodes[i];
        int global_size;
        MPI_Allreduce(&local_size, &global_size, 1, MPI_INT, MPI_SUM, comm);
        if(global_size==0) MAXDEPTH=i;
        if(!myrank) std::cout<<global_size<<' ';
      }
      if(!myrank) std::cout<<'\n';

      if(!myrank) std::cout<<"Leaf Nodes: ";
      for(int i=0;i<MAXDEPTH;i++){
        int local_size=leaf_nodes[i];
        int global_size;
        MPI_Allreduce(&local_size, &global_size, 1, MPI_INT, MPI_SUM, comm);
        if(!myrank) std::cout<<global_size<<' ';
      }
      if(!myrank) std::cout<<'\n';
    }

  }

  size_t n_coeff3=(cheb_deg+1)*(cheb_deg+2)*(cheb_deg+3)/6;
  PetscInt m=0,n=0,M=0,N=0;
  { // Get local and global size
    long long loc_size=0, glb_size=0;
    std::vector<FMMNode_t*> nlist=tree->GetNodeList();
    for(size_t i=0;i<nlist.size();i++)
      if(nlist[i]->IsLeaf() && !nlist[i]->IsGhost())
        loc_size+=n_coeff3; //nlist[i]->ChebData().Dim();
    MPI_Allreduce(&loc_size, &glb_size, 1, MPI_LONG_LONG , MPI_SUM, comm);
    n=loc_size*kernel->ker_dim[0];
    N=glb_size*kernel->ker_dim[0];
    m=loc_size*kernel->ker_dim[1];
    M=glb_size*kernel->ker_dim[1];
    num_oct=glb_size/n_coeff3;
  }
  if(TREE_ONLY) return 0;

  //Initialize FMM_Mat.
  fmm_mat->Initialize(mult_order,cheb_deg,comm,kernel);

  fmm_data->kernel =kernel ;
  fmm_data->fmm_mat=fmm_mat;
  fmm_data->tree   =tree   ;
  fmm_data->bndry  =bndry  ;
  fmm_data->m=m;
  fmm_data->n=n;
  fmm_data->M=M;
  fmm_data->N=N;

  return 0;
}

int FMMCreateShell(FMMData *fmm_data, Mat *A){
  FMM_Tree_t*   tree=fmm_data->tree   ;
  const MPI_Comm& comm=*tree->Comm();
  PetscInt m,n,M,N;
  m=fmm_data->m;
  n=fmm_data->n;
  M=fmm_data->M;
  N=fmm_data->N;

  { // Evaluate rho at Chebyshev node points.
    std::vector<double>& rho_vec=fmm_data->rho;

    std::vector<FMMNode_t*> nlist;
    { // Get non-ghost, leaf nodes.
      std::vector<FMMNode_t*>& nlist_=tree->GetNodeList();
      for(size_t i=0;i<nlist_.size();i++){
        if(nlist_[i]->IsLeaf() && !nlist_[i]->IsGhost()){
          nlist.push_back(nlist_[i]);
        }
      }
    }

    int cheb_deg=fmm_data->fmm_mat->ChebDeg();
    size_t n_nodes3=(cheb_deg+1)*(cheb_deg+1)*(cheb_deg+1);
    rho_vec.resize(n_nodes3*nlist.size());

    std::vector<double> cheb_node_coord3=pvfmm::cheb_nodes<double>(cheb_deg, 3);
    int omp_p=omp_get_max_threads();
#pragma omp parallel for
    for(size_t tid=0;tid<omp_p;tid++){
      size_t i_start=(nlist.size()* tid   )/omp_p;
      size_t i_end  =(nlist.size()*(tid+1))/omp_p;

      std::vector<double> cheb_node_coord3_(n_nodes3*3);
      std::vector<double> rho_val(n_nodes3);
      for(size_t i=i_start;i<i_end;i++){
        // Shift Cheb node points and evaluate rho
        double* coord=nlist[i]->Coord();
        double s=pow(0.5,nlist[i]->Depth());
        for(size_t j=0;j<n_nodes3;j++){
          cheb_node_coord3_[j*3+0]=cheb_node_coord3[j*3+0]*s+coord[0];
          cheb_node_coord3_[j*3+1]=cheb_node_coord3[j*3+1]*s+coord[1];
          cheb_node_coord3_[j*3+2]=cheb_node_coord3[j*3+2]*s+coord[2];
        }
        rho(&cheb_node_coord3_[0], n_nodes3, &rho_val[0]);

        size_t vec_offset=i*n_nodes3;
        for(size_t j=0;j<n_nodes3;j++){
          rho_vec[vec_offset+j]=rho_val[j];
        }
      }
    }
  }
  { // Evaluate rho at Chebyshev node points.
    std::vector<double>& rho_vec=fmm_data->rho1;
    std::vector<double>& u_ref_vec=fmm_data->u_ref;

    std::vector<FMMNode_t*> nlist;
    { // Get non-ghost, leaf nodes.
      std::vector<FMMNode_t*>& nlist_=tree->GetNodeList();
      for(size_t i=0;i<nlist_.size();i++){
        if(nlist_[i]->IsLeaf() && !nlist_[i]->IsGhost()){
          nlist.push_back(nlist_[i]);
        }
      }
    }

    int cheb_deg=fmm_data->fmm_mat->ChebDeg()*ERROR_RESOLUTION;
    size_t n_nodes3=(cheb_deg+1)*(cheb_deg+1)*(cheb_deg+1);
    rho_vec.resize(n_nodes3*nlist.size());
    u_ref_vec.resize(n_nodes3*nlist.size()*3);

    std::vector<double> cheb_node_coord3=pvfmm::cheb_nodes<double>(cheb_deg, 3);
    int omp_p=omp_get_max_threads();
#pragma omp parallel for
    for(size_t tid=0;tid<omp_p;tid++){
      size_t i_start=(nlist.size()* tid   )/omp_p;
      size_t i_end  =(nlist.size()*(tid+1))/omp_p;

      std::vector<double> cheb_node_coord3_(n_nodes3*3);
      std::vector<double> rho_val(n_nodes3);
      std::vector<double> u_ref_val(n_nodes3*3);
      for(size_t i=i_start;i<i_end;i++){
        // Shift Cheb node points and evaluate rho
        double* coord=nlist[i]->Coord();
        double s=pow(0.5,nlist[i]->Depth());
        for(size_t j=0;j<n_nodes3;j++){
          cheb_node_coord3_[j*3+0]=cheb_node_coord3[j*3+0]*s+coord[0];
          cheb_node_coord3_[j*3+1]=cheb_node_coord3[j*3+1]*s+coord[1];
          cheb_node_coord3_[j*3+2]=cheb_node_coord3[j*3+2]*s+coord[2];
        }
        rho(&cheb_node_coord3_[0], n_nodes3, &rho_val[0]);
        u_ref(&cheb_node_coord3_[0], n_nodes3, &u_ref_val[0]);

        size_t vec_offset=i*n_nodes3;
        for(size_t j=0;j<n_nodes3;j++){
          rho_vec[vec_offset+j]=rho_val[j];
          u_ref_vec[(vec_offset+j)*3+0]=u_ref_val[j*3+0];
          u_ref_vec[(vec_offset+j)*3+1]=u_ref_val[j*3+1];
          u_ref_vec[(vec_offset+j)*3+2]=u_ref_val[j*3+2];
        }
      }
    }
  }

  return MatCreateShell(comm,m,n,M,N,fmm_data,A);
}

int FMMDestroy(FMMData *fmm_data){
  delete fmm_data->fmm_mat;
  delete fmm_data->tree;
  return 1;
}

////////////////////////////////////////////////////////////////////////////////

int stokes_err(Mat M, Vec U){
  PetscErrorCode ierr;
  FMMData* fmm_data=NULL;
  MatShellGetContext(M, &fmm_data);
  FMM_Tree_t* tree=fmm_data->tree;
  const MPI_Comm* comm=tree->Comm();

  int myrank, np;
  MPI_Comm_rank(*comm, &myrank);
  MPI_Comm_size(*comm,&np);

  int cheb_deg=fmm_data->fmm_mat->ChebDeg()*ERROR_RESOLUTION;
  std::vector<double>& rho_vec=fmm_data->rho1;
  std::vector<double>& u_ref_vec=fmm_data->u_ref;
  int omp_p=omp_get_max_threads();
  pvfmm::Profile::Tic("StokesErr",comm,true);

  std::vector<FMMNode_t*> nlist;
  { // Get non-ghost, leaf nodes.
    std::vector<FMMNode_t*>& nlist_=tree->GetNodeList();
    for(size_t i=0;i<nlist_.size();i++){
      if(nlist_[i]->IsLeaf() && !nlist_[i]->IsGhost()){
        nlist.push_back(nlist_[i]);
      }
    }
  }
  assert(nlist.size()>0);

  // Cheb node points
  size_t n_coeff3=(cheb_deg/ERROR_RESOLUTION+1)*(cheb_deg/ERROR_RESOLUTION+2)*(cheb_deg/ERROR_RESOLUTION+3)/6;
  size_t n_nodes3=(cheb_deg+1)*(cheb_deg+1)*(cheb_deg+1);
  std::vector<double> cheb_node_coord1=pvfmm::cheb_nodes<double>(cheb_deg, 1);
#pragma omp parallel for
  for(size_t i=0;i<cheb_node_coord1.size();i++){
    cheb_node_coord1[i]=cheb_node_coord1[i]*2.0-1.0;
  }

  // Input Vector ( \rho * U )
  double vel_ss_inf   =0;
  double vel_ss_l2    =0;
  double vel_inf      =0;
  double vel_l2       =0;
  double vel_shift_inf=0;
  double vel_shift_l2 =0;
  double vel_rho_inf  =0;
  double vel_rho_l2   =0;
  double rho_inf      =0;
  double rho_l2       =0;
  {
    PetscInt U_size;
    ierr = VecGetLocalSize(U, &U_size);
    int data_dof=U_size/(n_coeff3*nlist.size());
    assert(data_dof*n_coeff3*nlist.size()==U_size);

    PetscScalar *U_ptr;
    ierr = VecGetArray(U, &U_ptr);
    //#pragma omp parallel for
    for(size_t tid=0;tid<omp_p;tid++){
      size_t i_start=(nlist.size()* tid   )/omp_p;
      size_t i_end  =(nlist.size()*(tid+1))/omp_p;
      pvfmm::Vector<double> coeff_vec(n_coeff3*data_dof);
      pvfmm::Vector<double> vel_vec(n_nodes3*data_dof);
      for(size_t i=i_start;i<i_end;i++){
        double s=std::pow(2.0,COORD_DIM*nlist[i]->Depth()*0.5*SCAL_EXP);

        { // coeff_vec: Cheb coeff data for this node
          size_t U_offset=i*n_coeff3*data_dof;
          for(size_t j=0;j<n_coeff3*data_dof;j++) coeff_vec[j]=PetscRealPart(U_ptr[j+U_offset])*s;
        }

        // vel_vec: Evaluate coeff_vec at Chebyshev node points
        cheb_eval(coeff_vec, cheb_deg/ERROR_RESOLUTION, cheb_node_coord1, cheb_node_coord1, cheb_node_coord1, vel_vec);

        {// Compute error
          double* rho=&rho_vec[i*n_nodes3];
          double s=std::pow(0.5,COORD_DIM*nlist[i]->Depth()); // octant volume
          double* ss=&u_ref_vec[i*n_nodes3*3];
          for(size_t j0=0;j0<data_dof;j0++){
            double* vel=&vel_vec[j0*n_nodes3];
            double shift=0.0;
            if(j0==0) shift=-1.0;
            for(size_t j1=0;j1<n_nodes3;j1++){
              double n_n0=rho[j1]/rho_;
              double omn_n0=1-n_n0;
              if(TEST_CASE==0 || TEST_CASE==1){
                n_n0=0.0;
                omn_n0=1.0;
                shift=0.0;
              }

              vel_ss_inf    =std::max(vel_ss_inf,fabs((vel[j1]+ss[j1*3+j0])*omn_n0)); // Analytical error
              vel_ss_l2    += (vel[j1]+ss[j1*3+j0])*omn_n0 * (vel[j1]+ss[j1*3+j0])*omn_n0 *s/n_nodes3/3;

              vel_inf       =std::max(vel_inf,fabs(vel[j1]*omn_n0)); // Velocity outsize
              vel_l2       += (vel[j1]*omn_n0) * (vel[j1]*omn_n0) *s/n_nodes3/3;

              vel_shift_inf =std::max(vel_shift_inf,fabs((vel[j1]+shift)*omn_n0)); // Velocity outsize
              vel_shift_l2 += (vel[j1]+shift)*omn_n0 * (vel[j1]+shift)*omn_n0 *s/n_nodes3/3;

              vel_rho_inf   =std::max(vel_rho_inf,fabs((vel[j1]+shift)*n_n0)); // Velocity inside
              vel_rho_l2   += ((vel[j1]+shift)*n_n0) * ((vel[j1]+shift)*n_n0) *s/n_nodes3/3;

              rho_inf       =std::max(rho_inf,fabs(n_n0)); // rho norm
              rho_l2       += (n_n0) * (n_n0) *s/n_nodes3/3;
            }
          }
        }
      }
    }
  }

  double glb_vel_ss_inf   =0;
  double glb_vel_ss_l2    =0;
  double glb_vel_inf      =0;
  double glb_vel_l2       =0;
  double glb_vel_shift_inf=0;
  double glb_vel_shift_l2 =0;
  double glb_vel_rho_inf  =0;
  double glb_vel_rho_l2   =0;
  double glb_rho_inf      =0;
  double glb_rho_l2       =0;
  MPI_Allreduce(&vel_ss_inf   , &glb_vel_ss_inf   , 1, MPI_DOUBLE , MPI_MAX, *comm);
  MPI_Allreduce(&vel_ss_l2    , &glb_vel_ss_l2    , 1, MPI_DOUBLE , MPI_SUM, *comm);
  MPI_Allreduce(&vel_inf      , &glb_vel_inf      , 1, MPI_DOUBLE , MPI_MAX, *comm);
  MPI_Allreduce(&vel_l2       , &glb_vel_l2       , 1, MPI_DOUBLE , MPI_SUM, *comm);
  MPI_Allreduce(&vel_shift_inf, &glb_vel_shift_inf, 1, MPI_DOUBLE , MPI_MAX, *comm);
  MPI_Allreduce(&vel_shift_l2 , &glb_vel_shift_l2 , 1, MPI_DOUBLE , MPI_SUM, *comm);
  MPI_Allreduce(&vel_rho_inf  , &glb_vel_rho_inf  , 1, MPI_DOUBLE , MPI_MAX, *comm);
  MPI_Allreduce(&vel_rho_l2   , &glb_vel_rho_l2   , 1, MPI_DOUBLE , MPI_SUM, *comm);
  MPI_Allreduce(&rho_inf      , &glb_rho_inf      , 1, MPI_DOUBLE , MPI_MAX, *comm);
  MPI_Allreduce(&rho_l2       , &glb_rho_l2       , 1, MPI_DOUBLE , MPI_SUM, *comm);

  glb_vel_ss_inf   =     glb_vel_ss_inf   ;
  glb_vel_ss_l2    =sqrt(glb_vel_ss_l2)   ;
  glb_vel_inf      =     glb_vel_inf      ;
  glb_vel_l2       =sqrt(glb_vel_l2)      ;
  glb_vel_shift_inf=     glb_vel_shift_inf;
  glb_vel_shift_l2 =sqrt(glb_vel_shift_l2);
  glb_vel_rho_inf  =     glb_vel_rho_inf  ;
  glb_vel_rho_l2   =sqrt(glb_vel_rho_l2)  ;
  glb_rho_inf      =     glb_rho_inf      ;
  glb_rho_l2       =sqrt(glb_rho_l2)      ;

  BG_INF_ERR =glb_vel_rho_inf/glb_rho_inf; // Velocity inside normalized by norm of rho
  BG_L2_ERR  =glb_vel_rho_l2 /glb_rho_l2 ;
  AVG_INF_ERR=glb_vel_rho_inf/glb_vel_shift_inf; // Velocity inside normalized by norm of velocity outside
  AVG_L2_ERR =glb_vel_rho_l2 /glb_vel_shift_l2 ;
  SS_INF_ERR =glb_vel_ss_inf /glb_vel_inf;
  SS_L2_ERR  =glb_vel_ss_l2  /glb_vel_l2 ; // Analytical error norm
  if(!myrank){
    std::cout<<"ERROR:  |u-s|_inf        :"<<glb_vel_ss_inf   <<"\n";
    std::cout<<"ERROR:  |u-s|_l2         :"<<glb_vel_ss_l2    <<"\n";
    std::cout<<"ERROR:  |u-u0|_inf       :"<<glb_vel_inf      <<"\n";
    std::cout<<"ERROR:  |u-u0|_l2        :"<<glb_vel_l2       <<"\n";
    std::cout<<"ERROR:  |u|_inf          :"<<glb_vel_shift_inf<<"\n";
    std::cout<<"ERROR:  |u|_l2           :"<<glb_vel_shift_l2 <<"\n";
    std::cout<<"ERROR:  |rho/rho0|_inf   :"<<glb_rho_inf      <<"\n";
    std::cout<<"ERROR:  |rho/rho0|_l2    :"<<glb_rho_l2       <<"\n";
    std::cout<<"ERROR:  |u*rho/rho0|_inf :"<<glb_vel_rho_inf  <<"\n";
    std::cout<<"ERROR:  |u*rho/rho0|_l2  :"<<glb_vel_rho_l2   <<"\n";
    std::cout<<'\n';
    std::cout<<"ERROR:  |u*rho|_inf / |rho|_inf :"<<BG_INF_ERR <<"\n";
    std::cout<<"ERROR:  |u*rho|_l2  / |rho|_l2  :"<<BG_L2_ERR  <<"\n";
    std::cout<<"ERROR:  |u*rho|_inf / |u|_inf   :"<<AVG_INF_ERR<<"\n";
    std::cout<<"ERROR:  |u*rho|_l2  / |u|_l2    :"<<AVG_L2_ERR <<"\n";
    std::cout<<"ERROR:  |u-s|_inf / |u-u0|_inf  :"<<SS_INF_ERR <<"\n";
    std::cout<<"ERROR:  |u-s|_l2  / |u-u0|_l2   :"<<SS_L2_ERR  <<"\n";
  }

  pvfmm::Profile::Toc();
  return 0;
}

int mult(Mat M, Vec U, Vec Y){
  PetscErrorCode ierr;
  FMMData* fmm_data=NULL;
  MatShellGetContext(M, &fmm_data);
  FMM_Tree_t* tree=fmm_data->tree;
  const MPI_Comm* comm=tree->Comm();
  int cheb_deg=fmm_data->fmm_mat->ChebDeg();
  std::vector<double>& rho_vec=fmm_data->rho;
  int omp_p=omp_get_max_threads();
  pvfmm::Profile::Tic("FMM_Mul",comm,true);

  std::vector<FMMNode_t*> nlist;
  { // Get non-ghost, leaf nodes.
    std::vector<FMMNode_t*>& nlist_=tree->GetNodeList();
    for(size_t i=0;i<nlist_.size();i++){
      if(nlist_[i]->IsLeaf() && !nlist_[i]->IsGhost()){
        nlist.push_back(nlist_[i]);
      }
    }
  }
  assert(nlist.size()>0);

  // Cheb node points
  size_t n_coeff3=(cheb_deg+1)*(cheb_deg+2)*(cheb_deg+3)/6;
  size_t n_nodes3=(cheb_deg+1)*(cheb_deg+1)*(cheb_deg+1);
  std::vector<double> cheb_node_coord1=pvfmm::cheb_nodes<double>(cheb_deg, 1);
#pragma omp parallel for
  for(size_t i=0;i<cheb_node_coord1.size();i++){
    cheb_node_coord1[i]=cheb_node_coord1[i]*2.0-1.0;
  }

  // Input Vector ( \rho * U )
  pvfmm::Profile::Tic("FMM_Input",comm,true);
  {
    PetscInt U_size;
    ierr = VecGetLocalSize(U, &U_size);
    int data_dof=U_size/(n_coeff3*nlist.size());
    assert(data_dof*n_coeff3*nlist.size()==U_size);

    PetscScalar *U_ptr;
    ierr = VecGetArray(U, &U_ptr);
#pragma omp parallel for
    for(size_t tid=0;tid<omp_p;tid++){
      size_t i_start=(nlist.size()* tid   )/omp_p;
      size_t i_end  =(nlist.size()*(tid+1))/omp_p;
      pvfmm::Vector<double> coeff_vec(n_coeff3*data_dof);
      pvfmm::Vector<double> val_vec(n_nodes3*data_dof);
      for(size_t i=i_start;i<i_end;i++){
        double s=std::pow(2.0,COORD_DIM*nlist[i]->Depth()*0.5*SCAL_EXP);

        { // coeff_vec: Cheb coeff data for this node
          size_t U_offset=i*n_coeff3*data_dof;
          for(size_t j=0;j<n_coeff3*data_dof;j++) coeff_vec[j]=PetscRealPart(U_ptr[j+U_offset])*s;
        }

        // val_vec: Evaluate coeff_vec at Chebyshev node points
        cheb_eval(coeff_vec, cheb_deg, cheb_node_coord1, cheb_node_coord1, cheb_node_coord1, val_vec);

        {// rho*val_vec
          double* rho=&rho_vec[i*n_nodes3];
          for(size_t j0=0;j0<data_dof;j0++){
            double* vec=&val_vec[j0*n_nodes3];
            for(size_t j1=0;j1<n_nodes3;j1++) vec[j1]*=rho[j1];
          }
        }

        { // Compute Chebyshev approx
          pvfmm::Vector<double>& coeff_vec=nlist[i]->ChebData();
          if(coeff_vec.Dim()!=(data_dof*(cheb_deg+1)*(cheb_deg+2)*(cheb_deg+3))/6){
            coeff_vec.ReInit((data_dof*(cheb_deg+1)*(cheb_deg+2)*(cheb_deg+3))/6);
          }
          pvfmm::cheb_approx<double,double>(&val_vec[0], cheb_deg, data_dof, &coeff_vec[0]);
          nlist[i]->DataDOF()=data_dof;
        }
      }
    }
  }
  pvfmm::Profile::Toc();

  // Run FMM ( Compute: G[ \rho * u ] )
  tree->ClearFMMData();
  tree->RunFMM();

  // Copy data from tree to Y
  pvfmm::Profile::Tic("tree2vec",comm,true);
  {
    PetscInt Y_size;
    ierr = VecGetLocalSize(Y, &Y_size);
    int data_dof=Y_size/(n_coeff3*nlist.size());

    PetscScalar *Y_ptr;
    ierr = VecGetArray(Y, &Y_ptr);

#pragma omp parallel for
    for(size_t tid=0;tid<omp_p;tid++){
      size_t i_start=(nlist.size()* tid   )/omp_p;
      size_t i_end  =(nlist.size()*(tid+1))/omp_p;
      for(size_t i=i_start;i<i_end;i++){
        pvfmm::Vector<double>& coeff_vec=((FMM_Mat_t::FMMData*) nlist[i]->FMMData())->cheb_out;
        double s=std::pow(0.5,COORD_DIM*nlist[i]->Depth()*0.5*SCAL_EXP);

        size_t Y_offset=i*n_coeff3*data_dof;
        for(size_t j=0;j<n_coeff3*data_dof;j++) Y_ptr[j+Y_offset]=coeff_vec[j]*s;
      }
    }

    ierr = VecRestoreArray(Y, &Y_ptr);
  }
  pvfmm::Profile::Toc();

  // Output Vector ( Compute:  U + G[ \rho * U ] )
  pvfmm::Profile::Tic("FMM_Output",comm,true);
  ierr = VecAXPY(Y,1,U);CHKERRQ(ierr);
  pvfmm::Profile::Toc();

  pvfmm::Profile::Toc();
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
int disp_result(){
  MPI_Comm comm=MPI_COMM_WORLD;
  PetscMPIInt    rank,size;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);

  if(rank==0){
    int width=10;
    std::cout<<std::setprecision(3)<<std::setiosflags(std::ios::right);
    std::cout<<"       ";
    std::cout<<"  "<<std::setw(width)<<"Device";
    std::cout<<"  "<<std::setw(width)<<"proc";
    std::cout<<"  "<<std::setw(5)<<"thrd";
    std::cout<<"  "<<std::setw(width)<<"ref_tol";
    std::cout<<"  "<<std::setw(width)<<"rho";
    std::cout<<"  "<<std::setw(5)<<"depth";
    std::cout<<"  "<<std::setw(width)<<"oct";
    std::cout<<"  |";
    std::cout<<"  "<<std::setw(4)<<"q";
    std::cout<<"  "<<std::setw(4)<<"m";
    std::cout<<"  "<<std::setw(width)<<"gmres_tol";
    std::cout<<"  |";
    std::cout<<"  "<<std::setw(   23)<<"|u-s|_inf/|u-u0|_inf";
    std::cout<<"  "<<std::setw(   23)<<"|u-s|_l2 /|u-u0|_l2 ";
    std::cout<<"  "<<std::setw(width)<<"time_ksp";
    std::cout<< "("<<std::setw(    4)<<"iter"<<")";
    std::cout<<"  |";
    std::cout<<"\n";

    std::cout<<"Result:";
#if defined(__INTEL_OFFLOAD)
    std::cout<<"  "<<                 std::setw(width)<<"Phi";
#elif defined(PVFMM_HAVE_CUDA)
    std::cout<<"  "<<                 std::setw(width)<<"GPU";
#else
    std::cout<<"  "<<                 std::setw(width)<<"CPU";
#endif
    std::cout<<"  "<<                 std::setw(width)<<size;
    std::cout<<"  "<<                 std::setw(5)<<omp_get_max_threads();
    std::cout<<"  "<<std::scientific<<std::setw(width)<<TOL;
    std::cout<<"  "<<std::scientific<<std::setw(width)<<rho_;
    std::cout<<"  "<<                 std::setw(5)<<MAXDEPTH;
    std::cout<<"  "<<                 std::setw(width)<<num_oct;
    std::cout<<"  |";
    std::cout<<"  "<<                 std::setw(4)<<CHEB_DEG;
    std::cout<<"  "<<                 std::setw(4)<<MUL_ORDER;
    std::cout<<"  "<<std::scientific<<std::setw(width)<<GMRES_TOL;
    std::cout<<"  |";
    std::cout<<"  "<<std::scientific<<std::setw(   23)<<SS_INF_ERR;
    std::cout<<"  "<<std::scientific<<std::setw(   23)<<SS_L2_ERR;
    std::cout<<"  "<<     std::fixed<<std::setw(width)<<time_ksp;
    std::cout<< "("<<                 std::setw(    4)<<iter_ksp<<")";
    std::cout<<"  |";
    std::cout<<"\n";
  }
  return 0;
}

int SolveStokes(int argc, char **args, FMMData& fmm_data) {

  PetscErrorCode ierr;
  PetscInitialize(&argc,&args,0,help);

  MPI_Comm comm=MPI_COMM_WORLD;
  PetscMPIInt    rank,size;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);

  // -------------------------------------------------------------------
  PetscOptionsGetInt (NULL,  "-test_case",&TEST_CASE  ,NULL);
  PetscOptionsGetInt (NULL,  "-vtk_order",&VTK_ORDER  ,NULL);
  PetscOptionsGetInt (NULL,        "-dof",&INPUT_DOF  ,NULL);
  PetscOptionsGetReal(NULL,       "-scal",& SCAL_EXP  ,NULL);
  PetscOptionsGetBool(NULL,   "-periodic",& PERIODIC  ,NULL);
  PetscOptionsGetBool(NULL,       "-tree",& TREE_ONLY ,NULL);

  PetscOptionsGetInt (NULL, "-max_depth" ,&MAXDEPTH   ,NULL);
  PetscOptionsGetInt (NULL, "-min_depth" ,&MINDEPTH   ,NULL);
  PetscOptionsGetReal(NULL,   "-ref_tol" ,&      TOL  ,NULL);
  PetscOptionsGetReal(NULL, "-gmres_tol" ,&GMRES_TOL  ,NULL);

  PetscOptionsGetInt (NULL,   "-fmm_q"   ,& CHEB_DEG  ,NULL);
  PetscOptionsGetInt (NULL,   "-fmm_m"   ,&MUL_ORDER  ,NULL);

  PetscOptionsGetInt (NULL, "-gmres_iter",& MAX_ITER  ,NULL);

  PetscOptionsGetInt (NULL,    "-pt_cnt" ,&  PT_CNT   ,NULL);
  PetscOptionsGetReal(NULL,    "-pt_rad" ,&  PT_RAD   ,NULL);
  PetscOptionsGetReal(NULL,"-jump_width" ,&JUMP_WIDTH ,NULL);

  PetscOptionsGetReal(NULL,       "-rho" ,&    rho_   ,NULL);
  PetscOptionsGetInt (NULL,        "-tn" ,&NTSTEP     ,NULL);

  // -------------------------------------------------------------------
  if(PT_CNT==1){
    pt_coord.push_back(0.5);
    pt_coord.push_back(0.5);
    pt_coord.push_back(0.5);
  }else if(PT_CNT>1){
    long long pt_cnt=0;
    std::vector<double> pt_coord_(PT_CNT*COORD_DIM,-1);
    for(size_t i=0;i<PT_CNT*COORD_DIM;i+=COORD_DIM){
      int max_iter=100000;
      while(max_iter>0){
        if(TEST_CASE==3){
          pt_coord_[i+0]=drand48()*(1+6*PT_RAD)-5*PT_RAD;
          pt_coord_[i+1]=drand48()*(1+4*PT_RAD)-2*PT_RAD;
          pt_coord_[i+2]=drand48()*(1+4*PT_RAD)-2*PT_RAD;
        }else{
          pt_coord_[i+0]=drand48()*(1-3*PT_RAD)+PT_RAD;
          pt_coord_[i+1]=drand48()*(1-2*PT_RAD)+PT_RAD;
          pt_coord_[i+2]=drand48()*(1-2*PT_RAD)+PT_RAD;
        }
        bool break_cond=true;
        for(size_t i_=0;i_<i;i_+=COORD_DIM){
          double r=0;
          r+=(pt_coord_[i_+0]-pt_coord_[i+0])*(pt_coord_[i_+0]-pt_coord_[i+0]);
          r+=(pt_coord_[i_+1]-pt_coord_[i+1])*(pt_coord_[i_+1]-pt_coord_[i+1]);
          r+=(pt_coord_[i_+2]-pt_coord_[i+2])*(pt_coord_[i_+2]-pt_coord_[i+2]);
          r=sqrt(r);
          if(r<2.1*PT_RAD) break_cond=false; // Overlap
        }
        if(break_cond) break;
        max_iter--;
      }
      if(max_iter==0) break;

      if(TEST_CASE==3){ // Pack
        for(int max_iter=0;max_iter<100000;max_iter++){
          double c[3]={pt_coord_[i+0],pt_coord_[i+1],pt_coord_[i+2]};
          pt_coord_[i+0]-=(drand48()-0.1)*0.001;
          pt_coord_[i+1]-=(drand48()-0.5)*0.001;
          pt_coord_[i+2]-=(drand48()-0.5)*0.001;
          bool overlap=false;
          if(pt_coord_[i+0]<-4*PT_RAD) overlap=true;
          if(pt_coord_[i+1]<-2*PT_RAD) overlap=true;
          if(pt_coord_[i+2]<-2*PT_RAD) overlap=true;
          if(pt_coord_[i+0]>1.0+2*PT_RAD) overlap=true;
          if(pt_coord_[i+1]>1.0+2*PT_RAD) overlap=true;
          if(pt_coord_[i+2]>1.0+2*PT_RAD) overlap=true;
          if(!overlap)
            for(size_t i_=0;i_<i;i_+=COORD_DIM){
              double r=0;
              r+=(pt_coord_[i_+0]-pt_coord_[i+0])*(pt_coord_[i_+0]-pt_coord_[i+0]);
              r+=(pt_coord_[i_+1]-pt_coord_[i+1])*(pt_coord_[i_+1]-pt_coord_[i+1]);
              r+=(pt_coord_[i_+2]-pt_coord_[i+2])*(pt_coord_[i_+2]-pt_coord_[i+2]);
              r=sqrt(r);
              if(r<2.1*PT_RAD) overlap=true; // Overlap
            }
          if(overlap){
            pt_coord_[i+0]=c[0];
            pt_coord_[i+1]=c[1];
            pt_coord_[i+2]=c[2];
          }
        }
      }
      pt_cnt++;
    }
    MPI_Bcast(&pt_coord_[0], PT_CNT*COORD_DIM, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&pt_cnt, 1, MPI_LONG_LONG, 0, comm);

    for(size_t i=0;i<pt_cnt*COORD_DIM;i++){
      pt_coord.push_back(pt_coord_[i]);
    }
    PT_CNT=pt_coord.size()/COORD_DIM;
    if(rank==0) std::cout<<"num_sph = "<<PT_CNT<<'\n';
  }

  { //f_max
    switch (TEST_CASE)
    {
      case 0: // Constant coefficient
        {
          double L=125;
          double a=sqrt((2.0-sqrt(11.0/4.0)));
          f_max=fabs(-4*a*exp(-a*a)*(2*a*a - 5)*sqrt(L*L*L));
        }
        break;
      case 1: // Synthetic variable coefficient
        {
          double L=500;
          double a=sqrt((2.0-sqrt(11.0/4.0)));
          f_max=fabs(-4*a*exp(-a*a)*(2*a*a - 5)*sqrt(L*L*L));
          f_max+=fabs(exp(-0.5)*sqrt(2*L)*rho_);
        }
        break;
      case 2: // Sphere
        {
          f_max=std::max(fabs(rho_),1.0);
        }
        break;
      case 3: // Porous media
        {
          f_max=std::max(fabs(rho_),1.0);
        }
        break;
      default:
        break;
    }
  }

  // Initialize FMM
  FMM_Init(comm, &fmm_data);

  {
    /* -------------------------------------------------------------------
       Compute the matrix and right-hand-side vector that define
       the linear system, Ax = b.
       ------------------------------------------------------------------- */
    if(TREE_ONLY){
      pvfmm::Profile::print(&comm);
      disp_result();
      ierr = PetscFinalize();
      return 0;
    }

    PetscInt       m,n, M,N;
    m=fmm_data.m; // local rows
    n=fmm_data.n; // local columns
    M=fmm_data.M; // global rows
    N=fmm_data.N; // global columns

    pvfmm::Profile::Tic("FMMCreateShell",&comm,true);
    Mat A;
    { // Create Matrix. A
      FMMCreateShell(&fmm_data, &A);
      MatShellSetOperation(A,MATOP_MULT,(void(*)(void))mult);
    }
    pvfmm::Profile::Toc();

    Vec f,x,b;
    { // Create vectors
      VecCreateMPI(comm,n,PETSC_DETERMINE,&f);
      VecCreateMPI(comm,m,PETSC_DETERMINE,&b); // b=G[f]
      VecCreateMPI(comm,n,PETSC_DETERMINE,&x); // Ax=b
    }

    pvfmm::Profile::Tic("Input_Vector_f",&comm,true);
    { // Create Input Vector. f
      tree2vec(fmm_data,f);
    }
    pvfmm::Profile::Toc();

    pvfmm::Profile::Tic("Input_Vector_b",&comm,true);
    { // Create Input Vector. b = G[ f ]
      fmm_data.tree->SetupFMM(fmm_data.fmm_mat);
      fmm_data.tree->RunFMM(); // Compute: G[ f ]
      fmm_data.tree->Copy_FMMOutput();
      tree2vec(fmm_data,b);
    }
    pvfmm::Profile::Toc();

    // Create solution vector
    pvfmm::Profile::Tic("Initial_Vector_x",&comm,true);
    ierr = VecDuplicate(f,&x);CHKERRQ(ierr);
    pvfmm::Profile::Toc();

    // Create linear solver context
    pvfmm::Profile::Tic("KSPCreate",&comm,true);
    KSP ksp  ; ierr = KSPCreate(PETSC_COMM_WORLD,&ksp  );CHKERRQ(ierr);
    pvfmm::Profile::Toc();

    // Set operators. Here the matrix that defines the linear system
    // also serves as the preconditioning matrix.
    pvfmm::Profile::Tic("KSPSetOperators",&comm,true);
    ierr = KSPSetOperators(ksp  ,A  ,A  ,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
    pvfmm::Profile::Toc();

    // Set runtime options
    KSPSetType(ksp  ,KSPGMRES);
    KSPSetNormType(ksp  , KSP_NORM_UNPRECONDITIONED);
    KSPSetTolerances(ksp  ,GMRES_TOL  ,PETSC_DEFAULT,PETSC_DEFAULT,MAX_ITER  );
    KSPGMRESSetRestart(ksp  , MAX_ITER  );
    ierr = KSPSetFromOptions(ksp  );CHKERRQ(ierr);

    // -------------------------------------------------------------------
    // Solve the linear system
    // -------------------------------------------------------------------
    pvfmm::Profile::Tic("KSPSolve",&comm,true);
    time_ksp=-omp_get_wtime();
    ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
    MPI_Barrier(comm);
    time_ksp+=omp_get_wtime();
    pvfmm::Profile::Toc();

    // View info about the solver
    KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

    // -------------------------------------------------------------------
    // Check solution and clean up
    // -------------------------------------------------------------------

    // Iterations
    PetscInt       its;
    ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Iterations %D\n",its);CHKERRQ(ierr);
    iter_ksp=its;

    // Compute error
    stokes_err(A,x);

    vec2tree(x, fmm_data);

    // Free work space.  All PETSc objects should be destroyed when they
    // are no longer needed.
    ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
    ierr = VecDestroy(&x);CHKERRQ(ierr);
    ierr = VecDestroy(&b);CHKERRQ(ierr);
    ierr = MatDestroy(&A);CHKERRQ(ierr);
    disp_result();
    ierr = PetscFinalize();
  }
  return 0;
}

int main(int argc,char **args){
  MPI_Init(&argc,&args);

  MPI_Comm comm=MPI_COMM_WORLD;
  int np, myrank;
  MPI_Comm_size(comm, &np);
  MPI_Comm_rank(comm, &myrank);

  // =========================================================================
  // SIMULATION PARAMETERS
  // =========================================================================
  parse_command_line_options(argc, args);

  int test = strtoul(commandline_option(argc, args, "-test",     "1", false,
                                        "-test <int> = (1)    : 1) Gaussian profile 2) Zalesak disk"),NULL,10);

  int merge = strtoul(commandline_option(argc, args, "-merge",     "1", false,
                                         "-merge <int> = (1)    : 1) no merge 2) complete merge 3) Semi-Merge"),NULL,10);

  tbslas::SimConfig* sim_config = tbslas::SimConfigSingleton::Instance();
  pvfmm::Profile::Enable(sim_config->profile);

  // =========================================================================
  // PRINT METADATA
  // =========================================================================
  if (!myrank) {
    MetaData_t::Print();
  }

  FMMData fmm_data;
  SolveStokes(argc, args, fmm_data);

  //=====================================================================
  // PREPARE THE VELOCITY FIELD
  //=====================================================================
  FMM_Tree_t* tvel_curr = fmm_data.tree;
  // tbslas::Profile<double>::Tic("convert_vel",false,5);
  switch(TEST_CASE) {
    case 0:
    case 1:
    case 2:
    case 3:
      // convert velocity tree
      {
        FMM_Mat_t::FMMNode_t* n_curr = tvel_curr->PostorderFirst();
        while (n_curr != NULL) {
          if(!n_curr->IsGhost() && n_curr->IsLeaf()) break;
          n_curr = tvel_curr->PostorderNxt(n_curr);
        }
        while (n_curr != NULL) {
          if (n_curr->IsLeaf() && !n_curr->IsGhost()) {
            n_curr->ChebData()[0] = n_curr->ChebData()[0] - 1;
          }
          n_curr = tvel_curr->PostorderNxt(n_curr);
        }
      } break;
    default:
      printf("Undefined test case.\n");
      return 1;
  }
  //pvfmm::Profile<double>::Toc();

  char out_name_buffer[300];
  snprintf(out_name_buffer, sizeof(out_name_buffer),
           "%s/stokes_vel_%d_", tbslas::get_result_dir().c_str(), 0);
  tvel_curr->Write2File(out_name_buffer, VTK_ORDER);

  // ======================================================================
  // TBSLAS
  // ======================================================================
   {
    // =========================================================================
    // TEST CASE
    // =========================================================================
    switch(test) {
      case 1:
        fn_con = get_gaussian_field_3d_wrapper_01<double,3>;
        break;
      case 2:
        fn_con = get_gaussian_field_cylinder_atT<double,3>;
        break;
      case 3:
        fn_con = get_slotted_cylinder_atT<double,3>;
        break;
      case 4:
        fn_con = get_gaussian_field_cylinder_atT<double,3>;
        break;
    }

    // =========================================================================
    // SIMULATION PARAMETERS
    // =========================================================================
    // sim_config->vtk_filename_prefix   = "advection";
    sim_config->vtk_filename_variable = "conc01";
    sim_config->bc                    = fmm_data.bndry;

    // =========================================================================
    // INIT THE CONCENTRATION TREE
    // =========================================================================
    tcurr = 0;
    FMM_Tree_t tconc_curr(comm);
    tbslas::ConstructTree<FMM_Tree_t>(sim_config->tree_num_point_sources,
                                      sim_config->tree_num_points_per_octanct,
                                      sim_config->tree_chebyshev_order,
                                      sim_config->tree_max_depth,
                                      sim_config->tree_adap,
                                      sim_config->tree_tolerance,
                                      comm,
                                      fn_con,
                                      1,
                                      tconc_curr);
    // char out_name_buffer[300];
    if (sim_config->vtk_save_rate) {
      snprintf(out_name_buffer,
               sizeof(out_name_buffer),
               sim_config->vtk_filename_format.c_str(),
               tbslas::get_result_dir().c_str(),
               sim_config->vtk_filename_variable.c_str(),
               0);
      tconc_curr.Write2File(out_name_buffer, sim_config->vtk_order);
    }

    // =========================================================================
    // RUN
    // =========================================================================
    // set the input_fn to NULL -> needed for adaptive refinement
    std::vector<FMMNode_t*>  ncurr_list = tconc_curr.GetNodeList();
    for(int i = 0; i < ncurr_list.size(); i++) {
      ncurr_list[i]->input_fn = (void (*)(const double* , int , double*))NULL;
    }

    switch(merge) {
      case 2:
        pvfmm::Profile::Tic("CMerge", &sim_config->comm, false, 5);
        tbslas::MergeTree(*tvel_curr, tconc_curr);
        pvfmm::Profile::Toc();
        break;
      case 3:
        pvfmm::Profile::Tic("SMerge", &sim_config->comm, false, 5);
        tbslas::SemiMergeTree(*tvel_curr, tconc_curr);
        pvfmm::Profile::Toc();
        break;
    }

    int timestep = 1;
    for (; timestep < sim_config->total_num_timestep+1; timestep++) {

      // (SEMI) MERGE TO FIX IMBALANCE
      switch(merge) {
        case 2:
          pvfmm::Profile::Tic("CMerge", &sim_config->comm, false, 5);
          tbslas::MergeTree(*tvel_curr, tconc_curr);
          pvfmm::Profile::Toc();
          break;
        case 3:
          pvfmm::Profile::Tic("SMerge", &sim_config->comm, false, 5);
          tbslas::SemiMergeTree(*tvel_curr, tconc_curr);
          pvfmm::Profile::Toc();
          break;
      }

      pvfmm::Profile::Tic("SL", &sim_config->comm, false, 5);
      tbslas::SolveSemilagInSitu(*tvel_curr,
                                 tconc_curr,
                                 timestep,
                                 sim_config->dt,
                                 sim_config->num_rk_step);
      pvfmm::Profile::Toc();

      // refine the tree according to the computed values
      pvfmm::Profile::Tic("RefineTree", &sim_config->comm, false, 5);
      tconc_curr.RefineTree();
      pvfmm::Profile::Toc();

      //Write2File
      if (sim_config->vtk_save_rate) {
        if (timestep % sim_config->vtk_save_rate == 0)
        tconc_curr.Write2File(tbslas::GetVTKFileName(timestep, sim_config->vtk_filename_variable).c_str(), sim_config->vtk_order);
      }

    }  // end for

    //Output Profiling results.
    pvfmm::Profile::print(&comm);
  }


   if (test == 1) {
    // =========================================================================
    // TEST CASE
    // =========================================================================
     fn_con = get_gaussian_field_3d_wrapper_02<double,3>;

    // =========================================================================
    // SIMULATION PARAMETERS
    // =========================================================================
    // sim_config->vtk_filename_prefix   = "advection";
    sim_config->vtk_filename_variable = "conc02";
    sim_config->bc                    = fmm_data.bndry;

    // =========================================================================
    // INIT THE CONCENTRATION TREE
    // =========================================================================
    tcurr = 0;
    FMM_Tree_t tconc_curr(comm);
    tbslas::ConstructTree<FMM_Tree_t>(sim_config->tree_num_point_sources,
                                      sim_config->tree_num_points_per_octanct,
                                      sim_config->tree_chebyshev_order,
                                      sim_config->tree_max_depth,
                                      sim_config->tree_adap,
                                      sim_config->tree_tolerance,
                                      comm,
                                      fn_con,
                                      1,
                                      tconc_curr);
    // char out_name_buffer[300];
    if (sim_config->vtk_save_rate) {
      snprintf(out_name_buffer,
               sizeof(out_name_buffer),
               sim_config->vtk_filename_format.c_str(),
               tbslas::get_result_dir().c_str(),
               sim_config->vtk_filename_variable.c_str(),
               0);
      tconc_curr.Write2File(out_name_buffer, sim_config->vtk_order);
    }

    // =========================================================================
    // RUN
    // =========================================================================
    // set the input_fn to NULL -> needed for adaptive refinement
    std::vector<FMMNode_t*>  ncurr_list = tconc_curr.GetNodeList();
    for(int i = 0; i < ncurr_list.size(); i++) {
      ncurr_list[i]->input_fn = (void (*)(const double* , int , double*))NULL;
    }

    switch(merge) {
      case 2:
        pvfmm::Profile::Tic("CMerge", &sim_config->comm, false, 5);
        tbslas::MergeTree(*tvel_curr, tconc_curr);
        pvfmm::Profile::Toc();
        break;
      case 3:
        pvfmm::Profile::Tic("SMerge", &sim_config->comm, false, 5);
        tbslas::SemiMergeTree(*tvel_curr, tconc_curr);
        pvfmm::Profile::Toc();
        break;
    }

    int timestep = 1;
    for (; timestep < sim_config->total_num_timestep+1; timestep++) {

      // (SEMI) MERGE TO FIX IMBALANCE
      switch(merge) {
        case 2:
          pvfmm::Profile::Tic("CMerge", &sim_config->comm, false, 5);
          tbslas::MergeTree(*tvel_curr, tconc_curr);
          pvfmm::Profile::Toc();
          break;
        case 3:
          pvfmm::Profile::Tic("SMerge", &sim_config->comm, false, 5);
          tbslas::SemiMergeTree(*tvel_curr, tconc_curr);
          pvfmm::Profile::Toc();
          break;
      }

      pvfmm::Profile::Tic("SL", &sim_config->comm, false, 5);
      tbslas::SolveSemilagInSitu(*tvel_curr,
                                 tconc_curr,
                                 timestep,
                                 sim_config->dt,
                                 sim_config->num_rk_step);
      pvfmm::Profile::Toc();

      // refine the tree according to the computed values
      pvfmm::Profile::Tic("RefineTree", &sim_config->comm, false, 5);
      tconc_curr.RefineTree();
      pvfmm::Profile::Toc();

      //Write2File
      if (sim_config->vtk_save_rate) {
        if (timestep % sim_config->vtk_save_rate == 0)
          tconc_curr.Write2File(tbslas::GetVTKFileName(timestep, sim_config->vtk_filename_variable).c_str(),
                                sim_config->vtk_order);
      }

    }  // end for

    //Output Profiling results.
    pvfmm::Profile::print(&comm);
   }

   if (test == 1) {
     // =========================================================================
     // TEST CASE
     // =========================================================================
     fn_con = get_gaussian_field_3d_wrapper_03<double,3>;

     // =========================================================================
     // SIMULATION PARAMETERS
     // =========================================================================
     // sim_config->vtk_filename_prefix   = "advection";
     sim_config->vtk_filename_variable = "conc03";
     sim_config->bc                    = fmm_data.bndry;

     // =========================================================================
     // INIT THE CONCENTRATION TREE
     // =========================================================================
     tcurr = 0;
     FMM_Tree_t tconc_curr(comm);
     tbslas::ConstructTree<FMM_Tree_t>(sim_config->tree_num_point_sources,
                                       sim_config->tree_num_points_per_octanct,
                                       sim_config->tree_chebyshev_order,
                                       sim_config->tree_max_depth,
                                       sim_config->tree_adap,
                                       sim_config->tree_tolerance,
                                       comm,
                                       fn_con,
                                       1,
                                       tconc_curr);

     if (sim_config->vtk_save_rate) {
       snprintf(out_name_buffer,
                sizeof(out_name_buffer),
                sim_config->vtk_filename_format.c_str(),
                tbslas::get_result_dir().c_str(),
                sim_config->vtk_filename_variable.c_str(),
                0);
       tconc_curr.Write2File(out_name_buffer, sim_config->vtk_order);
     }

     // =========================================================================
     // RUN
     // =========================================================================
     // set the input_fn to NULL -> needed for adaptive refinement
     std::vector<FMMNode_t*>  ncurr_list = tconc_curr.GetNodeList();
     for(int i = 0; i < ncurr_list.size(); i++) {
       ncurr_list[i]->input_fn = (void (*)(const double* , int , double*))NULL;
     }

     switch(merge) {
       case 2:
         pvfmm::Profile::Tic("CMerge", &sim_config->comm, false, 5);
         tbslas::MergeTree(*tvel_curr, tconc_curr);
         pvfmm::Profile::Toc();
         break;
       case 3:
         pvfmm::Profile::Tic("SMerge", &sim_config->comm, false, 5);
         tbslas::SemiMergeTree(*tvel_curr, tconc_curr);
         pvfmm::Profile::Toc();
         break;
     }

     int timestep = 1;
     for (; timestep < sim_config->total_num_timestep+1; timestep++) {

       // (SEMI) MERGE TO FIX IMBALANCE
       switch(merge) {
         case 2:
           pvfmm::Profile::Tic("CMerge", &sim_config->comm, false, 5);
           tbslas::MergeTree(*tvel_curr, tconc_curr);
           pvfmm::Profile::Toc();
           break;
         case 3:
           pvfmm::Profile::Tic("SMerge", &sim_config->comm, false, 5);
           tbslas::SemiMergeTree(*tvel_curr, tconc_curr);
           pvfmm::Profile::Toc();
           break;
       }

       pvfmm::Profile::Tic("SL", &sim_config->comm, false, 5);
       tbslas::SolveSemilagInSitu(*tvel_curr,
                                  tconc_curr,
                                  timestep,
                                  sim_config->dt,
                                  sim_config->num_rk_step);
       pvfmm::Profile::Toc();

       // refine the tree according to the computed values
       pvfmm::Profile::Tic("RefineTree", &sim_config->comm, false, 5);
       tconc_curr.RefineTree();
       pvfmm::Profile::Toc();

       //Write2File
       if (sim_config->vtk_save_rate) {
         if (timestep % sim_config->vtk_save_rate == 0)
           tconc_curr.Write2File(tbslas::GetVTKFileName(timestep, sim_config->vtk_filename_variable).c_str(), sim_config->vtk_order);
       }

     }  // end for

     //Output Profiling results.
     pvfmm::Profile::print(&comm);
   }

  // Delete fmm data
  FMMDestroy(&fmm_data);
  MPI_Finalize();
  return 0;
}

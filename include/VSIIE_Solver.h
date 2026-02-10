/*
This is the header file of the C++ class "VSIIE_Solver" which implements the Variable Stepsize 
Implicit-Implicit-Explicit solvers with support for
orders 1 to 4. The solvers apply adaptive time-stepping techniques to optimize the computational
process where the accuracy is ensured by error checking. The following public functions are highlighted 
in this class:

1)  "Time_Step": describes the computations necessary to perform a time step of a
kth-order VSIIE method. This function has an argument which is a pointer to a IVP ODE
object.

2)  "Variable_Time_Step" function implements one coarse time step and two fine time
steps followed by the approximation of the local truncation error.

3)  "Adaptive_dt_Integrate" implements the entire adaptive time-stepping method 
for a single time step of k-step VSIIE method.

VSIIE Copyright (C) 2025 Jose Miguel Mantas Ruiz (jmmantas@ugr.es) and Raed Ali Mara'Beh (raedmaraabeh@gmail.com)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/



#ifndef VSIIE_SOLVER_H
#define VSIIE_SOLVER_H


#include "IVP_ODE.h"

using namespace std;
#include <assert.h>



//*****************************************************************
// Class for the IVP Solver VSIIE-order 
//*****************************************************************
class VSIIE_Solver {

//*****************************************************
// PRIVATE VARIABLES AND FUNCTIONS
private:
  #define max_order 4
  double coef_Y[max_order], coef_F1[max_order+1], 
        coef_F2[max_order+1], coef_F3[max_order];    
  unsigned order; // order of the VSIIE IVP solver
  IVP_ODE* IVP; // pointer to the particular IVP Object
  

  // Coefficients for evaluations of F1 and F2 for the next time point
  double coef_idx_F1, coef_idx_F2, h_coef_idx_F1, h_coef_idx_F2; 
  int neqn; // number of equations of the IVP
  LIS_MATRIX As,As2; // Auxiliary sparse matrices to store Jacobian matrices
  LIS_SOLVER solver; // LIS Iterative solver
  LIS_VECTOR b, x; // Auxiliary LIS vectors
  double * R; //Auxiliary vector to store the residual of the modified NR iteration
  int idx; //leading index
  int nnz; // Number of non-zero entries in Jacobian matrix

  //**************************************************************************
  // Compute the variable term (in a time step) of Right Hand Side at each Newton Iteration
  //**************************************************************************
  void compute_RHS(const double t, const double h,  double* Y1,double* R);

  //**************************************************************
  // Check important errors between sparse matrices A1 and A2
  void compare_matrix_csr(LIS_MATRIX A1, LIS_MATRIX A2);
  //**************************************************************

  //**************************************************************
  // Print values of a sparse matrix A
  void print_matrix_csr(LIS_MATRIX A);
  //**************************************************************
  //**************************************************************************
  // Compute Y=Y+a*F1(t,Y) 
  //**************************************************************************
  void  Add_aF1Y(const double t, const double a,double* Y);
  //**************************************************************************
  // Compute Y=Y+a*F2(t,Y) 
  //**************************************************************************
  void  Add_aF2Y(const double t, const double a,double* Y);
  //**************************************************************************
  // Compute Y=Y+a*F3(t,Y) 
  //**************************************************************************
  void  Add_aF3Y(const double t, const double a,double* Y);

  //***************************************************
  // Init CSR Sparse neqn x neqn LIS matrix
  void init_CSR_LIS_matrix(LIS_MATRIX *A);
  //***************************************************

  //**************************************************************************
  // Compute the constant term (in a time step) of Right Hand Side 
  // at each Newton Iteration
  //**************************************************************************
  void compute_RHS0(const double t, const double * h_vector, double** Y);

  //***************************************************
  // Init every CSR Sparse LIS matrix
  void init_matrices();
  //***************************************************

  //******************************************************
  // Compute the values for every CSR Sparse LIS matrix
  //******************************************************
  void compute_matrix(const double t, const double h, double ** Y);
  //***************************************************


  //******************************************************
  // Update coefficients of fixed stepsize-order scheme
  //******************************************************
  void Fixed_update_coefs(const unsigned order); 
  //******************************************************

  //******************************************************
  // Update coefficients of VSIIE-order scheme
  //******************************************************
  void Update_coefs(const unsigned order,const double * h_vector);
  //******************************************************


  //******************************************************
  // Update coefficients for the  SBDF-4 scheme
  //******************************************************
  void Update_coefs_SBDF4(const double * h_vector);
  //******************************************************

  //******************************************************
  // Compute LTE estimate from:
  // the coarse solution Yc_sol and fine solution Yf_sol
  // and the step size coarse vector h_vector
  //******************************************************
  double compute_LTE(const double * h_vector, 
                      double * Yc_sol, double * Yf_sol, double * LTE);
  //******************************************************

  //**************************************************************************
  // Update intermediate vectors Yf and Y using the current vector solution Y1 
  // for the adaptive time stepper before the next integration step
  //**************************************************************************
  void Update_adaptive_intermediate_vectors(double *h_vector_half, double *h_vector_half2, 
                                            double *h_vector, double ** Yf, double ** Y, 
                                            double *Y1);

  //**************************************************************************
  // Update the coarse vector Y1_c from the fine vector Y1_f 
  // by using Richardson extrapolation 
  //**************************************************************************
  void Extrapolate(const double * h_vector, const double *Y1_f, double *Y1_c);
   



//******************************************************
// PUBLIC METHODS
public:
//******************************************************

  //******************************************************
  // Constructor of the class VSIIE_Solver
  //******************************************************
  VSIIE_Solver (const int _order, IVP_ODE* _IVP);

  //*******************************************************************
  // Function implementing one time step 
  //*******************************************************************
  void  Time_Step(const double t, const double *h_vector, double** Y, double* Y1, 
const bool variable_tstep, int & newton_iters);

  //**************************************************************************
  // Update intermediate vectors before the next integration step
  //**************************************************************************
  void Update_intermediate_vectors(double** Y, double* Y1);
  //******************************************************

  //******************************************************
  // Destructor of the class VSIIE_Solver
  //******************************************************
  ~VSIIE_Solver();

  //***************************************************
  // Function implementing the order 1-4 VSIIE Time Integrator
  // It assumes a constant time step h
  //***************************************************
  void Const_dt_Integrate(const double t0, const double tf,
                               const double h, double** Y_, 
                               double* Y1);


  //***************************************************
  // Function implementing the 4th-order SBDF Time Integrator
  // It assumes a constant time step h
  //***************************************************
  void Const_dt_Integrate_SBDF4(const double t0, const double tf,
    const double h, double** Y_, double* Y1);
  //***************************************************
                                 

  //***************************************************
  // Coarse and Fine Time steps using Variable VSIIE scheme
  //***************************************************
  void Variable_Time_Step(const double t, double * h_vector, 
                                     double * h_vector_half, double * h_vector_half2, double ** Y, 
                                     double ** Yf, double * Y1, double * LTE, double *epsilon_c, 
                                     int& total_iters);
  //***************************************************

  //***************************************************
// Function implementing the order 1-4 VSIIE Time Integrator
// It assumes an adaptive time step
//***************************************************
void Adaptive_dt_Integrate(const double t0, const double tf,
    const double h, double** Y_init, double ** Yf_init, double* Y1, 
    const double tol, int *nsteps, int * n_isteps, double & av_iters,
    double alpha=0.8, double eta_min= 0.5, 
    double eta_max = 4);        
//***************************************************


  //***************************************************
  // Store the results of an experiment in a data file
  //***************************************************
  void Store_result(const double tf, const double* Y1, 
                    const double tol);
  //***************************************************
   
};


#endif
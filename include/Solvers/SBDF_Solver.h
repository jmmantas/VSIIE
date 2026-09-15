/*
This is the header file of the C++ class "SBDF_Solver" which implements the VSSBDF solvers with support for
orders 1 to 4. The solvers apply adaptive time-stepping techniques to optimize the computational
process where the accuracy is ensured by error checking. The following public functions are highlighted 
in this class:

1)  "Time_Step": describes the computations necessary to perform a time step of a
kth-order SBDF method. This function has an argument which is a pointer to a IVP ODE
object.

2)  "Variable_Time_Step" function implements one coarse time step and two fine time
steps followed by the approximation of the local truncation error.

3)  "Adaptive_dt_SBDF_Integrate0" implements the entire adaptive time-stepping method 
for a single time step of k-step VSSBDF method.


VSSBDF Copyright (C) 2024 Jose Miguel Mantas Ruiz (jmmantas@ugr.es) and Raed Ali Mara'Beh (raedmaraabeh@gmail.com)

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


#ifndef SBDF_SOLVER_H
#define SBDF_SOLVER_H


#include "IVP_ODE.h"

using namespace std;
#include <assert.h>





//*****************************************************************
// Class for the IVP Solver SBDF-order 
//*****************************************************************
class SBDF_Solver {

//*****************************************************
// PRIVATE VARIABLES AND FUNCTIONS
private:
//*****************************************************

  unsigned order; // order of the SBDF IVP solver
  IVP_ODE* IVP; // pointer to the particular IVP Object
  // Coefficients of SBDF-order method 
  double coef_G[4]={1.0, 2.0/3, 6.0/11, 12.0/25};
  double coef_F[4][4]={    
    {1.0, 0.0, 0.0, 0.0 },
    {-2.0/3, 4.0/3, 0.0, 0.0 },
    { 6.0/11,  -18.0/11, 18.0/11, 0.0 },
    {-12.0/25, 48.0/25, -72.0/25, 48.0/25 }
  };
  double coef_Y[4][4]={    
    {1.0, 0.0, 0.0, 0.0 },
    {-1.0/3, 4.0/3, 0.0, 0.0 },
    {2.0/11, -9.0/11, 18.0/11, 0.0 },    
    {-3.0/25, 16.0/25, -36.0/25, 48.0/25 }
  };

  double stability_factor[4]= {5.0 ,2.414,1.501,1.101};
 // double stability_factor[4] = { 1.101 ,1.101,1.101,1.101 };


  double coef_idx;
  int neqn; // number of equations of the IVP
  LIS_MATRIX As; // Auxiliary sparse matrix to store Jacobian
  LIS_SOLVER solver; // LIS Iterative solver
  LIS_VECTOR b, x; // Auxiliary LIS vectors
  double * R; //Auxiliary vector to store the residual of the modified NR iteration
  int idx;
  int nnz; // Number of non-zero entries in Jacobian matrix

  //**************************************************************************
  // Compute the variable term (in a time step) of Right Hand Side at each Newton Iteration
  //**************************************************************************
  void  compute_RHS_SBDF(const double t, const double h, const double* R0, double* Y1,double* RHS);

  //**************************************************************************
  // Compute R=R+a*F(t,Y) 
  //**************************************************************************
  void  Add_aFY(const double t, const double a,double* Y);

  //***************************************************
  // Init CSR Sparse neqn x neqn LIS matrix
  void init_CSR_LIS_matrix(LIS_MATRIX *A);
  //***************************************************

  //**************************************************************************
  // Compute the constant term (in a time step) of Right Hand Side 
  // at each Newton Iteration
  //**************************************************************************
  void  compute_RHS0_SBDF(const double t, const double * h_vector, double** Y);


  //******************************************************
  // Compute the values for every CSR Sparse LIS matrix
  //******************************************************
  void compute_matrices(const double t, const double h, double ** Y);
  //***************************************************

  //******************************************************
  // Update coefficients of VSSBDF-order scheme
  //******************************************************
  void Update_coefs(const int order,const double * h_vector);
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
  // Constructor of the class SBDF
  //******************************************************
  SBDF_Solver (const int _order, IVP_ODE* _IVP);

  //*******************************************************************
// Function implementing one time step using the 
// SBDF-order Time Integrator (order=1,..,4)
//*******************************************************************
void  Time_Step(const double t, const double *h_vector, double** Y, double* Y1, 
const bool variable_tstep, int & newton_iters);

  //**************************************************************************
  // Update intermediate vectors before the next integration step
  //**************************************************************************
  void Update_intermediate_vectors(double** Y, double* Y1);
  //******************************************************

  //******************************************************
  // Destructor of the class SBDF
  //******************************************************
  ~SBDF_Solver();

  //***************************************************
  // Function implementing the order 1-4 SBDF Time Integrator
  // It assumes a constant time step h
  //***************************************************
  void Const_dt_SBDF_Integrate(const double t0, const double tf,
                               const double h, double** Y_, 
                               double* Y1);

  //***************************************************
  // Coarse and Fine Time steps using Variable SBDF scheme
  //***************************************************
  //void Variable_Time_Step(const double t, double * h_vector, 
  //                double * h_vector_half, double * h_vector_half2, double ** Y, 
  //                double ** Yf, double * Y1, double * LTE, double *epsilon_c);
  //***************************************************

  //***************************************************
  // Function implementing the order 1-4 SBDF Time Integrator
  // It assumes a adaptive time step taking into account a stability factor
  //***************************************************
  //void Adaptive_dt_Integrate(const double t0, const double tf,
  //                                const double h, double** Y_init, 
  //                                double ** Yf_init, double* Y1, 
  //                                const double tol, int *nsteps, 
  //                                int * n_isteps);        
  //***************************************************

//***************************************************
  // Coarse and Fine Time steps using Variable SBDF scheme
  //***************************************************
  void Variable_Time_Step(const double t, double * h_vector, 
                                     double * h_vector_half, double * h_vector_half2, double ** Y, 
                                     double ** Yf, double * Y1, double * LTE, double *epsilon_c, 
                                     int& total_iters);

  ///***************************************************
// Function implementing the order 1-4 SBDF Time Integrator
// It assumes a adaptive time step
//***************************************************
void Adaptive_dt_Integrate(const double t0, const double tf,
    const double h, double** Y_init, double ** Yf_init, double* Y1, 
    const double tol, int *nsteps, int * n_isteps, int * n_iters,
    double alpha=0.8, double eta_min= 0.5, 
    double eta_max = 4); 


//***************************************************
// Store the results of an experiment in a data file
//***************************************************
void Store_result(const double tf, 
                              const double* Y1, 
                              const double tol);
//***************************************************


};


#endif
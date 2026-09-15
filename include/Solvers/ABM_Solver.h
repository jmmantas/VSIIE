/*
This is the header file of the C++ class "ABM_Solver" which implements the Variable stepsize Adams-Bashforth-Moulton (ABM) solvers 
with support for orders 1 to 4. The solvers apply adaptive time-stepping techniques to optimize the computational
process where the accuracy is ensured by error checking. The following public functions are highlighted 
in this class:

1)  "Time_Step": describes the computations necessary to perform a time step of a
kth-order ABM method. This function has an argument which is a pointer to a IVP ODE
object.

2)  "Variable_Time_Step" function implements the approximation of the local truncation error.

3)  "Adaptive_dt_Integrate" implements the entire adaptive time-stepping method 
for a single time step of k-step ABM method.


VSIIE Copyright (C) 2026 Jose Miguel Mantas Ruiz (jmmantas@ugr.es) and Raed Ali Mara'Beh (raedmaraabeh@gmail.com)

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


#ifndef ABM_SOLVER_H
#define ABM_SOLVER_H


#include "IVP_ODE.h"
using namespace std;
#include <assert.h>




//*****************************************************************
// Class for the IVP Solver ABM-order 
//*****************************************************************
class ABM_Solver {

//*****************************************************
// PRIVATE VARIABLES AND FUNCTIONS
private:
//*****************************************************

  unsigned order; // order of the ABM IVP solver
  IVP_ODE* IVP; // pointer to the particular IVP Object
  
  // Coefficients of ABM-order method 
  double coef_AB[4][4]={    
    {1.0, 0.0, 0.0, 0.0 },
    {-0.5, 1.5, 0.0, 0.0 },
    { 5.0/12, -16.0/12, 23.0/12,  0.0 },
    {-9.0/24 , 37.0/24, -59.0/24, 55.0/24  }
  };
  
  double coef_AM[4][4]={    
    {1.0    , 0.0   , 0.0   , 0.0 },
    {0.5    , 0.5   , 0.0   , 0.0 },
    {-1.0/12, 8.0/12, 5.0/12, 0.0 },    
    {1.0/24 , -5.0/24, 19.0/24, 9.0/24  }
  };
  // Coefficients for the LTE estimate
  double LTE_coef[4]={0.5, 1.0/6.0, 0.1, 19.0/270}; 

  int neqn; // number of equations of the IVP
  int idx;

  //******************************************************
  // Update coefficients of ABM-order scheme
  //******************************************************
  void Update_coefs(const int order,const double * h_vector);
  //******************************************************

  //******************************************************
  // Compute LTE estimate from the predictor solution Y_AB 
  // and the corrector solution Y_AM
  //******************************************************
  double compute_LTE(double * Y_AB, double * Y_AM);
  //******************************************************

 
//******************************************************
// PUBLIC METHODS
public:
//******************************************************

  //******************************************************
  // Constructor of the class ABM
  //******************************************************
  ABM_Solver (const int _order, IVP_ODE* _IVP);


  //**************************************************************************
  // Adams-Bashforth (AB) predictor with variable time step
  //**************************************************************************
  void AB_Predictor(const double t, double *h_vector, double ** Y, double *Y1); 

  //**************************************************************************
  // Adams-Moulton (AM) corrector with variable time step with variable time step
  // and predictor solution Y_AB and its evaluation FY_AB
  //**************************************************************************
  void  AM_corrector(const double t, const double * h_vector, double** Y, double *FY_AB, double *Y1);

  //*******************************************************************
  // Function implementing one time step using the 
  // ABM-order Time Integrator (order=1,..,4)
  // Generate the corrector solution Y1 and the 
  // error estimate of the local truncation error
  //*******************************************************************
  void  Time_Step(const double t, const double *h_vector, double** Y, double* Y1, 
                                   const bool variable_tstep, double & error_estimate);

  //**************************************************************************
  // Update intermediate vectors and stepsizes before the next integration step
  //**************************************************************************
  void Update_intermediate_vectors(double *h_vector, 
                                                      double** Y, double* Y1); 
  //***************************************************************************

  //******************************************************
  // Destructor of the class ABM
  //******************************************************
  ~ABM_Solver();

  //***************************************************
  // Function implementing the order 1-4 ABM Time Integrator
  // It assumes a constant time step h
  //***************************************************
  void Const_dt_Integrate(const double t0, const double tf,
                               const double h, double** Y_, 
                               double* Y1);

  
///***************************************************
  // Function implementing the order 1-4 ABM Time Integrator
  // It assumes a adaptive time step
  //**************************************************
  void Adaptive_dt_Integrate(const double t0, const double tf,
    const double h, double** Y_init, double* Y1, 
    const double tol, int *nsteps, int * n_isteps, 
    double alpha=0.9, double eta_min= 0.5, double eta_max = 5);

  //***************************************************
  // Store the results of an experiment in a data file
  //***************************************************
  void Store_result(const double tf, 
                              const double* Y1, 
                              const double tol);
  //***************************************************


};


#endif

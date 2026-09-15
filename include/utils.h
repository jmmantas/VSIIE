/*
These are auxiliary functions for the driver programs which use the solvers.

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


#ifndef Included_UTILS_H

#define Included_UTILS_H


#include "lis.h"
#include "IVP_ODE.h"
using namespace std;
using namespace std::chrono;
typedef std::chrono::time_point<high_resolution_clock, nanoseconds> time_ns;

//**************************************************************************
//AUXILIARY  FUNCTIONS
//**************************************************************************

//**************************************************************
// Check important errors between sparse matrices A1 and A2
void compare_matrix_csr(LIS_MATRIX A1, LIS_MATRIX A2);

//**************************************************************
// Print values of sparse matrix A
void print_matrix_csr(LIS_MATRIX A);


//**************************************************************
//1st, 2nd and 3rd order Runge-Kutta Time Integrator
//**************************************************************
void RK(const int nstages, IVP_ODE* IVP, const double t0,
    const double tf, const double h0, const double* Y0, double* Y1);
//**************************************************************

//***************************************************
//Print the norm of a vector Y of size neqn
//***************************************************
void print_vector(const int neqn, const double* Y, 
    const string name);
//***************************************************





//**************************************************************************
#endif
/*This file contains a C++ subclass of "IVP_ODE" called "IVP_ODE_stiff_brusselator" which represents  
a stiff variation of the one-dimensional Brusselator model.

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




#ifndef STIFF_BRUSSELATOR_H
#define STIFF_BRUSSELATOR_H


#include "IVP_ODE.h"

using namespace std;



//*****************************************************************
// Class for the IVP-ODE  which represents the stiff version of the 
// the Brusselator  model 
//*****************************************************************

class IVP_ODE_stiff_brusselator : public IVP_ODE {
private:
    int nx; // Number of grid points in the spatial domain
    double dtx; // Spatial step size
    double dtx2; //dtx*dtx
    double dtx_2; //2*dtx
    double DD; //alpha/dtx2
    double AA; //rho/dtx_2
    double alpha=0.01; // Diffusion coefficients for u, v, and w

    //*********************************************************
    double rho=0.01; // Advection coefficients for u, v, and w
    //*********************************************************

    double a = 0.6;  // Reaction parameters
    double b = 2.0; // Reaction parameters
    double eps=0.01; // Parameter to control stiffness of the reaction term for w
    double pi = 3.14159265358979323846;
 

    // Indexation function which maps 2D spatial coordinates (i,j)
    // to a 1D position in a vector
    inline int idx(const int i, const int j) { return i * 3 + j; }



public:
    //******************************************************
    // Constructor of the class IVP_ODE_cusp
    //******************************************************
    IVP_ODE_stiff_brusselator(const int nx_points);

    //******************************************************
    // Initialize stage vector Y0 with neqn components
    //******************************************************
    inline void init(double* Y0);


  //***************************************************
  //vector system function for the 1st stiff term DY=F1(t,Y)
  //***************************************************
  void F1 (const double t, const double *Y, double *DY);

  //***************************************************
  //vector system function for the 2nd stiff term DY=F2(t,Y)
  //***************************************************
  void F2 (const double t, const double *Y, double *DY);

  //***************************************************
  //vector system function for the nonstiff term DY=F3(t,Y)
  //***************************************************
  void F3 (const double t, const double *Y, double *DY);

  //******************************************************
  //Compute matrix I-c1*J1-c2*J2 where J1= Jacobian of F1 and J2 the Jacobian of F2
  // in the IVP_ODE_advdiff1d
  //******************************************************     
  void Compute_Matrix_Exact(const double t, double c1, double c2,
    double *Y, LIS_MATRIX As);

  //******************************************************
  //Compute matrix Jf= Jacobian of the function 
  // feval defining the IVP_ODE
  //******************************************************       
  void Compute_Feval_Jacobian_exact(const double t,double *Y, 
                                    LIS_MATRIX As);


//******************************************************
// Print the data corrersponding to a state vector Y
// on a file 
//****************************************************** 
void save_vector(const std::string& filename, const double *Y);

};


//*****************************************************
// Constructor of the class IVP_ODE_brusselator2D
//*****************************************************
IVP_ODE_stiff_brusselator::IVP_ODE_stiff_brusselator(const int nx_points) {
     IVP_name = "stiff_Brusselator";
    nx = nx_points;
    // Number of ODEs
    neqn = 3 * nx;

    // Number of non-zero elements in the Jacobian matrix
    //nnz = 9 * nx-12;
    nnz = 14 * nx - 12;
    dtx = 1.0 / (nx - 1); // Adjusted for boundary conditions
    dtx2 = dtx * dtx;
    dtx_2 = 2*dtx ;
    DD = alpha / dtx2;
    AA = -rho / dtx_2;
}


//******************************************************
    // Initialize stage vector Y0 with neqn components
//******************************************************
    inline void IVP_ODE_stiff_brusselator::init(double *Y0)
    {
      for (int i = 0; i < nx; i++) {
        double x = i * dtx;
        Y0[idx(i, 0)] = a + 0.1 * sin(pi * x); // Initial condition for u
        Y0[idx(i, 1)] = b / a + 0.1 * sin(pi * x); // Initial condition for v
        Y0[idx(i, 2)] = b + 0.1 * sin(pi * x); // Initial condition for w
      }
    }



//***************************************************
//vector system function for the nonstiff term DY=F3(t,Y) 
// corresponding to the advection term 
//***************************************************
void IVP_ODE_stiff_brusselator::F3 (const double t, const double *Y, double *DY)
{
  // Apply stationary boundary conditions at the first and last points 
  // explicitly for advection terms
  // Assuming no change due to advection at the boundaries
  DY[idx(0, 0)] = 0;  DY[idx(0, 1)] = 0;  DY[idx(0, 2)] = 0;
  DY[idx(nx - 1, 0)] = 0;   DY[idx(nx - 1, 1)] = 0;   DY[idx(nx - 1, 2)] = 0;

  for (int i = 1; i < (nx - 1); i++) {
    // Indexes for u, v, and w at grid point i
    const int i0 = idx(i, 0), i1 = i0 + 1, i2 = i0 + 2;
    // Calculate central differences for advection terms
    DY[i0] = AA*(Y[idx(i + 1, 0)] - Y[idx(i - 1, 0)]);
    DY[i1] = AA*(Y[idx(i + 1, 1)] - Y[idx(i - 1, 1)]);
    DY[i2] = AA*(Y[idx(i + 1, 2)] - Y[idx(i - 1, 2)]);
  }
}



//***************************************************
//vector system function for the stiff term DY=F1(t,Y)
//***************************************************
void IVP_ODE_stiff_brusselator::F1 (const double t, const double *Y, double *DY)
{
   // First, set the boundary conditions to zero derivative (stationary)
   DY[idx(0, 0)] = 0;  DY[idx(0, 1)] = 0;  DY[idx(0, 2)] = 0;
   DY[idx(nx - 1, 0)] = 0;   DY[idx(nx - 1, 1)] = 0;   DY[idx(nx - 1, 2)] = 0;
   
   // Apply the diffusion term for interior points
   for (int i = 1; i < (nx - 1); i++) {
    // Indexes for u, v, and w at grid point i
    const int i0 = idx(i, 0), i1 = i0 + 1, i2 = i0 + 2;
    // Apply the diffusion term using the Laplacian discretization for interior points
    DY[i0] = DD * (Y[idx(i - 1, 0)] - 2.0 * Y[i0] + Y[idx(i + 1, 0)]);
    DY[i1] = DD * (Y[idx(i - 1, 1)] - 2.0 * Y[i1] + Y[idx(i + 1, 1)]);
    DY[i2] = DD * (Y[idx(i - 1, 2)] - 2.0 * Y[i2] + Y[idx(i + 1, 2)]);
   } 

}
//***************************************************

//***************************************************
//vector system function for the 2nd stiff term DY=F2(t,Y)
//***************************************************
void IVP_ODE_stiff_brusselator::F2 (const double t, const double *Y, double *DY)
{
  // First, set the boundary conditions to zero derivative (stationary)
  DY[idx(0, 0)] = 0;  DY[idx(0, 1)] = 0;  DY[idx(0, 2)] = 0;
  DY[idx(nx - 1, 0)] = 0;   DY[idx(nx - 1, 1)] = 0;   DY[idx(nx - 1, 2)] = 0;
  
  // Reaction term for interior points
  for (int i = 1; i < nx-1; i++) {
    const int i0 = idx(i, 0), i1 = i0 + 1, i2 = i0 + 2;
    const double u = Y[i0], v = Y[i1], w = Y[i2];
    // Reaction terms
    DY[i0] = a - (w + 1) * u + u * u * v ; // Adding reaction term for u
    DY[i1] = w * u - u * u * v ; // Adding reaction term for v
    DY[i2] = (b - w) / eps - w * u; // Adding reaction term for w
  }

  

}
//***************************************************


//******************************************************
//Compute matrix I-c1*J1-c2*J2 where J1= Jacobian of the Diffusive Term
// and J2= Jacobian of the reactive term in the IVP_ODE_advdiff1d
//******************************************************       
void IVP_ODE_stiff_brusselator::Compute_Matrix_Exact(const double t, double c1, double c2,
                                                 double *Y, LIS_MATRIX As)     {
  LIS_INT k,row;

  LIS_INT *ptr=As->ptr;
  LIS_INT *index=As->index;
  LIS_SCALAR *value=As->value;
  
  // Coefficient for the diagonal elements
  const double Aii=1 + 2.0 * c1 * DD; 
  const double other = -c1*DD;    



  init_row_insertion(ptr, k, row);

  int i;
  double u,v,w;
  
  // row  (0,0) 
  i=0; u = Y[idx(i, 0)];v = Y[idx(i, 1)];w = Y[idx(i, 2)];
  double K2=c2*(2.0*u*v-w);
  double J21=K2-c2;  double J22=-c2*u*u; double J23=c2*u;
  new_entry(index, value, k, idx(i, 0), 1.0-J21);
  new_entry(index, value, k, idx(i, 1),  J22 );
  new_entry(index, value, k, idx(i, 2), J23    );
  next_row(ptr, k, row);
  // row  (0,1) 
  new_entry(index, value, k, idx(i, 0) , K2 );
  new_entry(index, value, k, idx(i, 1), 1.0-J22);
  new_entry(index, value, k, idx(i, 2), -J23);
  next_row(ptr, k, row);
  // row  (0,2)
  new_entry(index, value, k, idx(i, 0), c2*w);
  new_entry(index, value, k, idx(i, 2), 1.0-c2*(-1.0/eps-u));
  next_row(ptr, k, row);
  

  // row(i,j), i=1,...,nx-2
  for ( i = 1; i <= nx - 2; i++) {
    u = Y[idx(i, 0)]; v = Y[idx(i, 1)]; w = Y[idx(i, 2)];
    K2=c2*(2.0*u*v-w); J21=K2-c2;  J22=-c2*u*u; J23=c2*u;
    // row  (i,0)
    new_entry(index, value, k, idx(i-1, 0), other);
    new_entry(index, value, k, idx(i, 0), Aii-J21);
    new_entry(index, value, k, idx(i, 1),  J22 );
    new_entry(index, value, k, idx(i, 2), J23    );
    new_entry(index, value, k, idx(i+1, 0), other);
    next_row(ptr, k, row);
    // row  (i,1)
    new_entry(index, value, k, idx(i-1, 1), other); 
    new_entry(index, value, k, idx(i, 0) , K2 );
    new_entry(index, value, k, idx(i, 1), Aii-J22);
    new_entry(index, value, k, idx(i, 2), -J23);
    new_entry(index, value, k, idx(i+1, 1), other);
    next_row(ptr, k, row);
    // row  (i,2)
    new_entry(index, value, k, idx(i-1, 2), other);
    new_entry(index, value, k, idx(i, 0), c2*w);
    new_entry(index, value, k, idx(i, 2), Aii-c2*(-1.0/eps-u));
    new_entry(index, value, k, idx(i+1, 2), other);
    next_row(ptr, k, row);    
  }


  // row (nx-1,0)
  i=nx-1; u = Y[idx(i, 0)];v = Y[idx(i, 1)];w = Y[idx(i, 2)];
  K2=c2*(2.0*u*v-w); J21=K2-c2;  J22=-c2*u*u; J23=c2*u;
  new_entry(index, value, k, idx(i, 0), 1.0-J21);
  new_entry(index, value, k, idx(i, 1),  J22 );
  new_entry(index, value, k, idx(i, 2), J23    );
  next_row(ptr, k, row);
  // row  (nx-1,1) 
  new_entry(index, value, k, idx(i, 0) , K2 );
  new_entry(index, value, k, idx(i, 1), 1.0-J22);
  new_entry(index, value, k, idx(i, 2), -J23);
  next_row(ptr, k, row);
  // row  (nx-1,2)
  new_entry(index, value, k, idx(i, 0), c2*w);
  new_entry(index, value, k, idx(i, 2), 1.0-c2*(-1.0/eps-u));
  next_row(ptr, k, row);
 
}



//******************************************************
// Save the data corresponding to a state vector Y
// on a file 
//****************************************************** 
void IVP_ODE_stiff_brusselator::save_vector(const std::string& filename, const double *Y)
{
  std::ofstream file(filename);
  if (!file) {
    std::cerr << "Error opening the file: " << filename << std::endl;
    return;
  }
  for (int i = 0; i < nx; ++i){ 
    const double x_i=i*dtx;
    file << x_i<< "     " << Y[idx(i,0)] << std::endl; 
  }
  file.close();
}

#endif

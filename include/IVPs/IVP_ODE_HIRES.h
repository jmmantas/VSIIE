/*
This file contains a C++ subclass of "IVP_ODE" called "HIRES" which represents  
the HIRES (High Irradiance Response) benchmark model written in three-additive form.

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



#ifndef HIRES_H
#define HIRES_H


#include "IVP_ODE.h"

using namespace std;



//*****************************************************************
// Class for the IVP-ODE representing a 1D Advection-Diffusion model 
//*****************************************************************
class IVP_ODE_HIRES:public IVP_ODE {

private:
  


public:
  //*****************************************************
  // Constructor of the class IVP_ODE_advdiff1d
  IVP_ODE_HIRES(const int nx_points);
  //*****************************************************

  //******************************************************
  // Initialize stage vector Y0 with neqn components
  //******************************************************
  inline void init(double *Y0); 


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
  //Compute matrix I-a*J1-b*J2 where J1= Jacobian of F1 and J2 the Jacobian of F2
  // in the IVP_ODE_HIRES
  //******************************************************     
  void Compute_Matrix_Exact(const double t, double a, double b,
                            double *Y, LIS_MATRIX As);

  //******************************************************
  // Save the data corresponding to a state vector Y
  // on a file called filename 
  //****************************************************** 
  void save_vector(const std::string& filename, const double *Y);

};


//*****************************************************
// Constructor of the class IVP_ODE_HIRES
IVP_ODE_HIRES::IVP_ODE_HIRES(const int nx_points)
//*****************************************************
{ 
  IVP_name="HIRES";
  // Number of ODEs
  neqn=8;  
  // Number of non-zeros in Jacobian of I-c1*J1-c2*J2
  nnz=22;
  
}
  //******************************************************
  // Initialize stage vector Y0 with neqn components
  //******************************************************
  inline void IVP_ODE_HIRES::init(double *Y0) 
  {
    Y0[0] = 1.0;
    Y0[1] = 0.0;
    Y0[2] = 0.0;
    Y0[3] = 0.0;
    Y0[4] = 0.0;
    Y0[5] = 0.0;
    Y0[6] = 0.0;
    Y0[7] = 0.0057;
  }
//***************************************************

//***************************************************
//vector system function for the nonstiff term DY=F3(t,Y)
//***************************************************
void IVP_ODE_HIRES::F3 (const double t, const double *Y, double *DY)
{
  const double K = 280.0 * Y[5] * Y[7];
  DY[0] = 0.0;
  DY[1] = 0.0;
  DY[2] = 0.0;
  DY[3] = 0.0;
  DY[4] = 0.0;
  DY[5] = -K;
  DY[6] =  K;
  DY[7] = -K;  
}
//***************************************************

//***************************************************
//vector system function for the stiff term DY=F1(t,Y)
//***************************************************
void IVP_ODE_HIRES::F1 (const double t, const double *Y, double *DY)
{
  DY[0] = -1.71 * Y[0];
  DY[1] =  1.71 * Y[0] - 8.75 * Y[1];
  DY[2] = -10.03 * Y[2];
  DY[3] = -1.12 * Y[3];
  DY[4] = -1.745 * Y[4];
  DY[5] = -0.43 * Y[5];
  DY[6] = -1.81 * Y[6];
  DY[7] =  0.0;
}
//***************************************************

//***************************************************
//vector system function for the 2nd stiff term DY=F2(t,Y)
//***************************************************
void IVP_ODE_HIRES::F2 (const double t, const double *Y, double *DY)
{
  DY[0] = 0.43 * Y[1] + 8.32 * Y[2] + 0.0007;
  DY[1] = 0.0;
  DY[2] = 0.43 * Y[3] + 0.035 * Y[4];
  DY[3] = 8.32 * Y[1] + 1.71 * Y[2];
  DY[4] = 0.43 * Y[5] + 0.43 * Y[6];
  DY[5] = 0.69 * Y[3] + 1.71 * Y[4] + 0.69 * Y[6];
  DY[6] = 0.0;
  DY[7] = 1.81 * Y[6];
}
//***************************************************


//******************************************************
//Compute matrix I-a*J1-b*J2 where J1= Jacobian of the 1st stiff Term
// and J2= Jacobian of the 2nd stiff term in the IVP_ODE_HIRES
//******************************************************       
void IVP_ODE_HIRES::Compute_Matrix_Exact(const double t, double a, double b,
                                                 double *Y, LIS_MATRIX As)
  
  { LIS_INT k,row;

  LIS_INT *ptr=As->ptr;
  LIS_INT *index=As->index;
  LIS_SCALAR *value=As->value;
  
  init_row_insertion(ptr, k, row);

  // row 0
  new_entry(index, value, k, 0, 1.0+a*1.71);
  new_entry(index, value, k, 1, -b*0.43);
  new_entry(index, value, k, 2, -b*8.32);
  next_row(ptr, k, row);
  // row  1 
  new_entry(index, value, k, 0 , -a*1.71 );
  new_entry(index, value, k, 1, 1.0+a*8.75);
  next_row(ptr, k, row);
  // row  2
  new_entry(index, value, k, 2, 1.0+a*10.03);
  new_entry(index, value, k, 3, -b*0.43);
  new_entry(index, value, k, 4, -b*0.035);
  next_row(ptr, k, row);
  // row 3
  new_entry(index, value, k, 1, -b*8.32);
  new_entry(index, value, k, 2, -b*1.71);
  new_entry(index, value, k, 3, 1.0+a*1.12);
  next_row(ptr, k, row);
  // row 4
  new_entry(index, value, k, 4, 1.0+a*1.745);
  new_entry(index, value, k, 5, -b*0.43);
  new_entry(index, value, k, 6, -b*0.43);
  next_row(ptr, k, row);
  // row 5
  new_entry(index, value, k, 3, -b*0.69);
  new_entry(index, value, k, 4, -b*1.71);
  new_entry(index, value, k, 5, 1.0+a*0.43);
  new_entry(index, value, k, 6, -b*0.69);
  next_row(ptr, k, row);
  // row 6
  new_entry(index, value, k, 6, 1.0+a*1.81);
  next_row(ptr, k, row);
  // row 7
  new_entry(index, value, k, 6, -b*1.81);
  new_entry(index, value, k, 7, 1.0);
  next_row(ptr, k, row);  

}
//******************************************************


//******************************************************
// Save the data corresponding to a state vector Y
// on a file called filename
//****************************************************** 
void IVP_ODE_HIRES::save_vector(const std::string& filename, const double *Y)
{
  std::ofstream file(filename);
  if (!file) {
    std::cerr << "Error opening the file: " << filename << std::endl;
    return;
  }
  for (int i = 0; i < neqn; ++i)
  {
    double x_i=(double)(i+1);
    file << x_i<< "     " << Y[i] << std::endl; 
  }
  file.close();
}

#endif

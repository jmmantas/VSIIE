/*
This file contains a C++ subclass of "IVP_ODE" called "IVP_ODE_advdiff1d" which represents  
a  nonlinear one-dimensional Advection-Diffusion-Reaction model.

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



#ifndef ADVDIFF1D_H
#define ADVDIFF1D_H


#include "IVP_ODE.h"

using namespace std;



//*****************************************************************
// Class for the IVP-ODE representing a 1D Advection-Diffusion model 
//*****************************************************************
class IVP_ODE_advdiff1d:public IVP_ODE {

private:
  int nx; // number of grid points at each dimension
  double dtx; // Spatial step
  double dtx2; //dtx*dtx
  //constant scalars representing the strength of advection and diffusion
 
  //const double a=10.0,d=1.0;
  const double a=1.0,d=1.0, c=1;
  const double pi=3.14159265358979;
  double f(const double x, const double t);


public:
  //*****************************************************
  // Constructor of the class IVP_ODE_advdiff1d
  IVP_ODE_advdiff1d(const int nx_points);
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
  // in the IVP_ODE_advdiff1d
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
// Constructor of the class IVP_ODE_advdiff1d
IVP_ODE_advdiff1d::IVP_ODE_advdiff1d(const int nx_points)
//*****************************************************
{ 
  IVP_name="1D_Adv-Diff";
  nx=nx_points;
  // Number of ODEs
  neqn=nx;  
  nnz=3*neqn;
  // Compute Spatial step
  dtx=1.0/nx;
  dtx2=dtx*dtx;
  
}
  //******************************************************
  // Initialize stage vector Y0 with neqn components
  //******************************************************
  inline void IVP_ODE_advdiff1d::init(double *Y0) 
  {
    for (int i=0;i<neqn;i++) { 
      double x_i=(double)(i+1)*dtx;
      Y0[i]=sin (2.0*pi*x_i);
    }
  }
//***************************************************


//***************************************************
//Auxiliary function f(x,t) used to compute the 2nd stiff term  
//***************************************************
double IVP_ODE_advdiff1d::f(const double x, const double t)
{ 
  const double pi2xpt=2.0*pi*x + t;
  const double C=cos(pi2xpt);
  const double S=sin(pi2xpt); 
  return( C + 2*a*pi* S*C + 4*d*pi*pi*S - c*S);
}
//***************************************************


//***************************************************
//vector system function for the nonstiff term DY=F3(t,Y)
//***************************************************
void IVP_ODE_advdiff1d::F3 (const double t, const double *Y, double *DY)
{
  const double dtx_4=4.0*dtx;
  // Compute partially DY in inner points
  for(int i=1;i<nx-1;i++)
    {   
      DY[i]=-a*(Y[i+1]*Y[i+1]-Y[i-1]*Y[i-1])/dtx_4;
    }
  // Compute partially DY in boundary points (i=0 and i=nx-1)
  DY[0]   =-a*(Y[1]*Y[1]-Y[nx-1]*Y[nx-1])/dtx_4;
  DY[nx-1]=-a*(Y[0]*Y[0]-Y[nx-2]*Y[nx-2])/dtx_4;
  


}
//***************************************************

//***************************************************
//vector system function for the stiff term DY=F1(t,Y)
//***************************************************
void IVP_ODE_advdiff1d::F1 (const double t, const double *Y, double *DY)
{
  // Compute partially DY in inner points
  for(int i=1;i<nx-1;i++)
    DY[i]=d*(Y[i+1] -2*Y[i]+ Y[i-1])/dtx2;

  // Compute partially DY in boundary points (i=0 and i=nx-1)
  DY[0]   =d*(Y[1] -2*Y[0]   + Y[nx-1])/dtx2;
  DY[nx-1]=d*(Y[0] -2*Y[nx-1]+ Y[nx-2])/dtx2;

}
//***************************************************

//***************************************************
//vector system function for the 2nd stiff term DY=F2(t,Y)
//***************************************************
void IVP_ODE_advdiff1d::F2 (const double t, const double *Y, double *DY)
{
  for(int i=0;i<nx;i++)  { 
    DY[i]= c*Y[i]+f((i+1)*dtx,t);
  }


}
//***************************************************


//******************************************************
//Compute matrix I-a*J1-bJ2 where J1= Jacobian of the Diffusive Term
// and J2= Jacobian of the reactive term in the IVP_ODE_advdiff1d
//******************************************************       
void IVP_ODE_advdiff1d::Compute_Matrix_Exact(const double t, double a, double b,
                                                 double *Y, LIS_MATRIX As)
  
  { LIS_INT i,k;
    
    LIS_INT *ptr=As->ptr;
    LIS_INT *index=As->index;
    LIS_SCALAR *value=As->value;
    const double t2=d/dtx2;
    const double Aii=1.0+2*a*t2*d-b*c;
    const double other=-a*t2*d;

    ptr[0] = 0;
    index[0] = 0;      value[0] = Aii; //Jn(0,0)
    index[1] = 1;      value[1] = other; //Jn(0,1)
    index[2] = neqn-1; value[2] = other; //Jn(0,nx-1)

    ptr[1]=3;
    k=3;
    for(i=1;i<neqn-1;i++)
    {  
      index[k] = i-1; value[k] = other;//Jn(i,i-1) 
      k++;
      index[k] = i;   value[k] = Aii; //Jn(i,i)
      k++;
      index[k] = i+1; value[k] = other; //Jn(i,i)
      k++;
      ptr[i+1] = k;
    }
    index[k] = 0;      value[k] = other; //Jn(nx-1,0)
    k++;
    index[k] = neqn-2; value[k] = other; //Jn(nx-1,nx-2)
    k++;
    index[k] = neqn-1;      value[k] = Aii; //Jn(nx-1,nx-1)
    k++;
    ptr[neqn] = k;
  
}
//******************************************************



//******************************************************
// Save the data corrersponding to a state vector Y
// on a file called filename
//****************************************************** 
void IVP_ODE_advdiff1d::save_vector(const std::string& filename, const double *Y)
{
  std::ofstream file(filename);
  if (!file) {
    std::cerr << "Error opening the file: " << filename << std::endl;
    return;
  }
  for (int i = 0; i < neqn; ++i)
  {
    double x_i=(double)(i+1)*dtx;
    file << x_i<< "     " << Y[i] << std::endl; 
  }
  file.close();
}

#endif

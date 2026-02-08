/*This file contains a C++ subclass of "IVP_ODE" called "IVP_ODE_combustion" which represents  
several combustion models.

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


#ifndef COMBUSTION_H
#define COMBUSTION_H

#include "IVP_ODE.h"

using namespace std;
#include <algorithm>


//*****************************************************************
// Class for the IVP-ODE  which represents 
// combustion model 
//*****************************************************************
class IVP_ODE_combustion:public IVP_ODE {

private:
  int nx; // number of grid points at each dimension
  double dtx; // Spatial step
  double dtx2; //dtx*dtx
  double pi_dtx_div_L;
  const double a=1.0;
  const double d=1.0; //constant scalars representing the strength of advection and diffusion 
  const double pi=3.14159265358979;
  const double gamma=0.1;
  const double beta=1.0;
  const double GoB=gamma/beta;
 // const double alpha1 = 0.1;
  const double alpha1 = 0.1;
  const double cs=0.6;
  const double alpha2 = alpha1/((1-cs)*(1-cs));
  const double m = 10;
  const double alpha3 = alpha1 / (4 * (pow(m / (m + 1), m) * (1 - m / (m + 1))));
  const double x0=20.5;
  const double xi=10.0;
  const double xf=50.0;
  const double sigma=10;
  const double U0=0.5;//0.99, 0.75
  // Type of combustion reaction term---> 1:FKPP, 2:Ignition, 3:Fisher
  int Type; 
  

  const double L=1;
  //double f_ignition(const double y);
  //double f_Fisher(const double y);
  double reaction(const double y);
  double reaction_derivative(const double y);

  inline double int_pow(const double y, const unsigned m) {
    double result=y;
    for (unsigned i=1;i<m;i++){result*=y;}
    return result; 
  }

public:


    

  //******************************************************
  // Constructor of the class IVP_ODE
  // Argument Type defines the  type of reaction 
  // 1: for FKPP, 2: for Ignition, and 3: for Fisher
  //******************************************************
  IVP_ODE_combustion (const int nx_points, const int Type_combustion);

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

  
};     //*****************************************************


//*****************************************************
// Constructor of the class IVP_ODE_combustion
// Argument Type defines the  type of reaction 
//  1: for FKPP,  2: for Ignition, and   3: for Fisher
IVP_ODE_combustion::IVP_ODE_combustion(const int nx_points, const int Type_combustion)
//*****************************************************
{ 
  switch (  Type_combustion) {
    case 1:
        //Combustion model with FKPP reaction term
        IVP_name="FKPP";
        break;
    case 2:
        //Combustion model with Ignition reaction term
        IVP_name="Ignition";
        break;
    case 3:
        //Combustion model with Fisher reaction term
        IVP_name="Fisher";
        break;
    default:
        // Invalid Reaction Type. Defaulting to FKPP
        IVP_name="FKPP";
        break;
  }
  Type=Type_combustion;
  nx=nx_points;
  // Number of ODEs
  neqn=nx;
  // Number of non-zero elements in the Jacobian matrix
  nnz=3*nx;
  // Compute Spatial step
  dtx=(xf-xi)/nx;
  dtx2=dtx*dtx;
  pi_dtx_div_L=pi*dtx/L;
  
 }
//***************************************************

double IVP_ODE_combustion::reaction(const double y)
{ 
    if (Type == 1) { 
        return (alpha1 * (y) * (1 - y)); 
    }
    else if (Type == 2) {
        return (alpha2 * (1 - y) * std::max(y - cs, 0.0));
    }
    else if (Type == 3) {
        return(alpha3 * int_pow(y, m) * (1 - y));
    }
    else {
        std::cerr << "Invalid Reaction Type." << std::endl;
        return 0.0;
    }
}



 double IVP_ODE_combustion::reaction_derivative(const double y)
{

    if (Type == 1) {
        return (alpha1 *  (1 - 2*y));

    }
    else if (Type == 2) {
         double r= cs > y ? 0 : alpha2 * (cs - 2 * y + 1);
          //  alpha2 * (1 - y) * std::max(y - cs, 0.0));
         return r;
    }
    else if (Type == 3) {
        const double pow1= int_pow(y, m-1);
        return(alpha3 *(m* pow1 - (m+1) * pow1*y));
    }
    else {
        std::cerr << "Invalid Reaction Type." << std::endl;
        return 0.0;

    }

}

//double IVP_ODE_combustion::f_ignition(const double y)
//{
   
 //   return(alpha2*(1-y)* std::max(y-cs,0.0));
//}

//double IVP_ODE_combustion::f_Fisher(const double y)
//{

 //   return(alpha3 * pow(y,m) * (1-y));
//}



  //******************************************************
  // Initialize stage vector Y0 with neqn components
  //******************************************************
  inline void IVP_ODE_combustion::init(double *Y0)
  {
    for (int i=0;i<neqn;i++){ 
    double x_i=xi+(double)i*dtx;
    double dif=x_i-x0;
    Y0[i]=exp(-dif*dif/sigma);
    }
  }
//***************************************************


//***************************************************
// Vector system function for the 1st stiff term DY=F1(t,Y)
//***************************************************
void IVP_ODE_combustion::F1(const double t, const double *Y, double *DY)
{
  // Compute partially DY in inner points
  for(int i=1;i<nx-1;i++)  {   
    DY[i]=(1+U0*sin(i*pi_dtx_div_L))*GoB*(Y[i+1] -2*Y[i]+ Y[i-1])/dtx2;
  }
  // Compute partially DY in boundary points (i=0 and i=nx-1)
  DY[0]   =GoB*(Y[1] -2*Y[0]+ Y[nx-1])/dtx2;
  DY[nx-1]=(1+U0*sin((nx-1)*pi_dtx_div_L))*GoB*(Y[0] -2*Y[nx-1]+ Y[nx-2])/dtx2;
}
//***************************************************



//***************************************************
// Vector system function for the 1st stiff term DY=F2(t,Y)
//***************************************************
void IVP_ODE_combustion::F2(const double t, const double *Y, double *DY)
{
 // Compute partially DY in inner points
  for (int i = 0; i < nx; i++) {
    DY[i] = reaction(Y[i]);
  };
}
//***************************************************


//***************************************************
//vector system function for the nonstiff term DY=F3(t,Y)
//***************************************************
void IVP_ODE_combustion::F3 (const double t, const double *Y, double *DY)
{
  // Compute partially DY in inner points
  for (int i = 1; i < nx - 1; i++) {
      // DY[i]=-a*(1+U0*sin(i*pi_dtx_div_L))*(Y[i+1]-Y[i-1])/(2*dtx)      +     alpha1*(Y[i])*(1-Y[i]);
      DY[i] = -a * (1 + U0 * sin(i * pi_dtx_div_L)) * (Y[i + 1] - Y[i - 1]) / (2*dtx);
  };
  // Compute partially DY in boundary points (i=0 and i=nx-1)
  //DY[0]   =-a*(Y[1]-Y[nx-1])/ (2 * dtx) +  alpha1*(Y[0])*(1-Y[0]);
  DY[0] = -a * (Y[1] - Y[nx - 1]) / (2 * dtx);
  // DY[nx-1]=-a*(1+U0*sin((nx-1)*pi_dtx_div_L))*(Y[0]-Y[nx-2])/ (2 * dtx)
                       //         +  alpha1*(Y[nx-1])*(1-Y[nx-1]);
  DY[nx-1]=-a*(1+U0*sin((nx-1)*pi_dtx_div_L))*(Y[0]-Y[nx-2])/ (2 * dtx);
}
 

//******************************************************
//Compute matrix I-a*J1-b*J2 where J1= Jacobian of F1 and J2 the Jacobian of F2
// in the IVP_ODE_advdiff1d
//******************************************************     
void IVP_ODE_combustion::Compute_Matrix_Exact(const double t, double a, double b,
    double *Y, LIS_MATRIX As)
  { 
    LIS_INT k, row;
    LIS_INT *ptr=As->ptr;
    LIS_INT *index=As->index;
    LIS_SCALAR *value=As->value;


    const double At2=a/dtx2;
    const double mAt2_GoB=-At2*GoB;
    double other=mAt2_GoB;
    double F2_i; //b * derivative_of_F2
    init_row_insertion(ptr, k, row);

    ptr[0] = 0;
    F2_i = b*reaction_derivative(Y[0]);
    new_entry(index, value, k, 0   , 1.0-2*other-F2_i);//Jn(0,0)
    new_entry(index, value, k, 1   , other);//Jn(0,1)
    new_entry(index, value, k, nx-1, other);//Jn(0,nx-1)
    next_row(ptr, k, row);

    for(int i=1;i<nx-1;i++) { 
      other=mAt2_GoB*(1+U0*sin(i*pi_dtx_div_L));
      F2_i = b*reaction_derivative(Y[i]);
      new_entry(index, value, k, i-1   , other);//Jn(i,i-1)
      new_entry(index, value, k, i     , 1.0-2*other-F2_i);//Jn(i,i)
      new_entry(index, value, k, i+1   , other);//Jn(i,i+1)
      next_row(ptr, k, row);
    }
    other=mAt2_GoB*(1+U0*sin((nx-1)*pi_dtx_div_L));
    F2_i = b*reaction_derivative(Y[nx - 1]);
    new_entry(index, value, k, 0   , other);//Jn(nx-1,0)
    new_entry(index, value, k, nx-2 , other);//Jn(nx-1,nx-2)
    new_entry(index, value, k, nx-1 , 1.0-2*other-F2_i);//Jn(nx-1,nx-1)
    next_row(ptr, k, row);
    
  }


//******************************************************
// Save the data corrersponding to a state vector Y
// on a file called filename
//****************************************************** 
void IVP_ODE_combustion::save_vector(const std::string& filename, const double *Y)
{
  std::ofstream file(filename);
  if (!file) {
    std::cerr << "Error opening the file: " << filename << std::endl;
    return;
  }
  for (int i=0;i<neqn;i++){
    double x_i=xi+(double)i*dtx;
    file << x_i<< "     " << Y[i] << std::endl; 
  }
  file.close();
}


//***************************************************





#endif

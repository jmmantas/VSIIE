/*
This file contains the C++ abstract class "IVP_ODE" which defines the virtual functions to be implemented
in the subclasses (representing particular IVPs for ODEs given with three components). It includes public
functions which are common to all the IVPs as: the evaluation of the nonstiff component F3(t, y),  
the evaluation of the two stiff components F1(t, y) and F2(t,y), to get the number of ODEs, 
to initialize the sparse matrices, etc.).

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

#ifndef Included_IVP_ODE_H

#define Included_IVP_ODE_H


#include <iostream>
#include <fstream>
#include <math.h>
#include <cstring>
#include <chrono>
#include "lis.h"
#include "lis_config.h"
#include <vector>
#include <cblas.h>
#include <sstream>
#include <iomanip>

using namespace std;

//*****************************************************************
// Class for the IVP-ODE representing a model 
//*****************************************************************
class IVP_ODE {

protected:
  int neqn; // number of equations
  string IVP_name; // Name of the IVP model
  int nnz; //number of non-zero elements in the Jacobian for the stiff term G+H
 

  //***************************************************
  // Introduce a new element "entry" in a CSR sparse matrix, 
  // in the k-th position of the index vector on the column "col"
  //*************************************************** 
  inline void new_entry(LIS_INT* index, LIS_SCALAR* value,  
                      LIS_INT & k, const int col, const double entry){
    index[k] = col;
    value[k] = entry;
    k++;
  }  
  //***************************************************

  //***************************************************
  // Insert a new row in a CSR sparse matrix with pointer vector ptr
  // this new row start in position k of the index vector
  //*************************************************** 
  inline void next_row(LIS_INT* ptr, const LIS_INT k, int & row){
    row++;
    ptr[row] = k;
  }  
  //***************************************************

  //***************************************************
  // Init row 0 in CSR sparse matrix with pointer vector ptr
  // this new row start in position k=0 of the index vector
  //*************************************************** 
  inline void init_row_insertion(LIS_INT* ptr, LIS_INT & k, int & row){
    row = 0;
    ptr[row] = 0;
    k = 0;
  }  
  //***************************************************



//***************************************************
// PUBLIC FUNCTIONS
public:
//***************************************************
  
  //******************************************************
  // Get the number of ODEs
  //******************************************************
  inline int get_num_ODEs() {return neqn;};
  //******************************************************

  //***************************************************
  // Get the name of the IVP model
  inline string get_name() {return IVP_name;}
  //***************************************************                                              


  //***************************************************
  // Get the number of nonzero elements in Jacobian matrix
  inline int get_nnz() {return nnz;}
  //***************************************************                                              


  //***************************************************
  // Init a CSR Sparse neqn x neqn LIS matrix
  //******************************************************
  inline void init_CSR_LIS_matrix(const int splitting_type, LIS_MATRIX *A) 
  {
    LIS_INT * ptr_A, * index_A;
    LIS_SCALAR* value_A;
    lis_matrix_create(0, A);
    lis_matrix_set_size(*A, 0, neqn);
    lis_matrix_malloc_csr(neqn, nnz, &ptr_A, &index_A, &value_A);
    lis_matrix_set_csr(nnz,ptr_A,index_A,value_A,*A);
    lis_matrix_assemble(*A);
  }
  //***************************************************                                              

 
  //******************************************************
  // Initialize the stage vector Y0 with neqn components
  //******************************************************
  virtual void init(double *Y0)=0; 


  //*************************************************************
  // ODE vector system function feval:
  // It denotes the adddition of the three components F1, F2 and F3.  
  // DY=feval(t,Y)=F1(t,Y)+F2(t,Y)+F3(t,Y)
  //*************************************************************
  inline void feval (const double t, const double *Y, double *DY)
  { double DY_tmp[neqn];
    F1 (t, Y, DY);
    F2 (t, Y, DY_tmp);
    cblas_daxpy(neqn,1.0,DY_tmp,1,DY,1);
    F3 (t, Y, DY_tmp);
    cblas_daxpy(neqn,1.0,DY_tmp,1,DY,1);

  }

  //******************************************************
  //vector system function for the 1st stiff term DY=F1(t,Y)
  //******************************************************
  virtual void F1 (const double t, const double *Y, double *DY)=0;
  //******************************************************
  //******************************************************
  //vector system function for the 2nd stiff term DY=F2(t,Y)
  //******************************************************
  virtual void F2 (const double t, const double *Y, double *DY)=0;
  //******************************************************
  //vector system function for the non-stiff term DY=F3(t,Y)
  //******************************************************
  virtual void F3 (const double t, const double *Y, double *DY)=0;
  //******************************************************



  //******************************************************
  // Compute matrix I-a*J1-b*J2 where J1= Jacobian of the 1st stiff Term 
  // J2= Jacobian of the 2nd stiff Term in the IVP_ODE
  //****************************************************** 
  virtual void Compute_Matrix_Exact(const double t, double a, double b,
                                    double *Y, LIS_MATRIX As)=0;
  //******************************************************
  

  ///******************************************************
  // Compute matrix  I-a*J1-b*J2, where J1 and J2 are the 
  // Forward Difference Approximation of Jacobianw of F1 and F2
  //******************************************************       
  void Compute_Matrix_FD(const double t, double a, double b, double *Y, LIS_MATRIX As)
  { 
    double * DY0_F1=new double[neqn];
    double * DY0_F2=new double[neqn];
    double * DY_F1=new double[neqn];
    double * DY_F2=new double[neqn];
    LIS_INT *ptr=As->ptr;
    LIS_INT *index=As->index;
    LIS_SCALAR *value=As->value;
    LIS_INT  row=0, k;
    init_row_insertion(ptr, k, row);
    double Aij;
    const double uround=1.0e-15;
    double ysafe[neqn];
    F1(t, Y, DY0_F1);
    F2(t, Y, DY0_F2);
    cblas_dcopy(neqn, Y,1, ysafe, 1);
    for(int i=0;i<neqn;i++){
      for(int j=0;j<neqn;j++){
        const double delta=sqrt(uround*max(1.e-5,abs(ysafe[j])));
        Y[j]=ysafe[j]+delta;
        F1(t, Y, DY_F1);
        F2(t, Y, DY_F2);
        Y[j]=ysafe[j];
        Aij=-a*(DY_F1[i]-DY0_F1[i])/delta-b*(DY_F2[i]-DY0_F2[i])/delta;
        if (i==j) Aij+=1.0;
        if (fabs(Aij)>uround){  
          new_entry(index, value, k, j, Aij);
          //nonzeros++;
        } 
      } 
      next_row(ptr, k, row);
    }  
    //cout<<"nonzeros="<<nonzeros<<endl;  
    delete[] DY0_F1;
    delete[] DY0_F2;
    delete[] DY_F1;
    delete[] DY_F2;
  }

//******************************************************



//******************************************************
// Save the data corresponding to a state vector Y
// on a file called filename 
//****************************************************** 
virtual void save_vector(const std::string& filename, const double *Y)=0; 
//******************************************************



};


#endif

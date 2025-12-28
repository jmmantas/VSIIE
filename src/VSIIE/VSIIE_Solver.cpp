/*
This is the implementation of the C++ class "VSIIE_Solver" which implements the VSIIE solvers with support for
orders 1 to 4. The solvers apply adaptive time-stepping techniques to optimize the computational
process where the accuracy is ensured by error checking. 

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

#include "IVP_ODE.h"
#include "VSIIE_Solver.h"

using namespace std;
#include <assert.h>

//*******************************************************************
//*******************************************************************
// IMPLEMENTATION OF THE PRIVATE FUNCTIONS
//*******************************************************************


//******************************************************
// Update coefficients of fixed stepsize-order scheme
//******************************************************
void VSIIE_Solver::Fixed_update_coefs(const unsigned order) {
  double *a=new double[order+1], *b1=new double[order+1], 
        *b2=new double[order+1], *b3=new double[order+1];
  switch (order) 
  { 
    case 1:
    {
      a[0]=-1.0;  a[1]=1.0;
      b1[0]=0.5;  b1[1]=0.5;
      b2[0]=-0.5;  b2[1]=1.5;
      b3[0]=1.0;
    }
    break;
    case 2:
    {
      // Definition of Parameters of the method IIE-2
      a[0]= -1.0;  a[1]= 0.0;   a[2]= 1.0;
      b1[0]=1.0;   b1[1]=0.0;   b1[2]= 1.0;
      b2[0]= 2.0;  b2[1]= -2.0; b2[2]= 2.0; 
      b3[0]= 0.0;  b3[1]= 2.0;
    }
    break;
    case 3:
    {  
      // Definition of Parameters of the method IIE-3
      a[0] =-2.0/11; a[1] =9.0/11;     a[2] = -18.0/11; a[3] = 1.0;
      b1[0]= 0.0;    b1[1] = 0.0;      b1[2] = 0.0;b1[3] = 6.0/11;
      b2[0] = 1.0/22;b2[1] = -3.0/22;  b2[2] = 3.0/22; b2[3] = 1.0/2;
      b3[0] =6.0/11; b3[1] = -18.0/11; b3[2] = 18.0/11;
    }  
    break;
    case 4:
    { 
      // Definition of Parameters of the method IIE-4
      a[0] = 3.0/25; a[1] = -16.0/25; a[2] = 36.0/25; a[3] = -48.0/25; a[4] = 1.0 ;
      b1[0]= b1[1] = b1[2] = b1[3]= 0.0;     b1[4]= 12.0/25;            
      //b2[0] = -6.0/25; b2[1] = 24.0/25; b2[2] = -36.0/25; b2[3] = 24.0/25; b2[4]=6.0/25;
      b2[0] = -24.0/25; b2[1] = 96.0/25; b2[2] = -144.0/25; b2[3] = 96.0/25; b2[4]=-12.0/25;
      b3[0] = -12.0/25; b3[1] = 48.0/25; b3[2] = -72.0/25;  b3[3] = 48.0/25;

    }
    break;
  }


 // Divide all by a[order] to obtain coef_Y, coef_F1, coef_F2 and coef_F3
  const double a_order=a[order];
  for(unsigned i=0;i<=order;i++){ 
    coef_F1[i]= b1[i]/a_order;
    coef_F2[i]= b2[i]/a_order;
    if (i<order)
      {
        coef_Y[i]= -a[i]/a_order; 
        coef_F3[i]= b3[i]/a_order;
      }  
  }
    coef_idx_F1=coef_F1[order];
    coef_idx_F2=coef_F2[order];

    // Free dynamically allocated memory
    delete[] a;
    delete[] b1;
    delete[] b2;
    delete[] b3;
}







//******************************************************
// Update coefficients of VSIIE-order scheme
//******************************************************
void VSIIE_Solver::Update_coefs(const unsigned order, const double * h_vector) {
  double w[idx];
  double *a=new double[order+1], *b1=new double[order+1], 
        *b2=new double[order+1], *b3=new double[order+1];
  for (int i=0;i<=idx;i++){
    w[i]=h_vector[i+1]/h_vector[i];
  } 
  switch (order) 
  { 


    case 1:
    {
      a[0]=-1.0;  a[1]=1.0;
      b1[0]=0.0; b1[1]=1.0;
      b2[0]=0.5; b2[1]=0.5;
      b3[0]=1.0;
    }
    break;
    case 2:
    {
      // Definition of Parameters of the method VSIIE-2  
      /*const double alpha=0.5,  beta=0.5,  gamma=0.5;
      const double w_n1=w[0];
      a[0]= ((2.0 * gamma - 1) * w_n1 * w_n1) / (1.0 + w_n1);
      a[1]= (1.0 - 2 * gamma) * w_n1 - 1;
      a[2]= (1.0 + 2 * gamma * w_n1) / (1 + w_n1);
      b1[0]=beta / 2.0;
      b1[1]=1.0 - gamma - ( (1.0 + 1.0 / w_n1) * (beta / 2.0) );
      b1[2]= gamma + (beta / (2 * w_n1));
      b2[0]= (alpha - gamma) * w_n1;
      b2[1]= (gamma - alpha) * (w_n1 + 1.0)-alpha + 1.0;
      b2[2]= alpha; 
      b3[0]= -gamma * w_n1;
      b3[1]= 1.0 + gamma * w_n1;*/

      const double alpha=0.5;
      const double w_n1=w[0];
      a[0]= w_n1 * w_n1 / (1.0 + w_n1);
      a[1]= -(1+w_n1);
      a[2]= (1.0 + 2 * w_n1) / (1 + w_n1);
      b1[0]=0.0;
      b1[1]=0.0;
      b1[2]= 1.0;
      b2[0]= -0.5 * w_n1;
      b2[1]=  0.5 * (1.0 + w_n1);
      b2[2]= alpha; 
      b3[0]= -w_n1;
      b3[1]= 1.0 + w_n1;

    }
    break;
    case 3:
    {  
      // Definition of Parameters of the method VSIIE-3
      /*
      const double alpha=0.0, beta=1.0,  gamma=1.0, theta=0.0;
      const double w_n1=w[0], w_n2=w[1];
      const double G1=1.0 + w_n1, G2=1.0 + w_n2, A1= 1.0 + w_n1*G2;
      
      
      const double K1= theta/3.0 + gamma*gamma -beta;
      const double K2= theta/2.0 + gamma       -beta;
      const double K3= theta     + 1.0         -beta;

      const double w_n1_sqr=w_n1*w_n1, w_n1_cube=w_n1*w_n1_sqr, w_n2_sqr=w_n2*w_n2;
      const double G1m6=6.0 * G1;
      const double w_n2m2p3=3.0 + 2.0 * w_n2;
      const double gamma_wn_2p1=gamma * w_n2 + 1.0; 

      // Calculation of the specific parameter values as a function of w_n1 and w_n2
      a[0] = -( w_n1_cube * w_n2_sqr * (3.0 * gamma*gamma* w_n2 + 2.0 * gamma * (1 - w_n2) - 1.0) ) / (G1 * A1);
      a[1] = -( w_n2_sqr * (gamma * w_n1 * w_n2 * (2.0 - 3.0 * gamma) + (1.0 - 2.0 * gamma) * G1) ) / G2;
      a[2] = -(  w_n1 * gamma_wn_2p1 * (1.0 + w_n2 * (3.0 * gamma - 2.0)) + (w_n2 * (2 * gamma - 1.0) + 1.0) )/G1 - theta;
      a[3] = (1.0 + 2.0 * gamma * w_n2 + w_n1 * gamma_wn_2p1 * (3.0 * gamma * w_n2 + 1.0) ) / (G2 * A1)   + theta;

      b1[0]= (theta * w_n1_sqr *w_n2* w_n2m2p3) / G1m6 - alpha;
      b1[1] = (alpha * G1 * A1 - w_n1_sqr * w_n2_sqr * gamma * (1.0 - gamma)) / (w_n1_sqr * G2)
            - (1.0 / 6.0) * w_n2 * theta * ( 3.0 + w_n1 * w_n2m2p3 );
      b1[2] = (w_n1_sqr * w_n2 * (1.0 - gamma) * gamma_wn_2p1 - alpha * A1) / (w_n1_sqr * w_n2)
            + theta * ( 1.0 + 0.5 * w_n2 + (w_n1 * w_n2 *w_n2m2p3) / G1m6 );
      b1[3] = ( w_n1_sqr * w_n2 * gamma * gamma_wn_2p1 + alpha * G1 ) / (w_n1_sqr * w_n2 * G2);


      b2[0] = (K1 * w_n2 + K2) * w_n1_sqr * w_n2 / G1;
      b2[1] = -(K1 * w_n1 * w_n2 + K2 * G1) * w_n2;
      b2[2] = (K1 * w_n1 * w_n2_sqr + K2 * w_n2 * (2 * w_n1 + 1.0)) / G1 + K3;
      b2[3] = beta;
    
      b3[0] =  gamma * w_n1_sqr * w_n2 * gamma_wn_2p1 / G1   +  (theta * w_n1_sqr * w_n2 * w_n2m2p3) / G1m6;
      b3[1] = -gamma * w_n2 * (1.0 + w_n1 * gamma_wn_2p1)- w_n2 * theta * (3.0 + w_n1 * w_n2m2p3)/6.0;
      b3[2] = (gamma_wn_2p1 * (1.0 + w_n1 * gamma_wn_2p1)) / G1 + 
                         theta * (1.0 + 0.5 * w_n2 + (w_n1 * w_n2 * w_n2m2p3) / G1m6);
     */
    
      const double alpha=0.5;
      const double w_n1=w[0], w_n2=w[1];
      const double G1=1.0 + w_n1, G2=1.0 + w_n2, A1= 1.0 + w_n1*G2;
      const double w_n1_sqr=w_n1*w_n1, w_n1_cube=w_n1*w_n1_sqr, w_n2_sqr=w_n2*w_n2;
     
        // Calculation of the specific parameter values as a function of w_n1 and w_n2
      a[0] = -( w_n1_cube * w_n2_sqr * G2) / (G1 * A1);
      a[1] = ( w_n2_sqr * A1) / G2;
      a[2] = -G2*A1/G1;
      a[3] = 1.0 + w_n2/G2 + (w_n1 * w_n2)/A1;

      b1[0]=  0.0;
      b1[1] = 0.0;
      b1[2] = 0.0;
      b1[3] = 1.0;


      b2[0] = 0.5*w_n1_sqr * w_n2*G2 / G1;
      b2[1] = -0.5* w_n2 * A1;
      b2[2] = 0.5*G2*A1 / G1;
      b2[3] = alpha;
    
      b3[0] =  w_n1_sqr * w_n2 * G2 / G1;
      b3[1] = -w_n2 * A1;
      b3[2] = G2*A1 / G1;


    }  
    break;
    case 4:
    { 
      const double alpha=0.5; // the new parameter
      const double w_n1=w[0], w_n2=w[1], w_n3=w[2];
      const double G1=1.0 + w_n1, G2=1.0 + w_n2, G3=1.0 + w_n3;
      const double A1= 1.0 + w_n1*G2, A2= 1.0 + w_n2*G3, A3=1.0 + w_n1*A2;
      const double w_n1_sqr=w_n1*w_n1;
      const double w_n1_p4=w_n1_sqr*w_n1_sqr, w_n1_cube=w_n1*w_n1_sqr;
      const double w_n2_sqr=w_n2*w_n2;
      const double w_n2_cube=w_n2*w_n2_sqr;
      const double w_n3_sqr=w_n3*w_n3;
      const double alpha_1=alpha-1.0, H = alpha_1*A2 *(G3*w_n1 * w_n2 + G1);
    
      a[0] = (G3 * A2 * w_n1_p4 * w_n2_cube * w_n3_sqr) / (G1 * A1 * A3);
      a[1] = -(w_n2_cube * w_n3_sqr * G3 * A3) / (G2 * A2);
      a[2] = w_n3 * (w_n3 / G3 + (w_n2 * w_n3 * (A3 + w_n1)) / G1);
      a[3] = -G3 - (w_n3 * w_n2 * G3 * (A1 + w_n1 * A2)) / (A1 * G2);
      a[4] = 1.0 + w_n3/G3 + w_n2*w_n3/A2 + w_n1*w_n2*w_n3/A3;

      b1[0]= b1[1]=b1[2]=b1[3]= 0.0;     b1[4]= 1.0;
    
      b2[0] = A2*w_n3*w_n1_cube*G3*w_n2_sqr*alpha_1/(A1*G1);
      b2[1] = -(w_n3*alpha_1*G3*w_n2_sqr*A3)/G2;
      b2[2] = w_n3*H/G1;
      b2[3] = -G3*H/(G2*A1);
      b2[4]=alpha;

      b3[0] = -(w_n1_cube * w_n2_sqr*w_n3*G3*A2)/(G1*A1);
      b3[1] = (w_n2_sqr*w_n3*G3*A3)/G2;
      b3[2] = -(A2*A3*w_n3)/G1;
      b3[3] = (w_n2*G3*G3*(A3 + w_n1) + G1*G3)/(G2 * A1);

    }
    break;
  }

 // Divide all by a[order] to obtain coef_Y, coef_F1, coef_F2 and coef_F3
  const double a_order=a[order];
  for(unsigned i=0;i<=order;i++){ 
    coef_F1[i]= b1[i]/a_order;
    coef_F2[i]= b2[i]/a_order;
    if (i<order)
      {
        coef_Y[i]= -a[i]/a_order; 
        coef_F3[i]= b3[i]/a_order;
      }  
  }

  coef_idx_F1=coef_F1[order];
  coef_idx_F2=coef_F2[order];
  
  // Free dynamically allocated memory
  delete[] a; delete[] b1; delete[] b2; delete[] b3;

}




//******************************************************
// Update coefficients for the  SBDF-4 scheme
//******************************************************
void VSIIE_Solver::Update_coefs_SBDF4(const double * h_vector) {
  const int order=4;
  const int idx4= order -1;
  double w[idx4];
  double *a=new double[order+1], *b1=new double[order+1], 
         *b2=new double[order+1], *b3=new double[order+1];
  for (int i=0;i<=idx4;i++){
    w[i]=h_vector[i+1]/h_vector[i];
  }  
  const double alpha=0.5; // the new parameter
  const double w_n1=w[0], w_n2=w[1], w_n3=w[2];
  const double G1=1.0 + w_n1, G2=1.0 + w_n2, G3=1.0 + w_n3;
  const double A1= 1.0 + w_n1*G2, A2= 1.0 + w_n2*G3, A3=1.0 + w_n1*A2;
  const double w_n1_sqr=w_n1*w_n1;
  const double w_n1_p4=w_n1_sqr*w_n1_sqr, w_n1_cube=w_n1*w_n1_sqr;
  const double w_n2_sqr=w_n2*w_n2;
  const double w_n2_cube=w_n2*w_n2_sqr;
  const double w_n3_sqr=w_n3*w_n3;
  const double alpha_1=alpha-1.0, H = alpha_1*A2 *(G3*w_n1 * w_n2 + G1);
    
  a[0] = (G3 * A2 * w_n1_p4 * w_n2_cube * w_n3_sqr) / (G1 * A1 * A3);
  a[1] = -(w_n2_cube * w_n3_sqr * G3 * A3) / (G2 * A2);
  a[2] = w_n3 * (w_n3 / G3 + (w_n2 * w_n3 * (A3 + w_n1)) / G1);
  a[3] = -G3 - (w_n3 * w_n2 * G3 * (A1 + w_n1 * A2)) / (A1 * G2);
  a[4] = 1.0 + w_n3/G3 + w_n2*w_n3/A2 + w_n1*w_n2*w_n3/A3;

  /* b2[0] = A2*w_n3*w_n1_cube*G3*w_n2_sqr*alpha_1/(A1*G1);
  b2[1] = -(w_n3*alpha_1*G3*w_n2_sqr*A3)/G2;
  b2[2] = w_n3*H/G1;
  b2[3] = -G3*H/(G2*A1);
  b2[4]=alpha; */


  for (unsigned i=0;i<order;i++){
    b1[i]= b2[i]=0.0;
  }  
  b1[4]=b2[4]=1.0;

  b3[0] = -(w_n1_cube * w_n2_sqr*w_n3*G3*A2)/(G1*A1);
  b3[1] = (w_n2_sqr*w_n3*G3*A3)/G2;
  b3[2] = -(A2*A3*w_n3)/G1;
  b3[3] = (w_n2*G3*G3*(A3 + w_n1) + G1*G3)/(G2 * A1);


  // Divide all by a[order] to obtain coef_Y, coef_F1, coef_F2 and coef_F3
  const double a_order=a[order];
  for(unsigned i=0;i<=order;i++){ 
    coef_F1[i]= b1[i]/a_order;
    coef_F2[i]= b2[i]/a_order;
    if (i<order)
      {
        coef_Y[i]= -a[i]/a_order; 
        coef_F3[i]= b3[i]/a_order;
      }  
  }
  coef_idx_F1=coef_F1[order];
  coef_idx_F2=coef_F2[order];
  
  // Free dynamically allocated memory
  delete[] a; delete[] b1; delete[] b2; delete[] b3;

}




//**************************************************************
// Print values of sparse matrix A
void VSIIE_Solver::print_matrix_csr(LIS_MATRIX A)
//**************************************************************
{
  LIS_INT	n;
	LIS_INT	i,j,jj;
	n   = A->n;
  LIS_SCALAR value1;
  for(i=0;i<n;i++){
		for(j=A->ptr[i];j<A->ptr[i+1];j++){
			jj = A->index[j]; 		value1 = A->value[j];
      cout<<"ROW  "<<i<<" ....  COL "<<jj<<" Value in MATRIX is  "<<value1<<endl;
    } 	
  }
} 





//**************************************************************
// Check important errors between sparse matrices A1 and A2
void VSIIE_Solver::compare_matrix_csr(LIS_MATRIX A1, LIS_MATRIX A2)
//**************************************************************
{
  LIS_INT	n,nnz;
	LIS_INT	i,j,jj1,jj2;
	n   = A1->n;
	nnz = A1->nnz;
    LIS_SCALAR value1, value2;
    int errors=0;
    const double tol=1.0e-2;
    double max_value=A1->value[0];
    max_value=0;
    for(i=1;i<nnz;i++) 
       max_value=max(max_value, fabs(A1->value[i]));

    for(i=0;i<n;i++)
	{
		for(j=A1->ptr[i];j<A1->ptr[i+1];j++)
		{
			jj1 = A1->index[j];
            jj2 = A2->index[j];
			value1 = A1->value[j];
            value2 = A2->value[j];
            double diff=fabs((value2-value1)/max_value);
            if (jj1!=jj2) {
                cout<<"ROW  "<<i<<" ....  COL "<<jj1<<" in A1 is COL "<<jj2<<" in A2"<<endl;
                errors++;
            }  
            else if (diff>tol) {
                cout<<"ROW  "<<i<<" ....  COL "<<jj1<<" Value in A1 is  "<<value1
                                    <<"    but is   "<<value2<<"  in A2"<<endl;
                errors++;
            } 
          
		}
	}
    if (errors>0) {
        cout<<endl<<endl<<"*****************  "<<errors<< 
                       " ERRORS  ***********************" <<endl; 
       // exit(-1);
    } 
}


  //**************************************************************************
  // Compute R=R+a*F1(t,Y) 
  //**************************************************************************
  void  VSIIE_Solver::Add_aF1Y(const double t, const double a, double* Y) {   
    double* DY = new double[neqn];
    IVP->F1(t, Y, DY);
    cblas_daxpy(neqn, a, DY, 1, R, 1);
    delete [] DY;
  }
  //**************************************************************************
  // Compute R=R+a*F2(t,Y) 
  //**************************************************************************
  void  VSIIE_Solver::Add_aF2Y(const double t, const double a, double* Y) {   
    double* DY = new double[neqn];
    IVP->F2(t, Y, DY);
    cblas_daxpy(neqn, a, DY, 1, R, 1);
    delete [] DY;
  }
//**************************************************************************
  // Compute R=R+a*F3(t,Y) 
  //**************************************************************************
  void  VSIIE_Solver::Add_aF3Y(const double t, const double a, double* Y) {   
    double* DY = new double[neqn];
    IVP->F3(t, Y, DY);
    cblas_daxpy(neqn, a, DY, 1, R, 1);
    delete [] DY;
  }




//***************************************************
  // Init CSR Sparse neqn x neqn LIS matrix A
  void VSIIE_Solver::init_CSR_LIS_matrix(LIS_MATRIX *A) {
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
// Compute the values for every CSR Sparse LIS matrix
//******************************************************

void VSIIE_Solver::compute_matrix(const double t, const double h, double ** Y) {
      //Compute matrix As=Id-a*J1-b*J2
      //init_CSR_LIS_matrix(&As2);
      //IVP->Compute_Matrix_FD(t, h_coef_idx_F1, h_coef_idx_F2, Y[idx], As2);
      //print_matrix_csr( As);
      IVP->Compute_Matrix_Exact(t, h_coef_idx_F1, h_coef_idx_F2, Y[idx], As);
      //compare_matrix_csr(As, As2);
  }
//***************************************************


///**************************************************************************
// Compute the variable term (in a time step) of Right Hand Side R at each Newton Iteration
//**************************************************************************
void  VSIIE_Solver::compute_RHS(const double t, const double h, double* Y1, double * RHS){ 
  double* DY = new double[neqn];
  //RHS=Y1-R0 
  cblas_dcopy(neqn, Y1, 1, RHS, 1);
  cblas_daxpy(neqn, -1.0, R, 1, RHS, 1);

  
  //RHS=RHS-h*coef_idx_F1*F1(t+h,Y1)
  IVP->F1(t + h, Y1, DY);
  cblas_daxpy(neqn, -h_coef_idx_F1, DY, 1, RHS, 1);

  //RHS=RHS-h*coef_idx_F2*F2(t+h,Y1)
  IVP->F2(t + h, Y1, DY);
  cblas_daxpy(neqn, -h_coef_idx_F2, DY, 1, RHS, 1);
  
  delete[] DY;
}


//**************************************************************************

//**************************************************************************
// Compute the constant term (in a time step) of Right Hand Side (R) at each Newton Iteration
//**************************************************************************
void  VSIIE_Solver::compute_RHS0(const double t, const double * h_vector, double** Y)
{
  const double h=h_vector[idx];
  double ti[4];
  ti[idx]=t;
  for (int i = idx-1; i >= 0; i--){
    ti[i]=ti[i+1]-h_vector[i];
  }

  //R=Y[0] 
  cblas_dcopy(neqn, Y[0], 1, R, 1);
  // R=a0*Y[0]
  cblas_dscal(neqn, coef_Y[0],      R, 1);

  for (int i = 1; i <= idx; i++) {
    cblas_daxpy(neqn, coef_Y[i],      Y[i], 1, R, 1);
  }

  if (order<4) 
    for (int i = 0; i <= idx; i++) {
      Add_aF1Y(ti[i],      coef_F1[i]*h,    Y[i]);
    }

  for (int i = 0; i <= idx; i++) {
    Add_aF2Y(ti[i],      coef_F2[i]*h,    Y[i]);
  }

  for (int i = 0; i <=idx; i++) {
    Add_aF3Y(ti[i],      coef_F3[i]*h,    Y[i]);
  }
}  
 













  //******************************************************
  // Compute LTE estimate from:
  // the coarse solution Yc_sol and fine solution Yf_sol
  // and the step size coarse vector h_vector
  //******************************************************
  double VSIIE_Solver::compute_LTE(const double * h_vector, 
                      double * Yc_sol, 
                      double * Yf_sol, double * LTE)
  //******************************************************
  { 
    double factor;
    const double h_np1=h_vector[idx]; // h_{n+1}
    const double h_n  =h_vector[idx-1]; //h_n
    const double h_n_sq  =h_n*h_n; 
    const double h_np1_sq=h_np1*h_np1;
 
    switch(order) {
      case 1:
        factor=2.0;  break;

      case 2:
      { 
        factor=h_np1_sq*(h_np1 + h_n); //h_{n+1}^2 *(h_{n+1}+h_n)
        const double denom=(5.0*h_np1+7.0*h_n)*h_np1_sq/8.0; 
        factor=factor/denom;
        
        
        //const double denom = (7.0*h_n+5.0*h_np1);
        //factor=8.0*(h_n+h_np1)/denom;

        //cout <<"FACTOR2="<<factor<<endl;


        break;
      }

      case 3:
      { const double h_nm1=h_vector[idx-2];
        const double denom=14.0*h_n_sq +16.0*h_n*h_nm1 
                            + 27.0* h_n*h_np1 + 16.0*h_np1*h_nm1 + 11.0*h_np1_sq;
        factor=16.0*(h_n+h_np1)*(h_np1 + h_n + h_nm1)/denom;


      //const double h_nm1=h_vector[idx-2];
      //const double denom=3.0* h_nm1*h_n + 2.0* h_nm1*h_np1 + h_n_sq
      //                 -2.0*h_n*h_np1-4.0*h_np1_sq;
      //factor=4.0*(h_n + h_np1)*(h_nm1 + h_n + h_np1)/denom;               

      break;
      }  

      case 4:  
      {  
        const double h_nm1  =h_vector[idx-2]; //h_{n-1}
        const double h_nm2  =h_vector[idx-3]; //h_{n-2}
        const double h_nm1_sq=h_nm1*h_nm1;

        const double C_1 = h_n_sq*(28.0*h_n+32.0*h_nm2+62.0*h_nm1+84.0*h_np1) 
                           + 32.0*h_n *(h_nm2*h_nm1+2.0*h_nm2*h_np1+h_nm1_sq);
                           
        const double C_2 = 125.0 * h_n * h_nm1 * h_np1 + 79.0 *h_n*h_np1_sq + 32.0*h_nm2*h_nm1*h_np1
            + 32.0 * h_nm2 * h_np1_sq + 32.0 * h_nm1_sq * h_np1 + 63.0 * h_nm1 * h_np1_sq
            + 23.0 * h_np1_sq * h_np1;
        //double C_1 = h_n_sq*(28.0*h_n+32.0*h_nm2+62.0*h_nm1+84.0*h_np1) 
        //  + 32.0*h_n *(h_nm2*h_nm1+2.0*h_nm2*h_np1+h_nm1_sq);
        //const double C_2 = 125.0 * h_n * h_nm1 * h_np1 + 79.0 *h_n*h_np1_sq + 32.0*h_nm2*h_nm1*h_np1
          //  + 32.0 * h_nm2 * h_np1_sq + 32.0 * h_nm1_sq * h_np1 + 63.0 * h_nm1 * h_np1_sq
         //   + 23.0 * h_np1_sq * h_np1;
        const double denom = C_1 + C_2;

        double term = (h_np1 + h_n);
        term= term*(term + h_nm1) * (term + h_nm1 + h_nm2);
        factor= (32.0 * term)/denom; 
       
        

        //const double C1=32.0*h_nm1_sq + h_nm1*(32.0*h_nm2 + 62*h_n+ 63.0*h_np1);
        //const double C2=28.0*h_n_sq + h_n*(32.0*h_nm2+56.0*h_np1) + 32.0*h_nm2*h_np1 + 23.0*h_np1_sq;
        //const double denom= C1+C2;
        //factor=32.0*(h_nm2+ h_nm1 + h_n + h_np1)*(h_nm1+h_n + h_np1);
        //factor=factor/denom;


        
        break;

      }
      default: 
        cout <<"order should belong to {1,2,3,4} "<<endl; exit(-1);
    }  
    // Compute difference LTE=factor*(Yc-Yf)
    cblas_dcopy(neqn, Yc_sol, 1, LTE, 1);
    cblas_daxpy(neqn, -1.0, Yf_sol, 1, LTE, 1);
    cblas_dscal(neqn, factor,      LTE, 1);
    double error = cblas_dnrm2(neqn, LTE, 1);// / pow(neqn, 0.5);
    return(error);
  }
  //******************************************************




  //**************************************************************************
  // Update intermediate vectors Yf and Y using the current vector solution Y1 
  // for the adaptive time stepper before the next integration step
  //**************************************************************************
  void VSIIE_Solver::Update_adaptive_intermediate_vectors(double *h_vector_half, 
                                                         double *h_vector_half2, 
                                                         double *h_vector, 
                                                         double ** Yf, double ** Y, 
                                                         double *Y1)
  {  
    //Update h_vectors for the coarse and fine approximation  
    for (int i=0;i<idx;i++) {
      h_vector_half[i]=h_vector_half2[i+1];
      h_vector[i]=h_vector[i+1];
    }

    // Update coarse grain approximation vectors Y
    Update_intermediate_vectors(Y, Y1);

    // Update fine grain approximation vectors Yf from coarse vectors Y
    switch(order) {
      case 1:
        cblas_dcopy (neqn, Y[idx],1,    Yf[0],1);   break;
      case 2:
      { cblas_dcopy (neqn, Yf[order],1, Yf[0],1);
        cblas_dcopy (neqn, Y[idx],1,    Yf[1],1);   break;
      }
      case 3:
      { cblas_dcopy (neqn, Y[1],1,      Yf[0],1);
        cblas_dcopy (neqn, Yf[order],1, Yf[1],1);
        cblas_dcopy (neqn, Y[idx],1,    Yf[2],1);   break;
      }
      case 4:
      { cblas_dcopy (neqn, Yf[2],1, Yf[0],1);
        cblas_dcopy (neqn, Y[2],1,        Yf[1],1);
        cblas_dcopy (neqn, Yf[4],1,   Yf[2],1);
        cblas_dcopy (neqn, Y[idx],1,      Yf[3],1); break;
      }
    }
  }    
  
  
   
//**************************************************************************
// Update intermediate vectors before the next integration step
//**************************************************************************
/*void Update_intermediate_vectors(const int order, const int neqn, double** Y, double* Y1)
{   double * Ytmp=Y[0];
    for (int i = 1; i < order; i++) {
        Y[i-1]=Y[i];
    }    
    Y[order-1]=Ytmp;
    cblas_dcopy(neqn, Y1, 1, Y[order-1], 1);
}*/

  //**************************************************************************
  // Update the coarse vector Y1_c from the fine vector Y1_f 
  // by using Richardson extrapolation 
  //**************************************************************************
  void VSIIE_Solver::Extrapolate(const double * h_vector, 
                                const double *Y1_f, double * Y1_c){
    double alpha, beta;
    if (order==1) {
      alpha=-1.0;
      beta=2.0;
    } 
    else if (order==2) {
      const double denom=7.0*h_vector[0]+5.0*h_vector[1];
      alpha=-(h_vector[0]+3.0*h_vector[1])/denom;
      beta=8.0*(h_vector[0]+h_vector[1])/denom;
    }  
    else if (order==3) {
      const double h1_sq=h_vector[1]*h_vector[1];
      const double h2_sq=h_vector[2]*h_vector[2];
      const double h1h0=h_vector[1]*h_vector[0];
      const double h2h1=h_vector[2]*h_vector[1];
      const double h2h0=h_vector[2]*h_vector[0];
      const double h2ph1=h_vector[2]+h_vector[1];
      const double denom=15.0*h1_sq +14.0*h1h0 + 30.0*h2h1+10.0*h2h0 + 11.0*h2_sq;
      alpha=-(h1_sq +2.0*(h1h0+h2h1+3.0*h2h0)+5.0*h2_sq)/denom;
      beta=16.0*(h2ph1*(h2ph1+h_vector[0]))  /denom;
    }
    else //(order==4) 
    {

      const double h_np1=h_vector[idx]; // h_{n+1}
      const double h_n  =h_vector[idx-1]; //h_n
      const double h_nm1  =h_vector[idx-2]; //h_{n-1}
      const double h_nm2  =h_vector[idx-3]; //h_{n-2}

      const double C_1 = 31.0 * pow(h_n, 3) + 30.0 * pow(h_n, 2) * h_nm2 + 60.0 * pow(h_n, 2) * h_nm1
            + 93.0 * pow(h_n, 2) * h_np1 + 28.0 * h_n * h_nm2 * h_nm1 + 60.0 * h_n * h_nm2 * h_np1
            + 28.0 * h_n * pow(h_nm1, 2);
      const double C_2 = 120.0 * h_n * h_nm1 * h_np1 + 93.0 * h_n * pow(h_np1, 2) + 20.0 * h_nm2 * h_nm1 * h_np1
            + 22.0 * h_nm2 * pow(h_np1, 2) + 20.0 * pow(h_nm1, 2) * h_np1 + 44.0 * h_nm1 * pow(h_np1, 2)
            + 23.0 * pow(h_np1, 3);
      const double term = (h_np1 + h_n) * (h_np1 + h_n + h_nm1) * (h_np1 + h_n + h_nm1 + h_nm2);
      const double denom = C_1 + C_2;
      alpha = ((denom-32)*term) / denom;
      beta = (32 *term) / denom;
    }  

    cblas_dscal (neqn, alpha, Y1_c , 1);
    cblas_daxpy (neqn, beta, Y1_f , 1, Y1_c, 1);    

  }



//*******************************************************************
//*******************************************************************
// IMPLEMENTATION OF THE PUBLIC FUNCTIONS
//*******************************************************************

//******************************************************
// Constructor of the class SBDF
//******************************************************
VSIIE_Solver::VSIIE_Solver (const int _order, IVP_ODE* _IVP)
//*****************************************************
{ 
  order=_order;
  idx=order-1; 
  IVP=_IVP;
  neqn = IVP->get_num_ODEs();
  nnz=IVP->get_nnz();
  // Create LIS vectors
  lis_vector_create(0, &b);
  lis_vector_create(0, &x);
  lis_vector_set_size(b, 0, neqn);
  lis_vector_set_size(x, 0, neqn);

  // Create Intermediate Residual Vector R
  R = new double[neqn];
  // Create and initialize LIS solver with the suitable switches
  lis_solver_create(&solver);
  char opts[] = "-i bicgstab -p ilut  -tol 1.0e-15";
  //char opts[] = "-i icr  -p jacobi -tol 1.0e-15";
  //char opts[] = "-i fgmres -restart 20 -p ilut  -tol 1.0e-12";
  lis_solver_set_option(opts, solver);
  init_CSR_LIS_matrix(&As);
 }

 //******************************************************
// Destructor of the class SBDF
//******************************************************
VSIIE_Solver::~VSIIE_Solver (){

  // Destroy LIS vectors b ad x
  lis_vector_destroy(x);
  lis_vector_destroy(b);

  // Destroy LIS matrices As and As2
  lis_matrix_destroy(As);
  lis_matrix_destroy(As2); 

  // Destroy the LIS solver
  lis_solver_destroy(solver);
  delete[] R;
  
}

//**************************************************************************
// Update intermediate vectors before the next integration step
//**************************************************************************
void VSIIE_Solver::Update_intermediate_vectors(double** Y, double* Y1){   
  double * Ytmp=Y[0];
  for (unsigned i = 1; i < order; i++) {
      Y[i-1]=Y[i];
  }    
  Y[idx]=Ytmp;
  cblas_dcopy(neqn, Y1, 1, Y[idx], 1);
}




//***************************************************
// Store the results of an experiment in a data file
//***************************************************
void VSIIE_Solver::Store_result(const double tf, 
                              const double* Y1, 
                              const double tol)
//***************************************************

{
  ostringstream sstr;
  sstr << IVP->get_name()<<"_" << neqn<<"_tf-"<<tf<<"_tol-"<<tol << ".dat";
  ofstream str;


  str.open( sstr.str().c_str(),ios_base::out );

  if (!str.is_open()) {
        std::cerr << "Failed to open: " << sstr.str().c_str()<< std::endl;
        return;
  }      
  str << scientific << setprecision(16);
  str << IVP->get_name() << endl;
  str << neqn << endl;
  for( int i=0; i<neqn; ++i )
    {
	    str << Y1[i] << endl;
     }

     str.close();

     return;
}


//*******************************************************************
// Function implementing one time step using the 
// VSIIE-order Time Integrator (order=1,..,4)
//*******************************************************************
void  VSIIE_Solver::Time_Step(const double t, const double *h_vector, double** Y, double* Y1, 
const bool variable_tstep, int & newton_iters){

  if (variable_tstep) {
    Update_coefs(order, h_vector); 
  }
  newton_iters=0;
  const double h=h_vector[idx];
  h_coef_idx_F1=coef_idx_F1*h;
  h_coef_idx_F2=coef_idx_F2*h;

  compute_matrix(t,h,Y);
  // Get initial approximation of Y1=Yidx+h*Yidx
  cblas_dcopy(neqn, Y[idx], 1, Y1, 1); 
  cblas_daxpy(neqn, h, Y[idx], 1, Y1, 1);
  

  // Modified Newton Iteration to approximate Y1 at the current time step
  double norm = 1.0e23;
  int it = 0;
  bool convergence = false;
  const double EPSTOL = 1.0e-12;
  const int max_iterations = 200;
  // Compute the constant term in a time step of vector R
  compute_RHS0(t, h_vector, Y);
  //cout<< "Initial RHS norm = " << cblas_dnrm2(neqn, R, 1)/ sqrt((double)neqn) << endl;
  
  while (!convergence) {
    newton_iters++;
    // Complete the  RHS vector, adding the variable terms  
    compute_RHS(t, h, Y1, b->value);
    
    //cout << "Iteration " << it << " with residual norm =" << cblas_dnrm2(neqn, b->value, 1)/ sqrt((double)neqn) << endl;


    // Initialize vector x of the linear solver and b=R
    // to solve iteratively the system A*x=b to approximate the residual vector x
    lis_vector_set_all(0, x);
    int error = lis_solve(As, b, x, solver);
    if (error != 0) {
      cout << "######## ############# ERROR IN LINEAR SYSTEM SOLUTION" << endl; exit(0);
    }
    else {
      double lis_time; int iters; 
      lis_solver_get_time(solver,&lis_time); lis_solver_get_iter(solver,&iters);
      //cout<<"............. Iterative Linear System Solution:  ITERS= "<<iters
      //                                 <<"       TIME = " <<lis_time<<"  seconds"<<endl;
    }

    
    //Update Y1=Y1-x
    cblas_daxpy(neqn, -1.0, x->value, 1, Y1, 1);
    
    // Compute the norm of the residual vector x
    norm = cblas_dnrm2(neqn, x->value, 1)/ sqrt((double)neqn);
    convergence = (norm < EPSTOL) || (it > max_iterations);
  }  // End of modified Newton Iteration
  
  //cout<<"........................................<< Newton Iteration= "<<newton_iters<<"  >>"<<endl;
  
  if (newton_iters >= max_iterations) {
    cerr << ".................. ##################### Max iterations achieved!!!!" << endl;
    exit(0);
  }

}


//***************************************************
// Function implementing the order 1-4 VSIIE Time Integrator
// It assumes a constant time step h
//***************************************************
void VSIIE_Solver::Const_dt_Integrate(const double t0, const double tf,
    const double h, double** Y_, double* Y1)  {        
    //***************************************************
    const double EPSTOL=1.0e-20;
    const int neqn = IVP->get_num_ODEs();

    // Initialize previous step vectors Y[0],..., Y[order-1]
    double *Y[order];
    for (unsigned i = 0; i < order; i++) {
        Y[i] = new double[neqn];
        cblas_dcopy(neqn, Y_[i], 1, Y[i], 1); 
    }
       
    double t = t0;  
    t = t + idx * h;

    // Set h_vector for a constant stepsize strategy
    double h_vector[4]={h,h,h,h};
    Update_coefs(order, h_vector);

    //Fixed_update_coefs(order);
    for (unsigned i = 0; i < order; i++) h_vector[i] = h;
    // Set the constant stepsize strategy
    bool variable_tstep=false;
    int N_steps=0;
    double h_now=h;
    bool end=false;
    //*********  TIME LOOP  ****************
    while (!end) {
    //**************************************
      // Compute the distance to the final time tf
      const double distance = tf - t;  
      // If (distance<h), Use variable stepsize strategy with h=distance
      if (distance<h) {    
        h_now=distance;
        h_vector[idx]=h_now;
        variable_tstep=true;  
      } 
      // Perform a new integration step with stepsizes in h_vector
      //cout<<"Perform a new integration VSIIE-"<<order<<" step with stepsize "<<h_vector[0]<< "with t0="<<t<<endl;
      int  newton_iters;
      Time_Step(t, h_vector, Y, Y1, variable_tstep, newton_iters);
      // Update previous step vector Y
      Update_intermediate_vectors(Y, Y1);
      // Update time t
      t+=h_now;  
      // Increase the number of time steps
      N_steps++;
      //Check the end of the time integration
      if (fabs(tf-t)<EPSTOL) {end=true;}
      // cout<<"VSIIE Solution obtained in tf="<<t<<endl;
    } // End of time stepping
    // Free memory for the intermediate vectors Y[0],..., Y[order-1] 
    for (unsigned i = 0; i < order; i++) 
        delete[] Y[i];
}


//***************************************************
// Function implementing the 4th-order SBDF Time Integrator
// It assumes a constant time step h
//***************************************************
void VSIIE_Solver::Const_dt_Integrate_SBDF4(const double t0, const double tf,
    const double h, double** Y_, double* Y1)  {        
    //***************************************************
    const double EPSTOL=1.0e-20;
    const int neqn = IVP->get_num_ODEs();
    const int order=4;
    // Initialize previous step vectors Y[0],..., Y[order-1]
    double *Y[order];
    for (unsigned i = 0; i < order; i++) {
        Y[i] = new double[neqn];
        cblas_dcopy(neqn, Y_[i], 1, Y[i], 1); 
    }
       
    double t = t0;  
    t = t + idx * h;

    // Set h_vector for a constant stepsize strategy
    double h_vector[4]={h,h,h,h};
    Update_coefs_SBDF4(h_vector);

    //Fixed_update_coefs(order);
    for (unsigned i = 0; i < order; i++) h_vector[i] = h;
    // Set the constant stepsize strategy
    bool variable_tstep=false;
    int N_steps=0;
    double h_now=h;
    bool end=false;
    //*********  TIME LOOP  ****************
    while (!end) {
    //**************************************
      // Compute the distance to the final time tf
      const double distance = tf - t;  
      // If (distance<h), Use variable stepsize strategy with h=distance
      if (distance<h) {    
        h_now=distance;
        h_vector[idx]=h_now;
        variable_tstep=true;  
      } 
      // Perform a new integration step with stepsizes in h_vector
      //cout<<"Perform a new integration VSIIE-"<<order<<" step with stepsize "<<h_vector[0]<< "with t0="<<t<<endl;
      int  newton_iters;
      Time_Step(t, h_vector, Y, Y1, variable_tstep, newton_iters);
      // Update previous step vector Y
      Update_intermediate_vectors(Y, Y1);
      // Update time t
      t+=h_now;  
      // Increase the number of time steps
      N_steps++;
      //Check the end of the time integration
      if (fabs(tf-t)<EPSTOL) {end=true;}
    } // End of time stepping
    // Free memory for the intermediate vectors Y[0],..., Y[order-1] 
    for (unsigned i = 0; i < order; i++) 
        delete[] Y[i];    
}








//********************************************************************************************
// Compute Coarse and Fine Time steps using Variable VSIIE scheme as well as the LTE
//********************************************************************************************
void VSIIE_Solver::Variable_Time_Step(const double t, double * h_vector, 
                                     double * h_vector_half, double * h_vector_half2, double ** Y, 
                                     double ** Yf, double * Y1, double * LTE, double *epsilon_c, 
                                     int& total_iters)
{ const bool  variable_tstep=true;

  // COARSE STEP
  int newton_iters; 
  // Reset the total number of Newton iterations
  total_iters=0;
  Time_Step(t, h_vector, Y, Y1, variable_tstep, newton_iters);
  total_iters+=newton_iters;
  h_vector_half[idx]=h_vector[idx]/2.0;  
  for (int i=0;i<idx;i++) {
    h_vector_half2[i]=h_vector_half[i+1];
  }
  h_vector_half2[idx]=h_vector[idx]/2.0;

  // FINE STEP 1
  Time_Step(t,h_vector_half, Yf, Yf[order], variable_tstep, newton_iters);
  total_iters+=newton_iters;
  //FINE STEP 2
  Time_Step(t+h_vector_half[idx], h_vector_half2, &(Yf[1]), Yf[order+1], variable_tstep, newton_iters);
  total_iters+=newton_iters;
  // Computation of the LTE vector and the scalar error
  *epsilon_c=compute_LTE(h_vector, Y1, Yf[order+1], LTE);

  //cout<<"LTE="<<*epsilon_c<<endl;
}


//***************************************************
// Function implementing the order 1-4 VSIIE Time Integrator
// It assumes a adaptive time step
//***************************************************
void VSIIE_Solver::Adaptive_dt_Integrate(const double t0, const double tf,
    const double h, double** Y_init, double ** Yf_init, double* Y1, 
    const double tol, int *nsteps, int * n_isteps, double & av_iters,
    double alpha=0.8, double eta_min= 0.5, 
    double eta_max = 4)  {        
//***************************************************
  
  const  double range_factor=0.2;
  const double range=tol*range_factor;
  const double  p_inv=1.0/(order+1);
  const double EPSTOL=1.0e-30;
  const int neqn = IVP->get_num_ODEs();
  // Initialize intermediate stage vectors Y, Tf and the LTE vector
  double *Y[order], *Yf[order+2];
  double *LTE=new double[neqn];
  for (unsigned int i = 0; i < (order+2); i++) {
    Yf[i] = new double[neqn]; 
  }
  for (unsigned int i = 0; i < order; i++) {
    Y[i] = new double[neqn];
    cblas_dcopy(neqn, Y_init[i], 1, Y[i], 1);
    cblas_dcopy(neqn, Yf_init[i], 1, Yf[i], 1);    
  }

  // Init the time counter
  double t = t0;  
  const int idx=order-1;
  t = t + idx * h;
  double h_vector[4]={h,h,h,h}; 
  double h_vector_half[4]={h/2,h/2,h/2,h/2};
  double h_vector_half2[4]; double next_dt=h;
  int N_steps=0; 
  int total_isteps=0;
  // Reset the total number of Newton iterations;
  av_iters=0.0;

  bool end=false;
  //*********  TIME LOOP  ****************
  while (!end) {
  //**************************************
    const double distance = tf - t;  
    bool accepted=false;
    // Local trucation error epsilon_c and factor to modify dt_new
    double epsilon_c,factor;
    // number of iteration for a single adaptive time step
    unsigned i_step=0;

    h_vector[idx]=next_dt;
    //**************************
    // INNER LOOP OF A TIME STEP
    //**************************
    while (! accepted)  {

      // If the hnew is greater than distance=tf-t, then new h=tf-t and last step is possible 
      if (h_vector[idx]>distance) {
        h_vector[idx]=distance; 
        //cout<<"** T= "<<t<<"   ***********COMPUTING FINAL DT="
        //    <<h_vector[idx]<<"  *****  h_const= "<<h<<endl<<flush;        
      }     

      // New coarse and fine integration step
      int total_iters;
      Variable_Time_Step(t, h_vector, h_vector_half, h_vector_half2, Y, Yf, Y1, LTE, 
        &epsilon_c, total_iters);
      av_iters+=total_iters;

      // Check the error condition
      double error_diff=epsilon_c-tol;
      accepted=(error_diff<=range);
      factor=min(max(alpha*pow(tol/epsilon_c,p_inv),eta_min  ),eta_max);

      if (accepted){
        if (fabs(error_diff)<=range) next_dt=h_vector[idx];
        else                         next_dt=h_vector[idx]*factor;
      } 
      else { 
       // cout<< "TOL="<<tol<<" ..... STEP NOT ACCEPTED    NEXT DT="<< next_dt
       //     <<"     ERROR= "<<epsilon_c<<"   FACTOR=  "<<factor<<endl<<flush;
        h_vector[idx]*=factor;
      }  
      i_step++;     
    //*************************************  
    }   // END OF INNER LOOP OF A TIME STEP
    //*************************************

    total_isteps+=i_step;
    // Update t using the old time step stored in h_vector[idx]
    t+=h_vector[idx]; N_steps++;

    // Check the end ot time integration 
    end=(fabs(tf-t)<EPSTOL); 

    // Richardson Extrapolation using LTE:   Y1=Y1-LTE
    cblas_daxpy(neqn, -1.0, LTE, 1, Y1, 1);

    if (!end) {
      // Update intermediate vectors and stepsizes
      Update_adaptive_intermediate_vectors(h_vector_half, h_vector_half2, h_vector, Yf, Y, Y1);
    }      
  //*********************
  } // End of time stepping
  //*********************
  *nsteps=N_steps;  
  *n_isteps=total_isteps;
  av_iters=av_iters/(3*total_isteps); 
  for (unsigned i = 0; i < order; i++)
    delete[] Y[i];
  for (unsigned i = 0; i < (order+2); i++)
    delete[] Yf[i];
  delete[] LTE;    
}

/*
This is the implementation of the C++ class "ABM_Solver" which implements the Variable stepsize Adams-Bashforth-Moulton (ABM) solvers 
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

#include "IVP_ODE.h"
#include "ABM_Solver.h"
#include "utils.h"
using namespace std;
#include <assert.h>



//*******************************************************************
//*******************************************************************
// IMPLEMENTATION OF THE PRIVATE FUNCTIONS
//*******************************************************************
//******************************************************
// Update coefficients of Variable stepsize ABM-order scheme
//******************************************************
void ABM_Solver::Update_coefs(const int order, const double * h_vector) {

 double w[idx];
 for (int i=0;i<=idx;i++){
   w[i]=h_vector[idx-i]/h_vector[idx];
 } 
 if (order==2){
    const double A1=w[0];
    const double denom=2.0*w[0];
    coef_AB[idx][1]=(1.0+denom)/denom;
    coef_AB[idx][0]=-1.0/denom; 

    const double LTE_denom=3.0*(1.0+A1);
    LTE_coef[idx]=1.0/LTE_denom;

  } 
  else if (order==3){ 
    const double A1=w[0];
    const double A2=w[0]+w[1];   
    coef_AB[idx][2]= (6.0*A1*A2+3.0*A1+3.0*A2+2.0)/(6.0*A1*A2);
    coef_AB[idx][1] = -(3.0*A2+2.0)/(6.0*A1*(A2-A1));
    coef_AB[idx][0] = (3.0*A1+2.0)/(6.0*A2*(A2-A1));

    coef_AM[idx][2]= (3.0*A1+2.0)/(6.0*(1.0+A1));
    coef_AM[idx][1] = (3.0*A1+1.0)/(6.0*A1);
    coef_AM[idx][0] = -1.0/(6.0*(1.0+A1));

    const double LTE_denom=2.0*(3.0*A1+2.0)*(1.0+A2); 
    LTE_coef[idx]=(2.0*A1+1.0)/LTE_denom;


  }
  else if (order==4){ 
    const double A1=w[0];
    const double A2=w[0]+w[1];
    const double A3=w[0]+w[1]+w[2];

    coef_AB[idx][3]=(3.0+4.0*(A1+A2+A3)+6.0*(A1*A2+A1*A3+A2*A3)+12.0*A1*A2*A3)/(12.0*A1*A2*A3);
    coef_AB[idx][2]=-(3.0+4.0*(A2+A3)+6.0*A2*A3)/(12.0*A1*(A2-A1)*(A3-A1));
    coef_AB[idx][1]=(3.0+4.0*(A1+A3)+6.0*A1*A3)/(12.0*A2*(A2-A1)*(A3-A2));
    coef_AB[idx][0]=-(3.0+4.0*(A1+A2)+6.0*A1*A2)/(12.0*A3*(A3-A1)*(A3-A2));

    const double A=A1, B=A2;
    coef_AM[idx][3]=(6.0*A*B+4.0*A+4.0*B+3.0)/(12.0*(1.0+A)*(1.0+B));
    coef_AM[idx][2]=(6.0*A*B+2.0*A+2.0*B+1.0)/(12.0*A*B);
    coef_AM[idx][1]=-(2.0*B+1.0)/(12.0*A*(1.0+A)*(B-A));
    coef_AM[idx][0]=(2.0*A+1.0)/(12.0*B*(1.0+B)*(B-A));


    const double LTE_num=(10.0*A1*A2 + 5.0*A1 + 5.0*A2 + 3.0); 
    const double LTE_denom=5.0 * (1.0+A3) *(6.0*A1*A2 + 4.0*A1 + 4.0*A2 + 3.0); 
    LTE_coef[idx]=LTE_num/LTE_denom;




  }

}



//******************************************************
// Compute LTE estimate from the predictor solution Y_AB 
// and the corrector solution Y_AM
//******************************************************
double ABM_Solver::compute_LTE(double * Y_AB, double * Y_AM)
  //******************************************************
  { double LTE[neqn];
    const double factor=LTE_coef[idx];
    // Compute difference LTE=factor*(Y_AM-Y_AB)
    cblas_dcopy(neqn, Y_AM, 1, LTE, 1);
    cblas_daxpy(neqn, -1.0, Y_AB, 1, LTE, 1);
    cblas_dscal(neqn, factor,      LTE, 1);
    return(cblas_dnrm2(neqn, LTE, 1));
  }
  //******************************************************


//*******************************************************************
//*******************************************************************
// IMPLEMENTATION OF THE PUBLIC FUNCTIONS
//*******************************************************************

//******************************************************
// Constructor of the class ABM
//******************************************************
ABM_Solver::ABM_Solver (const int _order, IVP_ODE* _IVP)
//*****************************************************
{ 
  order=_order;
  idx=order-1; 
  IVP=_IVP;
  neqn = IVP->get_num_ODEs(); 

 }

 //******************************************************
// Destructor of the class ABM
//******************************************************
ABM_Solver::~ABM_Solver (){
 
}

//**************************************************************************
// Update intermediate vectors and stepsizes before the next integration step
//**************************************************************************
void ABM_Solver::Update_intermediate_vectors(double *h_vector, 
                                             double** Y, double* Y1) { 
//***************************************************************************
  double * Ytmp=Y[0];
    //Update h_vectors 
    for (int i=0;i<idx;i++) {
      h_vector[i]=h_vector[i+1];
    }

  for (unsigned i = 1; i < order; i++) {
      Y[i-1]=Y[i];
  }    
  Y[idx]=Ytmp;
  cblas_dcopy(neqn, Y1, 1, Y[idx], 1);
}




//***************************************************
// Store the results of an experiment in a data file
//***************************************************
void ABM_Solver::Store_result(const double tf, 
            const double* Y1, const double tol)
//***************************************************
{
  ostringstream sstr;
  sstr << IVP->get_name()<<"_" << neqn<<"_" 
          <<"_tf-"<<tf<<"_tol-"<<tol << ".dat";
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


  //**************************************************************************
  // Adams-Bashforth (AB) predictor with variable time step
  //**************************************************************************
  void ABM_Solver::AB_Predictor(const double t, double *h_vector, 
                               double ** Y, double *Y1){
    double* DY = new double[neqn];
    const double h=h_vector[idx];
    double ti[4];
    ti[idx]=t;
    for (int i = idx-1; i >= 0; i--){
      ti[i]=ti[i+1]-h_vector[i];
    }

    //cout<<"AB_Predictor: T= "<<t<<"  h= "<<h_vector[idx]<<"    idx="<<idx<<endl<<flush;
    //Y1=Y[idx] 
    cblas_dcopy(neqn, Y[idx], 1, Y1, 1);

    for (int i = 0; i <=idx; i++) {
      IVP->feval(ti[i], Y[i], DY);
      const double coef=coef_AB[idx][i]*h;
      cblas_daxpy(neqn, coef, DY, 1, Y1, 1);
    }

            
    delete [] DY;
  } 



//**************************************************************************
// Adams-Moulton (AM) corrector with variable time step with variable time step
// and predictor solution Y_AB and its evaluation FY_AB
//**************************************************************************
void  ABM_Solver::AM_corrector(const double t, const double * h_vector, 
                                       double** Y, double *Y_AB, double *Y1) {
  double FY_AB[neqn];
  const double h=h_vector[idx];

  double ti[4];
  ti[idx]=t;
  for (int i = idx-1; i >= 1; i--){
    ti[i]=ti[i+1]-h_vector[i];
  }

  //Y1=Y_n
  cblas_dcopy(neqn, Y[idx], 1, Y1, 1);
  // Compute approximation  of the evaluation of feval in the 
  // next time step using the predictor solution Y_AB
  IVP->feval(t + h, Y_AB, FY_AB);  
  
 //Y1+= coef_AM[idx][idx]*h*feval(t+h,Y_AB) 
  cblas_daxpy(neqn, coef_AM[idx][idx]*h, FY_AB, 1, Y1, 1);

  for (int i = idx; i >=1; i-- ) {
      IVP->feval(ti[i], Y[i], FY_AB);
      cblas_daxpy(neqn, coef_AM[idx][i-1]*h, FY_AB, 1, Y1, 1);
  }          
} 







//*******************************************************************
// Function implementing one time step using the 
// ABM-order Time Integrator (order=1,..,4)
// Generate the corrector solution Y1 and the 
// error estimate of the local truncation error
//*******************************************************************
void  ABM_Solver::Time_Step(const double t, const double *h_vector, 
                            double** Y, double* Y1, 
                            const bool variable_tstep, 
                            double & error_estimate){
  double Y_AB[neqn];    
  double FY_AB[neqn];
     
  if (variable_tstep) {
    Update_coefs(order, h_vector); 
  }
  const double h=h_vector[idx];
  
  // Get initial approximation Y1 using 
  // the Adams-Bashforth predictor
  AB_Predictor(t, h_vector, Y, Y_AB);
  //print_vector(neqn, Y_AB, "YAB_PREDICTOR-AB");
  // Compute the corrector solution Y1 
  // using the Adams-Moulton corrector
  AM_corrector(t, h_vector, Y, Y_AB, Y1);
  //print_vector(neqn, Y1, "Y1_CORRECTOR-AM");
  error_estimate=compute_LTE(Y_AB, Y1);
  //cout << "TIME=" << t << "  LTE= " << error_estimate << endl;

                            
}





//***************************************************
// Function implementing the order 1-4 ABM Time Integrator
// It assumes a constant time step h
//***************************************************
void ABM_Solver::Const_dt_Integrate(const double t0, const double tf,
    const double h, double** Y_, double* Y1)  {        
    //***************************************************
    const double EPSTOL=1.5e-12;
    double error_estimate;
    const int neqn = IVP->get_num_ODEs();
   
    // Initialize previous step vectors Y[0],..., Y[order-1]
    double *Y[order];
    for (unsigned i = 0; i < order; i++) {
        Y[i] = new double[neqn];
        cblas_dcopy(neqn, Y_[i], 1, Y[i], 1); 
    }
       
    double t = t0;  
    const int idx=order-1;
    t = t + idx * h;
    // Set h_vector for a constant stepsize strategy
    double h_vector[4]={h,h,h,h};
    
    // Set the constrant stepsize strategy
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
        cout<<"** T= "<<t<<"   ***********COMPUTING FINAL DT IN ABM!!!="
            <<h_vector[idx]<<"  *****  h_const= "<<h<<endl<<flush; 
      } 
      // Perform a new integration step with stepsizes in h_vector
      int  newton_iters;
      Time_Step(t, h_vector, Y, Y1, variable_tstep, error_estimate);

      // Update previous step vector Y
      Update_intermediate_vectors(h_vector, Y, Y1);
      // Update time t
      t+=h_now;  
      // Increase the number of time steps
      N_steps++;
      //Check the end of the time integration
      if (fabs(tf-t)<EPSTOL) {end=true;}

      //cout << "ABM: T= " << t << "  h= " << h_now << "  #Steps= " << N_steps<<endl;
    } // End of time stepping
    // Free memory for the intermediate vectors Y[0],..., Y[order-1] 
    for (unsigned i = 0; i < order; i++) 
        delete[] Y[i];
}


///***************************************************
  // Function implementing the order 1-4 ABM Time Integrator
  // It assumes a adaptive time step
  //**************************************************
  void ABM_Solver::Adaptive_dt_Integrate(const double t0, const double tf,
    const double h, double** Y_init, double* Y1, 
    const double tol, int *nsteps, int * n_isteps, 
    double alpha=0.9, double eta_min= 0.5, double eta_max = 5)  {        
//***************************************************
  const double  p_inv=1.0/(order+1);
  const double EPSTOL=1.0e-30;
  const int neqn = IVP->get_num_ODEs();
  
  double epsilon;

  // Initialize intermediate stage vectors Y, Tf and the LTE vector
  double *Y[order];
  for (unsigned int i = 0; i < order; i++) {
    Y[i] = new double[neqn];
    cblas_dcopy(neqn, Y_init[i], 1, Y[i], 1);  
  }

  // Init the time counter
  double t = t0;  
  const int idx=order-1;
  t = t + idx * h;
  double h_vector[4]={h,h,h,h}; 
  int N_steps=0; 
  int total_isteps=0;
   double next_dt=h;  // Next time step

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

      Time_Step(t, h_vector, Y, Y1, true, epsilon);
      
      accepted=(epsilon<=tol);
      if  ( isnan(epsilon)  || isinf(epsilon) ) {

        cout<<"** T= "<<t<<"   ***********ERROR IN ABM: NAN or INF in LTE!!!="
            <<epsilon<<"  *****  h_const= "<<h<<endl<<flush;

        exit(EXIT_FAILURE);
      }
          
      factor=min(max(alpha*pow(tol/epsilon,p_inv),eta_min  ), eta_max);
      if (!accepted){
        h_vector[idx]*=factor;
      } 
      else { 
        next_dt=h_vector[idx]*factor;
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

    
    if (!end) {
      // Update intermediate vectors and stepsizes
      Update_intermediate_vectors(h_vector, Y, Y1);
    }      
  //*********************
  } // End of time stepping
  //*********************

  cout << "ABM: T= " << t << "  h= " << h_vector[idx] << "  #Steps= " << N_steps<<endl;
  *nsteps=N_steps;  
  *n_isteps=total_isteps;

  for (unsigned i = 0; i < order; i++)
    delete[] Y[i];
}

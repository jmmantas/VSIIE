/*
This C++ driver program run experiments associated with a given IVP model using a particular VSIIE method.  
It generates the numerical error with respect to a reference solution (previously generated in a text file) and the execution time. 

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


#include <fstream>
#include <iostream>
#include <iomanip>
#include <sys/stat.h>
#include <sstream>

#include  <omp.h>
#include "IVP_ODE_advdiff1d.h"
#include "IVP_ODE_stiff_brusselator.h"
#include "IVP_ODE_combustion.h"

#include "VSIIE_Solver.h"
#include "utils.h"
using namespace std;
using namespace std::chrono;

//******************************************************
//Set up condition Vectors for the 
// coarse approximation (with "order" items): Y_init[0]=Y[0] in t=t0,
// Y_init[1] in t=t0+h, ... Y_init[order-1] in t=t0+(order-1)*h
// or also for the fine approximation (with "order+2" items): 
// Yf_init[0] in t=t0, Yf_init[1] in t=t0+h/2,Yf_init[2] in t+h, ...
//******************************************************
void Setup_initial_vectors(const int num_items, 
                        const int neqn, double ** Y_init){
    for (int i = 0; i < num_items; i++) {
        Y_init[i] = new double[neqn];
    }
}


//**************************************************************************
// Free auxiliary vectors in dynamic memory 	
//**************************************************************************
void Free_initial_vectors(const int num_items, double **Y_init)
{
    for (int i = 0; i < num_items; i++) 
        delete[] Y_init[i];
}





//**************************************************************************
// Initialize  intermediate vectors (Y_init) using 4th order Runge-Kutta 
// to start the multistep integration. Assumes Y_init[0] is already set
//**************************************************************************
void Init_coarse_Vectors_RK(const int order, const int neqn,  const double t0, 
                  const double h,  IVP_ODE* IVP, double ** Y_init) 
{
  const int RK_order=4;
  // Approximation of the intermediate step initial values 
  // using RK4 and very small stepsize
	const double h_RK=1.0e-8;
  // Initialize initial condition vectors (Y_init[...]) for coarse time integration
  double t2=t0;           
  for (int i = 1; i < order; i++){  
    RK(RK_order, IVP, t2, t2+h, h_RK, Y_init[i-1], Y_init[i]);
    t2+=h;
  } 
}





//**************************************************************************
// Set the  initial coarse vectors Y_init using 4th order Runge-Kutta 
// to start the multistep integration
//**************************************************************************
void Init_coarse_Vectors_SBDF4(const int order, const int neqn,  const double t0, 
                  const double h,  IVP_ODE* IVP, double ** Y_init) 
{ 
  const int order_SBDF4=4;
  // Declare and  Setup initial condition vectors Y_init
  double *Y_init_SBDF4[order_SBDF4];
  // Setup coarse initial vectors Y_init_SBDF4
  Setup_initial_vectors(order_SBDF4, neqn, Y_init_SBDF4);

  IVP->init(Y_init[0]); //First  step initial value
  // Initialize initial condition vectors (Y_init_SBDF4[...]) for coarse time integration

  double t2=t0;
  double h_SBDF4=(h<=1.0e-6)?h/2.0:1.0e-6;
  // Init VSIIE solver 
  VSIIE_Solver * VSIIE_IVP=new VSIIE_Solver(order_SBDF4, IVP);           
  for (int i = 1; i < order; i++){
    cblas_dcopy (neqn, Y_init[i-1],1, Y_init_SBDF4[0] , 1);
    Init_coarse_Vectors_RK(order_SBDF4, neqn, t2, h_SBDF4, IVP, Y_init_SBDF4); 
    // Execution of the VSIIE Solver to approximate accurately the solution
    VSIIE_IVP->Const_dt_Integrate_SBDF4(t2, t2+h, h_SBDF4, Y_init_SBDF4, Y_init[i]);
    t2+=h;
  } 
  Free_initial_vectors(order_SBDF4, Y_init_SBDF4);
  delete VSIIE_IVP;
}

//**************************************************************************
// Set the fine initial approximation vectors Yf_init using 4th order Runge-Kutta 
// using the initial coarse approximation vector Y_init
//**************************************************************************
void Init_fine_Vectors_SBDF4(const int order, const int neqn,  const double t0, 
                  const double h,  IVP_ODE* IVP, double ** Y_init, double ** Yf_init) 
{
  const int order_SBDF4=4;
  // Declare and  Setup initial condition vectors Y_init
  double *Y_init_SBDF4[order_SBDF4];
  // Setup coarse initial vectors Y_init_SBDF4
  Setup_initial_vectors(order_SBDF4, neqn, Y_init_SBDF4);
  
  // Initialize initial condition vectors (Y_init_SBDF4[...]) for coarse time integration
  double h_SBDF4=(h<=1.0e-6)?h/4.0:0.5e-7;

  // Init VSIIE solver 
  VSIIE_Solver * VSIIE_IVP=new VSIIE_Solver(order_SBDF4, IVP);


  double t2;
  // Set initial condition vectors (Yf_init[...]) for 
  // fine time integration (half of the time step size)
  const double h_over2=h/2;
  if (order==2) {  // order 2
    t2=t0;
    // Execution of the VSIIE Solver to approximate accurately the solution
    cblas_dcopy (neqn, Y_init[0],1, Y_init_SBDF4[0] , 1);
    Init_coarse_Vectors_RK(order_SBDF4, neqn, t2, h_SBDF4, IVP, Y_init_SBDF4); 
    VSIIE_IVP->Const_dt_Integrate_SBDF4(t2, t2+h_over2, h_SBDF4, Y_init_SBDF4, Yf_init[0]);
  }
  else if (order==3) { // order 3
    cblas_dcopy (neqn, Y_init[1],1, Yf_init[0] , 1);
    t2=t0+h;
    cblas_dcopy (neqn, Yf_init[0],1, Y_init_SBDF4[0] , 1);
    Init_coarse_Vectors_RK(order_SBDF4, neqn, t2, h_SBDF4, IVP, Y_init_SBDF4); 
    VSIIE_IVP->Const_dt_Integrate_SBDF4(t2, t2+h_over2, h_SBDF4, Y_init_SBDF4, Yf_init[1]);            
  } 
  else if (order==4) { // order 4
    t2=t0+h;
    cblas_dcopy (neqn, Y_init[1],1, Y_init_SBDF4[0] , 1);
    Init_coarse_Vectors_RK(order_SBDF4, neqn, t2, h_SBDF4, IVP, Y_init_SBDF4);
    VSIIE_IVP->Const_dt_Integrate_SBDF4(t2, t2+h_over2, h_SBDF4, Y_init_SBDF4, Yf_init[0]);   
    cblas_dcopy (neqn, Y_init[2],1, Yf_init[1] , 1);

    t2=t0+2*h;
    cblas_dcopy (neqn, Y_init[2],1, Y_init_SBDF4[0] , 1);
    Init_coarse_Vectors_RK(order_SBDF4, neqn, t2, h_SBDF4, IVP, Y_init_SBDF4); 
    VSIIE_IVP->Const_dt_Integrate_SBDF4(t2, t2+h_over2, h_SBDF4, Y_init_SBDF4, Yf_init[2]);           
  }
  const int idx=order-1;
  cblas_dcopy (neqn, Y_init[idx],1, Yf_init[idx] , 1);
  
  delete VSIIE_IVP;
  Free_initial_vectors(order, Y_init_SBDF4);

}




//**************************************************************************
// Set the fine initial approximation vectors Yf_init using VSIIE-4 
// using the initial coarse approximation vector Y_init
//**************************************************************************
void Init_fine_Vectors_RK(const int order, const int neqn,  const double t0, 
                  const double h,  IVP_ODE* IVP, double ** Y_init, double ** Yf_init) 
{
  const int RK_order=4;
  double t2;
  // Approximation of the fine intermediate step initial values 
  // using RK4 and very small stepsize
  const double h_RK=1.0e-8; 
  // Set initial condition vectors (Yf_init[...]) for 
  // fine time integration (half of the time step size)
  const double h_over2=h/2;
  if (order==2) {
    t2=t0;
    RK(RK_order, IVP, t2, t2+h_over2, h_RK, Y_init[0], Yf_init[0]);
  }
  else if (order==3) {
    cblas_dcopy (neqn, Y_init[1],1, Yf_init[0] , 1);
    t2=t0+h;
    RK(RK_order, IVP, t2, t2+h_over2, h_RK, Yf_init[0], Yf_init[1]);            
  } else if (order==4) {
    t2=t0+h;
    RK(RK_order, IVP, t2, t2+h_over2, h_RK, Y_init[1], Yf_init[0]);
    cblas_dcopy (neqn, Y_init[2],1, Yf_init[1] , 1);
    t2=t0+2*h;
    RK(RK_order, IVP, t2, t2+h_over2, h_RK, Yf_init[1], Yf_init[2]);
  }
  const int idx=order-1;
  cblas_dcopy (neqn, Y_init[idx],1, Yf_init[idx] , 1);

}



//*************************************************************************
// Approximate numerical solution with the Explicit Runge-Kutta method
void Approximate_solution(IVP_ODE* IVP, const double t0, 
                                     const double tf, double * Y_solut){
//*************************************************************************


  const int neqn=IVP->get_num_ODEs(); 
  
  const int order_SBDF4=4;
  const double h_SBDF4=1.0e-6; // Stepsize for SBDF4 reference solution

  // Init VSIIE solver 
  VSIIE_Solver * VSIIE_IVP=new VSIIE_Solver(order_SBDF4, IVP);
  
   // Declare and  Setup initial condition vectors Y_init
  double *Y_init[order_SBDF4];
  time_ns start,end;
  start = high_resolution_clock::now();
  // Setup coarse initial vectors Y_init
  Setup_initial_vectors(order_SBDF4, neqn, Y_init);

  // Approximation of the intermediate step initial values 
  // using RK4 and very small stepsize
  Init_coarse_Vectors_SBDF4(order_SBDF4, neqn,  t0, h_SBDF4,  IVP, Y_init);
  //IVP->init(Y_init[0]); 
  //Init_coarse_Vectors_RK(order, neqn,  t0, h,  IVP, Y_init);   
  //const double h_RK=1.0e-8;
  //const double RK_order=4;
  //RK(RK_order, IVP, t0, tf, h_RK, Y_init[0], Y_solut);

  // Execution of the VSIIE Solver to approximate accurately the solution
 
  VSIIE_IVP->Const_dt_Integrate_SBDF4(t0, tf, h_SBDF4, Y_init, Y_solut);
  

  end = high_resolution_clock::now();
  double runtime = duration_cast<nanoseconds>(end - start).count() * 1e-9;
  cout << "................................................" << endl;
  cout << "Runtime  VSIIE Reference Solution=  " << runtime << endl;
  cout << "................................................" << endl;
  delete VSIIE_IVP;
  Free_initial_vectors(order_SBDF4, Y_init);
}


//******************************************************
// Obtain the reference solution and save it in a file 
// if it does not exist
// If a data file already exists, read the data from the file 
//****************************************************** 
void Obtain_reference_solution(const std::string & ref_sol_filename, 
           IVP_ODE* IVP, const double t0, const double tf, 
           double* Y1_RK) {

  // Get the size of the solution vector neqn
  const int neqn=IVP->get_num_ODEs(); 

  // Check if the reference solution file already exists
  //ifstream file(ref_sol_filename, ios::binary); 
  ifstream file(ref_sol_filename); 
  if (!file.is_open()) {
	  cout << "Reference solution file does not exist" << endl;
    cout << "... Creating reference solution file:  "<<ref_sol_filename<<endl;
    
    // Compute the reference solution using RK4 with a very small stepsize
    Approximate_solution(IVP, t0, tf, Y1_RK);

    // Save the reference solution in a file
    //ofstream output_file(ref_sol_filename, std::ios::binary); 
    ofstream output_file(ref_sol_filename);
    if (output_file.is_open()) {
      // Write the reference solution to the file
      //output_file.write(reinterpret_cast<char*>(Y1_RK), neqn*sizeof(double)); 
      //output_file.close();  // Close the file

      output_file << std::setprecision(17) << std::scientific;
      for (int i = 0; i < neqn; ++i) {
        output_file << Y1_RK[i] << '\n';
      }
      output_file.close(); // Close the file
    } 
    else {
      cerr << "Error creating file "<< ref_sol_filename << endl;
      return;
    } 
  }
  else {
    file.close();
    cout << "Reference solution file already exists." << endl;
  }

  //ifstream input_file(ref_sol_filename, ios::binary);
  ifstream input_file(ref_sol_filename);   
  if (input_file.is_open()) {
    cout << "Reading reference solution from file: "<< ref_sol_filename << endl;
    //input_file.read(reinterpret_cast<char*>(Y1_RK), neqn*sizeof(double));
    for (int i = 0; i < neqn; ++i) {
      input_file >> Y1_RK[i];
    }
    input_file.close();
  } 
  else {
    cerr << "Error while reading file: "<<ref_sol_filename << endl;
    return;
  }


}

//******************************************************
// Create a folder called "route" 
//****************************************************** 
void create_folder(const std::string & route) {
    if (mkdir(route.c_str(), 0777) == 0) { 
        cout << "The folder '" << route << "' created successfully." << endl;
    } else {
        // The folder already exists or there was an error
        if (errno == EEXIST) {
            cout << "The folder '" << route << "' already exists." << endl;
        } else {
            cerr << "Error creating folder '" << route << "': " << endl;
        }
    }
}



//******************************************************
// Save in a file called "filename" the data corresponding 
// to error, runtime and number of stepsizes for that experiment 
//****************************************************** 
void save_error_runtime_nsteps(const std::string& filename,  const double tol_h, 
      const double error, const double runtime, const int nsteps) {

  ofstream output_file( filename ,std::ios::app);        
  if (output_file.is_open()) {
    // Write a line and its endline
    output_file << tol_h << "     " << error<< "     " << runtime<< "     " << nsteps <<endl; 
    output_file.close(); // Close the file
  } 
  else {
        cerr << "File could not be opened." << endl;
  }

}


//******************************************************
// Run experiments for a particular IVP problem and a particular 
// VSIIE/IIE solver  obtaining information
// tol argument only makes sense when is_adaptive=true
//****************************************************** 
void perform_experiments(const int order,  IVP_ODE* IVP,  VSIIE_Solver * VSIIE_IVP, 
    const double t0, const double tf, const double h, const double * Y1_solut, 
    const double tol , const int n_repetitions, const bool is_adaptive, 
    const double alpha, const double eta_min, const double eta_max) {  

  string solver;
  solver=is_adaptive ? "VSIIE":"IIE";

  // Initial and final time instants;
  time_ns start,end;
   
  const int neqn = IVP->get_num_ODEs();
  // Declare vector Y1 to store the numerical approximation to the solution	 
  double * Y1 = new double[neqn];

  // n_steps: Number of valid stepsizes  
  int n_steps=0;
  // n_isteps: Total number of stepsizes (only for VSIIE)
  int n_isteps;

  // n_iters: Average Number of Newton iterations
  double  n_iters;
  // error in the numerical solution and runtime    
  double error=0.0,runtime=0.0;
    
  cout<< "_____________________________________________________________________"<<endl;
  cout<<"_________________ "<<solver<<"-"<<order<<" _________________"<<endl;

  // Declare and  Setup initial condition vectors Y_init and Yf_init
  double *Y_init[order];
  double *Yf_init[order+2];
  // Setup coarse initial vectors Y_init
  Setup_initial_vectors(order, neqn, Y_init);
  // Setup fine initial vectors Yf_init 
  Setup_initial_vectors(order+2, neqn, Yf_init);
  // Approximation of the intermediate step initial values 
  // using RK4 and very small stepsize
  //IVP->init(Y_init[0]); 
  //Init_coarse_Vectors_RK(order, neqn,  t0, h,  IVP, Y_init);    
  //if (is_adaptive){Init_fine_Vectors_RK(order, neqn,  t0, h,  IVP, Y_init, Yf_init);}

  Init_coarse_Vectors_SBDF4(order, neqn,  t0, h,  IVP, Y_init);   
  if (is_adaptive){Init_fine_Vectors_SBDF4(order, neqn,  t0, h,  IVP, Y_init, Yf_init);} 
  

  //******************************************************************************
  // Loop to obtain the minimum of  n_repetitions tests
  for (int j = 0; j < n_repetitions; j++){
  //******************************************************************************  
    double & ref_n_iters=n_iters;
    start = high_resolution_clock::now();

    if (is_adaptive){
      VSIIE_IVP->Adaptive_dt_Integrate(t0, tf, h, Y_init, Yf_init, Y1,
        tol, &n_steps, &n_isteps, ref_n_iters, alpha, eta_min, eta_max);
    }
    else { 
      VSIIE_IVP->Const_dt_Integrate(t0, tf, h, Y_init, Y1);
      n_steps=ceil((tf-t0)/h); 
    }
    end = high_resolution_clock::now();

    // Compute the difference between the numerical solution and the reference solution
    // Only for th 1st iteration
    if (j==0){
      cblas_daxpy(neqn, -1.0, Y1_solut, 1, Y1, 1); 
      error = cblas_dnrm2(neqn, Y1, 1)/pow(neqn, 0.5);
    }
    const double current_runtime=duration_cast<nanoseconds>(end - start).count() * 1e-9;
    runtime = (j==0) ?current_runtime: min(runtime,current_runtime);  
  }
	 
  cout <<endl<< IVP->get_name()<<" with Neqn = " <<neqn<<"......"<<solver<<"-"<<order<<endl;

  // Output results           

  if (is_adaptive){
    cout << endl<<"***** H0 = " << setw(4) <<h<<"******"; 
    cout <<".. N_STEPS = " << n_steps<< "...";
    cout <<".. N_TOT_STEP = " << n_isteps;  
    cout <<".. Av. Newton iters = " << n_iters<<endl; 
    cout <<"*****  TOL=" <<setw(4) <<tol<<"******";          
    cout <<"   ....ERROR= "  << setw(3)<< error;
    cout <<"   ....RUNTIME= "<<  setw(3)<<runtime;
    cout<<endl;
  }
  else {
    cout << endl<<"***** FIXED H=" <<setw(4) <<h<<"******"; 
    cout <<"   ....ERROR= "  << setw(3)<< error;
    cout <<"   ....RUNTIME= "  <<  setw(3)<<runtime;
    cout <<"   ....NSTEPS= "<< setw(3)<<n_steps<<endl;
    cout<<endl;
  }  

  cout << setw(3)<< error << "      " <<  setw(3)
                       <<runtime<< "      " <<  setw(3)<<n_steps<<endl;                                
  cout<<"____________________________________________________________________________________"<<endl;

  const double tol_h= is_adaptive ? tol : h;
  save_error_runtime_nsteps(solver+"-"+to_string(order)+"-"+IVP->get_name()+ "-" + to_string(neqn) + 
                                                                     ".txt",  tol_h, error, runtime, n_steps);
    
  // Free vectors in dynamic memory	
  Free_initial_vectors(order, Y_init);
  Free_initial_vectors(order+2, Yf_init);
  delete[] Y1;
}


//***************************************************
//************  MAIN PROGRAM  ***********************
int main(int argc, char** argv)
//***************************************************
{
  // Initializa LIS environment
  lis_initialize(&argc, &argv);
       
  // Number of arguments   
  const int num_args=argc;
  // Declaration of variables to store the program parameters
  int problem_id, order, Neqn;
  double tf, h, tol; 
  double alpha, eta_min, eta_max;
  int n_repetitions;
  bool error=false;
  string const_var;

  if (num_args<9) {
    error=true; cerr << "Invalid number of arguments" << endl<<endl;
  }
  else {
    // Stores the value of the program parameters in the corresponding variables
    const_var=argv[1]; // Variable Stepsize("var") or constant ("const")
    if (const_var!="var" && const_var!="const") {
      cerr << "Invalid const_var parameter. Please select 'var' or 'const'." << endl;
      error=true;
    }
    else { 
      if ((const_var=="const" && num_args!=9) || 
          (const_var=="var" && (num_args<9 || num_args>12)) ) {
        error=true; cerr <<endl<< "Invalid number of arguments!!" << endl;
      }  
      if (const_var=="var"){
        alpha=0.8; // Default value of alpha
        eta_min=0.5; // Default value of eta_min
        eta_max=4; // Default value of eta_max
      }
    }
    if (!error){
      problem_id=atoi(argv[2]); // Problem identifier
      if (problem_id>4 || problem_id<0){
        cerr << "Invalid problem_id parameter. Please select problem_id between 0 and 4." << endl;
        error=true;
      }
      order=atoi(argv[3]); // Convergence order to be tested
      if (order<1 || order>4){
        cerr << "Invalid order parameter. Please select order between 1 and 4." << endl;
        error=true;
      }
      Neqn = atoi(argv[4]); // Number of ODEs
      tf = atof(argv[5]);   // Final time
      h = atof(argv[6]); // Initial/Constant stepsize for VSIIE/IIE method
      tol = atof(argv[7]); // Tolerance
      n_repetitions=atof(argv[8]); // Number of repetitions
      if (const_var=="var") {
        if(num_args>=10){
          alpha = atof(argv[9]); // Alpha parameter for the adaptive time-stepping strategy
          if (alpha<0 || alpha>1){
            cerr << "Invalid alpha parameter. Please select alpha between 0 and 1." << endl;
            error=true;
          }
          if(num_args>=11) eta_min = atof(argv[10]); // eta_min parameter for the adaptive time-stepping strategy
          if(num_args>=12) eta_max = atof(argv[11]); // eta_max parameter for the adaptive time-stepping strategy        
        }  
      }    
    } 
  }

  // Check the number of parameters
  if (error) {
    // Tell the user how to run the program
    cerr << endl<<"Usage: " << argv[0] << 
      "<const_var> <problem id.> <order> <Neqn> <tf>  <stepsize> <tol> <num_reps> <alpha> <eta_min> <eta_max>" 
                      << endl<<endl;
    cerr << "<const_var>: " <<endl;
    cerr << "             const: constant stepsize (IIE method is used)" <<endl;
    cerr << "             var: variable stepsize (VSIIE method is used)" << endl;
    cerr << endl;

    cerr << "<problem.id>= " <<endl;
    cerr << "              0: 1D Advection-Diffusion.  |   1: Stiff Brusselator." <<endl;
    cerr << "              2: FKPP Comb.               |   3: Ignition Comb.     |     4: Fisher Comb." << endl;
    cerr << endl;

    cerr << "<Neqn>: Number of ODEs " <<endl;
    cerr<<endl;

    cerr << "<stepsize>: It denotes the constant stepsize if const_var=const and the initial stepsize if const_var=var" << endl;
    cerr << endl;

    cerr << "<tol>: Error tolerance which it is only used for VSIIE (const_var=var)" <<endl;
    cerr << endl;

    cerr << "<num_reps>:  Number of repetitions of the experiment (to obtain the minimum runtime)" <<endl;
    cerr << endl;

    cerr << "<alpha> <eta_min> <eta_max>(Optional):  Parameters alpha, eta_min and eta_max"<<endl; 
    cerr<<"of the adaptive time-stepping strategy (default: alpha=0.8, eta_min=0.5, eta_max=4)" <<endl;
    cerr << endl;

    return 1;
  }




  // Declaration of the C++ objects representing 
  // the IVP-ODE problems to be solved 
  IVP_ODE* IVP;
  IVP_ODE_advdiff1d   IVP_advdiff1d(Neqn);
  IVP_ODE_stiff_brusselator IVP_stiff_brusselator(Neqn/3);
  IVP_ODE_combustion IVP_combustion(Neqn, problem_id-1);
  if (problem_id==0) 
    IVP= &IVP_advdiff1d;
  else if (problem_id==1)
    IVP= &IVP_stiff_brusselator;
  else if (problem_id>=2 && problem_id<=4)
    IVP= &IVP_combustion;
  else {
    cerr << "Invalid problem id. Please select a valid problem id." << endl;
    return 1;
  }



  //Time-step conditions
  const double t0 = 0.0;

  // Get the number of ODEs in the test problem
  const int neqn = IVP->get_num_ODEs();
  cout << "................................................" << endl;
  cout<<endl<<    IVP->get_name()<<" IVP with Neqn = "<<neqn<<endl;
  cout << "................................................" << endl;

  cout.precision(10);

  //Declare vector Y1_RK  which stores the numerical approximation to the solution	 
  double* Y1_RK =new double[neqn]; 

  const string data_folder="data/";
  create_folder(data_folder);

  // Generate a string to store the final time tf using only 3 decimals
  stringstream ss;
  ss << fixed << setprecision(3) << tf;
  string s_tf = ss.str();
  // Generate the name of the file which  stores the reference solution data 
  const string ref_sol_filename=data_folder+"ref_sol-" + 
                                  IVP->get_name()+
                                 "-"+to_string(Neqn)+"-"+s_tf+".txt";
  //const string ref_sol_filename = data_folder + "r.txt";
  // Obtain the reference solution from a file or produce 
  // it and save it in a file if the file does not exist
  Obtain_reference_solution(ref_sol_filename, IVP, t0, tf, Y1_RK);

  // Init VSIIE solver 
  VSIIE_Solver * VSIIE_IVP=new VSIIE_Solver(order, IVP);

    
  // Run the experiments fot a particular VSIIE-order solver and IVP
  perform_experiments(order, IVP, VSIIE_IVP, t0, tf, 
                     h,  Y1_RK, tol, n_repetitions, const_var=="var", 
                     alpha, eta_min, eta_max);
       
  delete VSIIE_IVP;                     
  delete[] Y1_RK;
 
  //Finalize LIS environment 
  lis_finalize();
}

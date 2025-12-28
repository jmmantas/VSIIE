/*
This C++ driver program run several experiments associated with a given IVP model using a particular VSIIE method  
and generate the numerical errors in order to check the theoretical convergence order. 

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

#include  <omp.h>
#include "IVP_ODE_advdiff1d.h"
#include "IVP_ODE_stiff_brusselator.h"
#include "IVP_ODE_combustion.h"

#include "VSIIE_Solver.h"
#include "utils.h"
using namespace std;
using namespace std::chrono;

typedef std::chrono::time_point<high_resolution_clock, nanoseconds> time_ns;


//******************************************************
//Set up Initial condition Vectors for the fine and 
// coarse approximation: Y_init[0]=Y[0] in t=t0,
// Y_init[1] in t=t0+h, ... Y_init[order-1] in t=t0+(order-1)*h
// Yf_init[0] in t=t0, Yf_init[1] in t=t0+h/2, Yf_init[2] in t+h, ...
//******************************************************
void Setup_initial_vectors(const int order, const int neqn,  
                           double ** Y_init, double ** Yf_init){
    for (int i = 0; i < order; i++) {
        Y_init[i] = new double[neqn];
    }
    for (int i = 0; i < (order+2); i++) {
        Yf_init[i] = new double[neqn];
    }
}






//**************************************************************************
// Initialize  intermediate vectors (Y_init and Yf_init) using 4th order Runge-Kutta 
// to start the multistep integration
//**************************************************************************
void Init_Vectors(const int order, const int neqn,  const double t0, 
                  const double h,  IVP_ODE* IVP, double ** Y_init, double ** Yf_init) 
{
    const int RK_order=4, idx=order-1;
    IVP->init(Y_init[0]); //First  step initial value

    // Approximation of the intermediate step initial values 
    // using RK4 and very small stepsize
	const double h_RK=1.0e-8;
    // Initialize initial condition vectors (Y_init[...]) for coarse time integration
    double t2=t0;           
    for (int i = 1; i < order; i++){  
        RK(RK_order, IVP, t2, t2+h, h_RK, Y_init[i-1], Y_init[i]);
        t2+=h;
    } 
    // Set up initial condition vectors (Yf_init[...]) for 
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
    cblas_dcopy (neqn, Y_init[idx],1, Yf_init[idx] , 1);
}





//**************************************************************************
// Free auxiliary vectors in dynamic memory:
// Y_init(order*neqn) and Yf_init((order*neqn+2) 	
//**************************************************************************
void Free_initial_vectors(const int order, double **Y_init, double ** Yf_init)
{
    for (int i = 0; i < order; i++) 
        delete[] Y_init[i];
    for (int i = 0; i < (order+2); i++) 
        delete[] Yf_init[i];  
   
}




//***************************************************
// MAIN PROGRAM
int main(int argc, char** argv)
//***************************************************
{
      // Initializa LIS environment
    lis_initialize(&argc, &argv);
    
    // Order of the Explicit Runge-Kutta Solver 
    // used to obtain the reference solution
    const int RK_order=4; 
     
    
    // Number of arguments   
    const int num_args=8;
    // Check the number of parameters
    if (argc < num_args) {
        // Tell the user how to run the program
        cerr << "Usage: " << argv[0] 
        << " <problem id.> <conv. order> <Neqn(N)>  <tf(final time)> <RK stepsize> "
        << "<init. VSIIE stepsize>  <Num. experiments> " << endl<<endl;

        cerr << "<problem.id>= " <<endl;
        cerr << "             0: 1D Advection-Diffusion" <<endl;
        cerr << "             1: Stiff Brusselator" << endl;
        cerr << "             2: FKPP Combustion" << endl;
        cerr << "             3: Ignition Combustion" << endl;
        cerr << "             4: Fisher Combustion" << endl;
        cerr << endl;
        return 1;
    }

    // Stores the value of the program parameters in the corresponding variables
    int problem_id=atoi(argv[1]); // Problem identifier
    const unsigned order = atoi(argv[2]);
    const int Num_points = atoi(argv[3]);
    const double tf = atof(argv[4]);
    //Time integration step for RK4 and VSIIE methods
    double h_RK = atof(argv[5]); //Time integration stepsize for RK4 method
    double h_init = atof(argv[6]); //Initial Time integration stepsize for VSIIE method
    const int num_tests=atoi(argv[7]); // Number of experiments

    
    
    // Declaration of the C++ objects representing 
    // the IVP-ODE problems to be solved 
    IVP_ODE* IVP;
    if (problem_id==0)
    {
      IVP_ODE_advdiff1d   IVP_advdiff1d(Num_points);
      IVP= &IVP_advdiff1d;
    }
    else if (problem_id==1)
    {
      IVP_ODE_stiff_brusselator IVP_stiff_brusselator(Num_points);
      IVP= &IVP_stiff_brusselator;
    }
    else if (problem_id>=2 && problem_id<=4)
    {
      IVP_ODE_combustion IVP_combustion(Num_points,problem_id-1);
      IVP= &IVP_combustion;
    }
    else
    {
      cerr << "Invalid problem id. Please select a valid problem id." << endl;
      return 1;
    }


    //Time-step conditions
    const double t0 = 0.0;
    
    // Get the number of ODEs in the test problem
    const int neqn = IVP->get_num_ODEs();
    cout <<"NEQN="<<neqn<<endl;
    cout << "................................................" << endl;

    cout<<endl<<    IVP->get_name()<<" IVP with Neqn = "<<neqn<<endl;
    cout << "................................................" << endl;

    //Set up Initial condition Vectors: Y[0] in t=t0,
    // Y_init[1] in t=t0+h, ... Y_init[order-1] in t=t0+(order-1)*h
    double *Y_init[order];
    for (unsigned i = 0; i < order; i++) {
            Y_init[i] = new double[neqn];
    }
    //Set up Initial condition Vectors for the fine approximation: 
    // Yf_init[0] in t=t0, Yf_init[1] in t=t0+h/2, Yf_init[2] in t+h, ...
    // Setup initial condition vectors and final vectors, Yf_init, for the fine time steps
    double *Yf_init[order+2];
    for (unsigned int i = 0; i < (order+2); i++) {
        Yf_init[i] = new double[neqn];
    }
        
    //Declare vector Y1_variable to store the numerical approximation 
    // to the solution	 
    double* Y1_variable = new double[neqn];

    //Declare vector Y1 to store the numerical approximation 
    // to the solution with RK solver	 
    double* Y1_RK = new double[neqn];

    cout.precision(10);

    //**********************************************
    double runtime_RK;
    cout << "................................................" << endl;
    // Numerical solution with the Explicit Runge-Kutta method 
    cout << "<<<< Explicit Runge-Kutta order "<<RK_order
         <<" ...   Time Stepsize h=" << h_RK << " >>>> " << endl;
    cout << "<<<< N= " << neqn << "...t0=" << t0 << "...tf= " << tf << endl;

    //First  step initial value    
    IVP->init(Y_init[0]);

    // Execution of the RK Solver to approximate accurately the solution
    time_ns start,end;
    start = high_resolution_clock::now();
    
    RK(RK_order, IVP, t0, tf, h_RK, Y_init[0], Y1_RK);
 
    end = high_resolution_clock::now();
    runtime_RK = duration_cast<nanoseconds>(end - start).count() * 1e-9;
    cout << "................................................" << endl;
    cout << "Runtime RK=  " << runtime_RK << endl;
    cout << "................................................" << endl;
       
     double h0=h_init;
    // Numerical solution with the VSIIE-order Method
    cout << "<<<< VSIIE" << order << "...   Initial Time Step h=" << h0 << " >>>>" << endl;

    time_ns start_VSIIE, end_VSIIE;
    // Approximation of the intermediate step initial values 
    // using RK4 and very small stepsize
    Init_Vectors(order,neqn,  t0, h0,  IVP, Y_init, Yf_init);   
    
    // Init VSIIE solver 
    VSIIE_Solver * VSIIE_IVP;
    VSIIE_IVP=new VSIIE_Solver(order, IVP);
  
    // error in the numerical solution for each test    
    double error[num_tests];
    // runtime for each test   
    double runtime[num_tests];
    // convergence order for each pais of tests
    double conv_order[num_tests-1];
   

    //******************************************************************************
    // Perform experiments
    //******************************************************************************    
     for (int i = 0; i < num_tests; i++){   
    //******************************************************************************
        Init_Vectors(order,neqn,  t0, h0,  IVP, Y_init, Yf_init);   
        start = high_resolution_clock::now();
        VSIIE_IVP->Const_dt_Integrate(t0, tf, h0, Y_init, Y1_variable);
        end = high_resolution_clock::now();
        runtime[i] = duration_cast<nanoseconds>(end - start).count() * 1e-9;
        // Compute the difference between the numerical solution and the reference solution
        cblas_daxpy(neqn, -1.0, Y1_RK, 1, Y1_variable, 1);
        error[i] = cblas_dnrm2(neqn, Y1_variable, 1)/pow(neqn, 0.5);
        if (i>0) conv_order[i-1]=error[i-1]/error[i];
        cout << "CONSTANT H=  "<<h0<<".. RUNTIME= "<< runtime[i] 
                                        <<"... ERROR = " << error[i]<< "...........";
        if (i>0) cout << "ORDER = ERROR["<< i-1<<"]/ERROR["<<i<<"]= "<<conv_order[i-1];
        cout<<endl<<"____________________________________________________________________________________"<<endl<<endl;
        h0=h0/2.0;
    }

    cout << endl << endl;
    cout << IVP->get_name()<<" with Neqn = "
             <<neqn<<"......VSIIE"<<order<<" .....  ";
    for (int j = 0; j < num_args; j++){
        cout<<argv[j]<<" ";
    }     
    cout<<endl<<endl;
    
    for (int i = 0; i < num_tests; i++){ 
        cout << endl<<"***** H=" <<setw(4) <<h_init/pow(2.0,i)<<"******"; 
        cout <<"   ....ERROR= "  << setw(4)<< error[i];
        cout <<"   ....TIME= "<<  setw(4)<<runtime[i];
        if (i>0) cout <<"   ....EMPIRICAL CONV. ORDER= "  << setw(4)<< conv_order[i-1];  

    }
    cout<<endl;



    // Free vectors in dynamic memory	
    Free_initial_vectors(order, Y_init, Yf_init);
    delete[] Y1_RK;
    delete[] Y1_variable;
    
    delete VSIIE_IVP;

    //Finalize LIS environment 
    lis_finalize();

    cout << "Done!" << endl;
}
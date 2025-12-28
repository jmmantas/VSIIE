""""
This python module defines certain  functions which are used by other python scripts in order
to  generate numerical simulations and image files.

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
""""


import matplotlib.pyplot as plt
import numpy as np
import argparse
import subprocess
import sys
import os



#************************************************************************************
# Function to generate a plot from solver output files, 
#************************************************************************************
def view (model_name, neqn, VSIIE_min_order, VSIIE_max_order, IIE_min_order, IIE_max_order):
    x_vals = [[],[],[],[],[],[],[],[]]
    error_vals = [[],[],[],[],[],[],[],[]]
    n_steps_vals = [[],[],[],[],[],[],[],[]]
    filename=[]
    label_name=[]
    included_solver=[]
    for i in range(8):
        if (i<4): 
            order = i+1
            solver= "VSIIE" + "-" +str(order)
            included_solver.append ( (order >= VSIIE_min_order) and (order <= VSIIE_max_order) )
        else:
            order= i-3 
            solver  = "IIE"  + "-" +str(order)
            included_solver.append( (order >= IIE_min_order) and (order <= IIE_max_order) )
        label_name.append(solver)    
        filename.append(solver+"-"+model_name+"-"+str(neqn)+".txt")



    for i in range(8):
        if included_solver[i]:
            print(f"Processing File : {i} -> {filename[i]}")
            with open(filename[i], 'r') as file:
                for linea in file:
                    datos = linea.split()
                    if len(datos) == 4:
                        try:
                            x1 = float(datos[1])
                            e1 = float(datos[2])
                            n1= int(datos[3])
                            x_vals[i].append(x1)
                            error_vals[i].append(e1)
                            n_steps_vals[i].append(n1)
                        except ValueError:
                            print(f"line ignored due to incorrect format: {linea.strip()}")

            if not x_vals[i]:
                print("No valid data was found in the file.")
                return
    
      
    colors=['black','red','green','blue','black','red','green','b']
    markers=['o','s', '^','x',   '*','s','D','p','2', '>']
    plt.figure(figsize=(8, 6))

    
    for i in range(0,4):
        if included_solver[i]: 
            plt.plot(x_vals[i], error_vals[i], marker=markers[i], linestyle='-', color=colors[i], label=label_name[i])
    for i in range(4,8): 
        if included_solver[i]:
            plt.plot(x_vals[i], error_vals[i], marker=markers[i], linestyle='--', color=colors[i], label=label_name[i])

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Error')
    plt.ylabel('runtime')
    plt.title("Work-Precision. "+model_name+". N="+str(neqn))
    plt.legend()
    plt.grid(True)
    plot_file="plot_error" + ".png"
    plt.savefig(plot_file)
    
    
    plt.figure()
    for i in range(0,4): 
        if included_solver[i]:
            plt.plot(x_vals[i], n_steps_vals[i], marker=markers[i], linestyle='-', color=colors[i], label=label_name[i])
    for i in range(4,8): 
        if included_solver[i]:
            plt.plot(x_vals[i], n_steps_vals[i], marker=markers[i], linestyle='--', color=colors[i], label=label_name[i])
   
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Error')
    plt.ylabel('N_Steps')
    plt.title("Number of Time Steps. "+model_name+". N="+str(neqn))
    plt.legend()
    plt.grid(True)
    plot_file="plot_n_steps" + ".png"
    plt.savefig(plot_file)    
#************************************************************************************
    



#************************************************************************************
# Function to run an IIE solver with given parameters
#************************************************************************************
def run_IIE_solver(problem, order,  neqn, tf, stepsize, tol, n_repetitions):        
    command = [ '../IVP_Solver', "const",  str(problem), str(order), str(neqn),
        str(tf), str(stepsize), str(tol), str(n_repetitions)]
    try:
        print(f"Executing command: {' '.join(command)}")
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        print("Solver Output:")
        print(result.stdout)
    except FileNotFoundError:
        print("Error: executable was not found.")
    except subprocess.CalledProcessError as e:
        print(f"Error executing the program:")
        print(f"Return code: {e.returncode}")
        print(f"Error exit: {e.stderr}")
    except Exception as e:
        print(f"Unexpected error: {e}")
#************************************************************************************


#************************************************************************************
# Function to run a VSIIE solver with given parameters
#************************************************************************************
def run_VSIIE_solver(problem, order,  neqn, tf, stepsize, tol, n_repetitions, 
                     alpha, eta_min, eta_max):        
    command = [ '../IVP_Solver', "var",  str(problem), str(order), str(neqn),
        str(tf), str(stepsize), str(tol), str(n_repetitions), str(alpha), str(eta_min), str(eta_max)]
    try:
        print(f"Executing command: {' '.join(command)}")
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        print("Solver Output:")
        print(result.stdout)
    except FileNotFoundError:
        print("Error: executable was not found.")
    except subprocess.CalledProcessError as e:
        print(f"Error executing the program:")
        print(f"Return code: {e.returncode}")
        print(f"Error exit: {e.stderr}")
    except Exception as e:
        print(f"Unexpected error: {e}")
#************************************************************************************





#************************************************************************************
# Function to run an IIE solver with given parameters, including a list of stepsizes
#************************************************************************************
def run_IIE_instances(problem, order, neqn, tf, stepsize_array, tol, n_repetitions):        
    os_command=f"rm IIE-{order}*.txt"
    os.system(os_command)
    print(f"........................ Executing IIE-{order} Solvers \n")
    for stepsize in stepsize_array:
        run_IIE_solver(problem, order, neqn, tf, stepsize, tol, n_repetitions)
    print(f"........................ Finished IIE-{order} Solvers \n")
#************************************************************************************



#************************************************************************************
# Function to run an VSIIE solver with given parameters, 
# including a list of stepsizes and tolerances
#************************************************************************************
def run_VSIIE_instances(problem, order, neqn, tf, stepsize_array, tol_array, n_repetitions, 
                        alpha_array, eta_min, eta_max):        
    os_command=f"rm VSIIE-{order}*.txt"
    os.system(os_command)
    print(f"........................ Executing VSIIE-{order} Solvers \n")
    i=0
    for tol in tol_array:
        run_VSIIE_solver(problem, order, neqn, tf, stepsize_array[i], tol, n_repetitions, 
                         alpha_array[i], eta_min, eta_max)
        i+=1
    print(f"........................ Finished VSIIE-{order} Solvers \n")


#************************************************************************************
# Function to show the specific parameters, 
# including a list of stepsizes, tolerances and alphas
#************************************************************************************
def show_parameters (Problem, IVP_name, Neqn, Tf, stepsize_array, VSIIE_tol_array, alpha_array):
    print(f"\n ******************************************************************************** ")
    print(f"\n Problem id.  {Problem} : {IVP_name} with Neqn={Neqn}. Tf={Tf} ")
    # Fill VSIIE stepsize arrays based on IIE stepsize arrays
    for order in range(1,5):
        print(f"\n Parameters for Order {order} :")
        for i in range(len(stepsize_array[order-1])):
            step=stepsize_array[order-1][i]
            tol=VSIIE_tol_array[order-1][i]
            alpha=alpha_array[order-1][i]
            print(f" Order {order}: VSIIE Stepsize: {step}   Tol: {tol}  Alpha: {alpha}   ")   
    print(f"\n ******************************************************************************** ")


#************************************************************************************
# Function to run the corresponding experiments,
#************************************************************************************
def do_experiments(Problem, IVP_name, Neqn, Tf, FOLDER, stepsize_array, 
                   VSIIE_tol_array, alpha_array):
    correct_num_args=5
    print(f"Detected {len(sys.argv) - 1} arguments received.")
    if len(sys.argv) != correct_num_args + 1:
        print(f"Error: Expected {correct_num_args} arguments, but {len(sys.argv) - 1} received.")
        print(f"Usage: {sys.argv[0]}  VSIIE_min_order   VSIIE_max_order  IIE_min_order   IIE_max_order  N_Repetitions")
        sys.exit(1) # Exit with error code

    parser=argparse.ArgumentParser()
    parser.add_argument("VSIIE_min_order",type=int)
    parser.add_argument("VSIIE_max_order",type=int)
    parser.add_argument("IIE_min_order",type=int)
    parser.add_argument("IIE_max_order",type=int)
    parser.add_argument("N_Repetitions",type=int)
    args=parser.parse_args()
    print(f" CALL: {sys.argv[0]}  {args.VSIIE_min_order}  {args.VSIIE_max_order}  {args.IIE_min_order}  {args.IIE_max_order}  {args.N_Repetitions}")

    os_command=f"cd ..; make clean;make;cd {FOLDER}"
    os.system(os_command)


    #**********************************************************
    print (f"........................ Executing IIE Solvers\n")
    #**********************************************************
    tol=0.01
    for order in range(1, 5):
        if (args.IIE_min_order <= order) and (args.IIE_max_order >= order):
            run_IIE_instances(Problem, order, Neqn, Tf, stepsize_array[order-1], tol, args.N_Repetitions)
    #**********************************************************
    print (f"........................ Executing VSIIE Solvers\n")
    #**********************************************************
    eta_min=0.2
    eta_max=4.0
    range_factor=0.2
    for order in range(1, 5):
        if (args.VSIIE_min_order     <= order) and (args.VSIIE_max_order >= order):
            run_VSIIE_instances(Problem, order, Neqn, Tf, stepsize_array[order-1], VSIIE_tol_array[order-1], 
                            args.N_Repetitions, alpha_array[order-1], eta_min, eta_max)

    #****************************************************
    print(f"........................ Plotting Results \n")
    #****************************************************
      
    view(IVP_name, Neqn, args.VSIIE_min_order, args.VSIIE_max_order, args.IIE_min_order, args.IIE_max_order)



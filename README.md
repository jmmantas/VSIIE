# VSIIE:  Variable Stepsize Implicit-Implicit-Explicit (IIE)  Solvers in C++
**VSIIE** software package includes the C++ implementation of a family of Variable-Stepsize 3-additive Implicit-Implicit-Explicit methods to solve stiff IVP-ODEs derived from  Advection-Diffusion-Reaction Models. 
The Implicit-Implicit-Explicit (IIE) solvers  are a family of high-order time integration methods designed
to efficiently and accurately solve stiff systems of ordinary differential equations (ODEs) that often
arise from the spatial discretization of relevant models based on partial differential equations (PDEs)
such as diffusion-reaction-advection (DRA) models, which can be divided into three separate terms.
Here, the IIE solvers are extended to allow for variable stepsize, resulting in the VSIIE solvers.
In addition, an implementation of Variable Stepsize Semi-implicit Backward Differentiation Formula (VSSBDF) solvers 
is also included for the purposes of experimental comparison.


The IVP-ODE  is given by: 

dy/dt = F(t,y(t)) = f1(t,y(t)) + f2(t,y(t))  + f3(t,y(t)),   y(t_0)=y_0,

being f1(t,y(t)) and f2(t,y(t)) the stiff terms in F(t,y(t))  (usually f1(t,y(t)) represents the diffusion term and f2(t,y(t)) is the reaction term) 
and f3(t,y(t)) a nonstiff term (usually the advection term).


## Installation

Before compiling and link the software, you need to install the libraries LIS (https://www.ssisc.org/lis/index.en.html), and openblas (https://www.openblas.net/).

It is necessary to edit "Makefile" to indicate the path where the LIS library is installed by setting the variable LISROOT which is set by default to "usr/local".

Once the Makefile has been modified, you have to compile and link by issuing:

make

This will generate the executable file "IVP_Solver". 

## Usage
A User's manual, which describes how to use the solvers and several python scripts to perform automatically simulations, is
available in the file user-manual.pdf of the folder "doc".

## License

The code in this repository is released under GPLv3, 2025 Jose Miguel Mantas Ruiz (jmmantas@ugr.es) and Raed Ali Mara'Beh (raedmaraabeh@gmail.com).

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


## References

- Raed Ali Mara’Beh, J.M. Mantas, P. González, and Raymond J. Spiteri. Performance comparison
of variable-stepsize IMEX SBDF methods on Advection-Diffusion-Reaction models. Computers &
Mathematics with Applications, 190:41–56, 2025.

- Raed Ali Mara’Beh, J.M. Mantas, and Raymond J. Spiteri. Variable-stepsize 3-additive linear multistep
methods on Diffusion-Reaction-Advection models. Submitted to IEEE Access, 2026.

- Raed Ali Mara’Beh, Raymond J. Spiteri, P. González, and José M. Mantas. 3-additive linear multi-step
methods for Diffusion-Reaction-Advection models. Applied Numerical Mathematics, 183:15–38, 2023.

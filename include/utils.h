/*
These are auxiliary functions for the driver programs which use the VSIIE solvers.

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


#include "lis.h"
#include "IVP_ODE.h"
using namespace std;
using namespace std::chrono;
typedef std::chrono::time_point<high_resolution_clock, nanoseconds> time_ns;

//**************************************************************************
//AUXILIARY  FUNCTIONS
//**************************************************************************


//**************************************************************
// Check important errors between sparse matrices A1 and A2
void compare_matrix_csr(LIS_MATRIX A1, LIS_MATRIX A2)
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
//**************************************************************






//***************************************************
//1st, 2nd and 3rd order Runge-Kutta Time Integrator
//***************************************************
void RK(const int nstages, IVP_ODE* IVP, const double t0,
    const double tf, const double h0, const double* Y0, double* Y1)
//***************************************************
{
    const int neqn = IVP->get_num_ODEs();

    double* Y0_tmp = new double[neqn];
    double* DY1 = new double[neqn];
    double* DY2 = new double[neqn];
    double* DY3 = new double[neqn];
    double* DY4 = new double[neqn];

    double t = t0;
    bool end = false;
    double h = h0;
    int N = 0;
    cblas_dcopy(neqn, Y0, 1, Y0_tmp, 1);
    while (!end) {
        if ((tf-t) < 0.0) {cout << "Surpassing final time!!!!" << endl; return; }
        else if ((tf - t) <= h) { h = tf - t; end = true; }
        // DY1=F(t,Y0)  
        IVP->feval(t, Y0_tmp, DY1);
        switch (nstages)
        {
        case 1:
            cblas_daxpy(neqn, h, DY1, 1, Y0_tmp, 1);
            break;
        case 2:
            cblas_dscal(neqn, 0.5 * h, DY1, 1);
            cblas_daxpy(neqn, 1, Y0_tmp, 1, DY1, 1);
            IVP->feval(t + 0.5 * h, DY1, Y1);
            cblas_daxpy(neqn, h, Y1, 1, Y0_tmp, 1);
            break;
        case 3:
            cblas_dcopy(neqn, Y0_tmp, 1, Y1, 1);
            cblas_daxpy(neqn, 0.5 * h, DY1, 1, Y1, 1);
            IVP->feval(t + 0.5 * h, Y1, DY2);
 
            cblas_dcopy(neqn, Y0_tmp, 1, Y1, 1);
            cblas_daxpy(neqn, -h, DY1, 1, Y1, 1);
            cblas_daxpy(neqn, 2 * h, DY2, 1, Y1, 1);
            IVP->feval(t + h, Y1, DY3);

            cblas_daxpy(neqn, h / 6    , DY1, 1, Y0_tmp, 1);
            cblas_daxpy(neqn, 4 * h / 6, DY2, 1, Y0_tmp, 1);
            cblas_daxpy(neqn, h / 6    , DY3, 1, Y0_tmp, 1);
            break;
        // RK4 
        case 4:
            cblas_dcopy(neqn, Y0_tmp, 1, Y1, 1);
            cblas_daxpy(neqn, 0.5*h, DY1, 1, Y1, 1);
            IVP->feval(t + 0.5 * h, Y1, DY2);

            cblas_dcopy(neqn, Y0_tmp, 1, Y1, 1);
            cblas_daxpy(neqn, 0.5*h, DY2, 1, Y1, 1);
            IVP->feval(t + 0.5 * h, Y1, DY3);

            cblas_dcopy(neqn, Y0_tmp, 1, Y1, 1);
            cblas_daxpy(neqn, h, DY3, 1, Y1, 1);
            IVP->feval(t + h, Y1, DY4);

            cblas_daxpy(neqn, h / 6, DY1, 1, Y0_tmp, 1);
            cblas_daxpy(neqn, h / 3, DY2, 1, Y0_tmp, 1);
            cblas_daxpy(neqn, h / 3, DY3, 1, Y0_tmp, 1);
            cblas_daxpy(neqn, h / 6, DY4, 1, Y0_tmp, 1);
            break;
        }

        N++;
        t = t + h;
    }
    cblas_dcopy(neqn, Y0_tmp, 1, Y1, 1);
    delete[] DY1;
    delete[] DY2;
    delete[] DY3;
    delete[] DY4;
    delete[] Y0_tmp;
}
//***************************************************









//**************************************************************************

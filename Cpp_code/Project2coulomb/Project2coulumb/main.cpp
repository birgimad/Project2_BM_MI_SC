#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
#include <iomanip>
using std::setw;
#include "armadillo"
using namespace arma;

void find_max(mat &A, int &n, int &row_number, int &column_number)
// set row_number = 0 and column_number = 1, when running the code. These are the initial guesses for max(A(i,j))
{
    //Finding maximum matrix element (and its row and column numbers) of the non-diagonal A
    //we are only considering the upper part of the matrix when seaching for the largest element
    //since the matrix is symmetric
        double max = A(0,1);

        for (int i=0; i<n; i++) //Compare entrances in the upper part of the matrix. Choose largest one (abs-value).
        {
            for (int j=i+1; j<n; j++)
            {

            if (fabs(A(i,j)) > fabs(max))
            {
                max = A(i,j);
                row_number = i;
                column_number = j;
            }
            }
        }

        return;
}

void Jacobi (mat &A, int n, double epsilon)
{
    int row_number, column_number;
    row_number = 0;
    column_number = 1;
    find_max(A,n,row_number,column_number);
    double max, max_number_of_iterations;
    max = A(row_number, column_number);
    max_number_of_iterations = n*n*n;

    // Use simpler test:
    int m=0;
    while(pow(fabs(max),2)> epsilon && m<max_number_of_iterations)
    {
        // computing tau, tan, cos and sine
        double tau, c, s, t;
        double a_kk, a_ll, a_kl, a_ik, a_il;
        a_kk = A(row_number, row_number);
        a_ll = A(column_number, column_number);
        a_kl = A(row_number, column_number);
        tau = (a_ll-a_kk)/(2*a_kl);
        if (tau >= 0)
        {
            t = (-tau+sqrt(1+pow(tau,2)));
        }
        else
        {
            t = (-tau-sqrt(1+pow(tau,2)));
        }
        c = 1/sqrt(1+pow(t,2));
        s = t*c;

        // Computing the new matrix A
        for (int i = 0; i<n; i++)
        {
            if (i != row_number && i != column_number)  // determining A(i,k) for new matrix
            {
                a_ik = A(i,row_number);
                a_il = A(i,column_number);
                A(i,row_number) = a_ik*c - a_il*s;
                A(i,column_number) = a_il*c + a_ik*s;
                A(row_number,i) = A(i,row_number);
                A(column_number,i) = A(i,column_number);
            }
        }
        A(row_number, row_number) = a_kk*pow(c,2) - 2*a_kl*c*s+a_ll*pow(s,2);
        A(column_number, column_number) = a_ll*pow(c,2) + 2*a_kl*c*s+a_kk*pow(s,2);
        A(column_number, row_number) = 0.00;  //By choice of theta
        A(row_number, column_number) = 0.00;  //By choice of theta

        row_number = 0;
        column_number = 1;

        find_max(A,n,row_number,column_number);
        max = A(row_number,column_number);

        m += 1;
/*
        cout << "number of iterations = " << m << endl;
        cout << "new A = " << endl;
        cout << A << max << setw(10) << row_number << setw(10) << column_number << endl;
*/
    }
    cout << "number of iterations = " << m << endl;
    cout << "maximum number of iterations = " << max_number_of_iterations << endl;
}

int main()
{
    double rho_min = 0;
    double rho_max = 10;
    //double epsilon;
    //epsilon = 1.0e-8;
    int n;
    cout << "Please enter value of n:\n>";
    cin >> n;
    cout << "n = " << n << endl;
    cout << "rho_min = " << rho_min << endl;
    cout << "rho_max = " << rho_max << endl;

    double h = (rho_max - rho_min)/(n+1); //step length

    double w = 5.0;


    vec V(n);
    for (int i = 0; i < n; i++)
    {
        V(i) = pow(rho_min + (i+1)*h,2)*pow(w,2);    //V(0) = V_1 = (rho_min + 1*h)^2 etc + 1/(rho_min +(i+1)*h).
    }
    mat A(n,n);
    A.zeros();

    for (int i=0; i<n; i++)
    {
        A(i,i) = 2/pow(h,2) + V(i);
    }

    double off_diagonal;
    off_diagonal = -1/pow(h,2);

    for (int i=1; i<n; i++)
    {
        A(i,i-1) = off_diagonal;
    }
    for (int i=0; i<n-1; i++)
    {
        A(i,i+1) = off_diagonal;
    }

    /*
    Jacobi(A,n,epsilon);

    vec eigen_values(n);

    for (int i=0; i<n; i++)
    {
        eigen_values(i) = A(i,i);
    }

    //looking for the 1st, 2nd, and 3rd eigen values
    double min1;
    min1 = 100;
    for (int i=0; i<n; i++)
    {

    if (fabs(eigen_values(i)) < fabs(min1))
    {
        min1 = eigen_values(i);
    }
    }
    double min2;
    min2 = 100;
    for (int i=0; i<n; i++)
    {

    if (fabs(eigen_values(i)) < fabs(min2) && fabs(eigen_values(i)) > fabs(min1))
    {
        min2 = eigen_values(i);
    }
    }
    double min3;
    min3 = 100;
    for (int i=0; i<n; i++)
    {

    if (fabs(eigen_values(i)) < fabs(min3) && fabs(eigen_values(i)) > fabs(min2))
    {
        min3 = eigen_values(i);
    }
    }

    //Printing 1st, 2nd and 3rd eigenvalue
    cout << "min1 = " << min1 << endl;
    cout << "min2 = " << min2 << endl;
    cout << "min3 = " << min3 << endl;
    */



    vec eigval;
    mat eigvec;

    eig_sym(eigval, eigvec, A);

    ofstream myfile ("wave_omega_no_interaction_500.txt");
        if (myfile.is_open())
        {
            myfile << "omega = "<< w << setw(10)<< "n = "<< n << setw(20) << "eigen value = " << eigval(0) <<endl;
            myfile << "rho" << setw(20) << "eigvec" << endl;
            for(int i = 0; i<n ; i++){
            myfile << (i+1)*h << setw(20) << eigvec(i,0) << endl;
            }

        }
    cout << eigval(0) <<endl;
    cout << "rho" << setw(10) << "eigvec" << endl;
    for(int i = 0; i<n ; i++){
    cout << (i+1)*h << setw(20) << eigvec(i,0) << endl;
    }




    return 0;
}

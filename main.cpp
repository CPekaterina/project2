#include <src/UnitTest++.h>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <fstream>
#include <src/lib.h>
#include <main.h>
#include <time.h>
using namespace std;



TEST(WillFail) {
    CHECK(false);
}

int main()
{

    cout << "EIGENVALUES AND EIGENVECTORS" <<  endl
         << "SOLVER FOR 1-2 ELECTRONS IN AN OSCILLATOR POTENTIAL" << endl << endl;

    bool runagain;
    do
    {
    //read in the number n of steps
    int n_step;
    cout << "How many steps should be taken? n_step= " ;
    cin >> n_step;
    cout << endl;

    //size of the matrix is then:

    int n = n_step - 1;
    
    //read in the maximum radius (dimensionless) rho

    double rho_max;
    cout << "Up to what dimensionless radius should the solution be computed? rho_max= " ;
    cin >> rho_max;
    cout << endl;

    //compute step h

    double h = rho_max/double(n_step);
    cout << "The steplength is: " << h << endl;

    // write the matrix and the diagonal/nondiagonal arrays for Schrödingers equation for one electron:

        double *b; //used for tqli function
        double *nd; //used for tqli function
        b = new double [n];
        nd = new double [n];

    double **A;
    A = new double* [n];
    for (int i = 0; i < n; i++)
        A[i] = new double[n];
    for (int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            double rho=(double(i)+1.)*h;
            if(i==j){A[i][j]=double(2)/(h*h)+rho*rho;b[i]=A[i][j];} // d_i=2/h^2+rho_i^2 with rho_i=(i+1)*h
            else if(i==j+1 || i==j-1){A[i][j]=-1./(h*h);nd[i]=A[i][j];}                            // -e_i=-1/h^2
            else{
                A[i][j]=0;                                                         // nondiagonal elements
            }
        }
    }

    // write the matrix and the diagonal/nondiagonal arrays for two non interacting electrons:


    //read in the potential strength
    double omega;
    cout << "Type the oscillator strength omega_r=";
    cin >> omega;
    cout << endl;

    double **B;
    B = new double* [n];

    for (int i = 0; i < n; i++)
        B[i] = new double[n];
    for (int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            double rho=(double(i)+1.)*h;
            if(i==j){B[i][j]=double(2)/(h*h)+rho*rho*omega*omega+double(1)/rho;b[i]=B[i][j];} // d_i=2/h^2+rho_i^2*omega^2+1/rho_i with rho_i=(i+1)*h
            else if(i==j+1 || i==j-1){B[i][j]=-1./(h*h); nd[i]=B[i][j];}                         // -e_i=-1/h^2
            else{
                B[i][j]=0;                                                         // other nondiagonal elements
            }
        }
    }


    //set up eigenvector-matrix R

    double **R;
    R = new double* [n];
    for (int i = 0; i < n; i++)
        R[i] = new double[n];
    for (int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            if(i==j)R[i][j]=1;
            else R[i][j]=0;
        }
    }

    //index-array

    int*ind;
    ind = new int[n];
    for (int i=0;i<n;i++){ind[i]=i;}

    //eigenvalue-array

    double *v;
    v= new double[n];

    //main computations

    bool func;
    cout << "Should jacobi or tqli be used for the computations?(1/0)" << endl;
    cin >> func;

    if(func==1){

    jacobi(A,R,v,n); //find eigenvalues and -vectors
    bsort(v,ind,n); //sort the eigenvalues

    }
    else{
    //use of tqli function

    clock_t start, finish;
    start = clock();

    tqli(b,nd,n,R);
    bsort(b,ind,n);

    finish = clock();
    clock_t t = finish - start;

    cout << "Execution time tqli: " << ((double)t)/CLOCKS_PER_SEC << "sec." << endl;
    }

    //print the first few eigenvalues
    double sum=0;
    cout << endl << "The first few eigenvalues are:" << endl;
    for(int i=0;i<3;i++)
    {
        cout << v[i] << endl;
        sum+=v[i];
    }
    sum=(sum-21.);
    if (sum<3*10e-5){cout << "four leading digits reached for n=" << n_step << endl;}


    cout << endl << "relerr=" << sum;

    cout << endl;

    //print results into a .dat file

    double *z,*y;
    z= new double[n];
    y= new double[n];

    //of  which state of excitement?

    int state;
    cout << "Which state to you want to visualize? N=" << endl;
    cin >> state;
    int N=ind[state];

    for(int i=0;i<n;i++)
    {

        z[i]=(i+1.0)*h;
        y[i]=R[i][N]*R[i][N]/(h);
    }

    //read output filename

    char name[30] ;
    cout << "Name the output file: ";
    cin >> name;
    write(z,y,n,name);

    //repeat the procedure?

    cout << "Run again? (1/0)" << endl;
    cin >> runagain;
    cout << endl << endl;
    }while(runagain==true);

    return UnitTest::RunAllTests(); //default unit test, no functionality
}


// for testing purposes: print a matrix or an array

void printmatrix(double ** A, int n, int m)
{

    for(int i=0;i<n;i++)
    {
        cout << "| ";
        for(int j=0; j<m; j++)
        {
            cout << setw(6) << A[i][j] << " ";
        }
        cout << "|" << endl;
    }
}
void printarray(double *A, int n)
{
    cout << "| ";
    for(int i=0;i<n;i++)
    {
            cout << setw(6) << A[i] << " ";
    }
    cout << "|" << endl;
}
void printarray(int *A, int n)
{
    cout << "| ";
    for(int i=0;i<n;i++)
    {
            cout << setw(6) << A[i] << " ";
    }
    cout << "|" << endl;
}



//jacobi algorithm

void jacobi (double **A, double **R, double *v, int n)
{
    int k,l,z=0;
    int maxiter =10e6;

    clock_t start, finish;
    start = clock();


    while(max_nondig(A,n,&k,&l)==false && z<maxiter)     // while the max^2 of one element of lower-triangle is larger than epsilon
    {
        double t1, t2, t, tau, c, s;
        tau=(A[l][l]-A[k][k])/(double(2)*A[k][l]); // formulas to obtain cos() and sin()
        t1= -tau+sqrt(double(1)+tau*tau);
        t2= -tau-sqrt(double(1)+tau*tau);
        if(fabs(t1)>fabs(t2)){t=t2;}          //use the smaller of the roots
        else{t=t1;}
        c=double(1)/sqrt(double(1)+t*t);
        s=c*t;
        z++;                                //increase number of iterations
        if (z==maxiter-1)
        {
            cout << "Maximum number of iterations reached, I am tired and will stop working. Sorry." << endl;
            break;
        }
        jacobi_rot(s,c,k,l,n,A,R);            //rotate A with c and s to put A[k][l] to zero


    }

    finish = clock();
    clock_t t = finish - start;
    cout << "Execution time jacobi: " << ((double)t)/CLOCKS_PER_SEC << "sec." << endl;

    //print the eigenvalues into v and the number of iterations

    for(int i=0;i<n;i++)v[i]=A[i][i];

    cout << "Number of iterations: " << z << endl;


return;
}


//rotation function

void jacobi_rot (double s, double c,int k,int l,int n, double **A, double **R)
{
    for(int i=0;i<n;i++) //see also formulas from the lecture notes
    {
      if(i!=k && i!=l)
      {
          double tempik = A[i][k];
          double tempil = A[i][l];
          A[i][k]=tempik*c-tempil*s;
          A[i][l]=tempil*c+tempik*s;
          A[k][i]=A[i][k];
          A[l][i]=A[i][l];
      }
     }
      double tempkk=A[k][k];
      double templl=A[l][l];
      double tempkl=A[k][l];
      A[k][k]=tempkk*c*c-double(2)*tempkl*c*s+templl*s*s;
      A[l][l]=templl*c*c+double(2)*tempkl*c*s+tempkk*s*s;
      A[k][l]=A[l][k]=0; //hard coding of the result

    for(int i =0;i<n;i++)
    {
        //eigenvectors in columns of R

            double tempki=R[i][k];
            double templi=R[i][l];
            R[i][k]=tempki*c-templi*s;
            R[i][l]=templi*c+tempki*s;
    }


}


// finds the maximum of the non-diagonal matrix elements in the lower tri-diagonal and tests
// if it is larger than the tolerance E

bool max_nondig(double ** A, int n, int* k, int* l)
{
    *k=1;
    *l=0;
    double E = 1.e-12;                   //tolerance E
   double max = fabs(A[*k][*l]);

   for(int i=1;i<n;i++)                    //iteration over the non-diagonal matrix elements in the lower tri-diagonal
   {
       for(int j=0;j<i;j++)
       {
           if(fabs(A[i][j])>max)
           {
               max = fabs(A[i][j]);
               *k = i;                      // write row of the new maximum to k;
               *l = j;                      // write column of the new maximum to l;
           }
         if((max*max)<=E)A[i][j]=0.0; //set values below E to 0
       }
   }
   return((max*max)<E);
}

//bubblesort
void bsort(double*v,int *ind,int n)
{
    double min, temp;
    for(int j=0;j<n;j++)
    {
        int k=j;
        min=v[j];
        for(int i=j;i<n;i++)
        {
            if(min>v[i]){min=v[i];k=i;}

        }

        v[k]=v[j];
        v[j]=min;

       temp=ind[k];
       ind[k]=ind[j];
       ind[j]=temp;




    }

}


//writes two arrays in a .dat file

void write(double *z, double *y, int n, char *file)
{
    ofstream resout;
    resout.open(file);
    for (int i=0; i<n; i++)
    {
        resout << setprecision(15) << setw(19) << z[i] << " " << setprecision(15) << setw(19) << y[i] << endl;
    }
    resout.close();
}

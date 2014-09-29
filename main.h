#ifndef MAIN_H
#define MAIN_H


void printmatrix(double ** A, int n, int m);
/*prints a n x m matrix
 */

void printarray(double *A, int n);
/*prints an array
*/

void jacobi (double **A, double **R, double *v,int n);
/*jacobi diagonalization method. Diagonalizes a symmetric matrix A of size n x n,
 * saves eigenvalues in v, saves eigenvectors in columns of R
 */

void jacobi_rot (double s, double c,int k,int l,int n, double **A, double **R);
/*performs a similarity transformation on a n x n matrix A with a rotation matrix S with the values c (cos)
 *  and +/-s (sin) at S[k][k], S[l][l] and S[k][l], S[l][k] respectively. Rotates the eigenvector
 * matrix R
 */

bool max_nondig(double ** A, int n, int* k, int* l);
/*finds the maximum nondiagonal matrix element of a n x n matrix A and turns
 * values smaller than a minimal value to 0, sets A[k][l] to the maximum value,
 * if the maximum value squared is smaller than the minimum values returns 1
 */

void bsort(double*v,int *ind, int n);
/* bubblesorts array v of the dimension n and writes the adress changes into a pointer array ind.
 */


void write(double *z, double *y, int n, char *file);
/* writes the x-values from teh array z and the y-values from the array y of the dimenssion n into a file
 * named file.
 */

void printarray(int *A, int n);
void printarray(double *A, int n);
/* prints the values in an array A onto the screen
 */

#endif // MAIN_H

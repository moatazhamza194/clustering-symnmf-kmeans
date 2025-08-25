# ifndef SYMNMF_H_
# define SYMNMF_H_

int Dim(FILE *file);
int N(FILE *file);
void print_Matrix(double **M, int n,int k);
double Euclidean(double *a, double *b, int dim);
double **Matrix_Multi(double **A, double **B, int n,int m, int l);
double Frobenius(double **A, double **B, int n, int k);
double **Transpose(double **A, int n, int m);
void SQ_Matrix(double **D, int n);
void Matrix_Free(double **M, int n);
double **sym(double **Points, int n, int dim);
double **ddg(double **A, int n);
double **norm(double **A, double **D, int n);
double **symnmf(double **H, double **W, int n, int k);

# endif
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include <string.h>
#include "symnmf.h"

/*gets the Dimension of the Points in the input file */
int Dim(FILE *file){
  int c,dim=0;
  c=fgetc(file);

  while(c !='\n'){
    
    if(c==','){
      dim++;
    }
    c=fgetc(file);
  }
  dim++;

  return dim;
  
}

/*gets the Number of the Points in the input file*/
int N(FILE *file){
  int c,n=0;
  c=fgetc(file);

  while(c!=EOF){
    if(c=='\n'){
      n++;
    }
    c=fgetc(file);
  }
  return n;
  
  
}

/*Prints an (n,k) Matrix*/
void print_Matrix(double **M, int n,int k){
    int i,j;

    for(i=0;i<n;i++){
        for(j=0;j<k-1;j++){
            printf("%.4f,",M[i][j]);
        }
        printf("%.4f\n",M[i][k-1]);
        
    }

}

/*Squared Euclidean distance between 2 points*/

double Euclidean(double *a, double *b, int dim){
    int i;
    double result=0;

    for(i=0;i<dim;i++){
        result=result+ pow(a[i]-b[i],2);
    }
    return result;
}

/*Multiplying an (n,m) matrix with an (m,l) matrix*/

double **Matrix_Multi(double **A, double **B, int n,int m, int l){
    int i,j,t;
    double sum=0;
    double **result=(double**)calloc(n, sizeof(double*));
    if(result==NULL){
        printf("AN ERROR HAS OCCURED!");
            return NULL;
    }
    for(i=0;i<n;i++){
        result[i]=(double*)calloc(l,sizeof(double));
        if(result[i]==NULL){
        printf("AN ERROR HAS OCCURED!");
            return NULL;
    }
    }

    for(i=0;i<n;i++){
        for(j=0;j<l;j++){
            sum=0;
            for(t=0;t<m;t++){
                sum+=A[i][t]*B[t][j];
            }
            result[i][j]=sum;

        }

    }

    return result;

}


/*Squared Frobenius Norm for (n,n) matrices*/

double Frobenius(double **A, double **B, int n, int k){
    int i,j;
    double result=0;

    for(i=0;i<n;i++){
        for(j=0;j<k;j++){
            result+= pow(A[i][j]-B[i][j],2);

        }
    }
    return result;
}

/*Transpose of an (n,m) Matrix*/

double **Transpose(double **A, int n, int m){
    int i,j;
    double **result=(double**)calloc(m,sizeof(double*));
    if(result==NULL){
        printf("AN ERROR HAS OCCURED!");
            return NULL;
    }
    for(i=0;i<m;i++){
        result[i]=(double*)calloc(n,sizeof(double));
        if(result[i]==NULL){
        printf("AN ERROR HAS OCCURED!");
            return NULL;
    }
    }

    for(i=0;i<m;i++){
        for(j=0;j<n;j++){
            result[i][j]=A[j][i];
        }
    }
    return result;
}

/*Raising a Diagonal (n,n) Matrix to the power of -(1/2)*/

void SQ_Matrix(double **D, int n){
    int i;

    for(i=0;i<n;i++){
        D[i][i]=pow(D[i][i],-0.5);
    }
}

/*Freeing a (n,k) Matrix */

void Matrix_Free(double **M, int n){
    int i;
    for(i=0; i<n;i++){
        free(M[i]);
    }
    free(M);
}

double **sym(double **Points, int n, int dim){
    int i,j;
    double **A= (double**)calloc(n,sizeof(double*));
    if(A==NULL){
        printf("AN ERROR HAS OCCURED!");
            return NULL;
    }
    for(i=0;i<n;i++){
        A[i]=(double*)calloc(n,sizeof(double));
        if(A[i]==NULL){
        printf("AN ERROR HAS OCCURED!");
            return NULL;
    }
    }

    for(i=0; i<n; i++){
        for(j=0;j<n;j++){
            if(i!=j){
                A[i][j]=exp(-((Euclidean(Points[i],Points[j],dim))/2));
            }
            else{
                A[i][j]=0;
            }
            
        }
    }

    return A;
}

double **ddg(double **A, int n){
    int i,j;
    double sum;
    double **D;
    double *Row_Sum=(double*)calloc(n,sizeof(double));
    if(Row_Sum==NULL){
        printf("AN ERROR HAS OCCURED!");
            return NULL;
    }
    D= (double**)calloc(n,sizeof(double*));
    if(D==NULL){
        printf("AN ERROR HAS OCCURED!");
            return NULL;
    }
    for(i=0;i<n;i++){
        D[i]=(double*)calloc(n,sizeof(double));
        if(D[i]==NULL){
        printf("AN ERROR HAS OCCURED!");
            return NULL;
    }
    }

    for(i=0;i<n;i++){
        sum=0;
        for(j=0;j<n;j++){
            sum+=A[i][j];
        }
        
        Row_Sum[i]=sum;
        
    }

    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(i==j){
                D[i][i]=Row_Sum[i];

            }
            else{
                D[i][j]=0;

            }
        }
    }

    free(Row_Sum);

    return D;

}

double **norm(double **A, double **D, int n){
    double **W,**W_0;
    SQ_Matrix(D,n);
    W_0=Matrix_Multi(D,A,n,n,n);
    W=Matrix_Multi(W_0, D, n,n,n);
    Matrix_Free(W_0,n);
    return W;

}




int main(int argc, char **argv)
{
    char *goal;
    double **points,**result_sym,**result_ddg, **result_norm, num;
    int i,j,n,dim,check;
    FILE *file;
    

    /*Read Points from File into a Matrix*/
    file = fopen(argv[2], "r");
    if (file == NULL) {
        printf("AN ERROR HAS OCCURED!");
        return 1;
    }
    argc=argc;

    /*get the Dimensions of the Matrix*/
    n=N(file);
    fseek(file, 0, SEEK_SET );
    dim=Dim(file);
    fseek(file, 0, SEEK_SET );

    

    points=(double**)calloc(n,sizeof(double*));
    if(points==NULL){
        printf("AN ERROR HAS OCCURED!");
        return 1;
    }

    for ( i = 0; i < n; i++)
    {
        points[i]=(double*)calloc(dim,sizeof(double));
        if(points[i]==NULL){
        printf("AN ERROR HAS OCCURED!");
        return 1;
    }

    }

    /*FILL THE MATRIX*/

    for (i=0;i<n;i++){
        for(j=0;j<dim;j++){
            check=fscanf(file,"%lf",&num);
            points[i][j] = num;
            fgetc(file);
        }
    }

    check++;

    fclose(file);

    goal=argv[1];

    

    if(strcmp(goal,"sym")==0){
        result_sym=sym(points,n,dim);
        print_Matrix(result_sym,n,n);
        Matrix_Free(result_sym,n);


    }

    if(strcmp(goal,"ddg")==0){
        result_sym=sym(points,n,dim);
        result_ddg=ddg(result_sym,n);
        print_Matrix(result_ddg,n,n);
        Matrix_Free(result_sym,n);
        Matrix_Free(result_ddg,n);
    }

    if(strcmp(goal,"norm")==0){
        result_sym=sym(points,n,dim);
        result_ddg=ddg(result_sym,n);
        result_norm=norm(result_sym,result_ddg,n);
        print_Matrix(result_norm,n,n);
        Matrix_Free(result_sym,n);
        Matrix_Free(result_ddg,n);
        Matrix_Free(result_norm,n);
 
    }

    Matrix_Free(points,n);


return 0;
}



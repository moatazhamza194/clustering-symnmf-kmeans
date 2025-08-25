#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "symnmf.h"
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include <string.h>

double **symnmf(double **H, double **W, int n, int k){
    int iter=0,i,j;
    double **new_H, **WH,**HT,**HHTH,**HHT;
    double parm;

    while(iter<300){
        new_H=(double**)calloc(n,sizeof(double*));
        if(new_H==NULL){
        printf("AN ERROR HAS OCCURED!");
            return NULL;
    }
        for(i=0;i<n;i++){
            new_H[i]=(double*)calloc(k,sizeof(double));
            if(new_H[i]==NULL){
            printf("AN ERROR HAS OCCURED!");
            return NULL;
    }
        }

        WH=Matrix_Multi(W,H,n,n,k);
        HT=Transpose(H,n,k);
        HHT=Matrix_Multi(H,HT,n,k,n);
        HHTH=Matrix_Multi(HHT,H,n,n,k);

        for(i=0;i<n;i++){
            for(j=0;j<k;j++){
                parm=0.5 + 0.5*(WH[i][j]/HHTH[i][j]);
                new_H[i][j]=H[i][j]*parm;
                
            }
        }

        if(Frobenius(new_H,H,n,k)<0.0001){
            Matrix_Free(H,n);
            Matrix_Free(HT,k);
            Matrix_Free(HHT,n);
            Matrix_Free(HHTH,n);
            Matrix_Free(WH,n);
            break;
        }

        Matrix_Free(H,n);
        Matrix_Free(HT,k);
        Matrix_Free(HHT,n);
        Matrix_Free(HHTH,n);
        Matrix_Free(WH,n);
        H=new_H;
        iter++;

        
    }
    return new_H;


}

double **Build_Points(PyObject *points,int n, int dim){
    double **c_points;
    int i,j;
    PyObject *lst_in_lst;
    PyObject *item;

    /*creating the points*/
    c_points=(double **)calloc(n,sizeof(double *));
    if(c_points==NULL){
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        return NULL;
    }
    for(i=0;i<n;i++){
        c_points[i]=(double *)calloc(dim,sizeof(double));
        if(c_points[i]==NULL){
                PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
                return NULL;
            }

    }

    for(i=0;i<n;++i){
        lst_in_lst=PyList_GetItem(points,i);
        for(j=0;j<dim;++j){
            item=PyList_GetItem(lst_in_lst,j);
            c_points[i][j]=PyFloat_AsDouble(item);

        }
    }

    return c_points;

}

PyObject *C_to_Python(double **c_result,int n, int dim){
    int i,j;
    PyObject *res_lst;
    PyObject *res_lst_in_lst;
    PyObject *f_item;

    /*Converting the Result in C to Python*/

    res_lst=PyList_New(n);
    for(i=0;i<n;i++){
        res_lst_in_lst=PyList_New(dim);
        for(j=0;j<dim;j++){
            f_item=PyFloat_FromDouble(c_result[i][j]);
            PyList_SetItem(res_lst_in_lst,j,f_item);
        }
        PyList_SetItem(res_lst,i,res_lst_in_lst);

    }

    return res_lst;

}


static PyObject* sym_1(PyObject *self, PyObject *args)
{
    int i,n,dim;
    PyObject *points,*result;
    double **c_points;
    double **c_results;


    if(!PyArg_ParseTuple(args,"O", &points)) {
        return NULL; 
    }

    n=PyObject_Length(points);
    dim=PyObject_Length(PyList_GetItem(points,0));
    c_points=Build_Points(points,n,dim);
    

    /*Getting the Result in C*/

    c_results=sym(c_points,n,dim);
    if(c_results==NULL){
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        return NULL;
    }

    result=C_to_Python(c_results,n,n);

    /*Freeing Memory*/

    for(i=0; i<n ;i++){
        free(c_points[i]);
        free(c_results[i]);
    }

    free(c_points);
    free(c_results);

    return result;

}

static PyObject* ddg_1(PyObject *self, PyObject *args)
{
    int i,n,dim;
    PyObject *points,*result;
    double **c_points;
    double **c_results,**c_results_sym;


    if(!PyArg_ParseTuple(args,"O", &points)) {
        return NULL; 
    }

    n=PyObject_Length(points);
    dim=PyObject_Length(PyList_GetItem(points,0));
    c_points=Build_Points(points,n,dim);

    /*Getting the Result in C*/

    c_results_sym=sym(c_points,n,dim);
    if(c_results_sym==NULL){
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        return NULL;
    }
    c_results=ddg(c_results_sym,n);
    if(c_results==NULL){
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        return NULL;
    }

    result=C_to_Python(c_results,n,n);

    /*Freeing Memory*/

    for(i=0; i<n ;i++){
        free(c_points[i]);
        free(c_results[i]);
        free(c_results_sym[i]);
    }

    free(c_points);
    free(c_results);
    free(c_results_sym);

    return result;

}

static PyObject* norm_1(PyObject *self, PyObject *args)
{
    int i,n,dim;
    PyObject *points,*result;
    double **c_points;
    double **c_results,**c_results_sym,**c_results_ddg;


    if(!PyArg_ParseTuple(args,"O", &points)) {
        return NULL; 
    }

    n=PyObject_Length(points);
    dim=PyObject_Length(PyList_GetItem(points,0));
    c_points=Build_Points(points,n,dim);

    /*Getting the Result in C*/

    c_results_sym=sym(c_points,n,dim);
    if(c_results_sym==NULL){
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        return NULL;
    }
    c_results_ddg=ddg(c_results_sym,n);
    if(c_results_ddg==NULL){
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        return NULL;
    }
    c_results=norm(c_results_sym,c_results_ddg,n);
    if(c_results==NULL){
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        return NULL;
    }

    result=C_to_Python(c_results,n,n);

    /*Freeing Memory*/

    for(i=0; i<n ;i++){
        free(c_points[i]);
        free(c_results[i]);
        free(c_results_sym[i]);
        free(c_results_ddg[i]);
    }

    free(c_points);
    free(c_results);
    free(c_results_sym);
    free(c_results_ddg);

    return result;

}

static PyObject* symnmf_1(PyObject *self, PyObject *args)
{
    int i,n,k;
    PyObject *W,*H,*result;
    double **c_W,**c_H;
    double **c_results;


    if(!PyArg_ParseTuple(args,"OOi", &H,&W,&k)) {
        return NULL; 
    }
    n=PyObject_Length(W);
    c_W=Build_Points(W,n,n);
    c_H=Build_Points(H,n,k);

    /*Getting the Result in C*/

    c_results=symnmf(c_H,c_W,n,k);
    if(c_results==NULL){
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        return NULL;
    }

    result=C_to_Python(c_results,n,k);

    /*Freeing Memory*/


    for(i=0; i<n ;i++){
        free(c_W[i]);
        free(c_results[i]);
    }
    free(c_W);
    free(c_results);
    return result;

}

static PyMethodDef symnmfMethods[] = {
    {"sym",                   
      (PyCFunction) sym_1, 
      METH_VARARGS,         
      PyDoc_STR("The Function receives the datapoints X")},
      {"ddg",                   
      (PyCFunction) ddg_1, 
      METH_VARARGS,         
      PyDoc_STR("The Function receives the datapoints X")},
      {"norm",                   
      (PyCFunction) norm_1, 
      METH_VARARGS,         
      PyDoc_STR("The Function receives the datapoints X")},
      {"symnmf",                   
      (PyCFunction) symnmf_1, 
      METH_VARARGS,         
      PyDoc_STR("The Function receives the Matrix W,H and and Integer K")},

    {NULL, NULL, 0, NULL}     
};


static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "symnmf", 
    NULL, 
    -1,  
    symnmfMethods 
};

PyMODINIT_FUNC
PyInit_symnmf(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}

#include <Python.h>
#include <sys/types.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <float.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <setjmp.h>
#include "chooch.h"

static jmp_buf jump_buffer;
char *sElement;
char *sEdge;
char cScanTitle[TITLE]="";
int id1=0, id2=0;
int verbose, status, silent, kev;
double fpInfl, fppInfl, fpPeak, fppPeak, EInfl, EPeak;
double fE1=0.0, fE2=0.0, fE3=0.0, fE4=0.0;
double fEres=0.00014;

// forward decl.
PyObject* chooch_calc(double *fXraw, double *fYraw, int nDataPoints, const char* element, const char *edge, const char* outfile); 

static PyObject* PyChooch_calc(PyObject *self, PyObject *args) {
  const char *element;
  const char *edge;
  const char *filename;
  char *outfile;
  PyObject *chooch_data;
  PyObject *temp;
  int npoints, i;
  double* xdata, *ydata;
  PyObject* xi;
  PyObject* yi;
  PyObject* res;

  // check arguments
  outfile = (char *) NULL;
  if (!PyArg_ParseTuple(args, "Oss|s", &chooch_data, &element, &edge, &outfile)) 
    return NULL;

  if(!PySequence_Check(chooch_data)) {
    PyErr_SetString(PyExc_TypeError, "argument 1 must be a list of data points: ((X1, Y1), (X2, Y2), ..., (Xn, Yn))");
    return NULL;
  } 

  // get X data and Y data from chooch_data sequence
  npoints = PySequence_Size(chooch_data);
  xdata = malloc(npoints * sizeof(double));
  ydata = malloc(npoints * sizeof(double));
  for (i=0; i<npoints; i++) {
    temp = PySequence_GetItem(chooch_data, i);
    if (!PySequence_Check(temp)) {
      free(xdata);
      free(ydata);
      Py_DECREF(temp);
      Py_DECREF(chooch_data);
      PyErr_SetString(PyExc_ValueError, "data should be a list of data points: ((X1, Y1), (X2, Y2), ..., (Xn, Yn))");
      return NULL;
    }
    if (PySequence_Size(temp)!=2) {
      free(xdata);
      free(ydata);
      Py_DECREF(temp);
      Py_DECREF(chooch_data);
      PyErr_SetString(PyExc_ValueError, "data should be a list of data points: ((X1, Y1), (X2, Y2), ..., (Xn, Yn))");
      return NULL;
    }
    xi = PySequence_GetItem(temp, 0);
    if (! PyFloat_Check(xi)) {
      Py_DECREF(xi);
      free(xdata);
      free(ydata);
      Py_DECREF(temp);
      Py_DECREF(chooch_data);
      PyErr_SetString(PyExc_TypeError, "data point values should be float");
      return NULL;
    } else {
      xdata[i]=PyFloat_AS_DOUBLE(xi);
      Py_DECREF(xi);
    }
    yi = PySequence_GetItem(temp, 1);
    if (! PyFloat_Check(yi)) {
      Py_DECREF(yi);
      free(xdata);
      free(ydata);
      Py_DECREF(temp);
      Py_DECREF(chooch_data);
      PyErr_SetString(PyExc_TypeError, "data point values should be float");
      return NULL;
    } else {
      ydata[i]=PyFloat_AS_DOUBLE(yi);
      Py_DECREF(yi);
    }
    Py_DECREF(temp);
  }

  // do chooch calculation
  if (setjmp(jump_buffer) == 0) {
    res=chooch_calc(xdata, ydata, npoints, element, edge, outfile);
  } else {
    // Error occurred, handle it
    PyErr_SetString(PyExc_RuntimeError, "GSL error occurred");
    return NULL;
  }

  return res;
}


static PyMethodDef PyChoochMethods[] = {
    {"calc",  PyChooch_calc, METH_VARARGS,
     "Input arguments: raw data from energy scan in the form ((X1, Y1), ..., (Xn, Yn)) then Element (e.g 'Se') and Edge (e.g 'K')\nOutput: (Epeak, fppPeak, fpPeak, Einfl, fppInfl, fpInfl, ((X1, Yspline1, Yfp1),...,(Xn, Ysplinen, Yfpn)) )" },
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


void custom_err_handler(const char* reason, const char* file, int line, int gsl_errno) {
    fprintf(stderr, "GSL error: %s:%d: %s (error code: %d)\n", file, line, reason, gsl_errno);

    // Calling directly PyErr_SetString here didn't work
    // Jump back to the recovery point
    longjmp(jump_buffer, 1);
};

// Define the module
static struct PyModuleDef PyChoochModule = {
    PyModuleDef_HEAD_INIT,
    "PyChooch",   /* name of module */
    NULL,         /* module documentation, may be NULL */
    -1,           /* size of per-interpreter state of the module,
                     or -1 if the module keeps state in global variables. */
    PyChoochMethods
};

// Initialization function for the module
PyMODINIT_FUNC
PyInit_PyChooch(void)
{
    // Install our own error handler for the GSL
    gsl_set_error_handler(&custom_err_handler);

    // Create the module
    return PyModule_Create(&PyChoochModule);
}

static int my_efswrite(const char *filename, double *x, double *y1, double *y2, int n) {
  int    i;
  int    len;
  FILE   *ff;
  //
  if ((ff = fopen(filename, "w")) == NULL) {
     printf("Cannot open %s for write\n", filename);
     return EXIT_FAILURE;
  }

  //  printf("Title: %s\nNo. data points: %d\n", cScanTitle, *nDataPoints);
  if(!silent)printf("Writing anomalous scattering factors to %s\n", filename);
  for (i = 0; i < n; i++) {
     fprintf(ff, "%10.4f  %7.2f  %7.2f\n", x[i], y1[i], y2[i]);
     if(verbose>0)printf("%10.4f  %7.2f  %7.2f\n", x[i], y1[i], y2[i]);
     //    printf("%10.3f  %10.3f\n", x[i], y[i]);
  }
  fclose(ff);
  return EXIT_SUCCESS;
}

PyObject* chooch_calc(double *fXraw, double *fYraw, int nDataPoints, const char* element, const char *edge, const char *outfile) {
  PyObject *res;
  PyObject *graph;
  PyObject *graph_point;
  int i, j, err;
  float fXref, fYref, fXcur, fYcur;
  char  ch[1];
  char opt;
  //
  int nFit, nPoints, plotX=0, psplot=0, pngplot=0, display=0;
  int nSavWin;
  //
  double dE, fEdge;
  double fYfita[MAXSIZE], fYfitb[MAXSIZE];  // Output
  double fYspline[MAXSIZE], fXfpp[MAXSIZE];  // Output values of Energy and fpp
  double fYsmooth[MAXSIZE], fYnorm[MAXSIZE]; // Smoothed and normalised fluorescence spectra
  double fYfpp[MAXSIZE], fYfpps[MAXSIZE], fYfp[MAXSIZE]; // Raw f'', smoothed f''
  double fYDeriv1[MAXSIZE], fYDeriv2[MAXSIZE], fYDeriv3[MAXSIZE];
  double C[3], result, error, fMid;
  //
  double fC, fM;
  //
  verbose=kev=0;
  silent=1;

  sElement = element;
  sEdge = edge;
  /********************************
   * Start output and calculations
   ********************************/
  /*
   * Check input data for common errors
   */
  err=checks(nDataPoints, fXraw, fYraw, &dE);

  fMid=(fXraw[nDataPoints-1]+fXraw[0])/2.0;
  sEdge=get_Edge(sElement, fMid, &fEdge);
  if(!silent)printf("\nSpectrum over %s %s edge at theoretical energy of %8.2f eV\n", sElement, sEdge, fEdge);
  
  /**********************************
   * Determine Savitzky-Golay window
   **********************************/ 
  savwin(fEres, fEdge, dE, &nSavWin);

  /******************
   * Normalise data
   ******************/
  err=normalize(nDataPoints, fEdge, fXraw, fYraw, fYnorm, plotX, fYfita, fYfitb);

  /**************************************************
   * Impose on theoretical spectrum of f'' to obtain 
   * experimental equivalent
   **************************************************/
  if(verbose>0)printf(" Converting spectrum to f''\n");
  err=impose(nDataPoints, fEdge, fXraw, fYnorm, fYfpp);

  /*************************************************************************
   * Determine zeroth, first, second and third derivatives of smoothed data
   * and plot them on top of one another.
   *************************************************************************/
  err = smooth(nDataPoints, fYfpp, fYfpps, nSavWin, nSavWin, 4, 0);
  err = smooth(nDataPoints, fYfpp, fYDeriv1, nSavWin, nSavWin, 4, 1);
  err = smooth(nDataPoints, fYfpp, fYDeriv2, nSavWin, nSavWin, 4, 2);
  err = smooth(nDataPoints, fYfpp, fYDeriv3, nSavWin, nSavWin, 4, 3);

  if(verbose>2){
    for(i=0; i<nDataPoints; i++){
      printf("%f  %f  %f  %f \n", fYfpp[i], fYDeriv1[i], fYDeriv2[i], fYDeriv3[i]);
    }
  }

  /**********************************
   * Perform Kramer-Kronig transform
   **********************************/
  Integrate(nDataPoints, &nPoints, fEdge, fXraw, fXfpp, fYspline, fYfpps, fYDeriv1, fYDeriv2, fYDeriv3, fYfp);
  err=selwavel(nPoints, fXfpp, fYspline, fYfp);
  
  /*************************************
   * OUTPUT RESULTS (f' and f'' spectra)
   *************************************/
  if (outfile)
    // ASCII output
    err=my_efswrite(outfile, fXfpp, fYspline, fYfp, nPoints); 

  // create Python structure: ( Epeak, fppPeak, fpPeak, Einfl, fppInfl, fpInfl, ((X1, Yspline1, Yfp1),...,(Xn, Ysplinen, Yfpn)) )
  res = PyTuple_New(7);
  PyTuple_SetItem(res, 0, PyFloat_FromDouble(EPeak));
  PyTuple_SetItem(res, 1, PyFloat_FromDouble(fppPeak));
  PyTuple_SetItem(res, 2, PyFloat_FromDouble(fpPeak));
  PyTuple_SetItem(res, 3, PyFloat_FromDouble(EInfl));
  PyTuple_SetItem(res, 4, PyFloat_FromDouble(fppInfl));
  PyTuple_SetItem(res, 5, PyFloat_FromDouble(fpInfl));
  graph = PyTuple_New(nPoints);
  for (i=0; i<nPoints; i++) {
    graph_point = PyTuple_New(3);
    PyTuple_SetItem(graph_point, 0, PyFloat_FromDouble(fXfpp[i]));
    PyTuple_SetItem(graph_point, 1, PyFloat_FromDouble(fYspline[i]));
    PyTuple_SetItem(graph_point, 2, PyFloat_FromDouble(fYfp[i]));
    PyTuple_SetItem(graph, i, graph_point);
  }
  PyTuple_SetItem(res, 6, graph);

  /******************************
   * Print text table of results
   ******************************/
  if(!silent){
     printf("\n Table of results\n");
     printf("------------------------------------\n");
     printf("|      |  energy  |    f\'\' |   f\'  |\n");
     printf("| peak | %8.2f |  %5.2f | %5.2f |\n", EPeak, fppPeak, fpPeak);
     printf("| infl | %8.2f |  %5.2f | %5.2f |\n", EInfl, fppInfl, fpInfl);
     printf("------------------------------------\n");
  }

  return res;
}



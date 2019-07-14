#ifndef HEADER_H_INCLUDED
#define HEADER_H_INCLUDED

#include <iostream> // cout
#include <fstream>  // ofstream
#include <cmath>    // round
#include <time.h>

using namespace std;

void percentage(int a, int b);

void linspace(double a, double b, int N, double output[]);

bool test(double **yss, int ly, double **AT, int colat, double tol);

void RKF4(double IC[], int s, double params[], double t0, double t1, double h,
              void (*func)(double, double[], double[], double[]), double **y);

void pvi(double t, double y[], double param[], double ydot[]);

void fill_AT(double **y, int RowYStart, int RowYEnd, double **AT, int ColATStart, int ColATEnd);

double z1test(double x[], int size_x);

void cumsum(double x[], int size_x, double result[]);

double sum(double a[], int size_a);

double mean(double a[], int size_a);

double corr(double x[], double y[], int size_n);

double median(double x[], int size_x);

void my_sort(double x[], int size_x);
#endif // HEADER_H_INCLUDED

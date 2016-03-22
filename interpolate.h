#ifndef __INTERPOLATE_H__
#define __INTERPOLATE_H__


double spline(double *args, double *funcs, int n, double x, int &code);

double multi(double** matr, int sizeX, int sizeY, double x, int nX, double y,  int nY, int &code);

#endif
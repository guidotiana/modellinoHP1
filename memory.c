#include <stdlib.h>
#include <stdio.h>
#include "struct.h"
#include "io.h"

short *AlloShort(int n)
{
        short *x,i;

        x = (short *) calloc(n,sizeof(short));
        if (!x) Error("Cannot allocate int vector");

        for (i=0;i<n;++i) x[i]=0;

        return x;
}

int *AlloInt(int n)
{
        int *x,i;

        x = (int *) calloc(n,sizeof(int));
        if (!x) Error("Cannot allocate int vector");

        for (i=0;i<n;++i) x[i]=0;

        return x;
}

double *AlloDouble(int n)
{
        double *x;
        int i;

        x = (double *) calloc(n,sizeof(double));
        if (!x) Error("Cannot allocate double vector");

        for (i=0;i<n;++i) x[i]=0.;

        return x;
}

int **AlloIntMatrix(int l, int m)
{
  int **x;
  int i,j;

  x = (int **) calloc(l,sizeof(int *));
  if (!x) Error("Cannot allocate int matrix");

  for (i=0;i<l;++i)
    {
      *(x+i) = calloc(m,sizeof(int));
      if (!(*(x+i))) Error("Cannot allocate int matrix");
      for (j=0;j<m;++j) *(*(x+i)+j) = 0;
    }

  return x;
}

short **AlloShortMatrix(int l, int m)
{
  short **x;
  int i,j;

  x = (short **) calloc(l,sizeof(short *));
  if (!x) Error("Cannot allocate int matrix");

  for (i=0;i<l;++i)
    {
      *(x+i) = calloc(m,sizeof(short));
      if (!(*(x+i))) Error("Cannot allocate int matrix");
      for (j=0;j<m;++j) *(*(x+i)+j) = 0;
    }

  return x;
}

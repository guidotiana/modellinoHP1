#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

/*****************************************************************************
 Dumb generator of integer random numbers
 *****************************************************************************/
int irand(int r)
{
  int i;

  i=(int)(rand() / (((double)RAND_MAX + 1)/ r));
  return i;
}

/*****************************************************************************
 Dumb generator of double random numbers
 *****************************************************************************/
double frand()
{
  double q;

  q=(double) (rand()+1.0)/(RAND_MAX+1.0);
  return q;
}

/*****************************************************************************
 Initialize random number with seed (-1 to use computer time)
 *****************************************************************************/
long Randomize(int n)
{
  int i;

  if (n==-1)
  {
          n = (time(0) - 760650080) % 3600;
          fprintf(stderr,"Random seed = %d\n",n);
  }

  for (i=0;i<n;++i) irand(i);

  return n;
}

/*****************************************************************************
 Variate with power-law distribution
 *****************************************************************************/
double PowerLawVariate(double alpha, double min, double max)
{
  double r =  pow( ( pow(max,alpha+1.) - pow(min,alpha+1.) ) * frand() + pow(min,alpha+1.) , 1./(alpha+1.) );
  return r;
}

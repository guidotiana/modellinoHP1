#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <string.h>
#include "struct.h"
#include "io.h"
#include "rules.h"
#include "memory.h"
#include "random.h"



int main(int argc, char *argv[])
{

   int irun,istep,r,i,iStride,iMatStride, nContacts=0;
   short **neighbours;
   struct parm_s parms;
   short *methylation, *hp1, *nContactsSite;
   double *propensities;
   double time, tau;
   int stat[NREACTIONS];
   FILE *out,*out2,*fMatOut=NULL;

   fprintf(stderr, "pid = %ld\n", (long)getpid() );

   // Read parameters from file
   if (argc != 2) Error("Please specify input file");
   ReadParameters( argv[1], &parms );

   // Allocate stuff
   methylation = AlloShort( parms.length );
   hp1 = AlloShort( parms.length );
   nContactsSite = AlloShort( parms.length );
   propensities = AlloDouble( NREACTIONS );
   Randomize( parms.seed );
   neighbours = AlloShortMatrix(parms.length, parms.length);
   NameReactions(&parms);
   for (i=0;i<NREACTIONS;i++) stat[i]=0;

   // Output stuff
   if ( !strcmp(parms.fullOutput,"stdout") ) out = stdout;
   else if ( !strcmp(parms.fullOutput,"stderr") ) out = stderr;
   else out = fopen(parms.fullOutput,"w");

   if ( !strcmp(parms.fullOutput2,"stdout") ) out2 = stdout;
   else if ( !strcmp(parms.fullOutput2,"stderr") ) out2 = stderr;
   else out2 = fopen(parms.fullOutput2,"w");


   if ( strcmp(parms.nameMatOut,"") )
   {
     fMatOut = fopen(parms.nameMatOut,"w");
   }


   // Loop over independent runs
   for (irun=0;irun<parms.nrun;irun++)
   {
      InitialConditions( methylation, hp1, neighbours, parms, nContactsSite );

      time = 0.;
      iStride = 0;
      iMatStride = 0;

      // print stuff to output
      if (fMatOut!=NULL)
        PrintMatrix(neighbours, parms.length, fMatOut, irun, 0);


      // Time loop
      for (istep=0;istep<parms.nstep;istep++)
      {

         CalculatePropensities( propensities, methylation, hp1, &parms, neighbours, &nContacts, nContactsSite );

         tau = log( 1./ frand() ) / propensities[0];

         r = SelectReaction( propensities );
         stat[r]++;

         ApplyReaction( methylation, hp1, r, parms.length, neighbours,  parms, nContactsSite );
 
         time += tau;


         // print output 
         if (parms.debug) fprintf(stdout,"  tau = %lf \t r=%d \t time=%lf nContacts=%d \n",tau, r, time, CountContacts(neighbours,parms.length) );

         if (iStride >= parms.stride)
         {
           fprintf(stderr,"run=%3d step=%6d time=%lf nContacts=%d\n",irun,istep,time,nContacts);
           fprintf(out,"%3d %7.4lf ",irun,time);
           fprintf(out2,"%3d %7.4lf ",irun,time);
           for (i=0;i<parms.length;i++) fprintf(out,"%d",hp1[i]);
           for (i=0;i<parms.length;i++) fprintf(out2,"%d",methylation[i]);
           fprintf(out,"\n");
           fprintf(out2,"\n");
           iStride = 0;
           fflush(out);
           fflush(out2);
         }
         iStride++;

         if (iMatStride >= parms.matStride && fMatOut!=NULL)
         {
           PrintMatrix(neighbours, parms.length, fMatOut, irun, istep);
           fflush(fMatOut);
           iMatStride = 0;
         }
         iMatStride++;

         if (parms.debug) 
         {
             fprintf(stdout,"NEW methyl %3d %7.4lf ",irun,time);
             for (i=0;i<parms.length;i++) fprintf(stdout,"%d",methylation[i]);
             fprintf(stdout,"\n");
             fprintf(stdout,"NEW hp1 %3d %7.4lf ",irun,time);
             for (i=0;i<parms.length;i++) fprintf(stdout,"%d",hp1[i]);
             fprintf(stdout,"\n");

         }
      }

   }

   printf("Statistics reactions:\n");
   for (i=0;i<NREACTIONS;i++) printf("%3d %9d %s\n",i,stat[i],parms.reactionNames[i]);

   exit(0);
}

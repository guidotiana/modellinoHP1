#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "struct.h"
#include "io.h"
#include "random.h"
#include "memory.h"
#include "rules.h"

void RemoveSpaces(char* s) {
    char* d = s;
    do {
        while (*d == ' ') {
            ++d;
        }
    } while (*s++ == *d++);
}

void InitialConditions(short *c, short *h, short **neigh, struct parm_s p, short *nContactsSite)
{
   FILE *fin,*fin2;
   char aux[p.length+1];
   int i,j;

   for (i=0;i<p.length;i++) 
   {
       c[i]=0;
       h[i]=0;
       for (j=0;j<p.length;j++) neigh[i][j]=0;
       nContactsSite[i]=0;
   } 

   if ( strcmp(p.initialFile,"") )
   {
      if (p.debug) fprintf(stdout,"Reading initial condition met from file %s\n",p.initialFile);

      fin = fopen( p.initialFile, "r" );
      if (fin==NULL) Error("Cannot open initial condition file");
      if ( fscanf(fin,"%s",aux) != 1) Error("Cannot read initial-state file");
      for (i=0;i<p.length;i++) c[i]=(int) aux[i] - 48;
      fclose(fin);
   }

   if ( strcmp(p.initialFile2,"") )
   {
      if (p.debug) fprintf(stdout,"Reading initial condition hp1 from file %s\n",p.initialFile2);

      fin = fopen( p.initialFile2, "r" );
      if (fin==NULL) Error("Cannot open initial condition file");
      if ( fscanf(fin,"%s",aux) != 1) Error("Cannot read initial-state file");
      for (i=0;i<p.length;i++) h[i]=(int) aux[i] - 48;
      fclose(fin);
   }

   if ( !strcmp(p.initialMatFile,"RANDOM") )
      GenerateSpatialNeighbours( neigh, p.nLinks, p.length, p.type, p.alpha, p.debug );
   else if ( strcmp(p.initialMatFile,"") )
   {
      if (p.debug) fprintf(stdout,"Reading initial matrix from file %s\n",p.initialMatFile);

      fin2 = fopen( p.initialMatFile, "r" );
      for (j=0;j<p.length;j++)
      {
        if ( fscanf(fin2,"%s",aux) != 1) Error("Cannot read line in initial-matrix file");
        RemoveSpaces(aux);
        for (i=0;i<p.length;i++) neigh[j][i]=(int) aux[i] - 48;
      }

      fclose(fin2);
   }


  for (i=0;i<p.length;i++)
    for (j=0;j<p.length;j++)
      nContactsSite[i] += neigh[i][j];

}

void NameReactions(struct parm_s *p)
{
   strcpy(p->reactionNames[1],"randomMethylation");
   strcpy(p->reactionNames[2],"randomDemethylation");
   strcpy(p->reactionNames[3],"randomBindingHp1");
   strcpy(p->reactionNames[4],"randimUnbindingHp1");
   strcpy(p->reactionNames[5],"MethylationNeighbour");
   strcpy(p->reactionNames[6],"MethylationSpatial");
   strcpy(p->reactionNames[7],"ContactFormation");
   strcpy(p->reactionNames[8],"ContactBreakUn");
   strcpy(p->reactionNames[9],"ContactBreakMet");
   strcpy(p->reactionNames[10],"UnbindingHP1Hotspot");
   strcpy(p->reactionNames[11],"UnbindingHp1Methylated");
   strcpy(p->reactionNames[12],"UnbindingHp1MethylatedHot");

}

void CalculatePropensities( double *a, short *met, short *hp1, struct parm_s *p, short **m, int *nContacts, short *nContactsSite )
{
   int i,j, nMethylated=0, nUnmodified=0, nHp1Free=0, nHp1Occupied=0, nNeigh=0, nCanContact=0, nHp1OccupiedMet=0;
   int nSpatialNeigh=0, nContactBreakUn=0,  nContactBreakMethyl=0, nHp1Hotspots=0, isInContact=0, nMethylatedNoHp1=0, nMethHp1UnHotspot=0, cContact, cSpatial ;


   // Count modified states
   (*nContacts) = 0;
   for (i=0;i<p->length;i++)
   {    
      isInContact=0;
      cContact = 0;
      cSpatial = 0;

      for (j=0;j<p->length;j++) 
       if (i!=j)			// spatial neighbours that can be methylated
       {		        
           if ( m[i][j]>0 && met[i]==0 && met[j]>0 && hp1[j]>0 ) cSpatial = 1;
           if (nContactsSite[i] < p->valence && nContactsSite[j] < p->valence )
               if ( p->multipleContacts || m[i][j]==0 ) cContact = 1; 

           if ( m[i][j]>0 && ( met[i]==0 || met[j]==0 || hp1[i]==0 || hp1[j]==0) ) nContactBreakUn++;	// rule 8
           if ( m[i][j]>0 && met[i]>0 && met[j]>0 && hp1[i]>0 && hp1[j]>0) nContactBreakMethyl++;	// rule 9
           if  ( m[i][j]>0 ) { (*nContacts)++; isInContact=1; }
       }

      if ( cSpatial == 1 ) nSpatialNeigh++;				// rule 6
      if ( cContact == 1 ) nCanContact++;				// rule 7

      if ( met[i] > 0) nMethylated++;
      if ( met[i] > 0 && hp1[i]==0) nMethylatedNoHp1++;
      if ( met[i] == 0) nUnmodified++;
      if ( hp1[i] == 0) nHp1Free++;
      if ( hp1[i] > 0)
      { 
        if ( ( p->nHotspots==0 || p->hotspot[i]==0 ) && isInContact==0 ) // no hotspot or irrelevant
	{
           if (p->kMethHp1Unbind < SMALL || met[i]==0) nHp1Occupied++;	// rule 4
           else if ( met[i]>0 ) nHp1OccupiedMet ++;			// rule 11
	}
        else if ( p->hotspot[i]>0 && isInContact==0 ) 			// hotspot
        {
           if (p->kMethHp1Unbind==0 || met[i]==0) nHp1Hotspots++;	// rule 10
           else if ( met[i]>0 ) nMethHp1UnHotspot ++;			// rule 12
        }
      }
      if ( i>0 && i<p->length-1 ) if ( met[i]==0 && ( ( met[i-1]>0 && hp1[i-1]>0 ) ||  ( met[i+1]>0 && hp1[i+1]>0 ) ) ) nNeigh ++;
                                 
   }

   // (1) Random methylation
   a[1] = p->kRandMethyl * nUnmodified;

   // (2) Random demethylation
   if ( p->hp1CanDemeth==0 ) a[2] = p->kRandDemethyl * nMethylatedNoHp1;
   else a[2] = p->kRandDemethyl * nMethylated;

   // (3) Random binding of hp1
   if (p->maxHP1>0)
      a[3] = p->kRandHp1Bind * nHp1Free * ( 1. - nHp1Free / p->maxHP1 );
   else
      a[3] = p->kRandHp1Bind * nHp1Free;

   // (4) Random unbinding of hp1
   a[4] = p->kRandHp1Unbind * nHp1Occupied;

   // (5) Methylation of neighbour
   a[5] = p->kNeighMethyl *  nNeigh;

   // (6) Methylation of spatial neighbours
   a[6] = p->kSpatialNeighMethyl * nSpatialNeigh ;
 
   // (7) Contact formation
   a[7] = p->kContactFormation * nCanContact; 

   // (8) Contact breaking
   a[8] = p->kContactBreakingUn * nContactBreakUn ; 

   // (9) Contact breaking methylated
   a[9] = p->kContactBreakingMethyl * nContactBreakMethyl;

   // (10) Unbinding of hp1 from hotspots
   a[10] = p->kUnbindingHotspots * nHp1Hotspots;

   // (11) Unbinding of hp1 from methylated sites
   a[11] = p->kMethHp1Unbind * nHp1OccupiedMet;

   // (12) Unbinding of hp1 from methylated sites in hotspots
   a[12] = p->kMethHp1UnHotspot * nMethHp1UnHotspot;

   // (0) Sum of all propensities
   a[0] = 0;
   for (i=1;i<NREACTIONS;i++) a[0] += a[i];
 
 
   // print debug info
   if (p->debug) for (i=0;i<NREACTIONS;i++) 
   {
       fprintf(stdout,"  a[%d]=%12.6lf\t",i,a[i]);
       fprintf(stdout,"\n");
   }
}



int SelectReaction( double *propensities )
{
   int i;
   double r,aSum=0.;

   r = propensities[0] * frand();

   for (i=1;i<NREACTIONS;i++)
   {
     aSum += propensities[i];
     if ( r < aSum ) return i;
   }

   Error("All propensities equal to zero");
   return -1;
}



void ApplyReaction(  short *met, short *hp1, int r, int length, short  **m, struct parm_s parms, short *nContactsSite )
{
   int i,j,ok=0, chk=0, isInContact=0;


   while ( ok==0 )
   {
       i = irand( length );
       if (parms.debug) fprintf(stdout,"  ApplyReaction: i=%d r=%d ok=%d %s\n",i,r,ok,parms.reactionNames[r]);

       // (1) Random methylation
       if ( r==1 )
       {
         if ( met[i] == 0 ) { met[i]=1; ok=1; }
       }

       // (2) Random demethylation
       else if ( r==2 )
       {
         if (parms.hp1CanDemeth)
         { 
            if ( met[i] >0 ) { met[i]=0; ok=1; }
         }
         else
            if ( hp1[i]==0 && met[i] >0 ) { met[i]=0; ok=1; }
       }

       // (3)  Random binding of hp1
       else if ( r==3 )
       { 
         if ( hp1[i] ==0 ) { hp1[i]=1; ok=1; }
       }

       // (4)  Random unbinding of hp1
       else if ( r==4 )
       { 
         isInContact = 0;
         for (j=0;j<parms.length;j++) 
            if ( i!=j && m[i][j]>0 ) isInContact = 1;
         if ( hp1[i] >0 ) 
            if ( ( parms.nHotspots==0 || parms.hotspot[i]==0) && isInContact==0 ) 
                if (parms.kMethHp1Unbind==0 || met[i]==0) { hp1[i]=0; ok=1; }
       }

       // (5) Methylation of neighbour
       else if ( r== 5)
       { 
         if (i>0) if ( met[i] == 0 && met[i-1] >0 && hp1[i-1]>0 ) { met[i]=1; ok=1; }
         if (i<length-1)  if ( met[i] == 0 && met[i+1] >0 && hp1[i+1]>0 ) { met[i]=1; ok=1; }
       }

       // (6) Methylation of spatial neighbours
       else if ( r==6 )
       {
               if ( met[i]==0 ) for (j=0;j<length;j++) if ( i!=j && m[i][j]>0 && met[j]>0 && hp1[j]>0 )
               {
                  met[i]=1;
                  ok=1;
                  break;
               }
       }
       // (7) Contact formation
      else if ( r==7 )
      {
          if (nContactsSite[i] < parms.valence)			// methylated and below valence
          {
             j = GenerateContact(i, m, parms.alpha, length, parms.debug);          // choose binder from power-law on effective distance

             if (j>=0 && j<length && j != i)
             {
                if (nContactsSite[j]<parms.valence)
                  if  (parms.multipleContacts || m[i][j]==0 )
                  {
                    m[i][j] ++;
                    m[j][i] ++;
                    nContactsSite[i]++;
                    nContactsSite[j]++;
                    ok = 1;
                  }
             }
          } 
      }

      // (8) Contact breaking unmethylated
      else if ( r==8 )
      {
         j = irand(length);
         if ( m[i][j] >0 && ( met[i]==0 || met[j]==0 || hp1[i]==0 || hp1[j]==0) )
         {
              m[i][j] --;
              m[j][i] --;
              nContactsSite[i]--;
              nContactsSite[j]--;
              ok = 1;
         }
      }

     // (9) contact breaking methylated
     else if ( r==9 )
     {
         j = irand(length);
         if ( i!=j && met[i]>0 && met[j]>0 && m[i][j]>0 && hp1[i]>0 && hp1[j]>0 )
         {
              m[i][j] --;
              m[j][i] --;
              nContactsSite[i]--;
              nContactsSite[j]--;
              ok = 1;
         }
     }

     // (10) unbinding of Hp1 in hotspots
     else if ( r==10 )
     {  
         isInContact = 0;
         for (j=0;j<parms.length;j++) 
            if ( i!=j && m[i][j]>0 ) isInContact = 1;

         if ( hp1[i] >0 ) 
            if (parms.hotspot[i]>0 && isInContact==0 ) 
              if (parms.kMethHp1Unbind==0 || met[i]==0) { hp1[i]=0; ok=1; }
     }

     // (11) Unbinding of hp1 from methylated sites
     else if ( r==11 )
     {  
         isInContact = 0;
         for (j=0;j<parms.length;j++) 
            if ( i!=j && m[i][j]>0 ) isInContact = 1;

         if ( hp1[i] >0 && met[i]>0) 
            if ( (  parms.nHotspots==0  || parms.hotspot[i]==0) && isInContact==0 ) { hp1[i]=0; ok=1; }
     }

     // (12)  Random unbinding of hp1 from methylated hotspots
      else if ( r==12 )
      { 
         isInContact = 0;
         for (j=0;j<parms.length;j++) 
            if ( i!=j && m[i][j]>0 ) isInContact = 1;
         if ( hp1[i] >0 && parms.hotspot[i]>0 && met[i]>0 && isInContact==0) { hp1[i]=0; ok=1; }
      }

     else Error("Undefined reaction in ApplyReaction");

     // Avoid infinite loops
     chk++;
     //if (chk>9999*length) Error("Infinite loop in ApplyReaction");
     if (chk>999*length) 
     { 
        printf("Stuck in r=%d i=%d  %s\n",r,i,parms.reactionNames[r]);
        break;
     }
   }

   return;
}


void GenerateSpatialNeighbours( short **m, int nLinks, int length, char type[50], double alpha, int debug )
{
   int i,j,iLink=0;
   double sign;

   if (debug) fprintf(stdout,"Generate spatial neighbours:\n");

   if (!strcmp(type,"none")) return ;

   do
   {

      if (!strcmp(type,"flat"))
      {
         i = irand(length);
         j = irand(length);
      }

      else if (!strcmp(type,"powerlaw"))
      {
         i = irand(length);
         sign =  2. * (frand() - 0.5 ) ;
         if (sign>0) j = (int) (pow( ( pow(length-1,alpha+1.) - pow(i+1,alpha+1.) ) * frand() + pow(i+1,alpha+1.) , 1./(alpha+1.) ));
         else j = (int) ( pow( ( - pow(i,alpha+1) ) * frand() + pow(i,alpha+1) , 1./(alpha+1.) ));
      }

      else if (!strcmp(type,"restrFlat"))
      {
         if ( frand() < 0.5 )
           i = irand( (int) (length/4) ) + (int) (length/4);
         else
           i = irand( (int) (length/4) ) + (int) (3*length/4);
         j = i + irand( (int) (length/4) ) - (int) (length/8);
      }

      else Error("Type not defined in GenerateSpatialNeighbours");
       
      if (debug) fprintf(stdout,"i=%3d j=%3d iLink=%d\n",i,j,iLink);

      if ( i != j && i>0 && j>0 && i<length && j<length)
      {
         m[i][j] = 1;
         m[j][i] = 1;
         iLink ++;
      }

   } while (iLink<nLinks);


   return;
}

int GenerateContact(int i, short **m, double alpha, int length, int debug)
{
   int j,dl;
   double p,powerlaw;

   while (1)						// generate power-law variate by rejection
   {
     j = irand(length);
     dl = EffectiveDistance(i, j, m, length);			// calculate distance along the network
 
     p = frand();
     powerlaw = pow(dl, alpha);

     if (debug) fprintf(stdout,"  GenerateContact: i=%d j=%d dl=%d p=%lf powerlaw=%lf\n",i,j,dl,p,powerlaw);

     if (i != j && p < powerlaw) return j;

   }

}

// One would like to calculate the polymeric distance through the network of links (polymeric or through-space)
// This is lengthy, so the effective distance neglects links within loops
int EffectiveDistance(int i, int j, short **m, int length)
{
   int d=0,ok=0,next=i,k,n,dout;
 
   if (i>j) { k=i; i=j; j=k; next=i; }	// it must be i<j
                           
   do				// walk from i towards j
   {
     ok=0;
     for (k=j;k>next;k--)	// look for links between next and the closest to j
       if (m[next][k]==1) 
       {
          next = k;
          ok=1;
          break;
       }
     if (ok==0) next++;		// if there is no link, use polymeric link
       
     d++;
     
   } while (next<j);
  
   for (k=i;k>i-d;k--)		// look for outer loops
     for (n=j;n<j+d;n++)
       if (k>0 && n<length)
       {
         dout = i-k+n-j;
         if ( m[n][k]>0 && dout<d ) return dout;
       }
   return d;
}

int CountContacts( short **m, int length)
{
   int i,j,n=0;

   for (i=0;i<length;i++)
     for (j=0;j<length;j++) n += m[i][j];

  return n;
}

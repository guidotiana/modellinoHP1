#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "struct.h"
#include "io.h"
#include "memory.h"

void Error(char *text)
{
  fprintf(stderr,"* ERROR: %s\n",text);
  exit(1);
}

void ReadParS(char *s, char key[20], char *par)
{
        int l,r;

        l = strlen(key);
        if (!strncmp(s,key,l))
        {
                r = sscanf(s,"%*s %s",par);
                if (r<1) { fprintf(stderr,"ERROR: Cannot read %s in parameter file\n",key); exit(1); }
                fprintf(stderr,"%s = %s\n",key,par);
        }
}

void ReadParN(char *s, char key[20], int *par)
{
        int l;

        l = strlen(key);
        if (!strncmp(s,key,l))
        {
                *par=1;
                fprintf(stderr,"%s\n",key);
        }
}

void ReadParD(char *s, char key[20], int *par)
{
        int l,r;

        l = strlen(key);
        if (!strncmp(s,key,l))
        {
                r = sscanf(s,"%*s %d",par);
                if (r<1) { fprintf(stderr,"ERROR: Cannot read %s in parameter file (r=%d)\n",key,r); exit(1); }
                fprintf(stderr,"%s = %d\n",key,*par);
        }

}

void ReadParL(char *s, char key[20], long *par)
{
        int l,r;

        l = strlen(key);
        if (!strncmp(s,key,l))
        {
                r = sscanf(s,"%*s %ld",par);
                if (r<1) { fprintf(stderr,"ERROR: Cannot read %s in parameter file (r=%d)\n",key,r); exit(1); }
                fprintf(stderr,"%s = %ld\n",key,*par);
        }

}

void ReadParF(char *s, char key[20], double *par)
{
        int l,r;

        l = strlen(key);
        if (!strncmp(s,key,l))
        {
                r = sscanf(s,"%*s %lf",par);
                if (r<1) { fprintf(stderr,"ERROR: Cannot read %s in parameter file\n",key); exit(1); }
                fprintf(stderr,"%s = %lf\n",key,*par);
        }
}

void ReadParameters(char *filename, struct parm_s *parm)
{
  char aux[500],keyword[100];
  FILE *fp;
  int readHotspots=0,i;

  // open file for reading
  fp = fopen(filename,"r");
  if (fp==NULL)
    {
          fprintf(stderr,"\nCannot open file %s for reading.\n",filename);
          exit(1);
    }

  // defaults
  parm->nstep = 1;
  parm->nrun = 1;
  parm->length = 100;
  parm->seed = -1;
  parm->debug = 0;
  parm->kRandMethyl = 0.;
  parm->kRandDemethyl = 0.;
  parm->kRandHp1Bind = 0.;
  parm->kRandHp1Unbind = 0.;

  parm->kNeighMethyl = 0.;
  parm->kSpatialNeighMethyl = 0.;
  strcpy(parm->fullOutput, "stdout");
  strcpy(parm->fullOutput2, "stdout");
  parm->stride = 1;
  parm->matStride = 1;
  parm->nLinks = 0;
  parm->alpha = -1.5;
  strcpy(parm->type, "flat");
  parm->valence = 99999;
  parm->kContactFormation=0.;
  parm->kContactBreakingUn=0.;
  strcpy(parm->initialFile2,"");
  strcpy(parm->initialFile,"");
  strcpy(parm->initialMatFile,"");
  parm->nHotspots = 0;  
  parm->seed=-1;
  parm->kContactBreakingMethyl=0.;
  parm->maxHP1=-1;
  parm->multipleContacts=0;
  parm->kUnbindingHotspots=0;
  parm->hp1CanDemeth=0;
  parm->kMethHp1Unbind=0;
  parm->kMethHp1UnHotspot=0;
  parm->hotspot = NULL;



  // scan parameter file
  while(fgets(aux,500,fp)!=NULL)
  {
    ReadParD(aux,"nstep",&(parm->nstep));
    ReadParD(aux,"length",&(parm->length));
    ReadParD(aux,"nrun",&(parm->nrun));
    ReadParD(aux,"seed",&(parm->seed));
    ReadParN(aux,"debug",&(parm->debug));
    ReadParD(aux,"stride",&(parm->stride));
    ReadParD(aux,"matStride",&(parm->matStride));
    ReadParS(aux,"outputHp1",(parm->fullOutput));
    ReadParS(aux,"outputMet",(parm->fullOutput2));
    ReadParS(aux,"nameMatOut",(parm->nameMatOut));
    ReadParF(aux,"kRandMethyl",&(parm->kRandMethyl));
    ReadParF(aux,"kRandDemethyl",&(parm->kRandDemethyl));
    ReadParF(aux,"kRandHp1Bind",&(parm->kRandHp1Bind));
    ReadParF(aux,"kRandHp1Unbind",&(parm->kRandHp1Unbind));
    ReadParF(aux,"kMethHp1Unbind",&(parm->kMethHp1Unbind));
    ReadParF(aux,"kNeighMethyl",&(parm->kNeighMethyl));
    ReadParF(aux,"kSpatialNeighMethyl",&(parm->kSpatialNeighMethyl));
    ReadParF(aux,"kContactBreakingUn",&(parm->kContactBreakingUn));
    ReadParF(aux,"kContactBreakingMethyl",&(parm->kContactBreakingMethyl));
    ReadParS(aux,"type",(parm->type));
    ReadParD(aux,"nLinks",&(parm->nLinks));
    ReadParF(aux,"alpha",&(parm->alpha));
    ReadParD(aux,"valence",&(parm->valence));
    ReadParD(aux,"maxHP1",&(parm->maxHP1));
    ReadParF(aux,"kContactFormation",&(parm->kContactFormation));
    ReadParF(aux,"kContactBreakingUn",&(parm->kContactBreakingUn));
    ReadParF(aux,"kUnbindingHotspot",&(parm->kUnbindingHotspots));
    ReadParS(aux,"initialFileMet",(parm->initialFile));
    ReadParS(aux,"initialFileHp1",(parm->initialFile2));
    ReadParS(aux,"initialMatFile",(parm->initialMatFile));
    ReadParN(aux,"multipleContacts",&(parm->multipleContacts));
    ReadParN(aux,"hp1CanDemeth",&(parm->hp1CanDemeth));
    ReadParF(aux,"kMethHp1UnHotspot",&(parm->kMethHp1UnHotspot));
  
    // read hotspot list 
    if (FindKeyword(aux, keyword) && !strcmp(keyword,"endhotspots")) readHotspots=0;
    if (readHotspots==1)
    {
       if (parm->hotspot==NULL) parm->hotspot = AlloShort(parm->length);

       if ( sscanf(aux,"%d",&i)!=1)  { fprintf(stderr,"ERROR: Cannot read hotpost in parameter file\n"); exit(1); }
       if (i>=parm->length)  { fprintf(stderr,"ERROR: Hotspot not contained in the chain\n"); exit(1); }
       parm->hotspot[ i ] = 1;
       parm->nHotspots++;

    }
    if (FindKeyword(aux, keyword) && !strcmp(keyword,"hotspots") ) readHotspots=1;

  }

  if (parm->nHotspots==0) parm->hotspot = AlloShort(parm->length);

  fprintf(stderr,"Read %d hotspots.\n",parm->nHotspots);
}


void PrintMatrix(short **m, int l, FILE *fout, int iRun, int iStep)
{
   int i,j;

   fprintf(fout,"# %d %d\n",iRun,iStep);

   for (i=0;i<l;i++)
   {
     for (j=0;j<l;j++) fprintf(fout,"%d",m[i][j]);
     fprintf(fout,"\n"); 
   }
}

int FindKeyword(char *string, char *keyword)
{
   int i=0,key=0,k=0;
   char c;

   do {
                   c = string[i];
                   if (c=='[') key=1;
                   if (c==']') key=0;
                   if (key==1 && c!=' ' && c!='[')
                   {
                           keyword[k] = c;
                           ++k;
                   }
                   ++i;

           } while ( c != '\0' && i<500 );

        keyword[k] = '\0';

        if (!strcmp(keyword,"")) return 0;
        else return 1;
}

// Number of reactions, including no reactions (i.e., should be +1)
#define NREACTIONS 13

// Structure for input parameters
struct parm_s
{
   int length;
   int nrun;
   int nstep;
   int debug;
   int seed;
   int stride;
   int matStride;
   int nLinks;
   char type[50];		// type of distribution for spatial links
   double alpha;		// power-law exponent for spatial links
   char fullOutput[50];
   char fullOutput2[50];
   char nameMatOut[50];		// file to write contact map
   double kRandMethyl;
   double kRandDemethyl;
   double kRandHp1Bind;
   double kRandHp1Unbind;
   double kNeighMethyl;
   double kUnbindingHotspots;
   double kSpatialNeighMethyl;
   double kContactFormation;
   int valence;
   double kContactBreakingUn;
   char initialFile[50];
   char initialFile2[50];
   char initialMatFile[50];
   int nHotspots;
   short *hotspot;
   char reactionNames[NREACTIONS][50];
   double kContactBreakingMethyl;
   int maxHP1;
   int multipleContacts;
   int hp1CanDemeth;
   double kMethHp1Unbind;
   double kMethHp1UnHotspot;
};

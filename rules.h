#define SMALL 1E-14

void InitialConditions(short *c, short *h, short **, struct parm_s p, short *nContactsSite);
void CalculatePropensities( double *a, short *methyl, short *hp1, struct parm_s *p, short **m , int *nContacsti, short *nContactsSite);
int SelectReaction( double *propensities );
void ApplyReaction( short *met, short *hp1, int r, int length, short **m, struct parm_s parms, short *nContactsSite);
void GenerateSpatialNeighbours( short **m, int nLinks, int length, char *type, double alpha, int debug );
int GenerateContact(int i, short **m, double alpha, int length, int debug);
int EffectiveDistance(int i, int j, short **m, int length);
int CountContacts( short **m, int length);
void NameReactions(struct parm_s *);







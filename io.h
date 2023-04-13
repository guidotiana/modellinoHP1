void ReadParS(char *s, char key[20], char *par);
void ReadParN(char *s, char key[20], int *par);
void ReadParD(char *s, char key[20], int *par);
void ReadParL(char *s, char key[20], long *par);
void ReadParF(char *s, char key[20], double *par);
void Error(char *text);
void ReadParameters(char *filename, struct parm_s *parm);
void PrintMatrix(short **m, int l, FILE *fout, int iRun, int iStep);
int FindKeyword(char *string, char keyword[30]);




/* Includes
*/
#include <stdio.h>
#include <string.h>
#include <ctype.h>

/***************************************************************************/
/* Globals
*/
static int  mdm_score[25][25];
static char mdm_aalist[30];

/***************************************************************************/
main(int argc, char ** argv)
{
   int score;
   
   if(argc==3)
   {
      ReadMDM("amdata:mdm78.mat");
      score = calcscore(toupper(argv[1][0]),toupper(argv[2][0]));
      printf("Match %c %c scores %d\n",argv[1][0],argv[2][0],score);
   }
   else
   {
      printf("Usage: checkmat <seq> <seq>\n");
   }
}


/***************************************************************************/
/*>int ReadMDM(char *mdmfile)
   --------------------------
   Read mutation data matrix into static global arrays
   07.10.92 Original
*/
int ReadMDM(char *mdmfile)
{
   FILE *mdm = NULL;
   int  i, j;
   char ch;

   if((mdm=fopen(mdmfile,"r"))==NULL)
   {
      return(0);
   }

   /* Read in the mutation data matrix */
   for(i=0; i<25; i++)
   {
      for(j=0; j<25; j++)
      {
         fscanf(mdm,"%d",&(mdm_score[i][j]));
      }
   }
   for(i=0; i<25; i++)
   {
      while(ch=getc(mdm),(ch==10)||(ch==13)||(ch==' '));
      mdm_aalist[i]=ch;
   }
   
   fclose(mdm);
   
   return(1);
}



/***************************************************************************/
/*>int calcscore(char resa, char resb)
   -----------------------------------
   Calculate score from mutation data matrix
   07.10.92 Adapted from NIMR-written original
*/
static int calcscore(char resa,
                     char resb)
{
   int i,j;

   for(i=0;i<25;i++)
   {
      if(resa==mdm_aalist[i]) break;
   }
   if(i==25) 
   {
      printf("Residue %c not found in matrix\n",resa);
   }
   for(j=0;j<25;j++)
   {
      if(resb==mdm_aalist[j]) break;
   }
   if(j==25) 
   {
      printf("Residue %c not found in matrix\n",resb);
   }
   return(mdm_score[i][j]);
}                               

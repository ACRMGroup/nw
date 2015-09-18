/*************************************************************************

   Program:    nw
   File:       nw.c
   
   Version:    V3.14
   Date:       11.03.15
   Function:   Do Needleman & Wunsch sequence alignment
   
   Copyright:  (c) Dr. Andrew C. R. Martin / UCL 1990-2015
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   EMail:      andrew@bioinf.org.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! 

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============

   A simple Needleman & Wunsch Dynamic Programming alignment of 2
   sequences.  The sequences are read from PIR files specified on the
   command line.  By default, the Dayhoff mutation matrix is used to
   score the alignment, a gap penalty of 10 is used.  These defaults
   may be changed by the use of command line switches.  The program is
   limited to 1000 residue sequences (this may be changed by using the
   -l flag).  A window is not used so the program may be a bit slow on
   long sequences.

   The algorithm is described in:

   Needleman, Saul B. and Wunsch, Christian D. (1970)
   `A General Method Applicable to the Search for Similarities in the
   Amino Acid Sequence of Two Proteins'
   JMB, 48, 443--453

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  19.06.90   Original
   V2.0  07.10.92   Changed to use routines from the library
   V2.1  03.11.95   Changed for improved routines from library such as 
                    ReadPIR().
   V3.0  07.11.95   Complete rewrite with lots of new features.
   V3.1  20.11.95   Calc. of max chain length wasn't being done
                    if only one chain!
   V3.2  21.11.95   Fixed a few warnings under gcc
   V3.3  02.05.96   Reports normalised scores and percent identity
   V3.4  24.06.96   Now reports 3 types of % identity
   V3.5  01.07.96   Added 4th identity type
   V3.6  02.07.96   Fixed -q handling and added -qq
   V3.7  11.07.96   Added percentage homologies
   V3.8  01.11.96   Made the output easier to parse. Added % id against
                    first sequence length
   V3.9  04.03.97   Prints number of INDEL residues
   V3.10 21.01.98   -i flag was causing error messages in calculation of
                    max possible score for alignment. No longer prints
                    homologies when using an identity matrix
   V3.11 28.09.00   Modified to allow user specified extension penalty
                    and to call affinealine() rather than aline()
   V3.12 09.06.08   Improved error messages; checks for zero length seqs
                    Changed default gap penalties to 10/2 rather than 
                    5/0
   V3.13 23.08.10   Added -s - Merged in from home version
   V3.14 11.03.15   Added -x - show matches in alignment

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "bioplib/macros.h"
#include "bioplib/seq.h"
#include "bioplib/general.h"

/************************************************************************/
/* Defines and macros
*/
#define VERSION  "3.14"
#define MDMFILE  "mdm78.mat"
#define DEF_GAPPEN 10
#define DEF_EXTPEN 2
#define MAXCHAIN 16
#define MAXBUFF  160

#define ALIGN_NONE 0
#define ALIGN_PIR  1

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, BOOL *identity, int *GapPenalty, 
                  int *ExtPenalty, BOOL *verbose, char *mdmfile, 
                  char *infile1, char *infile2, int *quiet,
                  int *AlignStyle, char *AlignFile, BOOL *ScoreOnly,
                  BOOL *showMatches);
void Usage(void);
BOOL DoAlignment(FILE *in1, FILE *in2, BOOL identity, int GapPenalty, 
                 int ExtPenalty, char *mdmfile, BOOL verbose, int quiet,
                 int AlignStyle, FILE *out, BOOL showMatches);
void WritePIRAlignment(FILE *out, SEQINFO SeqInfo1, SEQINFO SeqInfo2,
                       STRINGLIST *Alignments1, STRINGLIST *Alignments2);
int CalcNumId(char *seq1, char *seq2, int length);
int CalcNumAligned(char *align1, char *align2, int align_len);
int FindAlnLenNoTail(char *align1, char *align2, int align_len);
int CalcIDScore(char *seq1, char *seq2, BOOL identity);
void GetIndelInfo(char *align1, char *align2, 
                  int *ndel, int *ndelg, int *nins, int *ninsg);
BOOL ReadAlignmentAndScore(FILE *in);
int FindSeqLenNoGaps(char *seq);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for N&W sequence alignment

   07.10.92 Original based on NIMR version. Extracted some routines.
            Various changes and improvements.
   08.10.92 Improved comments. Added check on identity when printing mdm
            in verbose mode.
   07.11.95 Completely rewritten.
   02.07.96 Changed quiet to int so when = 2 is completely silent
   09.06.08 Added "Error: " to error message
   09.06.08 Changed default gap penalties to something more sensible
   23.08.10 Added ScoreOnly
   11.03.15 Added showMatches
*/
int main(int argc, char **argv)
{
   BOOL identity,
        verbose,
        ScoreOnly = FALSE,
        showMatches = FALSE;
   int  GapPenalty,
        ExtPenalty,
        AlignStyle,
        quiet;
   char mdmfile[MAXBUFF],
        infile1[MAXBUFF],
        infile2[MAXBUFF],
        AlignFile[MAXBUFF];
   FILE *in1, *in2, *out;

   
   /* Set default values                                                */
   identity   = FALSE;
   quiet      = 0;
   verbose    = FALSE;
   GapPenalty = DEF_GAPPEN;
   ExtPenalty = DEF_EXTPEN;
   AlignStyle = ALIGN_NONE;
   strcpy(mdmfile,MDMFILE);

   if(ParseCmdLine(argc, argv, &identity, &GapPenalty, &ExtPenalty,
                   &verbose, mdmfile, infile1, infile2, &quiet, 
                   &AlignStyle, AlignFile, &ScoreOnly, &showMatches))
   {
      if(ScoreOnly)
      {
         if((in1=fopen(infile1,"r"))==NULL)
         {
            fprintf(stderr,"Unable to open input file %s.\n",infile1);
            return(1);
         }
         ReadAlignmentAndScore(in1);
      }
      else
      {
         /* Open the input PIR files                                    */
         if((in1=fopen(infile1,"r"))==NULL)
         {
            fprintf(stderr,"Error: Unable to open input file %s.\n",
                    infile1);
            return(1);
         }
         if((in2=fopen(infile2,"r"))==NULL)
         {
            fprintf(stderr,"Error: Unable to open input file %s.\n",
                    infile2);
            return(1);
         }
         if(AlignStyle != ALIGN_NONE)
         {
            if((out=fopen(AlignFile,"w"))==NULL)
            {
               fprintf(stderr,"Error: Unable to open output file %s.\n",
                       AlignFile);
               return(1);
            }
         }
         
         DoAlignment(in1, in2, identity, GapPenalty, ExtPenalty, mdmfile, 
                     verbose, quiet, AlignStyle, out, showMatches);
      }
   }
   else
   {
      Usage();
   }
   
   return(0);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, BOOL *identity, 
                     int *GapPenalty, int *ExtPenalty, BOOL *verbose, 
                     char *mdmfile, char *infile1, char *infile2, 
                     int *quiet, int *AlignStyle, char *AlignFile,
                     BOOL *ScoreOnly, BOOL *showMatches)
   ----------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  BOOL   *identity    Use an identity matrix
            int    *GapPenalty  Specified gap penalty
            int    *ExtPenalty  Specified extension penalty
            BOOL   *verbose     Switch on verbose mode in alignment
            char   *mdmfile     Mutation data matrix file
            char   *infile1     First sequence file
            char   *infile2     Second sequence file
            int    *quiet       Do not display alignment, only score
            int    *AlignStyle  ALIGN_PIR = Create a PIR alignment file
            char   *AlignFile   Filename for output alignment file
            BOOL   *ScoreOnly   Score an existing alignment
            BOOL   *showMatches Show matches in alignment
   Returns: BOOL                Success?

   Parse the command line
   
   07.11.94 Original    By: ACRM
   02.07.96 Changed quiet to int so when = 2 is completely silent
            Handles -qq
   28.09.00 Added ExtPenalty
   09.06.08 Changed handling of defaults for ExtPenalty for -i
   23.08.10 Added ScoreOnly
   11.03.15 Added showMatches
*/
BOOL ParseCmdLine(int argc, char **argv, BOOL *identity, int *GapPenalty, 
                  int *ExtPenalty, BOOL *verbose, char *mdmfile, 
                  char *infile1, char *infile2, int *quiet,
                  int *AlignStyle, char *AlignFile, BOOL *ScoreOnly,
                  BOOL *showMatches)
{
   BOOL UserGapPenalty = FALSE;
   BOOL UserExtPenalty = FALSE;
   
   argc--;
   argv++;

   if(argc < 2)
      return(FALSE);

   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'v':
            *verbose = TRUE;
            break;
         case 'i':
            if(!UserGapPenalty)
               *GapPenalty = 1;
            if(!UserExtPenalty)
               *ExtPenalty = 0;
            *identity = TRUE;
            break;
         case 'g':
            argc--;
            argv++;
            UserGapPenalty = TRUE;
            if(argc < 0)
               return(FALSE);
            if((sscanf(argv[0], "%d", GapPenalty)) == 0)
               return(FALSE);
            break;
         case 'x':
            argc--;
            argv++;
            UserExtPenalty = TRUE;
            if(argc < 0)
               return(FALSE);
            if((sscanf(argv[0], "%d", ExtPenalty)) == 0)
               return(FALSE);
            break;
         case 'm':
            argc--;
            argv++;
            if(argc < 0)
               return(FALSE);
            strcpy(mdmfile, argv[0]);
            break;
         case 'q':
            *quiet = 1;
            if(argv[0][2] == 'q')
               *quiet = 2;
            break;
         case 'p':
            *AlignStyle = ALIGN_PIR;
            argc--;
            argv++;
            if(argc < 0)
               return(FALSE);
            strcpy(AlignFile, argv[0]);
            break;
         case 's':
            *ScoreOnly = TRUE;
            break;
         case 'd':
            *showMatches = TRUE;
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are 2 arguments left                       */
         if(argc != 2)
            return(FALSE);
         
         /* Copy them to the output variables                           */
         strcpy(infile1, argv[0]);
         strcpy(infile2, argv[1]);
         
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
   Prints a usage message

   07.11.95 Original    By: ACRM
   21.11.95 Wasn't printing the MDMFILE variable
   02.07.96 Added -qq (silent mode)
   28.09.00 Added -x
   23.08.10 Added -s
   11.03.15 Added -
*/
void Usage(void)
{
   fprintf(stderr,"\nNW %s (c) 1990-2010 Dr. Andrew C.R. Martin, \
NIMR/SciTech Software/UCL/Reading\n", VERSION);

   fprintf(stderr,"\nUsage: nw [-g n][-x n][-i][-m <matrix>][-v]\
[-q[q]][-p <file>][-d] <file1> <file2>\n");
   fprintf(stderr," -or-  nw -s <aligned-file>\n");
   fprintf(stderr,"       -g <n>      Specify the gap penalty\n");
   fprintf(stderr,"                   [Default: %d for matrix or 1 for \
identity matrix]\n", DEF_GAPPEN);
   fprintf(stderr,"       -x <n>      Specify the gap extension \
penalty\n");
   fprintf(stderr,"                   [Default: %d for matrix or 0 for \
identity matrix]\n",DEF_EXTPEN);
   fprintf(stderr,"       -i          Use an identity matrix\n");
   fprintf(stderr,"       -m <matrix> Specify the mutation matrix \
[Default: %s]\n",MDMFILE);
   fprintf(stderr,"       -v          Turn on verbose mode\n");
   fprintf(stderr,"       -q          Quiet mode (gives only the \
score)\n");
   fprintf(stderr,"       -qq         Silent mode (no printed \
information; use with -p)\n");
   fprintf(stderr,"       -d          Display matches in alignment\n");
   fprintf(stderr,"       -p <file>   Write a PIR format alignment \
file\n");
   fprintf(stderr,"       <file1>     First PIR sequence file\n");
   fprintf(stderr,"       <file2>     Second PIR sequence file\n");
   fprintf(stderr,"       -s          Score an existing alignment \
(in PIR format)\n");
   
   fprintf(stderr,"\nNW is a simple Needleman and Wunsch sequence \
alignment program taking\n");
   fprintf(stderr,"PIR format input files. Only the first sequence in \
each file is aligned,\n");
   fprintf(stderr,"but multiple chains will be handled.\n\n");
   fprintf(stderr,"Note that the default gap penalties have changed in \
V3.12\n\n");
}


/************************************************************************/
/*>BOOL DoAlignment(FILE *in1, FILE *in2, BOOL identity, int GapPenalty, 
                    int ExtPenalty, char *mdmfile, BOOL verbose, 
                    int quiet, int AlignStyle, FILE *out, 
                    BOOL showMatches)
   ---------------------------------------------------------------------
   Main routine which does the alignment work. Reads the sequence files
   and calls the alignment code for each pair of chains in turn.
   Prints header and information as appropriate.

   07.11.95 Original    By: ACRM
   20.11.95 Moved initialisation of max lengths outside check on
            multiple chains!
   21.11.95 Corrected void return
   02.05.96 Added calculation and printing of normalised scores and
            percent identities
   24.06.96 Added % identities over shorter sequence and over aligned
            residues. Thanks to Alex May for pointing out there are
            multiple interpretations of % identity :-)
   02.07.96 Changed quiet to int so when = 2 is completely silent
            Fixed handling of quiet mode
   11.07.96 Calculates % 'homology'
   01.11.96 Made the output easier to parse. Added % id against
            first sequence length
   28.09.00 Changed to call affinealign(), ExtPenalty passed in as
            a parameter. Moved call to GetIndelInfo() so it works
            properly with -q. Changed 'homology' to 'similarity'
   09.06.08 Added "Error: " to error message
   11.03.15 Updated for BiopLib bl prefix
            Added showMatches code
*/
BOOL DoAlignment(FILE *in1, FILE *in2, BOOL identity, int GapPenalty, 
                 int ExtPenalty, char *mdmfile, BOOL verbose, int quiet,
                 int AlignStyle, FILE *out, BOOL showMatches)
{
   char       *seq1[MAXCHAIN],              /* Sequences to align       */
              *seq2[MAXCHAIN],
              *align1, *align2;             /* Aligned sequences        */
   BOOL       error, punct;
   SEQINFO    SeqInfo1,
              SeqInfo2;
   int        nchain1,
              nchain2,
              ai, aj,
              i,  j,
              len1, len2,
              chain,
              offset,
              score,
              align_len,
              maxlen1, maxlen2,
              count,
              TotalScore = 0,
              NumId,
              NumAligned,
              AlnLenNoTail,
              TotalNumId,
              TotalLength = 0,
              TotalMinLength = 0,
              TotalAlnLength = 0,
              TotalSeq1Length = 0,
              TotalAlnLenNoTail = 0,
              IDScore = 0,
              TotalIDScore = 0,
              nins, ninsg,
              ndel, ndelg,
              TotalNIns  = 0,
              TotalNInsG = 0,
              TotalNDel  = 0,
              TotalNDelG = 0;
   STRINGLIST *Alignments1 = NULL,
              *Alignments2 = NULL;
   
   /* Read sequences                                                    */
   nchain1 = blReadPIR(in1,FALSE,seq1,MAXCHAIN,&SeqInfo1,&punct,&error);
   if(error)
   {
      fprintf(stderr,"Error: Unable to read PIR sequence file. \
No memory\n");
      return(FALSE);
   }
   
   nchain2 = blReadPIR(in2,FALSE,seq2,MAXCHAIN,&SeqInfo2,&punct,&error);
   if(error)
   {
      fprintf(stderr,"Error: Unable to read PIR sequence file. \
No memory\n");
      return(FALSE);
   }

   maxlen1 = strlen(seq1[0]);
   maxlen2 = strlen(seq2[0]);

   /* Find the maximum length chain in each sequence file               */
   if(nchain1 > 1)
   {
      for(i=1; i<nchain1; i++)
      {
         if(strlen(seq1[i]) > maxlen1)
            maxlen1 = strlen(seq1[i]);
      }
   }
   if(nchain2 > 1)
   {
      for(i=1; i<nchain2; i++)
      {
         if(strlen(seq2[i]) > maxlen2)
            maxlen2 = strlen(seq2[i]);
      }
   }

   /* Allocate memory for the alignment results                         */
   align1 = (char *)malloc((maxlen1+maxlen2) * sizeof(char));
   align2 = (char *)malloc((maxlen1+maxlen2) * sizeof(char));
   if(align1==NULL || align2==NULL)
   {
      fprintf(stderr,"Error: No memory for alignment storage\n");
      return(FALSE);
   }


   /* Print header information                                          */
   if(!quiet)
   {
      printf("\nNeedleman and Wunsch Sequence Alignment Program   \
%s\n",VERSION);
      printf("====================================================\
==\n");
      printf("Copyright Andrew C.R. Martin, NIMR/SciTech Software/UCL \
1990-1998\n");
      
      printf("\nParameters for this run\n");
      printf("-----------------------\n");
      printf("Mutation Matrix:      %s\n",(identity) ?"Identity":mdmfile);
      printf("Gap Penalty:          %d\n",GapPenalty);
      printf("Extension Penalty:    %d\n",ExtPenalty);
      printf("Verbose mode:         %s\n",(verbose)  ?"ON"      :"OFF");

      if((nchain1 > 1 || nchain2 > 1) && (nchain1 != nchain2))
      {
         printf("\nNumber of chains differs; the first %d will be \
aligned\n\n", MIN(nchain1, nchain2));
      }
   
      /* Display sequences                                              */
      printf("\nSequence 1\n----------\n");
      printf("Maximum sequence length:  %d\n",maxlen1);
      printf("Sequence identifier:      %s\n",SeqInfo1.code);
      printf("Sequence name:            %s\n",SeqInfo1.name);
      printf("Sequence source:          %s\n\n",SeqInfo1.source);
      for(chain = 0; chain<nchain1; chain++)
      {
         count=1;
         for(i=0;i<strlen(seq1[chain]);i++,count++)
         {
            if(count>80)
            {
               count=1;
               printf("\n");
            }
            printf("%c",seq1[chain][i]);
         }
         printf("*\n");
      }
      
      printf("\nSequence 2\n----------\n");
      printf("Maximum sequence length:  %d\n",maxlen2);
      printf("Sequence identifier:      %s\n",SeqInfo2.code);
      printf("Sequence name:            %s\n",SeqInfo2.name);
      printf("Sequence source:          %s\n\n",SeqInfo2.source);
      for(chain = 0; chain<nchain2; chain++)
      {
         count=1;
         for(i=0;i<strlen(seq2[chain]);i++,count++)
         {
            if(count>80)
            {
               count=1;
               printf("\n");
            }
            printf("%c",seq2[chain][i]);
         }
         printf("*\n");
      }

      printf("\nAlignment...\n------------\n");
   }
   
   
   /* Read mutation data matrix                                         */
   if(!identity)
   {
      if(!blReadMDM(mdmfile))
      {
         fprintf(stderr,"Error: Unable to read mutation matrix: %s\n",
                 mdmfile);
         for(chain=0; chain<nchain1; chain++)
            free(seq1[chain]);
         for(chain=0; chain<nchain2; chain++)
            free(seq2[chain]);
         return(FALSE);
      }
   }

   /* Do the alignments                                                 */
   for(chain=0; chain<MIN(nchain1, nchain2); chain++)
   {
      len1 = strlen(seq1[chain]);
      len2 = strlen(seq2[chain]);

      if((len1==0) || (len2==0))
      {
         fprintf(stderr,"Error: %s sequence is of zero length\n",
                 ((len1==0)?"First":"Second"));
         return(FALSE);
      }
      
      score = blAffinealign(seq1[chain], len1, seq2[chain], len2, 
                            verbose, identity, GapPenalty, ExtPenalty,
                            align1, align2, &align_len);
      if(!score)
      {
         fprintf(stderr,"Error: No memory for alignment matrix\n");
         return(FALSE);
      }
      /* Calculate various scores                                       */
      TotalScore        += score;
      IDScore            = CalcIDScore(seq1[chain], seq2[chain],
                                       identity);
      TotalIDScore      += IDScore;
      TotalLength       += align_len;
      TotalMinLength    += MIN(len1, len2);
      TotalSeq1Length   += len1;
      NumId              = CalcNumId(align1, align2, align_len);
      TotalNumId        += NumId;
      NumAligned         = CalcNumAligned(align1, align2, align_len);
      TotalAlnLength    += NumAligned;
      AlnLenNoTail       = FindAlnLenNoTail(align1, align2, align_len);
      TotalAlnLenNoTail += AlnLenNoTail;

      align1[align_len] = align2[align_len] = '\0';

      /* Store the sequences for later                                  */
      if(AlignStyle != ALIGN_NONE)
      {
         Alignments1 = StoreString(Alignments1, align1);
         Alignments2 = StoreString(Alignments2, align2);
      }

      /* Display the alignment                                          */
      if(!quiet)
      {
         if(nchain1 > 1 || nchain2 > 1)
         {
            printf("Chain %d\n",chain+1);
         }
         
         offset = 0;
         /* This loop prints seqa                                       */
         for(i=0,ai=0,aj=0; ai<align_len; ai++)
         {
            i++;
            /* If we've printed 80 chars, we print the equiv section of 
               seqb  
            */
            if(i>80)
            {
               i=1;
               printf("\n");
               if(showMatches)
               {
                  for(j=offset; j<80+offset; j++)
                  {
                     if(align1[j] == align2[j])
                        fputc('|', stdout);
                     else if(blCalcMDMScore(align1[j], align2[j]) > 0)
                        fputc('.', stdout);
                     else
                        fputc(' ', stdout);
                  }
                  
                  printf("\n");
               }
               
               for(j=offset; j<80+offset; j++) printf("%c",align2[j]);
               printf("\n\n");
               offset += 80;
            }
            printf("%c",align1[ai]);
         }
         /* Now print the remains of seqb                               */
         printf("\n");
         if(showMatches)
         {
            for(j=offset; j<align_len; j++)
            {
               if(align1[j] == align2[j])
                  fputc('|', stdout);
               else if(blCalcMDMScore(align1[j], align2[j]) > 0)
                  fputc('.', stdout);
               else
                  fputc(' ', stdout);
            }
            
            printf("\n");
         }
         for(j=0+offset; j<align_len; j++) printf("%c",align2[j]);
         printf("\n\n");

      }  /* if(!quiet)                                                  */

      if(quiet < 2)
      {
         GetIndelInfo(align1, align2, &ndel, &ndelg, &nins, &ninsg);
         TotalNDel  += ndel;
         TotalNDelG += ndelg;
         TotalNIns  += nins;
         TotalNInsG += ninsg;

         if(quiet)
            printf("CHAIN %d\n",chain+1);
            
         /* align() in identity mode scores 2 for a match, so we divide 
            the score by 2
         */
         printf("Score                                    SCORE: \
%d\n",score/(identity?2:1));
         printf("Score (normalised by alignment length)  NSCORE: \
%.2f\n",       (REAL)score/((REAL)align_len * ((identity?2.0:1.0))));
         if(!identity)
            printf("Percentage similarity                    HOMOL: \
%.2f%%\n",     (REAL)100.0 * (REAL)score/(REAL)IDScore);
         printf("Identity over alignment length         IDALLEN: \
%.2f%%\n",     (REAL)100.0 * (REAL)NumId / (REAL)align_len);
         printf("Identity over alignment with no tails IDNOTAIL: \
%.2f%%\n",     (REAL)100.0 * (REAL)NumId / (REAL)AlnLenNoTail);
         printf("Identity over shorter sequence         IDSHORT: \
%.2f%%\n",     (REAL)100.0 * (REAL)NumId / (REAL)MIN(len1,len2));
         printf("Identity over aligned residues          IDLONG: \
%.2f%%\n",     (REAL)100.0 * (REAL)NumId / (REAL)NumAligned);
         printf("Identity over sequence 1               IDFIRST: \
%.2f%%\n",     (REAL)100.0 * (REAL)NumId / (REAL)len1);
         printf("Deletions in sequence 1: residues(sites) NDELS: \
%d (%d)\n", ndel, ndelg);
         printf("Insertions in sequence 1: residues(sites) NINS: \
%d (%d)\n", nins, ninsg);
      }
   }  /* For each chain                                                 */

   if(quiet < 2)
   {
      /* Display the total score                                        */
      if(nchain1 > 1 || nchain2 > 1)
      {
         printf("Total score                                  \
   SCORETOT: %d\n",         TotalScore/(identity?2:1));
         printf("Total score (normalised by alignment length) \
  NSCORETOT: %.2f\n",
                (REAL)TotalScore/
                ((REAL)TotalLength * ((identity?2.0:1.0))));
         if(!identity)
            printf("Total percentage similarity                  \
   HOMOLTOT: %.2f%%\n",
                (REAL)100.0 * (REAL)TotalScore/(REAL)TotalIDScore);
         printf("Total Identity over alignment length         \
 IDALLENTOT: %.2f%%\n",
                (REAL)100.0 * (REAL)TotalNumId / (REAL)TotalLength);
         printf("Identity over alignment with no tails        \
IDNOTAILTOT: %.2f%%\n",
                (REAL)100.0 * (REAL)TotalNumId / (REAL)TotalAlnLenNoTail);
         printf("Total Identity over shorther sequence        \
 IDSHORTTOT: %.2f%%\n",
                (REAL)100.0 * (REAL)TotalNumId / (REAL)TotalMinLength);
         printf("Total Identity over aligned residues         \
  IDLONGTOT: %.2f%%\n\n",
                (REAL)100.0 * (REAL)TotalNumId / (REAL)TotalAlnLength);
         printf("Total Identity over first sequence           \
 IDFIRSTTOT: %.2f%%\n\n\n",
                (REAL)100.0 * (REAL)TotalNumId / (REAL)TotalSeq1Length);
         printf("Total deletions in sequence 1: residues(sites)\
     NDELS: %d (%d)\n", TotalNDel, TotalNDelG);
         printf("Total insertions in sequence 1: residues(sites)\
     NINS: %d (%d)\n", TotalNIns, TotalNInsG);
      }
   }

   switch(AlignStyle)
   {
   case ALIGN_PIR:
      WritePIRAlignment(out,SeqInfo1,SeqInfo2,Alignments1,Alignments2);
      break;
   default:
      break;
   }
   
   /* Free up memory used to store alignment strings                    */
   if(AlignStyle)
   {
      FreeStringList(Alignments1);
      FreeStringList(Alignments2);
   }
    
   /* Free allocated memory                                             */
   for(chain = 0; chain<nchain1; chain++)
      free(seq1[chain]);
   for(chain = 0; chain<nchain2; chain++)
      free(seq2[chain]);
   free(align1);
   free(align2);

   return(TRUE);
}


/************************************************************************/
/*>void WritePIRAlignment(FILE *out, SEQINFO SeqInfo1, SEQINFO SeqInfo2,
                          STRINGLIST *Alignments1, 
                          STRINGLIST *Alignments2)
   ---------------------------------------------------------------------
   Writes the alignment in PIR sequence alignment format.

   07.11.95 Original    By: ACRM
   09.06.08 Added "Error: " to error message
*/
void WritePIRAlignment(FILE *out, SEQINFO SeqInfo1, SEQINFO SeqInfo2,
                       STRINGLIST *Alignments1, STRINGLIST *Alignments2)
{
   STRINGLIST *a;
   int        count,
              i,
              length;

   if(Alignments1 == NULL || Alignments2 == NULL)
   {
      fprintf(stderr,"Error: No alignments stored. Unable to write PIR \
alignment file.\n");
      return;
   }
   
   fprintf(out,">P1;%s\n",SeqInfo1.code);
   fprintf(out,"%s - %s\n",SeqInfo1.name, SeqInfo1.source);
   for(a=Alignments1; a!=NULL; NEXT(a))
   {
      length=strlen(a->string);
      count=1;
      for(i=0;i<length;i++,count++)
      {
         if(count>60)
         {
            count=1;
            fprintf(out,"\n");
         }
         fprintf(out,"%c",a->string[i]);
      }
      fprintf(out,"*\n");
   }
   fprintf(out,"\n");

   fprintf(out,">P1;%s\n",SeqInfo2.code);
   fprintf(out,"%s - %s\n",SeqInfo2.name, SeqInfo2.source);
   for(a=Alignments2; a!=NULL; NEXT(a))
   {
      length=strlen(a->string);
      count=1;
      for(i=0;i<length;i++,count++)
      {
         if(count>60)
         {
            count=1;
            fprintf(out,"\n");
         }
         fprintf(out,"%c",a->string[i]);
      }
      fprintf(out,"*\n");
   }
   fprintf(out,"\n");
}


/************************************************************************/
/*>int CalcNumId(char *seq1, char *seq2, int length)
   -------------------------------------------------
   Input:   char    *seq1     First sequence
            char    *seq2     Second sequence
            int     length    Alignment length
   Returns: int               Number of identities

   Calculates the number of identities between two sequences. Doesn't
   count ??, XX or -- as an identity

   02.05.96 Original   By: ACRM
*/
int CalcNumId(char *seq1, char *seq2, int length)
{
   int numid = 0,
       i;
   
   for(i=0; i<length; i++)
   {
      if(seq1[i] == seq2[i])
      {
         if(seq1[i] != '?' && seq1[i] != 'X' && seq1[i] != '-')
            numid++;
      }
   }
   
   return(numid);
}


/************************************************************************/
/*>int CalcNumAligned(char *align1, char *align2, int align_len)
   -------------------------------------------------------------
   Calculates the number of residues which were aligned. i.e. misses
   out all those which were aligned against a gap.

   24.06.96 Original   By: ACRM
*/
int CalcNumAligned(char *align1, char *align2, int align_len)
{
   int numaln = 0,
       i;
   
   for(i=0; i<align_len; i++)
   {
      if(align1[i] != '-' && align2[i] != '-')
         numaln++;
   }
   
   return(numaln);
}


/************************************************************************/
/*>int FindAlnLenNoTail(char *align1, char *align2, int align_len)
   ---------------------------------------------------------------
   Calculate the alignment length without any tails
   i.e. ---XXXXXX
        XXXXXXXXX
   gives a length of 6 rather than the input alignment length of 9

   01.07.96 Original   By: ACRM
*/
int FindAlnLenNoTail(char *align1, char *align2, int align_len)
{
   int NTail = 0,
       j;

   for(j=0, NTail=0; j<align_len; j++)
   {
      if(align1[j] != '-' && align2[j] != '-')
         break;
      NTail++;
   }

   for(j=align_len-1; j>=0; j--)
   {
      if(align1[j] != '-' && align2[j] != '-')
         break;
      NTail++;
   }

   return(align_len-NTail);
}


/************************************************************************/
/*>int CalcIDScore(char *seq1, char *seq2, BOOL identity)
   ------------------------------------------------------
   Calculates the maximum possible score resulting from the identical
   sequence

   11.07.96 Original   By: ACRM
   21.01.98 Added identity flag
   11.03.15 Updated for BiopLib bl prefix
*/
int CalcIDScore(char *seq1, char *seq2, BOOL identity)
{
   int score1 = 0,
       score2 = 0,
       i,
       seqlen1 = strlen(seq1),
       seqlen2 = strlen(seq2);

   if(identity)
   {
      score1 = seqlen1;
      score2 = seqlen2;
   }
   else
   {
      for(i=0; i<seqlen1; i++)
      {
         if(isalpha(seq1[i]))
            score1 += blCalcMDMScore(seq1[i], seq1[i]);
      }
      for(i=0; i<seqlen2; i++)
      {
         if(isalpha(seq2[i]))
            score2 += blCalcMDMScore(seq2[i], seq2[i]);
      }
   }
   
   return(MIN(score1, score2));
}


/************************************************************************/
/*>void GetIndelInfo(char *align1, char *align2, 
                     int *ndel, int *ndelg, int *nins, int *ninsg)
   ---------------------------------------------------------------
   Input:     char      *align1      Sequence for counting inserts
              char      *align2      Sequence for counting deletes
   Output:    int       *ndel        Number of deletes in align1
                        *ndelg       Number of delete groups in align1
                        *nins        Number of deletes in align2
                                     (inserts in align1)
                        *ninsg       Number of delete groups in align2
                                     (insert groups in align1)

   04.03.97 Original   By: ACRM
   06.03.97 Initialised ndelg and ninsg (Oops!)
   11.03.15 Updated for BiopLib bl prefix
*/
void GetIndelInfo(char *align1, char *align2, 
                  int *ndel, int *ndelg, int *nins, int *ninsg)
{
   int i;
   char prev;
   
   *ndel = blCountchar(align1,'-');
   *nins = blCountchar(align2,'-');
   *ndelg = 0;
   *ninsg = 0;

   for(i=0, prev = '*'; align1[i]; i++)
   {
      if(align1[i]=='-' && prev != '-')
         (*ndelg)++;
      prev = align1[i];
   }
   for(i=0, prev = '*'; align2[i]; i++)
   {
      if(align2[i]=='-' && prev != '-')
         (*ninsg)++;
      prev = align2[i];
   }
}

/************************************************************************/
/*>BOOL ReadAlignmentAndScore(FILE *in)
   ------------------------------------
   Input:    FILE        *in       PIR file pointer - file contains two
                                   aligned sequences with headers
   Returns:  BOOL                  Success?

   28.03.03 Original    By: ACRM
*/
BOOL ReadAlignmentAndScore(FILE *in)
{
   char       *align1[MAXCHAIN],*align2[MAXCHAIN];
   SEQINFO    SeqInfo;
   int        nchain1, nchain2, align_len, NumId, NumAligned,
              AlnLenNoTail, len1, len2;
   BOOL       punct, error;

   nchain1 = ReadPIR(in,TRUE,align1,MAXCHAIN,&SeqInfo,&punct,&error);
   nchain2 = ReadPIR(in,TRUE,align2,MAXCHAIN,&SeqInfo,&punct,&error);
   align_len = strlen(align1[0]);

   if((nchain1!=1) || (nchain2!=1) || (align_len!=strlen(align2[0])))
   {
      fprintf(stderr,"Error: number of chains in each sequence must be \
one and aligned\n");
      fprintf(stderr,"sequences must be of the same length\n");
      return(FALSE);
   }

   NumId              = CalcNumId(align1[0], align2[0], align_len);
   NumAligned         = CalcNumAligned(align1[0], align2[0], align_len);
   AlnLenNoTail       = FindAlnLenNoTail(align1[0], align2[0], align_len);
   len1               = FindSeqLenNoGaps(align1[0]);
   len2               = FindSeqLenNoGaps(align2[0]);

   printf("Identity over alignment length         IDALLEN: \
%.2f%%\n",     (REAL)100.0 * (REAL)NumId / (REAL)align_len);
   printf("Identity over alignment with no tails IDNOTAIL: \
%.2f%%\n",     (REAL)100.0 * (REAL)NumId / (REAL)AlnLenNoTail);
   printf("Identity over shorter sequence         IDSHORT: \
%.2f%%\n",     (REAL)100.0 * (REAL)NumId / (REAL)MIN(len1,len2));
   printf("Identity over aligned residues          IDLONG: \
%.2f%%\n",     (REAL)100.0 * (REAL)NumId / (REAL)NumAligned);
   printf("Identity over sequence 1               IDFIRST: \
%.2f%%\n",     (REAL)100.0 * (REAL)NumId / (REAL)len1);

   return(TRUE);
}

/************************************************************************/
/*>int FindSeqLenNoGaps(char *seq)
   -------------------------------
   Input:    char       *seq       Sequence containing - gap characters
   Returns:  int                   Length excluding gaps

   28.03.03 Original    By: ACRM
*/
int FindSeqLenNoGaps(char *seq)
{
   char *chp;
   int  len = 0;

   for(chp=seq; *chp; chp++)
   {
      if(*chp != '-') len++;
   }
   
   return(len);
}


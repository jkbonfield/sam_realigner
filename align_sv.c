/* A PACKAGE FOR SEQUENCE COMPARISON WITH AFFINE WEIGHTS */
/* Here we maximize the similarity score and won't penalize the first
   and last gaps */

/* Globally passed params and macros */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
// #include "align.h"
// #include "uascii.gbl"
// #include "align_lib_old.h"

typedef int FastInt;
#define xmalloc malloc
#define xfree free
#define vmessage printf
#define MAX(a,b) ((a)>(b)?(a):(b))

static FastInt (*w)[128];				/* w = W */
static FastInt g, h, m;				/* g = G, h = H, m = g+h */

#define gap(k)  ((k) <= 0 ? 0 : g+h*(k))	/* k-symbol indel cost */

static FastInt *sapp;				/* Current script append ptr */
static FastInt  last;				/* Last script op appended */
static FastInt sl;

char base_val[128];

void init_base_val(void) {
    int i;
    for (i=0; i<128; i++)
	base_val[i] = 5;
    
    base_val['A'] = 0;
    base_val['a'] = 0;
    base_val['C'] = 1;
    base_val['c'] = 1;
    base_val['G'] = 2;
    base_val['g'] = 2;
    base_val['T'] = 3;
    base_val['t'] = 3;
    base_val['U'] = 3;
    base_val['u'] = 3;
    base_val['*'] = 4;
}

/* Turn vector into single (most likely) nucleotide */
char consen_6(FastInt B[6]) {
  FastInt a = B[0], b = 0;

  if (B[1] > a)
    a = B[1], b = 1;
  if (B[2] > a)
    a = B[2], b = 2;
  if (B[3] > a)
    a = B[3], b = 3;
  if (B[4] > a)
    a = B[4], b = 4;
  if (B[5] > a)
    a = B[5], b = 5;

  return a == 0 ? '-' : "ACGT*-"[b];
}

						/* Append "Delete k" op */
#define DEL(k)				\
{ sl += k;				\
  if (last < 0)				\
    last = sapp[-1] -= (k);		\
  else					\
    last = *sapp++ = -(k);		\
}
						/* Append "Insert k" op */
#define INS(k)				\
{ sl += k;				\
  if (last > 0)				\
    last = sapp[-1] += (k);		\
  else					\
    last = *sapp++ = (k);		\
}

#define REP { ++sl; last = *sapp++ = 0; }	/* Append "Replace" op */

#define NMAX 6000

static FastInt (*CD)[2];	/* Forward cost-only vectors */
static FastInt (*RS)[2];	/* Reverse cost-only vectors */
static char     *A2;
static FastInt (*B2)[6];

/* align(A,B,M,N,tb,te) returns the cost of an optimum conversion between
   A[1..M] and B[1..N] that begins(ends) with a delete if tb(te) is zero
   and appends such a conversion to the current script.                   */
/* topr, botr, lc, rc are used to trace the boundary lines */

static FastInt align(A,B,M,N,tb,te,topr,botr,lc,rc)
    char *A;
    FastInt (*B)[6];
    FastInt M, N;
    FastInt tb, te; char topr, botr, lc, rc;

{        FastInt   midi, midj, type;	/* Midpoint, type, and cost */
         FastInt midc;

{ register FastInt   i, j;
  register FastInt c, e, d, s;
           FastInt t, wa;
/* Boundary cases: M <= 1 or N == 0 */

  if (N <= 0)
    { if (M > 0) DEL(M)
      if (topr || botr) return 0;
      else return -gap(M);
    }
  if (M <= 1) {
      if (M <= 0)
        { INS(N);
	  if (topr || botr) return 0;
          else return -gap(N);
        }
      if (topr) {
	 midc = rc ? 0 : te-h;
	 midj = 0;
	 wa = A[1];
	 for (j = 1; j <= N; j++) {
	     c = B[j][wa] - gap(N-j);
             if (c > midc) {
		 midc = c;
                 midj = j;
               }
           }
      } else if (botr) {
	 if (lc) midc = 0;
	 else midc = tb-h;
	 midj = 0;
	 wa = A[1];
	 for (j = 1; j <= N; j++)
	   { c = -gap(j-1) + B[j][wa];
             if (c > midc)
               { midc = c;
                 midj = j;
               }
           }
      } else {
         if (tb < te) tb = te;
	 if (lc || rc) midc = -gap(N);
         else midc = (tb-h) - gap(N);
         midj = 0;
         wa = A[1];
         for (j = 1; j <= N; j++)
           { c = -gap(j-1) + B[j][wa] - gap(N-j);
             if (c > midc)
               { midc = c;
                 midj = j;
               }
           }
      }
      if (midj == 0)
        { DEL(1) INS(N) }
      else
        { if (midj > 1) INS(midj-1)
          REP
          if (midj < N) INS(N-midj)
        }
      return midc;
    }

/* Divide: Find optimum midpoint (midi,midj) of cost midc */

  midi = M/2;			/* Forward phase:                          */
  CD[0][0] = 0;			/*   Compute C(M/2,k) & D(M/2,k) for all k */
  if (topr) {
	for (j = 1; j <= N; j++)
	  { CD[j][0] = 0;
	    CD[j][1] = -g;
	  }
  } else {
     	t = -g;
     	for (j = 1; j <= N; j++)
       	{ CD[j][0] = t = t-h;
       	  CD[j][1] = t-g;
        }
  }
  t = tb;
  for (i = 1; i <= midi; i++)
    { s = CD[0][0];
      if (lc) {
	CD[0][0] = c = 0;
	e = -g;
      } else {
        CD[0][0] = c = t = t-h;
        e = t-g;
      }
      wa = A[i];
      for (j = 1; j <= N; j++) {
	  if ((c   - m) > (e =   e   - h)) e = c - m;
	  if ((j == N) && rc) {
             if ((c = CD[j][0]) > (d = CD[j][1])) d = c;
	  } else {   
             if ((c = CD[j][0] - m) > (d = CD[j][1] - h)) d = c;
	  }
          c = s + B[j][wa];
          if (e > c) c = e;
          if (d > c) c = d;
          s = CD[j][0];
          CD[j][0] = c;
          CD[j][1] = d;
      }
    }
  CD[0][1] = CD[0][0];

  RS[N][0] = 0;			/* Reverse phase:                          */
  				/*   Compute R(M/2,k) & S(M/2,k) for all k */
  if (botr) {
	for (j = N-1; j >= 0; j--)
	  { RS[j][0] = 0;
	    RS[j][1] = -g;
          }
  } else {
  	t = -g;
  	for (j = N-1; j >= 0; j--)
    	{ RS[j][0] = t = t-h;
      	  RS[j][1] = t-g;
    	}
  }
  t = te;
  for (i = M-1; i >= midi; i--)
    { s = RS[N][0];
      if (rc) {
	RS[N][0] = c = 0;
	e = -g;
      } else {
      	RS[N][0] = c = t = t-h;
      	e = t-g;
      }
      wa = A[i+1];
      for (j = N-1; j >= 0; j--)
        { if ((c   - m) > (e =   e   - h)) e = c - m;
	  if ((j == 0) && lc) {
             if ((c = RS[j][0]) > (d = RS[j][1])) d = c;
	  } else {
             if ((c = RS[j][0] - m) > (d = RS[j][1] - h)) d = c;
	  }
          c = s + B[j+1][wa];
          if (e > c) c = e;
          if (d > c) c = d;
          s = RS[j][0];
          RS[j][0] = c;
          RS[j][1] = d;
        }
    }
  RS[N][1] = RS[N][0];

  midc = CD[0][0]+RS[0][0];		/* Find optimal midpoint */
  midj = 0;
  type = 1;
  for (j = 0; j <= N; j++)
    if ((c = CD[j][0] + RS[j][0]) >= midc)
      if (c > midc || CD[j][0] != CD[j][1] && RS[j][0] == RS[j][1])
        { midc = c;
          midj = j;
        }
  if (rc) {
    if ((c = CD[N][1] + RS[N][1]) > midc)
      { midc = c;
        midj = N;
        type = 2;
      }
  } else {
    if ((c = CD[N][1] + RS[N][1] + g) > midc)
      { midc = c;
        midj = N;
        type = 2;
      }
  }
  for (j = N-1; j > 0; j--)
    if ((c = CD[j][1] + RS[j][1] + g) > midc)
      { midc = c;
        midj = j;
        type = 2;
      }
  if (lc) {
    if ((c = CD[0][1] + RS[0][1]) > midc)
      { midc = c;
        midj = 0;
        type = 2;
      }
  } else {
    if ((c = CD[0][1] + RS[0][1] + g) > midc)
      { midc = c;
        midj = 0;
        type = 2;
      }
  }
}

/* Conquer: recursively around midpoint */

  if (midj == 0 || midj == N) {
     if (type == 1)
       { align(A,B,midi,midj,tb,-g,topr,0,lc,rc);
         align(A+midi,B+midj,M-midi,N-midj,-g,te,0,botr,lc,rc);
       }
     else
       { align(A,B,midi-1,midj,tb,0,topr,0,lc,rc);
         DEL(2);
         align(A+midi+1,B+midj,M-midi-1,N-midj,0,te,0,botr,lc,rc);
       }
  } else {
     if (type == 1)
       { align(A,B,midi,midj,tb,-g,topr,0,lc,0);
         align(A+midi,B+midj,M-midi,N-midj,-g,te,0,botr,0,rc);
       }
     else
       { align(A,B,midi-1,midj,tb,0,topr,0,lc,0);
         DEL(2);
         align(A+midi+1,B+midj,M-midi-1,N-midj,0,te,0,botr,0,rc);
       }
  }
  return midc;
}

#ifdef never_used

/* CHECK_SCORE - return the score of the alignment stored in S */

static FastInt CHECK_SCORE(A,B,M,N,S,EG)
    char A[];
    FastInt (*B)[6];
    FastInt M, N;
    FastInt S[];
    char EG;
{ 
  register FastInt   i,  j, op;
  FastInt score;

  score = i = j = op = 0;
  while (i < M || j < N) {
	op = *S++;
	if (EG == 1 && i == 0 && j == 0 && op != 0) {
		if (op > 0) j = j+op;
		else i = i-op;
	} else if (EG == 1 && (i == M || j == N)) {
		i = M;
		j = N;
	} else if (op == 0) 
		score = w[A[++i]][consen_6(B[++j])] + score;
	else if (op > 0) {
		score = score - (g+op*h);
		j = j+op;
	} else {
		score = score - (g-op*h);
		i = i-op;
	}
  }
  return(score);
}
#endif

/* Accuracy - in hundreths */

#define ACC (100)

#define MM -4
#define MC 4

static int tmp[6][6] = {
    {MC*ACC, MM*ACC, MM*ACC, MM*ACC, -1*ACC,  1*ACC},
    {MM*ACC, MC*ACC, MM*ACC, MM*ACC, -1*ACC,  1*ACC},
    {MM*ACC, MM*ACC, MC*ACC, MM*ACC, -1*ACC,  1*ACC},
    {MM*ACC, MM*ACC, MM*ACC, MC*ACC, -1*ACC,  1*ACC},
    {-1*ACC, -1*ACC, -1*ACC, -1*ACC, MC*ACC, -1*ACC},
    { 1*ACC,  1*ACC,  1*ACC,  1*ACC, -1*ACC, MC*ACC}};

static char *nuascii = "ACGT*-";

/* Interface and top level of comparator */
/* EG specifies if we want to penalize the end gaps */

FastInt align_sv(A,B,M,N,low,up,W,G,H,S,s1,s2,e1,e2)
    char A[];
    FastInt (*B)[6];
    FastInt M,N;
    FastInt W[][128],G,H;
    FastInt S[];
    FastInt low,up;
    FastInt s1,s2,e1,e2;
{ 
  FastInt c;
  /* FastInt t; */
  /* FastInt ck; */
  /* char EG = 1; */
  FastInt i,j;

  CD = xmalloc(sizeof(FastInt) * 2 * (N+1));
  RS = xmalloc(sizeof(FastInt) * 2 * (N+1));
  A2 = xmalloc(sizeof(char)    * 1 * (N+1));
  B2 = xmalloc(sizeof(FastInt) * 6 * (N+1));

  if (!CD || !RS || !A2 || !B2)
      return -1;

  A--; B--;

  w = W;			/* Setup global parameters */
  g = G*ACC;
  h = H*ACC;
  m = g+h;
  sapp = S;
  last = 0;
  sl = 0;

  /*
   * Transform ascii sequence A[] into sequence of base_val's suitable for
   * indexing into B[]
   */
  for (i = 1; i <= M; i++) {
      A2[i] = base_val[A[i]];
  }

  //#ifdef notdef
#if 1
  /*
   * Having got our sequence vector we now change it from count to
   * score based. Not sure what to do with 'N's, we currently simply
   * treat them as any other base.
   */
  for (i = 1; i <= N; i++) {
      int k;
      
      for (j = 0; j < 6; j++) {
	  B2[i][j] = 0;
	  for (k = 0; k < 6; k++) {
	      B2[i][j] += B[i][k] * W[nuascii[j]][nuascii[k]];
	  }
      }
  }
#else
#ifdef notdef
  for (i = 1; i <= N; i++) {
      int k;
      
      for (j = 0; j < 6; j++) {
	  B2[i][j] = 0;
	  for (k = 0; k < 6; k++) {
	      B2[i][j] += B[i][k] * tmp[j][k];
	  }
      }
  }
#else
  for (i = 1; i <= N; i++) {
      int k, sum;
      sum = B[i][0] + B[i][1] + B[i][2] + B[i][3] + B[i][4] + B[i][5];
      
      for (j = 0; j < 6; j++) {
	  B2[i][j] = 0;
	  for (k=0; k<6; k++) {
	      B2[i][j] += B[i][k] * tmp[j][k];
	  }
	  B2[i][j] /= sum;
      }
  }
#endif
#endif

  c = align(A2,B2,M,N,0,0,1,1,1,1);   /* OK, do it */

#ifdef notdef
  if (EG == 1) {
     c = align(A2,B2,M,N,0,0,1,1,1,1);   /* OK, do it */
  } else {
     c = align(A2,B2,M,N,-g,-g,0,0,0,0);   /* OK, do it */
  }

  c /= ACC;

  if (EG == 1) {
     ck = CHECK_SCORE(A,B2,M,N,S,1);
  } else {
     ck = CHECK_SCORE(A,B2,M,N,S,0);
  }
  if (c != ck) printf("Check_score error. c=%d, ck=%d\n",c,ck);
#endif

  xfree(CD);
  xfree(RS);
  xfree(A2);
  xfree(B2);

  return c;
}

/* Alignment display routine */

static char ALINE[51], CLINE[51];
static FastInt BLINE[51][6];

void display_sv(A,B,M,N,S,AP,BP)
    char A[];
    FastInt (*B)[6];
    FastInt M, N;
    FastInt S[], AP, BP;
{ register char *a, *c;
  FastInt (*b)[6];
  register FastInt   i,  j, op;
           FastInt   lines, ap, bp;

  i = j = op = lines = 0;
  A--;
  B--;
  ap = AP;
  bp = BP;
  a = ALINE;
  b = BLINE;
  c = CLINE;
  while (i < M || j < N)
    { if (op == 0 && *S == 0)
        { op = *S++;
          *a = A[++i];
          memcpy(*b, B[++j], 6*sizeof(FastInt));
          *c++ = (*a++ == consen_6(*b++)) ? '|' : ' ';
        }
      else
        { if (op == 0)
            op = *S++;
          if (op > 0)
            { *a++ = ' ';
              memcpy(*b++, B[++j], 6*sizeof(FastInt));
              op--;
            }
          else
            { *a++ = A[++i];
              memset(*b++, 0, 6*sizeof(FastInt));
              op++;
            }
          *c++ = '-';
        }
      if (a >= ALINE+50 || i >= M && j >= N) {
	  char *x;

	  *a = *c = '\0';
          vmessage("\n%5d ",50*lines++);
          for (x = ALINE+10; x <= a; x += 10)
            vmessage("    .    :");
          if (x <= a+5)
            vmessage("    .");
          vmessage("\n%5d %s\n      %s\n",ap,ALINE,CLINE);

	  {
	      int tmp;
	      int i, j;

	      do {
		  tmp = 0;
		  for (i=0; i<c-CLINE; i++) {
		      for (j=0; j<6; j++) {
			  if (BLINE[i][j]) {
			      if (!tmp)
				  vmessage("%5d ", bp);
			      putchar("ACGT*-"[j]);
			      BLINE[i][j]--;
			      tmp=1;
			      break;
			  }
		      }
		      if (j == 6)
			  putchar(' ');
		  }
		  putchar('\n');
	      } while (tmp == 1);
	  }

	  ap = AP + i;
	  bp = BP + j;
          a = ALINE;
          b = BLINE;
          c = CLINE;
        }
    }
}


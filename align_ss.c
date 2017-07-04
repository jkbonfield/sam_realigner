/* A PACKAGE FOR SEQUENCE COMPARISON WITH AFFINE WEIGHTS */
/* Here we maximize the similarity score and won't penalize the first
   and last gaps */

/* Globally passed params and macros */

#include <stdio.h>
#include <stdlib.h>
//#include "align.h"

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

static FastInt (*CD)[2];	/* Forward cost-only vectors */
static FastInt (*RS)[2];	/* Reverse cost-only vectors */

/* align(A,B,M,N,tb,te) returns the cost of an optimum conversion between
   A[1..M] and B[1..N] that begins(ends) with a delete if tb(te) is zero
   and appends such a conversion to the current script.                   */
/* topr, botr, lc, rc are used to trace the boundary lines */

static FastInt align(A,B,M,N,tb,te,topr,botr,lc,rc) unsigned char *A, *B; FastInt M, N;
FastInt tb, te; char topr, botr, lc, rc;

{        FastInt   midi, midj, type;	/* Midpoint, type, and cost */
         FastInt midc;

{ register FastInt   i, j;
  register FastInt c, e, d, s;
           FastInt t, *wa;
/* Boundary cases: M <= 1 or N == 0 */

  if (N <= 0)
    { if (M > 0) DEL(M)
      if (topr || botr) return 0;
      else return -gap(M);
    }
  if (M <= 1)
    { if (M <= 0)
        { INS(N);
	  if (topr || botr) return 0;
          else return -gap(N);
        }
      if (topr) {
	 if (rc) midc = 0;
	 else midc = te-h;
	 midj = 0;
	 wa = w[A[1]];
	 for (j = 1; j <= N; j++)
	   { c = wa[B[j]] - gap(N-j);
             if (c > midc)
               { midc = c;
                 midj = j;
               }
           }
      } else if (botr) {
	 if (lc) midc = 0;
	 else midc = tb-h;
	 midj = 0;
	 wa = w[A[1]];
	 for (j = 1; j <= N; j++)
	   { c = -gap(j-1) + wa[B[j]];
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
         wa = w[A[1]];
         for (j = 1; j <= N; j++)
           { c = -gap(j-1) + wa[B[j]] - gap(N-j);
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
      wa = w[A[i]];
      for (j = 1; j <= N; j++) {
	  if ((c   - m) > (e =   e   - h)) e = c - m;
	  if ((j == N) && rc) {
             if ((c = CD[j][0]) > (d = CD[j][1])) d = c;
	  } else {   
             if ((c = CD[j][0] - m) > (d = CD[j][1] - h)) d = c;
	  }
          c = s + wa[B[j]];
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
      wa = w[A[i+1]];
      for (j = N-1; j >= 0; j--)
        { if ((c   - m) > (e =   e   - h)) e = c - m;
	  if ((j == 0) && lc) {
             if ((c = RS[j][0]) > (d = RS[j][1])) d = c;
	  } else {
             if ((c = RS[j][0] - m) > (d = RS[j][1] - h)) d = c;
	  }
          c = s + wa[B[j+1]];
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

/* Interface and top level of comparator */
/* EG specifies if we want to penalize the end gaps */

FastInt align_ss(A,B,M,N,low,up,W,G,H,S,s1,s2,e1,e2)
char A[],B[]; FastInt M,N; FastInt W[][128],G,H; FastInt S[];
FastInt low,up;
FastInt s1,s2,e1,e2;
{ 
  FastInt c;

  CD = xmalloc(sizeof(FastInt) * 2 * (N+1));
  RS = xmalloc(sizeof(FastInt) * 2 * (N+1));

  if (!CD || !RS)
      return -1;

  A--;B--;

  w = W;			/* Setup global parameters */
  g = G;
  h = H;
  m = g+h;
  sapp = S;
  last = 0;
  sl = 0;

  c = align(A,B,M,N,-g,-g,s1,s2,e1,e2);   /* OK, do it */

  xfree(CD);
  xfree(RS);

  return c;
}

/* Alignment display routine */

static char ALINE[51], BLINE[51], CLINE[51];

void display_ss(A,B,M,N,S,AP,BP)
    char A[], B[]; FastInt M, N; FastInt S[], AP, BP;
{ register char *a, *b, *c;
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
          *b = B[++j];
          *c++ = (*a == *b) ? '|' : (toupper(*a) == toupper(*b) ? ':' : ' ');
          a++; b++;
        }
      else
        { if (op == 0)
            op = *S++;
          if (op > 0)
            { *a++ = ' ';
              *b++ = B[++j];
              op--;
            }
          else
            { *a++ = A[++i];
              *b++ = ' ';
              op++;
            }
          *c++ = '-';
        }
      if (a >= ALINE+50 || i >= M && j >= N)
        { *a = *b = *c = '\0';
          vmessage("\n%5d ",50*lines++);
          for (b = ALINE+10; b <= a; b += 10)
	      vmessage("    .    :");
          if (b <= a+5)
	      vmessage("    .");
          vmessage("\n%5d %s\n      %s\n%5d %s\n",ap,ALINE,CLINE,bp,BLINE);
	  ap = AP + i;
	  bp = BP + j;
          a = ALINE;
          b = BLINE;
          c = CLINE;
        }
    }
}

#ifdef never_used
/* CHECK_SCORE - return the score of the alignment stored in S */

static FastInt CHECK_SCORE(A,B,M,N,S,EG)
char A[], B[]; FastInt M, N; FastInt S[]; char EG;
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
		score = w[A[++i]][B[++j]] + score;
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

#ifdef TEST_MAIN
FastInt W128[128][128];

void init_W128(void) {
    int i;

    memset(&W128[0][0], -1, 128*128*sizeof(W128[0][0]));
    for (i = 0; i < 8; i++)
	W128["ACGTacgt"[i]]["ACGTacgt"[i]] = 2;
}

void init_W128_score(int mis, int mat) {
  int i, j;

  for (i = 0; i < 128; i++)
    for (j = 0; j < 128; j++)
      W128[i][j] = mis;
  for (i = 0; i < 16; i++)
    W128["ACGTacgtacgtACGT"[i]]["ACGTacgtACGTacgt"[i]] = mat;

  for (i = 0; i < 128; i++)
    W128['n'][i] = W128['n'][i] = W128[i]['N'] = W128[i]['n'] = 0;
}

int main(int argc, char **argv) {
    char *A = argv[1];
    char *B = argv[2];
    int s1 = atoi(argv[3]);
    int s2 = atoi(argv[4]);
    int e1 = atoi(argv[5]);
    int e2 = atoi(argv[6]);
    int score;
    FastInt S[10000];
    int G = atoi(argv[7]); // open
    int H = atoi(argv[8]); // extend
    int mis = atoi(argv[9]);
    int mat = atoi(argv[10]);

    init_W128_score(mis,mat);

    score = align_ss(A, B, strlen(A), strlen(B),
                     0, 0, W128, G, H, S, s1, s2, e1, e2);
    
    printf("Score=%d\n", score);
    display_ss(A, B, strlen(A), strlen(B), S, 0, 0);

    return 0;
}

/* Eg:
   align_ss GATCGAGTGGA GATCTGAGTG 0 0 0 0 4 1 -2 1
 */
#endif

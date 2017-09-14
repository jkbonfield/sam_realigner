// Copy of assem_bam3.c, but for use within bam_problem_regions instead of
// as a standalone tool.


// cc -g -I$HOME/work/samtools_master/htslib assem_bam.c hash_table.c pooled_alloc.c string_alloc.c align_ss.c -DKMER=70 -lm -L$HOME/work/samtools_master/htslib -lhts -lz -pthread

// ./a.out a.fasta ref.fasta > a.sam && sfdp -Gstart=24 -Gmaxiter=10 -GK=.5 -Gsplines=true -Efontsize=24 -Nshape=point  g.dot -Tpdf -o g.pdf; evince g.pdf

// 1. Hash seqs and correct
// 2. Build graph
// 3. Find linear strings in graphs and assign a seq-vector to each head.
// 4. Find bubbles and Collapse, forming new vectors.
// 5. Add ref and repeat step 4.

/*
The main benefit to this approach is that the hash keys can still be
tied to the location their last (or first) base matches in the graph,
but the key itself isn't used as part of the vector so that indels
don't cause frameshifts when computing the new vector seq.

It's still possible (trying it here first) to have one vector per node
with the vector being unrelated to the hash key, but instead to the
appropriate location within the string.
 */


/*
 * TODO
 *
 * - Merge bubbles in better order.  Most significant first?  Smallest first?
 *   See eg13.fa for an example fail.
 *
 * - Don't assume the first KMER bases are all match.  Only the last we know
 *   to match.  Instead align back to ref.
 *
 * - Fix eg4h/H.  We could try merging longest bubbles before shorter ones?
 *
 * - Sanitise edges.  Why in and out have separate edges structures?
 *   Would make sense for n1->n2 edge to have n1->out[?] == n2->in[?] == edge.
 *   In this case edge->n[0] == n1 and edge->n[1] == n2.  Right now however
 *   n[0] is always <this node> and n[1] is always <other node>, which forces
 *   having two edge structs.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <limits.h>

#include "hash_table.h"
#include "string_alloc.h"
#include "str_finder.h"

#include "htslib/sam.h"
#include "htslib/kstring.h"

#define IS_REF 1
#define IS_CON 2

//---------------------------------------------------------------------------

int W128[128][128];
extern char base_val[128];
void init_W128_score(int mis, int mat) {
  int i, j;

  for (i = 0; i < 128; i++)
    for (j = 0; j < 128; j++)
      W128[i][j] = mis;
  for (i = 0; i < 16; i++)
    W128["ACGTacgtacgtACGT"[i]]["ACGTacgtACGTacgt"[i]] = mat;

  for (i = 0; i < 128; i++)
    W128['n'][i] = W128['n'][i] = W128[i]['N'] = W128[i]['n'] = 0;

  extern init_base_val(void);
  init_base_val();
}

// Ambiguity codes include base-pad combinations.
// We set the scores so that given the choice of aligning
// A vs AC ambig or A* ambig, the AC scores higher.  Thus by
// elimination pads preferentially align to A* than AC.
int X128[128][128];
void init_X128_score(int mis, int mat) {
    memset(X128, 0, 128*128*sizeof(int)); // why 0 and not mis?
    int i, j;

    // Matches
    for (i = 0; i < 4; i++)
	X128["ACGT"[i]]["ACGT"[i]] = mat;

    // Ambiguity codes; partial match
    for (i = 0; i < 3; i++)
	X128['A']["MRW"[i]] = X128["MRW"[i]]['A'] = mat/2;
    for (i = 0; i < 3; i++)
	X128['C']["MSY"[i]] = X128["MSY"[i]]['C'] = mat/2;
    for (i = 0; i < 3; i++)
	X128['G']["RSK"[i]] = X128["RSK"[i]]['G'] = mat/2;
    for (i = 0; i < 3; i++)
	X128['T']["WYK"[i]] = X128["WYK"[i]]['T'] = mat/2;
    // Ambiguity codes; mismatch (B = not A; D = not C; etc)
    X128['A']['B'] = X128['B']['A'] = mis/2;
    X128['C']['D'] = X128['D']['C'] = mis/2;
    X128['G']['H'] = X128['H']['G'] = mis/2;
    X128['T']['V'] = X128['V']['T'] = mis/2;
    
    // Ambiguity vs pads.
    // All to start with, to set mismatch (A vs G*) baseline.
    for (i = 0; i < 11; i++) {
	for (j = 0; j < 11; j++) {
	    X128["AMRWCSYGKTN"[i]][tolower("AMRWCSYGKTN"[j])] = mis;	
	    X128[tolower("AMRWCSYGKTN"[i])]["AMRWCSYGKTN"[j]] = mis;
	    // A* vs C* is a partial match (pad wise). How to score?
	    X128[tolower("AMRWCSYGKTN"[i])][tolower("AMRWCSYGKTN"[j])] = mis/2;
	}
    }
    // Now fix up the partial matches (A vs A*).
    for (i = 0; i < 4; i++) {
	X128["ACGT"[i]][tolower("ACGT"[i])] = mat/2-1;
	X128[tolower("ACGT"[i])]["ACGT"[i]] = mat/2-1;
	X128[tolower("ACGT"[i])][tolower("ACGT"[i])] = mat/2-1;
    }

    for (i = 0; i < 3; i++) {
	X128['A'][tolower("MRW"[i])] = X128[tolower("MRW"[i])]['A'] = mis/2;
	X128['C'][tolower("MSY"[i])] = X128[tolower("MSY"[i])]['C'] = mis/2;
	X128['G'][tolower("RSK"[i])] = X128[tolower("RSK"[i])]['G'] = mis/2;
	X128['T'][tolower("WYK"[i])] = X128[tolower("WYK"[i])]['T'] = mis/2;
    }

    // Base vs pad is bad, worse than opening a new gap.
    for (i = 0; i < 11; i++) {
	X128["AMRWCSYGKTN"[i]]['-'] = X128['-'][tolower("AMRWCSYGKTN"[i])] = mis*2;
	X128[tolower("AMRWCSYGKTN"[i])]['-'] = X128['-']["AMRWCSYGKTN"[i]] = mis*2;
    }
}

// Convert a 6-way ACGT*N vector to X128 lookup value.
int vec2X(int V[6]) {
#if 0
    char b = "NTGKCYSBAWRDMHVN"[(!!V[0]<<3)+(!!V[1]<<2)+(!!V[2]<<1)+(!!V[3]<<0)];
    if (V[4]>0) b = tolower(b);
    return b;
#else
    int t=V[0]+V[1]+V[2]+V[3];
    int T=t>>2;
    char b = "NTGKCYSBAWRDMHVN"[((V[0]>T)<<3)+((V[1]>T)<<2)+((V[2]>T)<<1)+((V[3]>T)<<0)];
    if (V[4]>0) b = tolower(b);
    if (t < V[4]/4) b = '-';
    return b;
#endif
}

char *vec2seq(int (*v)[5], int len) {
    int i;
    char *s = malloc(len+1);
    if (!s) return NULL;

    for (i = 0; i < len; i++)
	s[i] = vec2X(v[i]);
    s[i] = 0;

    return s;
}

int align_vv(int (*v1)[5], int (*v2)[5], int len1, int len2,
	     int G, int H, int *S, int s1, int s2, int e1, int e2) {
    char *seq1 = vec2seq(v1, len1); //puts(seq1);
    char *seq2 = vec2seq(v2, len2); //puts(seq2);
    
    int score = align_ss(seq1, seq2, len1, len2, 0,0, X128, G,H, S, s1,s2,e1,e2);

//    printf("Score=%d\n", score);
//    display_ss(seq1, seq2, len1, len2, S, 0, 0);

    free(seq1);
    free(seq2);

    return score;
}

typedef struct haps {
    int   pos;
    char *seq;
    char *name;
    uint8_t *qual; // not our copy; do not free
} haps_t;

//#include "hap2s.h"
//#include "haps.h"

#ifndef KMER
#  define KMER 14
#endif

#ifndef MAX_KMER
#  define MAX_KMER 100
#endif

#ifndef MAX_SEQ
#  define MAX_SEQ 1024
#endif

// FIXME: needs to be adaptive
#ifndef MAX_EDGE
#  define MAX_EDGE 16
#endif

#define GAP_OPEN 3
#define GAP_EXTEND 1

struct edge;

typedef struct node_t {
    HashItem **hi;
    int n_hi;
    int bases[MAX_KMER][5];
    int count;
    struct edge *in[MAX_EDGE], *out[MAX_EDGE]; // incoming and outgoing edges
    int n_in, n_out;
    int pruned;
    int id;
    int incoming;
    int visited;
    int is_head, is_tail; // fixme: pruned, is_hair etc should be flags

    // True if node is part of the reference
    int ref;
    int pos;
    int ins; // Nth base of an insertion between two mapped pos.
    int *posa; // points to an array of positions

    // Path tracking for breadth first search.
    //int path_id;
    // in[0] will be swapped around to ensure it is always the path parent.
} node_t;

typedef struct edge {
    int n[2];   // n[0] = from, n[1] = to. n[0] only used in hash key.
    int count;
} edge_t;

typedef struct {
    int kmer;

    node_t **node;
    int nnodes, anodes;  // number used and allocated
    edge_t **edge;
    int nedges, aedges;

    HashTable *node_hash;  // "seq" -> node
    HashTable *edge_hash;  // id[2] -> edge

    string_alloc_t *spool;
} dgraph_t;

//---------------------------------------------------------------------------
// Paths
typedef struct path {
    node_t *n, *pn; // curr and prev node
    //node_t **p;
    int np, ap;
    struct path *next;
    int id;
} path_t;

path_t *path_create(int id) {
    path_t *p = calloc(1, sizeof(path_t));
    p->id = id;
    return p;
}

int path_add(path_t *p, node_t *n) {
//    if (p->np >= p->ap) {
//	p->ap = p->ap ? p->ap*2 : 8;
//	node_t **p2 = realloc(p->p, p->ap * sizeof(*p->p));
//	if (!p2)
//	    return -1;
//	p->p = p2;
//    }
//
//    p->p[p->np++] = n;
    p->n = n;

    if (n->visited)
	return n->visited;

    n->visited = p->id;
    return 0;
}

void path_destroy(path_t *p) {
    if (!p)
	return;
//    free(p->p);
    free(p);
}

//---------------------------------------------------------------------------
// Graphs

void graph_destroy(dgraph_t *g) {
    if (!g)
	return;

    if (g->node_hash)
	HashTableDestroy(g->node_hash, 0);

    if (g->edge_hash)
	HashTableDestroy(g->edge_hash, 0);

    if (g->spool)
	string_pool_destroy(g->spool);

    int i;
    for (i = 0; i < g->nnodes; i++) {
	free(g->node[i]->hi);
	free(g->node[i]);
    }
    free(g->node);

    for (i = 0; i < g->nedges; i++)
	free(g->edge[i]);
    free(g->edge);

    free(g);
}

dgraph_t *graph_create(int kmer) {
    dgraph_t *g = calloc(1, sizeof(*g));
    if (!g) return NULL;

    g->kmer = kmer;
    g->node_hash = HashTableCreate(8, HASH_DYNAMIC_SIZE | HASH_POOL_ITEMS);
    if (!g->node_hash) {
	graph_destroy(g);
	return NULL;
    }

    g->edge_hash = HashTableCreate(8, HASH_DYNAMIC_SIZE | HASH_POOL_ITEMS);
    if (!g->edge_hash) {
	graph_destroy(g);
	return NULL;
    }

    g->spool = string_pool_create(1024*1024);
    if (!g->spool) {
	graph_destroy(g);
	return NULL;
    }

    return g;
}

node_t *add_node(dgraph_t *g, char *seq, int len) {
    node_t *node = calloc(1, sizeof(*node));
    if (!node)
	return NULL;

    if (g->nnodes >= g->anodes) {
	int a = g->anodes ? g->anodes*2 : 1024;
	node_t **n = realloc(g->node, a * sizeof(*n));
	if (!n) {
	    free(node);
	    return NULL;
	}
	g->node = n;
	g->anodes = a;
    }

    g->node[g->nnodes] = node;
    if (seq) {
	//node->len = len;
	//node->seq = string_alloc(g->spool, node->len);
	//memcpy(node->seq, seq, node->len);
	node->pos = INT_MIN;
	int j;
	for (j = 0; j < len; j++)
	    node->bases[j][base_val[seq[j] & 0x7f]]++;
	node->count = 0;

	HashData hd;
	hd.p = node;
	HashItem *hi = HashTableAdd(g->node_hash, seq, len, hd, 0);

	node->hi = realloc(node->hi, (node->n_hi+1) * sizeof(*node->hi));
	node->hi[node->n_hi++] = hi;
//    } else {
//	node->seq = NULL;
//	node->len = 0;
    }
    
    // Debugging
    node->id = g->nnodes++;

    return node;
}

// FIXME: speed up via hash or similar
node_t *find_node(dgraph_t *g, char *seq, int len, int add) {
    HashItem *hi = HashTableSearch(g->node_hash, seq, len);
    if (hi)
	return (node_t *)hi->data.p;

    // try uppercasing
    char tmp[MAX_KMER];
    int i;
    for (i = 0; i < len && i < MAX_KMER; i++)
	tmp[i] = toupper(seq[i]);

    hi = HashTableSearch(g->node_hash, tmp, len);
    if (hi)
	return (node_t *)hi->data.p;

    return add ? add_node(g, seq, len) : NULL;
}

edge_t *add_edge(dgraph_t *g, node_t *n1, node_t *n2) {
    edge_t *edge = malloc(sizeof(*edge));
    if (!edge)
	return NULL;

    // For simplicity of deallocation
    if (g->nedges+1 >= g->aedges) {
	int a = g->aedges ? g->aedges*2 : 1024;
	edge_t **e = realloc(g->edge, a * sizeof(*e));
	if (!e) {
	    free(edge);
	    return NULL;
	}
	g->edge = e;
	g->aedges = a;
    }
    g->edge[g->nedges++] = edge;

    // Add edge itself
    if (n1->n_out >= MAX_EDGE)
	return NULL;
    n1->out[n1->n_out++] = edge;
    edge->n[0] = n1->id;
    edge->n[1] = n2->id;
    edge->count = 0;

    HashData hd;
    hd.p = edge;
    HashTableAdd(g->edge_hash, (char *)edge->n, 2*sizeof(*edge->n), hd, 0);

    // Reverse edge (incoming)
    edge_t *i_edge = malloc(sizeof(*edge));
    if (!i_edge)
	return NULL;
    g->edge[g->nedges++] = i_edge;
    if (n1->n_in >= MAX_EDGE)
	return NULL;
    n2->in[n2->n_in++] = i_edge;
    i_edge->n[0] = n2->id;
    i_edge->n[1] = n1->id;
    i_edge->count = 0;

    return edge;
}

// Move edge e with new source (n2->n3) to n1 (n1->n3)
void move_edge_in(dgraph_t *g, edge_t *e, node_t *n1, node_t *n2) {
    node_t *n3 = g->node[e->n[1]];
    int j;
    if (n1->n_out >= MAX_EDGE)
	return;
    n1->out[n1->n_out++] = e;
    for (j = 0; j < n3->n_in; j++)
	if (n3->in[j]->n[1] == n2->id)
	    n3->in[j]->n[1] = n1->id;

    // Correct edge_hash too.
    HashItem *hi = HashTableSearch(g->edge_hash, (char *)e->n, 2*sizeof(*e->n));
    assert(hi);
    HashData hd;
    hd.p = e;
    e->n[0] = n1->id;
    HashTableAdd(g->edge_hash, (char *)e->n, 2*sizeof(*e->n), hd, 0);
}

// Move edge (n1->l1) to new destination n1->n2
void move_edge_out(dgraph_t *g, node_t *n1, node_t *l1, node_t *n2) {
    int j;
    edge_t *e;

    for (j = 0; j < n1->n_out; j++)
	if (n1->out[j]->n[1] == l1->id)
	    break;
    assert(j != n1->n_out);

    // Correct edge_hash too.
    e = n1->out[j];
    HashItem *hi = HashTableSearch(g->edge_hash, (char *)e->n, 2*sizeof(*e->n));
    assert(hi);

    HashData hd;
    hd.p = e;
    e->n[1] = n2->id;
    HashTableAdd(g->edge_hash, (char *)e->n, 2*sizeof(*e->n), hd, 0);
}

edge_t *find_edge(dgraph_t *g, node_t *n1, node_t *n2) {
    int n[2] = {n1->id, n2->id};
    HashItem *hi = HashTableSearch(g->edge_hash, (char *)n, 2*sizeof(*n));
    if (hi)
	return (edge_t *)hi->data.p;

    return add_edge(g, n1, n2);
}

edge_t *incr_edge(dgraph_t *g, char *seq1, int len1, char *seq2, int len2, int ref) {
    node_t *n1 = find_node(g, seq1, len1, 1);
    node_t *n2 = find_node(g, seq2, len2, 1);

    if (!n1 || !n2)
	return NULL;

    edge_t *e = find_edge(g, n1, n2);
    if (!e)
	return NULL;

    e->count++;
    n2->count++;

    return e;
}

int decr_edge(dgraph_t *g, char *seq1, int len1, char *seq2, int len2, int ref) {
    node_t *n1 = find_node(g, seq1, len1, 1);
    node_t *n2 = find_node(g, seq2, len2, 1);

    if (!n1 || !n2)
	return -1;

    edge_t *e = find_edge(g, n1, n2);
    if (!e)
	return -1;

    e->count--;
    n2->count--;

    if (e->count <= 0) {
	HashTableRemove(g->edge_hash, (char *)e->n, 2*sizeof(*e->n), 0);
	int i, j;
	for (i = j = 0; i < n1->n_out; i++)
	    if (n1->out[i]->n[1] != n2->id)
		n1->out[j++] = n1->out[i];
	n1->n_out = j;
	for (i = j = 0; i < n2->n_in; i++)
	    if (n2->in[i]->n[1] != n1->id)
		n2->in[j++] = n2->in[i];
	n2->n_in = j;
    }

    return 0;
}

// Recurse through each node checking the visitor status.
// If it's set already then we've been here at least once, so bail out - either
// pass or fail depending on value.
// If it's not set, then set it, recurse, and "unset" (by decrementing).
//
// The decrement avoids too many calls when we start at multiple points in the graph.
// If we reset to 0 then every starting point would need validating.
static int loop_check_recurse(dgraph_t *g, node_t *n, int check_num, int level,
			      int *last_node, int *last_edge) {
    //fprintf(stderr, "n=%d level=%d visit=%d\n", n->id, level, n->visited);

    //fprintf(stderr, "    n=%d visit=%d/%d\n", n->id, n->visited, check_num);
    if (n->id >= g->nnodes || n->visited >= check_num)
	// loop
	return -1;

    if (n->visited > 0 && n->visited < check_num)
	// checked before
	return 0;

    n->visited = check_num;

    int i;
    for (i = 0; i < n->n_out; i++) {
	if (loop_check_recurse(g, g->node[n->out[i]->n[1]], check_num, level+1, last_node, last_edge) < 0) {
	    n->visited--;
	    // FIXME?: return node with pair of paths that link to it (known?) and let
	    // caller pick which one to break based on edge count?
	    // Current method: just assume the one that got there first was correct as it
	    // represents the shorter assembly.
	    if (last_node && last_edge && *last_node == INT_MIN) {
		*last_node = n->id;
		*last_edge = i;
	    }
	    return -1;
	}
    }

    n->visited--;

    return 0;
}

int loop_check(dgraph_t *g, int loop_break) {
    int i;

    // Preserve visited status so we can use this within other algorithms.
    int *v = malloc(g->nnodes*sizeof(*v));
    for (i = 0; i < g->nnodes; i++)
	v[i] = g->node[i]->visited;

    for (i = 0; i < g->nnodes; i++)
	g->node[i]->visited = 0;

    int check_num = INT_MAX/2;
    for (i = 0; i < g->nnodes; i++) {
	node_t *n = g->node[i];

	// Have we been here before? If so skip it.
	if (n->pruned || n->visited)
	    continue;
	
	// We have a head node, but it may link up with
	// previously scanned head nodes.  All we care
	// about is following this starting point doesn't
	// get back to a node we're visiting in this cycle.
	int last_node = INT_MIN;
	int last_edge = 0;
	if (loop_check_recurse(g, n, ++check_num, 0, &last_node, &last_edge) < 0) {
	    if (loop_break) {
		node_t *last = g->node[last_node];
		node_t *to = g->node[last->out[last_edge]->n[1]];
		fprintf(stderr, "last_node=%d -> %d\n", last_node, to->id);
		memmove(&last->out[last_edge], &last->out[last_edge+1], 
			(last->n_out - last_edge) * sizeof(*last->out));
		last->n_out--;
		int z,k;
		for (z = k = 0; z < to->n_in; z++) {
		    if (to->in[z]->n[1] != last->id)
			to->in[k++] = to->in[z];
		}
		assert(k == to->n_in-1);
		to->n_in--;
	    }

	    for (i = 0; i < g->nnodes; i++)
		g->node[i]->visited = v[i];
	    free(v);

	    return -1;
	}
    }


    for (i = 0; i < g->nnodes; i++)
	g->node[i]->visited = v[i];
    free(v);

    return 0;
}


int add_seq(dgraph_t *g, char *seq, int len, int ref) {
    int i, j;
    static int counter = 1;

    if (!len)
	len = strlen(seq);

    if (!(ref & IS_REF)) {
	// Prune trailing Ns.
	while (*seq == 'N' && len > 0)
	    seq++, len--;
	while (seq[len-1] == 'N' && len > 0)
	    len--;

	if (len == 0)
	    return 0;
    }

    counter++;

    char *s = malloc(len);
    for (i = j = 0; i < len; i++)
	if (seq[i] != '*')
	    s[j++] = toupper(seq[i]);
    len = j;

    // split into successive pairs of kmers, incrementing edges
    for (i = 0; i < len - g->kmer; i++) {
	edge_t *e;
	if (!(e = incr_edge(g, s+i, g->kmer, s+i+1, g->kmer, ref)))
	    return -1;

	if (g->node[e->n[1]]->visited == counter) {
	    // cull edge
	    free(s);
	    return -1; // loop
	}

	g->node[e->n[1]]->visited = counter;

	if (ref) {
	    node_t *n1 = g->node[e->n[0]];
	    node_t *n2 = g->node[e->n[1]];
	    // FIXME: detect possible loops here. Ref must be a clean pass through!
	    //n1->pos = i+g->kmer-1;
	    if (ref & IS_REF)
		n1->pos = i;
	    n1->ref |= ref;
	    
	    //n2->pos = i+g->kmer;
	    if (ref & IS_REF)
		n2->pos = i+1;
	    n2->ref |= ref;
	}
    }

    free(s);

    return 0;
}

int del_seq(dgraph_t *g, char *seq, int len, int ref) {
    int i, j;

    if (!len)
	len = strlen(seq);

    if (!(ref & IS_REF)) {
	// Prune trailing Ns.
	while (*seq == 'N' && len > 0)
	    seq++, len--;
	while (seq[len-1] == 'N' && len > 0)
	    len--;

	if (len == 0)
	    return 0;
    }

    char *s = malloc(len);
    for (i = j = 0; i < len; i++)
	if (seq[i] != '*')
	    s[j++] = toupper(seq[i]);
    len = j;

    // split into successive pairs of kmers, incrementing edges
    for (i = 0; i < len - g->kmer; i++) {
	decr_edge(g, s+i, g->kmer, s+i+1, g->kmer, ref);
    }

    free(s);

    return 0;
}

// Returns true if the edge is part of the reference
int ref_edge(dgraph_t *g, edge_t *e) {
    node_t *n1 = g->node[e->n[0]];
    node_t *n2 = g->node[e->n[0]];
    return (n1->ref & IS_REF) || (n2->ref & IS_REF);
}

int ref_node(dgraph_t *g, int id) {
    return g->node[id]->ref & IS_REF;
}

// Returns the best previous node, scored by incoming edge count
node_t *best_prev_node(dgraph_t *g, node_t *n) {
    if (!n)
	return NULL;

    int best_i = 0, best_count = INT_MIN;
    int i;
    for (i = 0; i < n->n_in; i++) {
	if (best_count < n->in[i]->count) {
	    best_count = n->in[i]->count;
	    best_i = i;
	}
    }
    return best_count != INT_MIN ? g->node[n->in[best_i]->n[1]] : NULL;
}

// TODO:
//
// Where a node has multiple edges, if one edge has a low count and
// the path to the end (or where it merges again) is continuously low,
// remove it and all nodes.  (Ideally we want to collapse those reads back
// in to the graph somehow; maybe realign to haplotype.)
//
// Instead: take local greedy approach.  Good enough?

void prune_tails(dgraph_t *g, int min_count) {
    int done_something;

    do {
	int i, j;
	done_something = 0;

	for (i = g->nnodes-1; i >= 0; i--) {
	    node_t *n = g->node[i];
	    if (n->pruned)
		continue;

	    //fprintf(stderr, "Node %d:\n", i);

	    int k;
	    for (j = k = 0; j < n->n_out; j++, k++) {
		n->out[k] = n->out[j];
		edge_t *e = n->out[j];
		if (e->count >= min_count || ref_node(g, e->n[1]))
		    continue;

		if (g->node[e->n[1]]->n_out)
		    continue;

		//fprintf(stderr, "Pruned %d\n", e->n[1]);
		g->node[e->n[1]]->pruned = 1;
		k--;
		done_something = 1;
	    }
	    n->n_out = k;
	}
    } while (done_something);
}

void prune_heads(dgraph_t *g, int min_count) {
    int done_something;

    do {
	int i, j;
	done_something = 0;

	for (i = 0; i < g->nnodes; i++) {
	    node_t *n = g->node[i];
	    if (n->pruned || n->n_in)
		continue;

	    for (j = 0; j < n->n_out; j++)
		if (n->out[j]->count >= min_count || ref_edge(g, n->out[j]))
		    break;

	    if (j != n->n_out)
		continue;

	    for (j = 0; j < n->n_out; j++)
		g->node[n->out[j]->n[1]]->n_in--;//-=n->out[j]->count;
	    n->pruned = 1;
	    done_something = 1;
	    //fprintf(stderr, "Pruned %d\n", n->id);
	}
    } while (done_something);
}


// Node n_end has parents p1 and p2 which meet up again at some common
// node n_start.  Find n_start.
void node_common_ancestor(dgraph_t *g, node_t *n_end, node_t *p1, node_t *p2) {
    node_t **path1 = malloc(g->nnodes * sizeof(node_t *));
    node_t **path2 = malloc(g->nnodes * sizeof(node_t *));
    int np1 = 0, np2 = 0, i, j;

//    printf("Node_common_ancestor %d<-%d %d<-%d\n",
//	   n_end->id, p1->id,
//	   n_end->id, p2->id);

    // Backtrack up p1 marking set.
    node_t *n = p1;
    while (n) {
	n->visited |= (1<<31);
	n = n->n_in ? g->node[n->in[0]->n[1]] : NULL;
    }

    // Backtrack up p2 until we find 'visited' set
    n = p2;
    while (n) {
	if (n->visited & (1<<31))
	    break;
	n = n->n_in ? g->node[n->in[0]->n[1]] : NULL;
    }
    int start = n->id;

    // Clear visited status
    n = p1;
    while (n) {
	n->visited &= ~(1<<31);
	n = n->n_in ? g->node[n->in[0]->n[1]] : NULL;
    }

    // Create (reversed) paths
    n = p1;
    while (n && n->id != start) {
	//printf("1[%d]: %d\n", np1, n->id);
	path1[np1++] = n;
	n = n->n_in ? g->node[n->in[0]->n[1]] : NULL;
    }

    n = p2;
    while (n && n->id != start) {
	//printf("2[%d]: %d\n", np2, n->id);
	path2[np2++] = n;
	n = n->n_in ? g->node[n->in[0]->n[1]] : NULL;
    }
    

//    // Report
//    n = p1;
//    printf(" => path 1: ");
//    while (n && n->id != start) {
//	printf(" %d", n->id);
//	n = n->n_in ? g->node[n->in[0]->n[1]] : NULL;
//    }
//    printf(" %d, len=%d\n", start, np1);
//
//    n = p2;
//    printf(" => path 2: ");
//    while (n && n->id != start) {
//	printf(" %d", n->id);
//	n = n->n_in ? g->node[n->in[0]->n[1]] : NULL;
//    }
//    printf(" %d, len=%d\n", start, np2);

    // Align vectors.
    // FIXME: These should be n->bases[][] * n->count?
    int (*v1)[5] = calloc(np1+g->kmer, sizeof(*v1));
    int (*v2)[5] = calloc(np2+g->kmer, sizeof(*v2));
    int (*vn)[5] = calloc(np1+np2+2*g->kmer, sizeof(*vn));
    n = p1;
    int len1 = 0;
    node_t *l = NULL;
    while (n && n->id != start) {
	memcpy(v1[len1++], n->bases[g->kmer-1], 5*sizeof(int));
	l = n;
	n = n->n_in ? g->node[n->in[0]->n[1]] : NULL;
    }
    for (i = g->kmer-2; i >= 0; i--)
	memcpy(v1[len1++], l->bases[i], 5*sizeof(int));

    n = p2;
    int len2 = 0;
    while (n && n->id != start) {
	memcpy(v2[len2++], n->bases[g->kmer-1], 5*sizeof(int));
	l = n;
	n = n->n_in ? g->node[n->in[0]->n[1]] : NULL;
    }
    for (i = g->kmer-2; i >= 0; i--)
	memcpy(v2[len2++], l->bases[i], 5*sizeof(int));

    int *S = malloc((len1+len2)*sizeof(int));
    align_vv(v1, v2, len1, len2, GAP_OPEN, GAP_EXTEND, S, 0,0,0,0);

    // Remove path2 from start/end nodes.
    n = g->node[start];
    for (i = j = 0; i < n->n_out; i++) {
	if (n->out[i]->n[1] == (np2 ? path2[np2-1]->id : n_end->id)) {
	    // n->out[i]->count += ?
	    continue;
	} else {
	    n->out[j++] = n->out[i];
	}
    }
    n->n_out = j;

    n = n_end;
    for (i = j = 0; i < n->n_in; i++) {
	if (n->in[i]->n[1] == (np2 ? path2[0]->id : start))
	    continue;
	else
	    n->in[j++] = n->in[i];
    }
    n->n_in = j;

    // Merge v1/v2 into vn.  FIXME v1/v2 into v2 is simpler
    {
	int first = 1, *S2 = S;
	int x1 = 0, x2 = 0, p = 0;

	while (x1 < len1 && x2 < len2) {
	    int op = *S2++;
	    if (op == 0) {
		// match
		memcpy(vn[p], v1[x1], sizeof(v1[x1]));
		vn[p][0] += v2[x2][0];
		vn[p][1] += v2[x2][1];
		vn[p][2] += v2[x2][2];
		vn[p][3] += v2[x2][3];
		vn[p][4] += v2[x2][4];
		x1++; x2++;
		p++;
	    } else if (op < 0) {
		// Already in path1 only, nothing to do with path2;
		while (op++ && x1 < np1) {
		    memcpy(vn[p], v1[x1++], sizeof(v1[p]));
		    vn[p][4]++;
		    p++;
		}
	    } else {
		// Insertion in path2, link into path1
		while (op--) {
		    memcpy(vn[p], v1[x2++], sizeof(v1[p]));
		    vn[p][4]++;
		    p++;
		}
	    }
	}

	//printf("x1=%d of %d, x2=%d of %d, p=%d\n", x1, len1, x2, len2, p);
	// gap at end
	while (x1 < len1) {
	    // Already in path1 only, nothing to do with path2
	    vn[p++][4]++;
	    x1++;
	}
	while (x2 < len2) {
	    // Already in path1 only, nothing to do with path2
	    vn[p++][4]++;
	    x2++;
	}

//	for (x1 = 0; x1 < p; x1++) {
//	    printf("vn[%02d]={%d,%d,%d,%d,%d}\n",
//		   x1, vn[x1][0], vn[x1][1], vn[x1][2], vn[x1][3], vn[x1][4]);
//	}
    }

    // Merge path2 into path1.  These are in reverse order.
    // Ie. path1[0] and path2[0] both have *a* shared end (n_end);
    // path1[np1-1] and path2[np2-1] both have *a* shared parent (start).
    // => path 1:  5 4 3
    // => path 2:  14 13 12

    if (1) {
	int first = 1, *S2 = S;
	int x1 = 0, x2 = 0, p = 0;

	node_t *l1 = n_end;
	node_t *l2 = n_end;
	while (x1 < np1 && x2 < np2) {
	    node_t *n1 = path1[x1];
	    node_t *n2 = path2[x2];
	    int op = *S2++;
	    if (op == 0) {
		// March along both paths, so merge n2 into n1.
		//printf("M %2d %2d\n", n1->id, n2->id);

		// Merge coordinates
		if (n1->pos == INT_MIN) {
		    n1->pos = n2->pos;
		    n1->posa = n2->posa;
		}

		// Base frequencies
		//memcpy(n1->bases, vn[p++], g->kmer*5*sizeof(int));
		for (j = 0; j < g->kmer; j++)
		    memcpy(n1->bases[j], vn[p+(g->kmer-j)-1], 5*sizeof(int));
		p++;
		//for (j = 0; j < g->kmer; j++) {
		//    for (i = 0; i < 5; i++)
		//	n1->bases[j][i] += n2->bases[j][i];
		//}

		// Fix node hash keys
		//fprintf(stderr, "%s (%d) vs %s (%d)\n", n1->hi[0]->key, n1->count, n2->hi[0]->key, n2->count);
		n1->hi = realloc(n1->hi, (n1->n_hi + n2->n_hi)*sizeof(*n1->hi));
		if (n1->count >= n2->count) {
		    for (i = 0; i < n2->n_hi; i++) {
			HashItem *hi = n2->hi[i];
			hi->data.p = n1;
			n1->hi[n1->n_hi++] = hi;
		    }
		} else {
		    memmove(&n1->hi[n2->n_hi], &n1->hi[0], n1->n_hi*sizeof(*n1->hi));
		    for (i = 0; i < n2->n_hi; i++) {
			HashItem *hi = n2->hi[i];
			hi->data.p = n1;
			n1->hi[i] = hi;
		    }
		    n1->n_hi += n2->n_hi;
		}

		n1->count += n2->count;

		// Merge incoming.
		for (i = 0; i < n2->n_in; i++) {
		    for (j = 0; j < n1->n_in; j++)
			if (n1->in[j]->n[1] == n2->in[i]->n[1] ||
			    (x2+1 < np2 && n2->in[i]->n[1] == path2[x2+1]->id))
			    break;
		    if (j == n1->n_in && n2->in[i]->n[1] != start) { // not already a parent
			if (n1->n_in >= MAX_EDGE)
			    continue;
			n1->in[n1->n_in++] = n2->in[i];
			node_t *n3 = g->node[n2->in[i]->n[1]];
			move_edge_out(g, n3, n2, n1);
			//for (j = 0; j < n3->n_out; j++)
			//    if (n3->out[j]->n[1] == n2->id)
			//	n3->out[j]->n[1] = n1->id;
		    } else if (j != n1->n_in && n2->in[i]->n[1] != start &&
			       ((x2+1 >= np2) ||
				(x2+1 < np2 && path2[x2+1]->id != n2->in[i]->n[1]))) {
			// shared node, parent(n1) is also parent(n2) and
			// not on path, so ensure we delete the edge linking
			//to n2 // although we don't need to add it to n1 as
			// it's already there.
			node_t *n3 = g->node[n2->in[i]->n[1]];
			// remove edge out
			int k, l;
			for (k = l = 0; k < n3->n_out; k++) {
			    if (n3->out[k]->n[1] == n2->id)
				continue;
			    n3->out[l++] = n3->out[k];
			}
			n3->n_out = l;
		    }
		}

		// Merge outgoing
		for (i = 0; i < n2->n_out; i++) {
		    // Migrate n2->out unless n2->out on // n1 path
		    for (j = 0; j < n1->n_out; j++)
			if (n1->out[j]->n[1] == n2->out[i]->n[1] ||
			    (x2 > 0 && path2[x2-1] && n2->out[i]->n[1] == path2[x2-1]->id))
			    break;
		    if (j == n1->n_out) { // not already a child
			if (n2->out[i]->n[1] != n_end->id)
			    move_edge_in(g, n2->out[i], n1, n2);
		    } else {
			n1->out[j]->count += n2->out[i]->count;
		    }
		}
		// Cull n1->out that are on n2 path
		int k;
		for (j = k = 0; j < n1->n_out; j++) {
		    for (i = 0; i < n2->n_out; i++)
			if (x2 > 0 && path2[x2-1] && n1->out[j]->n[1] == path2[x2-1]->id)
			    break;
		    if (i == n2->n_out)
			n1->out[k++] = n1->out[j];
		}
		n1->n_out = k;

		// Merge positions
		//printf("Merge pos %d with %d\n", n1->pos, n2->pos);

		n2->pruned = 1;
		for (j = 0; j < n2->n_out; j++) {
		    node_t *t = g->node[n2->out[j]->n[1]];
		    int l;
		    for (l = k = 0; k < t->n_in; k++) {
			if (t->in[k]->n[1] != n2->id)
			    t->in[l++] = t->in[k];
		    }
		    t->n_in = l;
		}
		first = 0;
		x1++, x2++;
		l1 = n1;
		l2 = n2;
	    } else if (op < 0) {
		// Already in path1 only, nothing to do with path2
		while (op++ && x1 < np1) {
		    //printf("D %2d  -\n", n1->id);
		    //n1->bases[g->kmer-1][4]++;
		    memcpy(n1->bases, vn[p++], 5*sizeof(int));
		    n1 = path1[x1++];
		}
		l1 = n1;
	    } else {
		// Insertion in path2, link into path1
		if (!n1)
		    n1 = g->node[start]; // if inserting to end of alignment
		int first_ins = 1;
		while (op-- && x2 < np2) {
		    //printf("I  - %2d\n", n2->id);

		    if (first_ins) {
			// Link l1 in to n2
			// Link n2 out to l1
			l1->in[0]->n[1] = n2->id; // FIXME: plus edge hash

			for (j = 0; j < n2->n_out; j++)
			    if (n2->out[j]->n[1] == l1->id)
				break;
			if (j == n2->n_out)
			    // Move out n2->l2 to n2->l1
			    move_edge_out(g, n2, l2, l1);

			memcpy(n2->bases, vn[p++], 5*sizeof(int));
			first_ins = 0;
		    }

		    if (op == 0 || x2 == np2-1) { // last
			// Link n2 in to n1
			// Link n1 out to n2
			// Cull parent(n2) out (to n2)
			node_t *p2 = g->node[n2->in[0]->n[1]];
			n2->in[0]->n[1] = n1->id;
			for (i = 0; i < n1->n_out; i++)
			    if (n1->out[i]->n[1] == n2->id)
				break;
			if (i == n1->n_out)
			    // Move n1->l1 to n1->n2
			    move_edge_out(g, n1, l1, n2);

			if (p2) {
			    for (i = 0; i < p2->n_out; i++) {
				if (p2->out[i]->n[1] == n2->id) {
				    memmove(&p2->out[i], &p2->out[i+1],
					    (p2->n_out - (i+1)) * sizeof(p2->out[0]));
				    p2->n_out--;
				    i--;
				}
			    }
			}

			p++;
		    }

		    //n2->bases[g->kmer-1][4]++;

		    path2[x2] = NULL;
		    l2 = n2;
		    n2 = path2[++x2];
		}

		l1 = l2;
		n2 = n1;
	    }
	}

	// gap at end
	// TODO: Consider just culling these nodes so cigar generation turns
	// into soft-clips.  
	while (x1 < np1) {
	    // Already in path1 only, nothing to do with path2
	    node_t *n1 = path1[x1];
	    //printf("d %2d  -\n", n1->id);
	    n1->bases[g->kmer-1][4]++;
	    n1 = path1[x1++];
	}

	// TODO: Consider just culling these nodes so cigar generation turns
	// into soft-clips.  
	node_t *n1 = g->node[start];
	node_t *n2 = path2[x2];
	if (x2 < np2) {
	    // Insertion of remainder of path2 (1 or more nodes)
	    // into path1 between l1 and n1.
	    // The path2 nodes already form a link, so we only need
	    // to concern ourselves with the head of the path
	    // path2[np2-1] and the current tail path2[x2] (n2).
	    // Note these may be the same node if there is only 1 left.

	    //printf("i  - %2d\n", n2->id);

	    // tail:
	    // link n2 out to l1 in
	    for (i = 0; i < n2->n_out; i++)
		if (n2->out[i]->n[1] == l1->id)
		    break;
	    if (i == n2->n_out)
		move_edge_out(g, n2, l2, l1);
	    for (i = 0; i < l1->n_in; i++)
		if (l1->in[i]->n[1] == n1->id)
		    l1->in[i]->n[1] = n2->id; // FIXME: plus edge hash?

	    x2 = np2-1;
	    n2 = path2[x2];
	    // tail:
	    // link n1 out to n2 in.
	    for (j = 0; j < n1->n_out; j++)
		if (n1->out[j]->n[1] == n2->id)
		    break;
	    if (j == n1->n_out)
		move_edge_out(g, n1, l1, n2);
	}
    }
    
    free(S);
    free(v1);
    free(v2);
    free(vn);
    free(path1);
    free(path2);

    if (0) {
	static int n=0;
	char buf[100];
	sprintf(buf, "_bub%d.dot", n++);
	graph2dot(g, buf, 0);
    }
}

// Any paths that include 'n' and are below 'n_end' have been
// merged into n and can be culled.  We also need to clear
// the visited status for such paths too so we don't detect fake
// bubbles after merging.
void rewind_paths(dgraph_t *g, path_t **phead, int n_id, node_t *n_end) {
    printf("Find redundant paths from %d including node %d\n",
	   n_end->id, n_id);

    path_t *p, *pnext;
    for (p = *phead; p; p = pnext) {
	pnext = p->next;

//	if (p->n)
//	    printf("alive p=%d", p->n->id);
//	else
//	    printf("dead  ");
	node_t *n = p->pn, *ln = NULL;
	while (n && n->id != n_id) {
//	    printf("->%d", n->id);
	    ln = n;
	    n = n->n_in ?  g->node[n->in[0]->n[1]] : NULL;
	}
//	if (n) printf("->%d", n->id);
//	printf("\n");

	if (n && n->id == n_id && ln && ln->id == n_end->id) {
	    printf("    Found path %d->%d .. %d\n", p->n ? p->n->id : -1, p->pn->id, n_id);
	    // Trundle from p->n/p->pn up to n_end clearing visitor flag
	    n = p->n ? p->n : p->pn;
	    while (n && n->id != n_end->id) {
		printf("Clear visited (was %d) on node %d\n", n->visited, n->id);
		n->visited = 0;
		n = n->n_in ?  g->node[n->in[0]->n[1]] : NULL;
	    }

	    // kill it. Better way needed
	    p->n = 0;
	}
    }
}

// Breadth first search starting from id.
// We keep a track of nodes and paths from every fork.
// As we progress one node at a time per path, we check 'visited'.
// If this is set, then we know we have a bubble and can tell
// which paths are involved.
int find_bubble_from2(dgraph_t *g, int id, int use_ref) {
    int p_id = 1;
    int i, v;
    path_t *head = path_create(p_id++);
    int active = 1;
    int ret = 1;

    // Initialise with one path
    //printf("find_bubble_from2 node %d\n", id);
    path_add(head, g->node[id]);
//    printf("A: Node %d visisted by %d\n", id, head->id);
    g->node[id]->visited = head->id;

    // Iterate
    while (head && active) {
	path_t *p, *lp = NULL, *llp, *pnext = NULL;

//	printf("Active paths @");
//	for (p = head; p; p = p->next)
//	    printf(" %d", p->n ? p->n->id : -1);
//	printf("\n");

	active = 0;
	// Step one, move down one layer, adding/culling if needed
	for (p = head; p; p = pnext) {
	    pnext = p->next;
	    node_t *n = p->n;

	    if (!p->n) continue;

	    active = 1;

	    // Leaf node; cull it
	    if (n->n_out == 0) {
		//printf("Cull path ending at %d\n", n->id);
//		if (lp)
//		    lp->next = p->next;
//		else
//		    head = p->next;
//		path_destroy(p);
		p->pn = n;
		p->n = 0;
		continue;
	    }

	    p->pn = n;
	    p->n = g->node[n->out[0]->n[1]];

	    // New fork
	    if (n->n_out > 1) {
		for (i = 1; i < n->n_out; i++) {
		    if (!use_ref && n->out[i]->count == 1 &&
			(n->ref & IS_REF) && (g->node[n->out[i]->n[1]]->ref & IS_REF)) {
			//printf("Skipping ref edge %d to %d\n", n->id, n->out[i]->n[1]);
			continue;
		    }
		    //printf("New path starting at %d\n", n->out[i]->n[1]);
		    path_t *p2 = path_create(p_id++);
		    p2->n = g->node[n->out[i]->n[1]];
		    p2->pn = n;

		    p2->next = p->next;
		    p->next = p2;
		    p = p2;
		}
	    }

	    lp = p;
	}

//	// Report current layer
//	printf("Nodes:");
//	for (p = head; p; p = p->next) {
//	    //node_t *n = p->p[p->np-1];
//	    node_t *n = p->n;
//	    if (!n) continue;
//	    printf(" %d(%d,%d)", n->id, n->in[0]->n[1], p->pn->id);
//	}
//	printf("\n");

	// Step two: check & mark all active nodes for visited status
	for (lp = NULL, p = head; p; lp = p, p = pnext) {
	    pnext = p->next;
	    node_t *n = p->n, *pn = p->pn;
	    
	    if (!n) continue;

	    if (n->visited) {
		int merge_id = n->in[0]->n[1];
//		printf("Bubble detected with node %d (%d, %d)\n",
//		       n->id, merge_id, pn->id);
		//graph2dot(g, "_before.dot", 0);

		//rewind_paths(g, &head, merge_id, n);

		// Note this could form a loop, but if so it gets broken for us.
		node_common_ancestor(g, n, pn, g->node[merge_id]);

		// FIXME: merge forks.
		// Test: just cull one instead.
		if (lp)
		    lp->next = pnext;
		else
		    head->next = pnext;
		if (p != head) path_destroy(p);

		if(0) {
		    path_t *p;
		    printf(" => merged nodes:");
		    for (p = head; p; p = p->next) {
			//node_t *n = p->p[p->np-1];
			node_t *n = p->n;
			if (!n) continue;
			printf(" %d(%d,%d)", n->id, n->in[0]->n[1], p->pn->id);
		    }
		    printf("\n");
		}

		//graph2dot(g, "_after.dot", 0);

		// FIXME: one of these paths may have progressed further beyond
		// the merge point.  We merge right into left, so right one must
		// be followed and clear visitors.  (This may even be several paths
		// stemming from right.)

		// Add p->root for first node in path and p->from so we can find
		// sub paths that forked off us (but only after 'n'?).

		// OR.... clear all visited flags after a merge.
		// Then when finished, if we did any merges start again for
		// a new round.  (Brute force approach.)

		// Taking the easy route - bail out and restart!
		//fprintf(stderr, "visited; err\n");
		goto err;
	    }

	    // Ensure parent node is always incoming[0].
	    if (n->in[0]->n[1] != pn->id) {
		int i;
		for (i = 0; i < n->n_in; i++)
		    if (n->in[i]->n[1] == pn->id)
			break;
		assert(i < n->n_in);
		struct edge *tmp = n->in[0];
		n->in[0] = n->in[i];
		n->in[i] = tmp;
	    }
	    
	    n->visited = p->id;
//	    printf("B: Node %d visisted by %d\n", n->id, p->id);
	}
    }

    ret = 0;
    path_t *p, *pnext;
 err:
    pnext = NULL;
    //fprintf(stderr, "head=%p, next=%p\n", head, head?head->next:NULL);
    for (p = head; p; p = pnext) {
	pnext = p->next;
	path_destroy(p);
    }

    return ret;
}

void find_bubbles(dgraph_t *g, int use_ref) {
    int i, found;

    do {
	found = 0;

	for (i = 0; i < g->nnodes; i++)
	    g->node[i]->visited = 0;

	for (i = 0; i < g->nnodes; i++) {
	    int j;
	    // Needed?
	    for (j = 0; j < g->nnodes; j++)
		g->node[j]->visited = 0;

	    if (g->node[i]->n_in == 0 && !g->node[i]->pruned) {
		//printf("Graph start at %d\n", i);

		int b = find_bubble_from2(g, i, use_ref);

		if (b) {
		    found += b;
		    loop_check(g, 1);
		}
		//break; // only really need main start point?  Unknown...
	    }
	}
    } while (found);
}

// Should be run on a collapsed graph without bubbles.
// Any forks should only have one route that follows the reference.
void find_insertions(dgraph_t *g) {
    int i, j;

    // Find ref boundaries
    int r_min = INT_MAX, r_max = INT_MIN;
    for (i = 0; i < g->nnodes; i++) {
	node_t *n = g->node[i];
	if (n->pruned || n->pos == INT_MIN)
	    continue;
	if (r_min > n->pos)
	    r_min = n->pos;
	if (r_max < n->pos)
	    r_max = n->pos;
    }
    if (r_min == INT_MAX)
	return;

    fprintf(stderr, "Ref covered %d to %d\n", r_min, r_max);

    for (i = 0; i < g->nnodes; i++)
	g->node[i]->visited = 0;

    for (i = 0; i < g->nnodes; i++) {
	if (/*g->node[i]->n_in != 0 ||*/ g->node[i]->pruned || g->node[i]->visited)
	    continue;

	// Found a starting point, so scan forward from here.
	node_t *n = g->node[i];
	while (n && n->pos != INT_MIN) {
	    n->visited = 1;
	    for (j = 0; j < n->n_out; j++) {
		if (g->node[n->out[j]->n[1]]->pos != INT_MIN)
		    break;
	    }
	    if (j == n->n_out)
		break;

	    n = g->node[n->out[j]->n[1]];
	}

	if (n->pos == INT_MIN)
	    continue;

	int pass, ins_count = 0, end_pos = 0;
	for (pass = 0; pass < 2; pass++) {
	    node_t *nr_1 = n;    // last-ref and 1st-insertion
	    node_t *ni_1 = NULL;
	    for (j = 0; j < n->n_out; j++) {
		node_t *t = g->node[n->out[j]->n[1]];
		if (t->ref & IS_CON)
		    ni_1 = t;
	    }

	    if (!ni_1)
		continue;

	    // Find the other end of the insertion, if it exists
	    node_t *ni_2 = ni_1, *nr_2 = NULL;
	    do {
		if (ni_2->visited)
		    break;

		if (pass > 0) {
		    ni_2->ins = ++ins_count;
		    ni_2->pos = end_pos;
		}

		node_t *tmp = NULL;
		for (j = 0; j < ni_2->n_out; j++) {
		    node_t *t = g->node[ni_2->out[j]->n[1]];
		    if (t->pos != INT_MIN) {
			nr_2 = t;
			break;
		    } else if ((t->ref & IS_CON) || !n) {
			tmp = t;
		    }
		}
		if (nr_2)
		    break;
		ni_2 = tmp;
	    } while (ni_2);

	    if (!nr_2)
		continue;

	    fprintf(stderr, "ins from %d(%d) to (%d)%d\n",
		    nr_1->id, ni_1->id, ni_2->id, nr_2->id);
	    end_pos = nr_2->pos;
	}
    }
}


int tag_tail_from(dgraph_t *g, int parent, int id) {
    printf("Check tail at node %d->%d\n", parent, id);

    node_t *n = g->node[id], *l = n;
    if (n->visited == 't') {
	printf("Already checked: tail=%d\n", n->is_tail);
	return n->is_tail;
    }

    while (n) {
	printf("\t%d\n", n->id);
	if (n->n_out == 0) {
	    int ret;
	    n->is_tail = 1;
	    n = g->node[id];
	    if (n->n_in > 1 && !n->is_head) {
		printf("TAIL at %d++ ->%d\n", id, n->id);
		printf("But %d failed incoming check\n", n->id);
		n = g->node[l->out[0]->n[1]];
		ret = 0;
	    } else {
		printf("TAIL at %d ->%d\n", id, n->id);
		n = l;
		ret = 1;
	    }

	    while (n) {
		n->is_tail = 1;
		if (n->n_out != 1)
		    break;
		n = g->node[n->out[0]->n[1]];
	    }

	    return ret;
	}
	if (n->visited == 't')
	    break;
	n->visited = 't';
	if (n->n_out > 1)
	    break;

	n = g->node[n->out[0]->n[1]];
    }

    int i, t = 1;
    for (i = 0; i < n->n_out; i++) {
	int x = tag_tail_from(g, n->id, n->out[i]->n[1]);
	printf("tag_tail_from %d,%d -> %d\n", n->id, n->out[i]->n[1], x);
	if (x == 0 || (n->n_in > 1 && !n->is_head)) {
	    printf("%d->%d: but failed incoming check\n", n->id, n->out[i]->n[1]);
	    t = 0;
	} else {
	    printf("%d->%d: possible tail\n", n->id, n->out[i]->n[1]);
	}
    }
    
    if (t)
	n->is_tail = 1;
    printf("%s at %d\n", t ? "tail" : "non-tail", n->id);
    return t;
}

void tag_head_from(dgraph_t *g, int parent, int is_head, int id) {
    node_t *n = g->node[id];
    node_t *f = n, *l = n;

    //printf(">>>%d,%d,%d\n", parent, is_head, id);

    // Walk graph
    while (n) {
	n->visited = (n->visited & ~0xff) | 'h';
	n->visited += 256;

	//printf("Check node %d, visited %d of %d times\n", n->id, n->visited >> 8, n->n_in);
	if (n->n_in > 1 && (n->visited >> 8) < n->n_in) {
	    //printf("Stop1 at node %d\n", n->id);
	    return;
	}
	if (n->n_out > 1) {
	    //printf("Stop2 at node %d\n", n->id);
	    break;
	}

	if (is_head) {
//	    printf("Node %d is hair\n", n->id);
	    n->is_head = 1;
	}

	if (n->n_out == 0)
	    break;

	l = n;
	n = g->node[n->out[0]->n[1]];
    }
}

// Recursively find hairs (end points) and tag them so we
// know which nodes lead to dead ends (ie have no bubbles).
//
// FIXME: not bullet proof
void tag_hairs(dgraph_t *g) {
    int i;
    for (i = 0; i < g->nnodes; i++) {
	if (g->node[i]->n_in == 0 && !g->node[i]->pruned) {
	    tag_head_from(g, -1, 1, i);
	}
    }    

    for (i = 0; i < g->nnodes; i++) {
	if (g->node[i]->n_in == 0 && !g->node[i]->pruned) {
	    printf("Hair start at %d\n", i);
	    tag_tail_from(g, -1, i);
	}
    }    
}


#if 0
// Merge nodes n1 onto n2, with relationships n1->c1 and n2->c2.
// NB: c2 is unrequired, so not specified as an argument.
// We are culling node n1, so c1 incoming is removed.
// (FIXME: see dup1 comments regarding potential duplicated code).
int merge_head_node(dgraph_t *g, node_t *n1, node_t *n2, node_t *c1) {
    int i, j;

    if (n2->n_in + n1->n_in > MAX_EDGE || n2->n_out + n1->n_out > MAX_EDGE)
	return -1;
    HashItem **hi_ = realloc(n2->hi, (n2->n_hi + n1->n_hi) * sizeof(*n2->hi));
    if (!hi_)
	return -1;
    n2->hi = hi_;

    // Discard n1, keeping n2
    n1->pruned = 1;

    // various counts
    for (j = 0; j < g->kmer; j++)
	for (i = 0; i < 5; i++)
	    n2->bases[j][i] += n1->bases[j][i];
    n2->count += n1->count;

    // Migrate hash keys pointing to n1 to n2
    for (i = 0, j = n2->n_hi; i < n1->n_hi; i++, j++) {
	n2->hi[j] = n1->hi[i];
	n2->hi[j]->data.p = n2;
    }

    // Merge in[] arrays; no loops exist, so a pure tree.
    for (i = 0; i < n1->n_in; i++) {
	int in1 = n1->in[i]->n[1];
	// Unnecessary loop - here for bug checking only.
	for (j = 0; j < n2->n_in; j++) {
	    if (in1 == n2->in[j]->n[1])
		break;
	}
	if (j != n2->n_in)
	    abort(); // NB: shouldn't happen; implies not a tree.

	// n1's parent out edge now links to n2 instead of n1
	move_edge_out(g, g->node[in1], n1, n2);

	// Add n2's in edge to n1.  n2 has it, but is defunct.
	n1->in[i]->n[0] = n2->id;
	n2->in[n2->n_in++] = n1->in[i];
    }

    // Cull c1's in edge from n1.
    for (i = j = 0; i < c1->n_in; i++) {
	if (c1->in[i]->n[1] == n1->id)
	    continue;
	c1->in[j++] = c1->in[i];
    }
    c1->n_in = j;

    return 0;
}

// N1 and n2 have a common child c.  Insert n1 between n2 and c,
// Ie we have c->n1->n2 instead of c->n1 and c->n2.
int ins_head_node(dgraph_t *g, node_t *n1, node_t *n2, node_t *c) {
    int i, j;

    if (n1->n_in+1 >= MAX_EDGE || n1->n_out >= MAX_EDGE)
	return -1;

    // Add pad based on c->n2 frequency
    int *d = n2->bases[g->kmer-1];
    int depth = d[0] + d[1] + d[2] + d[3] + d[4];
    for (j = 0; j < g->kmer; j++)
	n1->bases[j][4] += depth;

    // Move out edge c->n2 to n1->n2
    move_edge_out(g, n2, c, n1);
    
    // The above doesn't move the incoming edge linking
    // n2->c, so we do that manually.
    for (i = j = 0; i < c->n_in; i++) {
	if (c->in[i]->n[1] == n2->id) {
	    n1->in[n1->n_in++] = c->in[i];
	    c->in[i]->n[0] = n1->id;
	} else {
	    c->in[j++] = c->in[i];
	}
    }
    c->n_in = j;

    return 0;
}

// Merge n_x1 nodes from n1 into n2 according to alignment S.
int merge_head_tip(dgraph_t *g, node_t *n1, node_t *n2, node_t *c, int n_x1, int *S,
		   char *ref_r, char *hseq) { // ref_r/hseq for debugging only
    // Worth while doing a merge.
    int x1 = 0, x2 = 0;
    int i, ret = -1;

    // NB: n_x1 is potentially up to g->kmer-1 too high.
    while (n1 && n2 && x1 < n_x1) {
	int op = *S++;
	if (op == 0) {
	    // match
	    //printf("MERGE\t%d/%c %d/%c  %d/%d\n", x1, hseq[x1], x2, ref_r[x2], n1->id, n2->id);
	    if (merge_head_node(g, n1, n2, c) < 0)
		return -1;
	    x1++; x2++;
	    n1 = n1->n_in ? g->node[n1->in[0]->n[1]] : NULL;
	    c = n2;
	    node_t *new_n2 = NULL;
	    for (i = 0; i < n2->n_in; i++) {
		if (g->node[n2->in[i]->n[1]]->ref) {
		    new_n2 = g->node[n2->in[i]->n[1]];
		    break;
		}
	    }
	    if (!new_n2)
		break; // not in ref, so S alignment matrix is invalid from here on.
	    ret = 0;
	    n2 = new_n2;
	} else if (op < 0) {
	    // ins in cons
	    while (n1 && op++) {
		node_t *new_n1;
		new_n1 = n1->n_in ? g->node[n1->in[0]->n[1]] : NULL;
		//printf("INS\t%d/%c %d/-  %d/%d\n", x1, hseq[x1], x2, n1->id, n2->id);
		if (ins_head_node(g, n1, n2, c) < 0)
		    return -1;
		x1++;
		c = n1;
		n1 = new_n1;
		ret = 0;
	    }
	} else {
	    // del in cons
	    while (n2 && op--) {
		//printf("SKIP\t%d/- %d/%c  %d/%d\n", x1, x2, ref_r[x2], n1->id, n2->id);
		x2++;
		for (i = 0; i < n2->n_in; i++) {
		    if (g->node[n2->in[i]->n[1]]->ref) {
			n2 = g->node[n2->in[i]->n[1]];
			break;
		    }
		}
		if (i == n2->n_in)
		    n2 = NULL;
	    }
	}
    }

    return ret;
}

/*
 * Scan along reference finding nodes that are incoming.  As no
 * bubbles, all of these represent head-tips.  Scan back up that
 * produce the sequence, align that to the reference seq and then
 * fold the head tip back into the main branch.  (Anything remaining is
 * to be marked as a soft-clip.)
 */
int merge_head_tips(dgraph_t *g, char *ref, int len) {
    int merged = 0;
    int i, j;

    for (i = 0; i < len - g->kmer; i++) {
	node_t *h = find_node(g, ref+i, g->kmer, 0);
	if (!h || h->n_in <= 1)
	    continue;

	int ref_len = i+g->kmer-1; // prefix length to match against

	// More than one node => head tip
	//printf("Head tip at node %d\n", h->id);

	int hlen = i+10, hidx = 0;
	int (*head)[5] = malloc(hlen * sizeof(*head));
	if (!head)
	    return 0;

	node_t *rn = NULL, *n = NULL;
	for (j = 0; j < h->n_in; j++) {
	    if (g->node[h->in[j]->n[1]]->ref == 0 && !n) {
		n = g->node[h->in[j]->n[1]];
	    } else if (g->node[h->in[j]->n[1]]->ref) {
		rn = g->node[h->in[j]->n[1]];
	    }
	}

	if (!n || !rn)
	    continue;

	//printf("Check node prefix from %d and %d\n", n->id, rn->id);

	node_t *last = n;
	node_t *n1 = n;
	// Backtrack from here to unvisited non-ref branches.
	// We stop at the first branch and merge that in a subsequent round,
	// if at all.
	while (n && hidx < hlen) {
	    memcpy(&head[hidx++], n->bases[g->kmer-1], sizeof(*head));
	    last = n;
	    n = n->n_in == 1 ? g->node[n->in[0]->n[1]] : NULL;
	}
	for (j = g->kmer-2; j >= 0 && hidx < hlen; j--)
	    memcpy(&head[hidx++], last->bases[j], sizeof(*head));
	
//	// Reverse the list; shared code with assign_coords
//	int (*v1)[5] = head;
//	int I, J, K;
//	for (I = 0, J = hidx-1; I < J; I++, J--) {
//	    for (K = 0; K < 5; K++) {
//		int tmp = v1[I][K];
//		v1[I][K] = v1[J][K];
//		v1[J][K] = tmp;
//	    }
//	}

	// Reverse the 'ref' seq instead as we want to process alignment
	// from right to left anyway.
	char *ref_r = malloc(ref_len+1);
	int k;
	for (k = 0; k < ref_len; k++)
	    ref_r[ref_len-1-k] = ref[k];
	ref_r[ref_len] = 0;

	char *hseq = vec2seq(head, hidx);
	//printf("%.*s\n", ref_len, ref_r);
	//puts(hseq);

	int *S = malloc((ref_len + hidx) * sizeof(*S));

	// Penalise left edge but not right edge; 1 = free, 0 = costs.
	// Note this is left edge of reversed seq, so it would be the right
	// edge in the graph.
	int score = align_ss(hseq, ref_r, hidx, ref_len, 0, 0, X128, 3,1, S, 0,1,0,1);
	//printf("Score=%d\n", score);
	//display_ss(hseq, ref_r, hidx, ref_len, S, 0, 0);
	//fflush(stdout);

	// Scan along the graph finding the best scoring point.
	int x1 = 0, x2 = 0, *S2 = S;
	int best_score = 0;
	int best_x1 = 0;
	score = 0;
	while (x1 < hidx && x2 < ref_len) {
	    int op = *S2++;
	    if (op == 0) {
		// match
		score += X128[hseq[x1]][ref_r[x2]]-2;
		//printf("%d\t%d/%c %d/%c  %d\n", score, x1, hseq[x1], x2, ref_r[x2], X128[hseq[x1]][ref_r[x2]]);
		x1++; x2++;
	    } else if (op < 0) {
		// ins in cons
		while (op++) {
		    score += islower(hseq[x1]) ? 0 : -4;
		    //printf("%d\t%d/%c %d/-\n", score, x1, hseq[x1], x2);
		    x1++;
		}
	    } else {
		// del in cons
		while (op--) {
		    score += islower(ref_r[x2]) ? 0 : -4;
		    //printf("%d\t%d/- %d/%c\n", score, x1, x2, ref_r[x2]);
		    x2++;
		}
	    }
	    if (best_score <  score) {
		best_score = score;
		best_x1 = x1; // before this but not including this pos
	    }
	}
	//printf("Best x1 = %d / %d\n", best_x1, best_score);

	if (best_score > 0)
	    merged += (merge_head_tip(g, n1, rn, h, best_x1, S, ref_r, hseq) >= 0);

	free(hseq);
	free(ref_r);
	free(S);
    }

    return merged;
}

int merge_tail_tips(dgraph_t *g) {
    return 0;
}
#endif

void prune(dgraph_t *g, int min_count) {
    int done_something;

    //    break_weak(g, min_count); FIXME for ref
    //burst_bubbles(g, min_count);

    prune_tails(g, min_count);
    prune_heads(g, min_count);

    // Leave bubbles if we're using this to realign BAM?
    // Theory is if they link back into the main graph then pruning them will
    // make realignment back a complicated process.  We're better to only
    // prune hairs (head/tails) as they're soft-clip candidates, but we cannot
    // soft-clip the internals.

    //prune_bubbles(g, min_count);
}

void validate_graph(dgraph_t *g) {
    int i, j, k;

    for (i = 0; i < g->nnodes; i++) {
	node_t *n = g->node[i];
	if (n->pruned)
	    continue;

	for (j = 0; j < n->n_in; j++) {
	    //assert(n->in[j]->n[0] == n->id);
	    node_t *p = g->node[n->in[j]->n[1]];
	    for (k = 0; k < p->n_out; k++)
		if (p->out[k]->n[1] == n->id)
		    break;
	    assert(k < p->n_out);
	}


	for (j = 0; j < n->n_out; j++) {
	    //assert(n->out[j]->n[0] == n->id);
	    node_t *p = g->node[n->out[j]->n[1]];
	    for (k = 0; k < p->n_in; k++)
		if (p->in[k]->n[1] == n->id)
		    break;
	    assert(k < p->n_in);
	}
    }
}

#if 1
int graph2dot(dgraph_t *g, char *fn, int wrap) {
    return 0;
}
#else
int graph2dot(dgraph_t *g, char *fn, int wrap) {
    FILE *fp = stdout;

    if (fn) {
	fp = fopen(fn, "w");
	if (!fp) {
	    perror(fn);
	    return -1;
	}
    }

    int i, j;
    fprintf(fp, "digraph g {\n");
    //fprintf(fp, "  graph [overlap=scale]\n");
    //fprintf(fp, "  node [shape=plaintext]\n");
    //fprintf(fp, "  node [shape=point]\n");
    //printf("  node [style=invis]\n");
    for (i = 0; i < g->nnodes; i++) {
	node_t *n = g->node[i];
	if (n->pruned)
	    continue;

//	if (wrap) {
//	    fprintf(fp, "  n_%d [label=\"%d@%d/", n->id, n->id, n->pos);
//	    int j;
//	    int inc = sqrt(n->len)*2 < 20 ? 20 : sqrt(n->len)*2;
//	    for (j = 0; j < n->len; j += inc) {
//		fprintf(fp, "%.*s\\n", n->len - j < inc ? n->len-j : inc, n->seq + j);
//	    }
//	    fprintf(fp, "\"]\n");
//	} else {
//	    fprintf(fp, "  n_%d [label=\"%d@%d/%c\"]\n",
//		    n->id, n->id, n->pos, n->seq[n->len-1]);
//	}

//	char bases[5];
//	for (j = 0; j < 5; j++) bases[j] = n->bases[j] + '0';
//	fprintf(fp, "  n_%d [label=\"%d/%d/%.5s\", color=\"%s\", penwidth=3, group=%d];\n",
//		n->id, n->id, n->pos, bases,
//		n->is_head ? "green3" : n->is_tail ? "skyblue" : "black",
//		n->id);


	fprintf(fp, "  n_%d [label=\"%d @ %d x %d", n->id, n->id, n->pos, n->count);
	//if (n->ins) fprintf(fp, ".%d", n->ins);
	//fprintf(fp, "  n_%d [label=\"%d", n->id, n->id);
	//for (j = n->n_in ? g->kmer-1 : 0; j < g->kmer; j++) {
	if (n->posa) {
	    for (j = 0; j < g->kmer; j++) {
		int k;
		int m = 0, b = 0;
		for (k = 0; k < 5; k++)
		    if (m <= n->bases[j][k])
			m = n->bases[j][k], b=k;
		//if (g->kmer <= 12 || j%2==0)
		fprintf(fp, "\\n%c ", "ACGT*"[b]);
		//else
		//	fprintf(fp, "   %c ", "ACGT*"[b]);
		for (k = 0; k < 5; k++)
		    fputc(n->bases[j][k] + '0', fp);
		fprintf(fp, " @ %d", n->posa ? n->posa[j] : 1000 + n->pos + j - (g->kmer-1));
	    }
	}

//	for (j = 0; j < n->n_hi; j++)
//	    fprintf(fp, "\\n%.*s", n->hi[j]->key_len, n->hi[j]->key);
	fprintf(fp, "\", color=\"%s\", penwidth=3];\n",
		n->is_head ? "green3" : n->is_tail ? "skyblue" : "black");
	//fprintf(fp, "  n_%d [label=\"%d/%c\"];\n", n->id, n->id, n->seq[n->len-1]);
	//fprintf(fp, "  n_%d [label=\"%d/%d\"];\n", n->id, n->id, n->len);
	for (j = 0; j < n->n_out; j++) {
//	    if (n->out[j]->count == 0)
//		continue;
	    int is_ref = (n->ref && g->node[n->out[j]->n[1]]->ref);
	    fprintf(fp, "  n_%d -> n_%d [penwidth=%f,label=\"%d\",edges=%d,color=\"%s\"]\n",
		    n->id, n->out[j]->n[1],
		    //n->out[j]->count*1.0,
		    //(sqrt(n->out[j]->count)-1)*2+1,
		    log(n->out[j]->count+1)+2,
		    n->out[j]->count,
		    n->n_out,
		    is_ref ? n->out[j]->count == 1 ? "gold" : "red1" : "black");
//	    fprintf(fp, "  n_%d -> n_%d [penwidth=%f,label=\"%c\",edges=%d]\n",
//		    n->id, n->out[j]->n[1],
//		    0.2,//(sqrt(n->out[j]->count)-1)*2+1,
//		    n->seq[n->len-1],
//		    n->n_out);
	}

#if 0
	// incoming edges
	for (j = 0; j < n->n_in; j++) {
	    fprintf(fp, "  n_%d -> n_%d [color=\"blue\"]\n", n->id, n->in[j]->n[1]);
	}
#endif
    }

#ifndef MIN
#  define MIN(a,b) ((a)<(b)?(a):(b))
#endif

#define MAX_GRAPH_KEY 100

#if 1
	// Hash key -> node map
	HashIter *iter = HashTableIterCreate();
	HashItem *hi;
	while (hi = HashTableIterNext(g->node_hash, iter)) {
	    node_t *n = hi->data.p;
	    fprintf(fp, " n_%.*s [label=\"%.*s\", color=\"grey\", group=%d]\n",
	    	    hi->key_len, hi->key,
		    MIN(MAX_GRAPH_KEY, hi->key_len), hi->key, n->id);
	    fprintf(fp, " n_%.*s -> n_%d [color=\"grey\"]\n",
	    	    hi->key_len, hi->key, n->id);
	    //fprintf(fp, " n_%.5s [label=\"%.*s\", color=\"grey\", group=%d]\n",
	    //	    hi->key_len-5 + hi->key, hi->key_len, hi->key, n->id);
	    //fprintf(fp, " n_%.5s -> n_%d [color=\"grey\"]\n",
	    //	    hi->key_len-5 + hi->key, n->id);
	}

	HashTableIterDestroy(iter);
#endif

	fprintf(fp, "}\n");


	validate_graph(g);
	if (fn)
	    return fclose(fp);

	return 0;
}
#endif

typedef struct hseq {
    struct hseq *next;
    char *seq;  // ambig for alignments
    char *seq2; // real kmers, for matching to graph
    int len;
    int score;
} hseqs;

// Given a haplotype sequence, a reference sequence and the alignment between
// them, assign a ref coordinate to each kmer from that haplotype indicating
// the location (of the terminating base) in the reference.
//
// Where the haplotype has insertions, we use pos -1 to indicate ins.
// Where the haplotype has a deletion, we just skip that pos.
//
// The plan is we can then go through the debruijn construction process again
// with the reads, but instead of constructing the graph we use the updated
// graph to construct a new CIGAR string for each read.
//
// TODO: This is just all hap vs ref one at a time.  This may not give
// sensible hap to hap alignment though.  For that we need a true multiple
// sequence alignment (including the reference).
void pad_ref_hap(dgraph_t *g, hseqs *h, haps_t *ref, int *S) {
    int i = 0, j = 0;
    char *B = h->seq2, *A = ref->seq;
    int N = h->len, M = strlen(ref->seq);

    int *pos = calloc(N + g->kmer, sizeof(int));

    // Pos is hooked to last base in KMER.
    while (i < M || j < N) {
	int op = *S++;
	if (op == 0) {
	    // match/mismatch
	    pos[j] = i;
	    if (j >= g->kmer-1) {
		node_t *n = find_node(g, B + j-(g->kmer-1), g->kmer, 0);
		if (n) {
		    printf("%3d %3d  %c (%d)\n", j, pos[j], B[j], n->id);
		    //putchar('M');
		    n->pos = i;
		    // FIXME: what if n->posa is already set?
		    // We may want to prioritise the highest count route.
		    n->posa = &pos[j-(g->kmer-1)];
		} else {
		    fprintf(stderr, "%3d  %c: Failed to find tail? node %.*s\n", i, B[j],
			    g->kmer, B+j-(g->kmer-1));
		}
	    } else {
		printf("%3d %3d  %c (before 1st KMER)\n", j, pos[j], B[j]);
	    }
	    i++, j++;
	} else if (op > 0) {
	    // ins to read
	    //j += op;
	    int ilen = op;
	    while (op) {
		pos[j] = op-ilen-1;
		printf("%3d %3d +%.*s\n", j, pos[j], op, &B[j]);
		//pos[j] = -1;
		node_t *n = j >= g->kmer-1 ? find_node(g, B + j-(g->kmer-1), g->kmer, 0) : NULL;
		if (n) {
		    n->pos = -1;
		    n->posa = &pos[j-(g->kmer-1)];
		}
		op--;
		j++;
		//putchar('I');
	    }
	} else if (op < 0) {
	    // del in read.
	    printf("%3d %3d -%.*s\n", j, pos[j], -op, &A[i]);
	    i -= op;
	}
    }
    //printf("\n");
}


#define ADD_CIGAR(op,len)						\
    do {								\
        if (cig_op != op) {						\
	    if (cig_op != 999) {					\
		cig_a[cig_ind++] = cig_op | (cig_len << BAM_CIGAR_SHIFT); \
		if (cig_ind > MAX_CIGAR) return -1;			\
	    }								\
            cig_len = 0;						\
            cig_op = op;						\
        }								\
        cig_len+=len;							\
    } while(0)


#define MAX_CIGAR 65535

// Taken from samtools/padding.c.
// FIXME: even though it's a direct copy from samtools, it needs fixing to
// return a value and check mallocs.
static void replace_cigar(bam1_t *b, int n, uint32_t *cigar)
{
    if (n != b->core.n_cigar) {
        int o = b->core.l_qname + b->core.n_cigar * 4;
        if (b->l_data + (n - b->core.n_cigar) * 4 > b->m_data) {
            b->m_data = b->l_data + (n - b->core.n_cigar) * 4;
            kroundup32(b->m_data);
            b->data = (uint8_t*)realloc(b->data, b->m_data);
        }
        memmove(b->data + b->core.l_qname + n * 4, b->data + o, b->l_data - o);
        memcpy(b->data + b->core.l_qname, cigar, n * 4);
        b->l_data += (n - b->core.n_cigar) * 4;
        b->core.n_cigar = n;
    } else memcpy(b->data + b->core.l_qname, cigar, n * 4);
}


// seq2cigar based on the newer find_bubbles and common_ancestor output.
int seq2cigar_new(dgraph_t *g, char *ref, int shift, bam1_t *b, char *seq, int *new_pos) {
    int i;
    node_t *n = NULL, *last = NULL;
    int cig_op = 999, cig_len = 0;
    int first_go = 1;
    int seq_start = 0;
    char *orig_seq = seq;
    int len = strlen(seq);

    char *sub = malloc(g->kmer + len + 1), *sub_k = sub + g->kmer;
    memcpy(sub_k, seq, len);
    sub_k[len] = 0;

    int cig_a[MAX_CIGAR];
    int cig_ind = 0;

    //fprintf(stderr, "Orig seq = %.*s\n", len, seq);

    int last_ins = 0;
    int since_mat = 0;
 try_again:
    // First node
    for (i = seq_start; i < len-g->kmer; i++) {
	if ((n = find_node(g, seq+i, g->kmer, 0)) && n->pos >= 0 && n->pruned == 0)
	    break;
    }

    if (!n) {
	fprintf(stderr, "No match found for seq\n");
	goto unmapped;
    }

//    if (n->pos < 0) {
//	// No mapped kmer, so start with the first unmapped one instead
//	for (i = seq_start; i < len-g->kmer; i++) {
//	    if ((n = find_node(g, seq+i, g->kmer, 0)) && n->pruned == 0)
//		break;
//	}
//	if (n->pruned) {
//	    goto unmapped;
//	}
//    }

    int pos = 0, j, padded = 0;

    if (n->pos >= 0 /*&& (first_go || n->ins == 0)*/ /*|| (n->ref & IS_CON)*/) {
//	int np_dist = 0;
//	if (n->pos < 0) {
//	    // last base of 1st kmer is in an insertion.  Backtrack to see
//	    // if the sequence keeps matching the consensus (todo - we need
//	    // to check the actual sequence does!).
//	    node_t *np = n;
//	    while (np->n_in > 0 /*&& np->pos < 0*/ && (np->ref & IS_CON)) {
//		np = g->node[np->in[0]->n[1]];
//		np_dist++;
//	    }
//
//	    // Did we come out of the insertion and into data mapped to a ref?
//	    // If so we have xMyI CIGAR.
//	    if (np->pos > 0) {
//		pos = np->pos-1;
//		ADD_CIGAR(...
//	    }
//	    
//	} else {
	    pos = n->pos-1;
//	}

	if (first_go) {
	    memset(sub, 'N', g->kmer);
	    memset(sub_k, 'N', i);
	    sub_k += i;
	
	    // First node found has kmer bases, but we only know the position of
	    // the last base.  Greedy match backwards to determine the correct
	    // path for earlier bases.  We simply scan upwards on our (now linear)
	    // graph to find the precursor, replacing Ns with the base found.
	    // This then permits us to start the process again, but this time
	    // using the known prefix sequence to anchor all bases rather than
	    // just the last base in the kmer from our hash hit.
	    //printf("\nFound %.*s => %p, subset of %.*s\n", g->kmer, seq+i, n, g->kmer*2, sub);
	    for (j = 1; j < g->kmer; j++) {
		int k;
		node_t *sn = NULL;
	    up_one:
		//printf("Search for %.*s", g->kmer, sub_k-j);
		for (k = 0; k < n->n_in; k++) {
		    node_t *s2 = g->node[n->in[k]->n[1]];
		    int h;
		    for (h = 0; h < s2->n_hi; h++) {
			if (memcmp(sub_k-j+1, s2->hi[h]->key+1, g->kmer-1) == 0) {
			    // Matches bar first N, fix base
			    *(sub_k-j) = *s2->hi[h]->key;
			    sn = s2;
			    break;
			}
		    }
		    if (sn)
			break;
		}
		if (!sn) {
		    //printf(" => not found, try one up\n");
		    n = best_prev_node(g, n);
		    if (n)
			goto up_one;
		} else {
		    n = sn;
		}

		if (!n)
		    break;
	    }

	    seq_start = -(j-1)+i;
	    seq = sub_k-i;
	    //printf("New seq = %.*s\n", len - seq_start, seq + seq_start);
	    first_go = 0;
	    goto try_again;
	}

	// Start of read is unaligned.  Align back to ref or just soft-clip.
	//b->core.pos = pos+shift;
	*new_pos = pos+shift;
	if (i + g->kmer-1 > 0) {
	    int sc = i + g->kmer-1;
	    ADD_CIGAR(BAM_CSOFT_CLIP, sc);
	}

	// Trace path for each successive node.
	char *s1 = seq+i;
	char s2[MAX_KMER*2];
	for (; i <= len - g->kmer; i++) {
	    last = n;
	    if (i == len-g->kmer) {
		memcpy(s2, seq+i, g->kmer);
		memcpy(s2+g->kmer, ref+pos+2, g->kmer);
		s1 = s2;
	    }
	    if (!(n = find_node(g, s1++, g->kmer, 0)) || n->pruned)
		break;

	    if (n->pos == pos+1 && n->ins == 0) {
		ADD_CIGAR(BAM_CMATCH, 1);
		pos++;
		last_ins = 0;
	    } else if (n->pos > pos+1){
		ADD_CIGAR(BAM_CDEL, n->pos - (pos+1));
		if (n->ins) {
		    if (n->ins > 1)
			ADD_CIGAR(BAM_CPAD, n->ins-1);
		    ADD_CIGAR(BAM_CINS, 1);
		    last_ins = n->ins;
		    pos = n->pos-1;
		} else {
		    ADD_CIGAR(BAM_CMATCH, 1); // Is this correct?  Not in multiple alignment?
		    last_ins = 0;
		    pos = n->pos;
		}
		//} else if (n->pos <= -1 && n->pos != INT_MIN) {
	    } else if (n->pos == pos+1 && n->ins) {
		// An unmapped base in reference; insertion, but possibly
		// also a node prior to this that we didn't observe in our
		// kmer stepping? (ie deletion)

		// bubble up to check if D or P before I
		node_t *np = NULL, *nn = n;
		int path[MAX_SEQ] = {0}, path_ind = 0;
		if (n->ins > last_ins+1)
		    ADD_CIGAR(BAM_CPAD, n->ins - (last_ins+1));
	    
		ADD_CIGAR(BAM_CINS, 1);
		last_ins = n->ins;
		// FIXME: mismatches at the end need to be soft-clips and
		// not insertions.  See eg4c.{dat,ref}
		//
		// Also.. need to merge tails in so they get folded in to
		// the main branch and become matches as far as they can.
		// (Needed for long kmers.)

	    } else if (n->pos == INT_MIN) {
		// An unmapped base in reference; insertion, but possibly
		// also a node prior to this that we didn't observe in our
		// kmer stepping? (ie deletion)

		// bubble up to check if D or P before I
		node_t *np = NULL, *nn = n;
		int path[MAX_SEQ] = {0}, path_ind = 0;
		// Recheck all of this.
		// It dates back to before we added 2D coordinates of
		// Nth base inserted at Mth ref pos.
		for (;last && last != np && path_ind < MAX_SEQ; nn=np) {
		    np = nn->n_in ? g->node[nn->in[0]->n[1]] : NULL;
		    if (!np)
			break;
		
		    // Track path of which out[x] leads from last to n.
		    int i;
		    for (i = 0; i < np->n_out; i++)
			if (np->out[i]->n[1] == nn->id)
			    break;
		    path[path_ind++] = i;
		}

		if (path_ind == MAX_SEQ)
		    goto fail;

		// Now replay the path in order from last->n
		np = last;
		if (--path_ind > 0 && np->n_out > path[path_ind])
		    np = g->node[np->out[path[path_ind]]->n[1]];
		while (--path_ind >= 0) {
		    //printf("Path %d\n", np->out[path[path_ind]]->n[1]);
		    if (np->pos >= 0 && !np->ins) {
			ADD_CIGAR(BAM_CDEL, 1);
			pos = np->pos;
		    } else {
			ADD_CIGAR(BAM_CPAD, 1);
		    }
		    if (path[path_ind] >= np->n_out)
			goto fail; // under what scenario does this happen?
		    np = g->node[np->out[path[path_ind]]->n[1]];
		}

		ADD_CIGAR(BAM_CINS, 1);
	    } else {
		// soft clip from here on.
		// OR... an internal mismatch that was pruned.
		// We need to realign pruned reads back to the graph?
		break;
	    }
	}

	if (cig_op == BAM_CINS && n && n->ins == 0) {
	    fprintf(stderr, "ending on ins, n->ins=%d last_ins=%d cig_len=%d\n", n->ins, last_ins, cig_len);
	    // partial match to known insertion plus partial mismatch => softclip?
	    //
	    // Eg should be 5M 10I, but want 5M 8I 2S due to seq error
	    // in 2nd to last base.
	    if (last_ins > 1 && last_ins < cig_len) {
		int sc_len = cig_len - last_ins;
		cig_len = last_ins;
		ADD_CIGAR(BAM_CSOFT_CLIP, sc_len);
	    } else if (last_ins <= 1) {
		cig_op = BAM_CSOFT_CLIP;
	    }
	    //fprintf(stderr, "%d,%d: ins at %d, last %d, last_ins %d\n", n->id, last->id, n->ins, last->ins, last_ins);
	} else if (cig_op == BAM_CINS) {
	    // insertion and no node means we ran out part way.  Change to soft-clip
	    // as it's likely not an insertion at all, just a region with no
	    // alignment to ref.
	    cig_op = BAM_CSOFT_CLIP;
	}

    fail:
	if (i+g->kmer-1 < len) {
	    // trailing unaligned. FIXME: align it or soft-clip as appropriate.
	    ADD_CIGAR(BAM_CSOFT_CLIP, len-(i+g->kmer-1));
	}

	ADD_CIGAR(999, 0); // flush
	b->core.flag &= ~BAM_FUNMAP;
    } else {
    unmapped:
	// Unmapped
	b->core.flag |= BAM_FUNMAP;
    }

    // TODO: don't replace cigar if alignment fails.
    // Instead keep original, as likely still valid (it's against reference still)
    replace_cigar(b, cig_ind, cig_a);
    
    free(sub);
    return 0;
}

// Scan consensus to look for short tandem repeats.
// Clip cigar based on trailing ends overlapping but not spanning an STR.
//
// Note this is STR specific.
// For example (NB: for real data we'll only clip when het indel is present,
// and ideally only then when the indel is as long or longer than the repeat
// unit).
//
// BEFORE
// Cons: ATGCTAGTGTGTGTGTCTCTCTCTCAGCTAG
// Seq   ATGCTAGTGTGTG
// Seq   ATGCTAGTGTGTGTGTCTCTC
// Seq      CTAGTGTGTGTGTCTCTCTCTCAGCTAG
// Seq          TGTGTGTGTCTCTCTCTCAGCTAG
// Seq                    TCTCTCTCAGCTAG
//
// AFTER
// STR:  ------1111111111---------------
// STR:  ---------------2222222222------
// Cons: ATGCTAGTGTGTGTGTCTCTCTCTCAGCTAG
// Seq   ATGCTAgtgtgtg                    doesn't span 1
// Seq   ATGCTAGTGTGTGTGtctctc            spans 1, but not 2
// Seq      CTAGTGTGTGTGTCTCTCTCTCAGCTAG  spans
// Seq          tgtgtgtgtCTCTCTCTCAGCTAG  spans 2, but not 1
// Seq                    tctctctcAGCTAG  doesn't span 2
int trim_cigar_STR(char *ref, int start, bam1_t **bams, int nbam, int *new_pos) {
    // Compute short tandem repeats.
    uint32_t *str = NULL, i;
    int len = strlen(ref);

    // Find STR and mark nodes as belonging to specific STR numbers.
    rep_ele *reps, *elt, *tmp;
    int str_num = 0;

    reps = find_STR(ref, len, 0);
    str = calloc(len, sizeof(*str));

    DL_FOREACH_SAFE(reps, elt, tmp) {
	// Compute markers for nodes.
	// FIXME: only do this if seq between start/end contains
	// lowercase letters (het ins).  Otherwise STR doesn't matter.
	if (elt->start < len && elt->start > 0 && ref[elt->start] != 'N') {
	    for (i = elt->start; i < elt->end && i < len; i++)
		str[i] |= (1<<str_num);
	    str_num = (str_num+1)&31;

	    //fprintf(stderr, "STR: %2d .. %2d %.*s\n", elt->start, elt->end,
	    //	    elt->end - elt->start+1, &ref[elt->start]);
	}
	DL_DELETE(reps, elt);
	free(elt);
    }

    // left edge;
    char *cig_str = NULL;
    int cig_str_len = 0, cig_str_ind;
    for (i = 0; i < nbam; i++) {
	int left_trim = 0;
	int left_shift = 0;
	bam1_t *b = bams[i];
	uint32_t STR = str[new_pos[i]-start];
	if (new_pos[i] > start)
	    STR |= str[new_pos[i]-start-1]; // incase we start in an insertion
	int sp, rp; // seq & ref pos
	//fprintf(stderr, "Name: %s\t", bam_get_qname(b));
	uint8_t *seq = bam_get_seq(b);
	int op_len = 0, op, cig_ind = 0;
	uint32_t *cig = bam_get_cigar(b);
	int adjacent_STR = 0;

	if (cig_str_len < b->core.l_qseq + len) {
	    cig_str_len = b->core.l_qseq + len;
	    cig_str = realloc(cig_str, cig_str_len);
	}
	cig_str_ind = 0;
	//for (sp = 0, rp = b->core.pos; sp < b->core.l_qseq; ) {
	for (sp = 0, rp = new_pos[i]; sp < b->core.l_qseq; ) {
	    char sbase = "=ACMGRSVTWYHKDBN"[bam_seqi(seq, sp)];
	    char rbase = rp >= start && rp < start+len ? ref[rp-start] : 'N';
	    if (op_len == 0) {
		if (cig_ind < b->core.n_cigar) {
		    op = bam_cigar_op(cig[cig_ind]);
		    op_len = bam_cigar_oplen(cig[cig_ind++]);
		} else {
		    op = BAM_CSOFT_CLIP;
		    op_len = INT_MAX;
		}
	    }
	    cig_str[cig_str_ind] = op;
	    //printf("%c/%d\t%c/%d\t%c\t%08x\t%x\n", sbase, sp, rbase, rp, "MIDNSHP=XB"[op], STR, str[rp-start]);
	    if (STR) {
		// check backwards for pads
		int j = cig_str_ind-1;
		while (j >= 0 && cig_str[j] == BAM_CPAD)
		    j--;
		cig_str_ind = j+1;

		if (bam_cigar_type(op) & 2) //ref
		    left_shift++;
		if (!(bam_cigar_type(op) & 1) && op != BAM_CHARD_CLIP) {
		    // del, skip, pad; just delete them from cigar
		    cig_str_ind--;
		} else if (bam_cigar_type(op) & 1) { //seq
		    cig_str[cig_str_ind] = BAM_CSOFT_CLIP;
		    left_trim++;
		}
		adjacent_STR = 1;
	    } else if (adjacent_STR) {
		if (!(bam_cigar_type(op) & 2)) { //!ref, so pad, ins, etc
		    if (op != BAM_CHARD_CLIP) {
			left_trim++;
			if (bam_cigar_type(op) & 1) // seq
			    cig_str[cig_str_ind] = BAM_CSOFT_CLIP;
			else
			    cig_str_ind--; // trim - eg pad
		    }
		} else {
		    adjacent_STR = 0;
		}
	    }

	    if ((bam_cigar_type(op) & 2) && rp >= start && rp < start+len)
		// and off existing STRs that finish, but don't acquire
		// new ones that we didn't start within.
		STR &= ~(STR ^ str[rp-start]);
	    sp += bam_cigar_type(op) & 1;
	    rp += (bam_cigar_type(op) & 2) ? 1 : 0;
	    op_len--;
	    cig_str_ind++;
	}
	// FIXME: trailing hard-clip at end of cigar; not in seq so exit loop early.

	if (0) {
	    fprintf(stderr, "trim %d, shift %d\t", left_trim, left_shift);
	    int k;
	    for (k = 0; k < cig_str_ind; k++)
		fputc("MIDNSHP=XB"[cig_str[k]], stderr);
	    fputc('\n', stderr);
	}

	if (left_trim || left_shift) {
	    uint32_t *cig2 = malloc(b->core.l_qseq * sizeof(*cig2)), cig2_ind = 0;
	    int k = 0;
	    while (k < cig_str_ind) {
		int j = k;
		while (j < cig_str_ind && cig_str[k] == cig_str[j])
		    j++;
		cig2[cig2_ind++] = bam_cigar_gen(j-k, cig_str[k]);
		k = j;
	    }
	    new_pos[i] += left_shift;
	    replace_cigar(b, cig2_ind, cig2);
	    free(cig2);
	}
    }

    free(cig_str);
    free(str);
}

int int64_compar(const void *vp1, const void *vp2) {
    return ((*(const uint64_t *)vp2)>>32) - ((*(const uint64_t *)vp1)>>32);
}


// Starting from a single node, we recurse down (most) paths making
// sure we cover all nodes that can be reached.  Each route is allocated
// a new haplotype structure.
//
// To avoid combinatorial explosion, we don't try every route and prefer
// to stick to the busiest ones instead when we call this multiple
// times with different head nodes.
hseqs *follow_graph(dgraph_t *g, int x, char *prefix, char *prefix2, int len, hseqs *h, int *nh) {
    int seq_idx = 0, seq_alloc = len + g->kmer*2;
    char *seq = calloc(1, seq_alloc);
    char *seq2 = calloc(1, seq_alloc);

    //fprintf(stderr, "Follow graph n=%d with prefix '%.*s'\n", x, len, prefix);

    memcpy(seq, prefix, len);
    memcpy(seq2, prefix2, len);
    seq_idx += len;

    node_t *n = g->node[x], *norig = n;
    if (!prefix) {
	// Note: currently uses first tagged seq only per node, but this is
	// valid as we should only use follow_graph before bubble collapsing.
	memcpy(seq + len, n->hi[0]->key, g->kmer);
	memcpy(seq2 + len, n->hi[0]->key, g->kmer);
	seq_idx += g->kmer;
    } else {
	seq2[seq_idx] = n->hi[0]->key[g->kmer-1];
	seq[seq_idx++] = vec2X(n->bases[g->kmer-1]);
    }

    while (n->n_out == 1 /*&& n->visited != 'h'*/) {
	n->visited = 'h';
	n = g->node[n->out[0]->n[1]];
	if (seq_idx == seq_alloc) {
	    seq_alloc *= 2;
	    seq = realloc(seq, seq_alloc);
	    seq2 = realloc(seq2, seq_alloc);
	}
	seq2[seq_idx] = n->hi[0]->key[g->kmer-1];
	seq[seq_idx++] = vec2X(n->bases[g->kmer-1]);
    }

    if (n->n_out > 1 /*&& n->visited != 'h'*/) {
	int i;
	n->visited == 'h';

	// Recurse down best route first
	uint64_t *counts = malloc(n->n_out * sizeof(*counts));
	for (i = 0; i < n->n_out; i++)
	    counts[i] = (((uint64_t)n->out[i]->count)<<32) | i;
	qsort(counts, n->n_out, sizeof(*counts), int64_compar);

	// If we have any node unvisited, then follow it.
	int followed = 0;
	for (i = 0; i < n->n_out; i++) {
	    int j = counts[i]&0xffffffff;
	    node_t *x = g->node[n->out[j]->n[1]];
	    if (x->visited != 'h') {
		h = follow_graph(g, n->out[j]->n[1], seq, seq2, seq_idx, h, nh);
		followed = 1;
	    }
	}
	// If we didn't follow any routes because all have been visited,
	// still follow one (the best one) so our haplotypes are full length.
	if (!followed)
	    h = follow_graph(g, n->out[counts[0]&0xffffffff]->n[1], seq, seq2, seq_idx, h, nh);
	free(counts);
	free(seq);
	free(seq2);
    } else {
	//fprintf(stderr, "haplotype: %d..%d %.*s\n", norig->id, n->id, seq_idx, seq);
	*nh++;
	hseqs *h2 = calloc(1, sizeof(*h2));
	h2->next = h;
	h2->seq = seq;
	h2->seq2 = seq2;
	h2->len = seq_idx;
	h = h2;
    }

    return h;
}

#ifndef ERRK
#define ERRK 30
#endif

// FIXME: min_count needs to be depth based.  Find mean count and
// use this to cap min_count?  So low coverage would reduce,
// min_count, but high coverage or lots of low complexity data won't
// increase it.
HashTable *kmer_hash = NULL, *neighbours = NULL;
int correct_errors(haps_t *h, int n, int errk, int min_count, int min_qual) {
    //HashTable *hash, *neighbours;
    HashItem *hi;
    int i, counth = 0, countw = 0;

    HashTable *hash = kmer_hash;
    if (hash) HashTableDestroy(hash, 0);
    if (neighbours) HashTableDestroy(neighbours, 0);
    
    hash = HashTableCreate(8, HASH_DYNAMIC_SIZE | HASH_POOL_ITEMS);
    kmer_hash = hash;

    // Hash words
    for (i = 0; i < n; i++) {
	char *seq = h[i].seq;
	int len = strlen(seq), j;
	for (j = 0; j < len-errk; j++) {
	    HashData hd;
	    HashItem *hi;
	    int nw, k;

	    hd.i = 0;
	    hi = HashTableAdd(hash, seq+j, errk, hd, &nw);
	    hi->data.i++;
	    counth++;
	    countw+=nw;
	}
    }

    fprintf(stderr, "%d unique words, %d total words\n", countw, counth);

    // Find common words and produce neighbourhoods
    int avg = 2*ceil((double)countw / counth); // FIXME: improve this
    neighbours = HashTableCreate(8, HASH_DYNAMIC_SIZE | HASH_POOL_ITEMS);

    HashIter *hiter = HashTableIterCreate();
    while ((hi = HashTableIterNext(hash, hiter))) {
	//printf("%.*s %d\n", errk, hi->key, hi->data.i);
	if (hi->data.i >= avg) {
	    int j;
	    char *s = hi->key;
	    char N[errk];
	    memcpy(N, hi->key, errk);
	    for (j = 0; j < errk; j++) {
		int nw, k;
		HashData hd;
		hd.p = hi->key;
		int base = N[j];
		for (k = 0; k < 4; k++) {
		    HashItem *hi2;
		    if ("ACGT"[k] == base) continue;
		    N[j] = "ACGT"[k];
		    hi2 = HashTableAdd(neighbours, N, errk, hd, &nw);
		    if (!nw) hi2->data.p = NULL; // mark the clash
		}
		N[j] = base;
	    }
	}
    }

//    // Debug: find rare words and report if they have a common neighbour
//    HashTableIterReset(hiter);
//    while ((hi = HashTableIterNext(hash, hiter))) {
//	//if (hi->data.i >= min_count) continue;
//	HashItem *hi_n = HashTableSearch(neighbours, hi->key,  errk);
//	printf("%.*s %d -> %.*s\n", errk, hi->key, (int)hi->data.i,
//	       errk, hi_n ? (char *)hi_n->data.p : 0);
//	if (!hi_n) continue;
//    }

    HashTableIterDestroy(hiter);


    // Auto-correct rare words
    int nc = 0;
    for (i = 0; i < n; i++) {
	char *seq = h[i].seq;
	char *qual = h[i].qual;
	char *s2 = strdup(seq);
	int len = strlen(seq), j;
#define EDGE_DIST 3
	for (j = EDGE_DIST; j < len-errk-EDGE_DIST; j++) {
	    HashItem *hi, *hi2;

	    hi = HashTableSearch(hash, seq+j, errk);
	    if (!hi || hi->data.i >= min_count)
		continue;

	    hi2 = HashTableSearch(neighbours, hi->key, errk);
	    if (!hi2 || !hi2->data.p) {
		//fprintf(stderr, "No correction for %.*s %d\n", errk, hi->key, hi->data.i);
		continue;
	    }

	    //memcpy(s2+j, hi2->data.p, errk);

	    int k;
	    for (k = 0; k < errk; k++)
		if (s2[j+k] != ((char *)hi2->data.p)[k])
		    break;
	    if (k == errk)
		continue;

	    if (qual && min_qual && qual[k] >= min_qual)
		continue;

	    s2[j+k] = ((char *)hi2->data.p)[k];
	    nc++;
	    //fprintf(stderr, "Correct %.*s %d -> %.*s\n", errk, hi->key, hi->data.i, errk, hi2->data.p);
	}
	free(seq);
	h[i].seq = s2;
    }
    fprintf(stderr, "Corrected %d (%5.2f%% %5.2f%%)\n", nc, 100.0*nc/countw, 100.0*nc/counth);


//    HashTableDestroy(hash, 0);
//    HashTableDestroy(neighbours, 0);

    return nc;
}

static unsigned char complementary_base[256] = {
    0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
    16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
    32, '!', '"', '#', '$', '%', '&', '\'','(', ')', '*', '+', ',', '-', '.', '/',
    '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', ':', ';', '<', '=', '>', '?',
    '@', 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
    'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z', '[', '\\',']', '^', '_',
    '`', 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
    'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', '{', '|', '}', '~', 127,
    128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
    144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
    160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
    176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
    192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
    208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
    224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
    240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255,
};

/* Reverse complements a sequence in-place */
void complement_seq ( char *seq, int seq_len ) {
    int i, middle, j;
    unsigned char temp;

    middle = seq_len/2;
    for ( i = 0, j = seq_len-1; i < middle; i++, j--) {
        temp = (unsigned char) seq[i];
        seq[i] = complementary_base [ (unsigned char) seq[j] ];
        seq[j] = complementary_base [ temp ];
    }

    if ( seq_len % 2 )
	seq[middle] = complementary_base [ (unsigned char) seq[middle] ];
}


// Need ADAPTER
#define ADAPTER_KMER 14
#define ADAPTER_COUNT 10

int trim_adapters(haps_t *h, int n, char *fn, int kmer, int min_qual) {
    int i,j;
    HashTable *hash;
    hash = HashTableCreate(1024, HASH_DYNAMIC_SIZE | HASH_POOL_ITEMS);

    // expected fasta of short seqs
    FILE *fp = fopen(fn, "r");
    if (!fp) {
	perror(fn);
	return -1;
    }

    // Hash adapter file
    char line[1024], line2[1024];
    while (fgets(line, 1024, fp)) {
	if (*line == '>')
	    continue;
	line[1023]=0;
	int len = strlen(line)-1;
	line[len] = 0; // trim \n
	memcpy(line2, line, len);
	int dir;
	for (dir = 0; dir <= 1; dir++) {
	    for (i = 0; i < len-kmer; i++) {
		HashData hd;
		hd.i = dir;
		HashTableAdd(hash, line+i, kmer, hd, NULL);

		// Neighbours
		for (j = 0; j < kmer; j++) {
		    line2[i+j] = 'A'; HashTableAdd(hash, line2+i, kmer, hd, NULL);
		    line2[i+j] = 'C'; HashTableAdd(hash, line2+i, kmer, hd, NULL);
		    line2[i+j] = 'G'; HashTableAdd(hash, line2+i, kmer, hd, NULL);
		    line2[i+j] = 'T'; HashTableAdd(hash, line2+i, kmer, hd, NULL);
		    line2[i+j] = line[i+j];
		}
	    }

	    complement_seq(line, len);
	}
    }
    fclose(fp);

    // Trim seqs
    for (i = 0; i < n; i++) {
	char *seq = h[i].seq;
	int len = strlen(seq);

	int left_end = INT_MAX;
	int right_end = INT_MIN;
	int left_count = 0;
	int right_count = 0;

	for (j = 0; j < len-kmer; j++) {
	    HashItem *hi;
	    if ((hi = HashTableSearch(hash, seq+j, kmer))) {
		switch (hi->data.i == 0) {
		case 0: // left
		    left_end = j+kmer-1;
		    left_count++;
		    break;

		case 1: // right
		    if (right_end == INT_MIN)
			right_end = j;
		    right_count++;
		    break;
		}
		    
		//printf("%s: %d/%d %.*s %d\n", h[i].name, j, len-kmer, kmer, seq+j, (int)hi->data.i);
	    }
	}

	if (left_count >= ADAPTER_COUNT) {
	    fprintf(stderr, "%s: trim left at %d: %.*s\n", h[i].name, left_end, left_end, h[i].seq);
	    memset(h[i].seq, 'N', left_end);
	}
	if (right_count >= ADAPTER_COUNT) {
	    fprintf(stderr, "%s: trim right at %d: %.*s\n", h[i].name, right_end, len-right_end, h[i].seq+right_end);
	    memset(h[i].seq+right_end, 'N', len-right_end);
	}

	for (j = 0; j < len; j++) {
	    if (h[i].qual[j] < min_qual)
		seq[j] = 'N';
	    //else if (j == len-1 || h[i].qual[j+1] >= min_qual)
	    else
		break;
	}

	for (j = len-1; j >= 0; j--) {
	    if (h[i].qual[j] < min_qual)
		seq[j] = 'N';
	    //else if (j == 0 || h[i].qual[j-1] >= min_qual)
	    else
		break;
	}
    }

    HashTableDestroy(hash, 0);
    fflush(stdout);

    return 0;
}

// // Also makes them writeable, but leaks the memory later on if we don't free.
// void fix_seqs(haps_t *h, int n) {
//     int i;
//     for (i = 0; i < n; i++) {
// 	char *s = strdup(h[i].seq);
// 	int j, k, l = strlen(s);
// 	for (j = k = 0; j < l; j++) {
// 	    if (s[j] != '*')
// 		s[k++] = toupper(s[j]);
// 	}
// 	s[k] = 0;
// 	h[i].seq = s; // leads to memleak
//     }
// }

// Computes a consensus via a greedy approach.
// Find the node with the highest count.
// Backtrack up from here and recurse down from here following
// the highest scoring path in all cases.
//
// Add kmer-1 Ns to front. This helps in the CIGAR generation later.
haps_t *compute_consensus(dgraph_t *g) {
    int i, j;

    // Update edge 'in' counts
    for (i = 0; i < g->nnodes; i++) {
	node_t *n = g->node[i];
	if (n->pruned)
	    continue;

	for (j = 0; j < n->n_out; j++) {
	    int k;
	    node_t *n2 = g->node[n->out[j]->n[1]];
	    for (k = 0; k < n2->n_in; k++) {
		if (n2->in[k]->n[1] == n->id) {
		    n2->in[k]->count = n->out[j]->count;
		    break;
		}
	    }
	}
    }

    // Find the best scoring node.
    int best_n = -1, best_score = 0;
    for (i = 0; i < g->nnodes; i++) {
	node_t *n = g->node[i];
	if (n->pruned)
	    continue;

	int cnt = 0;
	for (j = 0; j < n->n_in; j++)
	    cnt += n->in[j]->count;
	for (j = 0; j < n->n_out; j++)
	    cnt += n->out[j]->count;

	if (best_score < cnt) {
	    best_score = cnt;
	    best_n = i;
	}
    }
    //printf("Best node = %d (cnt %d)\n", best_n, best_score);

    if (best_n < 0)
	return NULL;

    // Scan backwards from here
    kstring_t seq = {0};

    node_t *n = g->node[best_n];
    while (n) {
	int best_i = -1;
	best_score = 0;
	for (i = 0; i < n->n_in; i++) {
	    if (best_score < n->in[i]->count) {
		best_score = n->in[i]->count;
		best_i = i;
	    }
	}
	if (best_i >= 0) {
	    //printf("Best prev = %d\n", n->in[best_i]->n[1]);
	    n = g->node[n->in[best_i]->n[1]];
	    if (islower(vec2X(n->bases[g->kmer-1])))
		kputc(tolower(n->hi[0]->key[g->kmer-1]),&seq);
	    else
		kputc(n->hi[0]->key[g->kmer-1],&seq);
	} else {
	    for (j = g->kmer-1; j > 0; j--)
		kputc(n->hi[0]->key[j-1], &seq);
	    n = NULL;
	}
    }
    for (i = 0; i < g->kmer-1; i++)
	kputc('N', &seq);
    for (i = 0, j = seq.l-1; i < j; i++, j--) {
	char tmp = seq.s[i];
	seq.s[i] = seq.s[j];
	seq.s[j] = tmp;
    }

    // Scan downwards
    n = g->node[best_n];
    while (n) {
	if (islower(vec2X(n->bases[g->kmer-1])))
	    kputc(tolower(n->hi[0]->key[g->kmer-1]), &seq);
	else
	    kputc(n->hi[0]->key[g->kmer-1], &seq);

	int best_i = -1;
	best_score = 0;
	for (i = 0; i < n->n_out; i++) {
	    if (best_score < n->out[i]->count) {
		best_score = n->out[i]->count;
		best_i = i;
	    }
	}
	if (best_i >= 0) {
	    //printf("Best next = %d\n", n->out[best_i]->n[1]);
	    n = g->node[n->out[best_i]->n[1]];
	} else {
	    n = NULL;
	}
    }

    haps_t *h = malloc(sizeof(*h));
    if (!h)
	return NULL;

    h->name = "cons";
    h->seq = seq.s;
    h->pos = 0;
    h->qual = NULL;

    //printf("Cons = %s\n", seq.s);

    return h;
}

haps_t *bam2haps(bam1_t **bams, int nrecs) {
    int i;
    haps_t *haps = malloc(sizeof(*haps) * nrecs);
    if (!haps)
	return NULL;

    for (i = 0; i < nrecs; i++) {
	bam1_t *b = bams[i];
	uint8_t *seq = bam_get_seq(b);
	int j;

	haps[i].pos  = 0;
	haps[i].name = bam_get_qname(b);
	haps[i].seq  = malloc(b->core.l_qseq+1);
	if (!haps[i].seq)
	    return NULL; // leak on failure
	haps[i].qual = bam_get_qual(b);

	for (j = 0; j < b->core.l_qseq; j++)
	    haps[i].seq[j] = "=ACMGRSVTWYHKDBN"[bam_seqi(seq, j)];
	haps[i].seq[b->core.l_qseq] = 0;
    }

    return haps;
}

void free_haps(haps_t *h, int n) {
    if (!h)
	return;

    int i;

    for (i = 0; i < n; i++)
	free(h[i].seq);
    free(h);
}

static void dump_input(bam_hdr_t *hdr, bam1_t **bams, int nbams, char *ref, int ref_len, int ref_start) {
#if 1
    samFile *fp;
    int i;

    if (!(fp = sam_open("_tmp.sam", "w")))
	abort();
    if (sam_hdr_write(fp, hdr) < 0)
	abort();
    for (i = 0; i < nbams; i++) {
	if (sam_write1(fp, hdr, bams[i]) < 0)
	    abort();
    }
    if (sam_close(fp) < 0)
	abort();

    FILE *f = fopen("_tmp.ref", "w");
    fprintf(f, ">%d\n%.*s\n", ref_start, ref_len, ref);
    fclose(f);
#endif
}

int bam_realign(bam_hdr_t *hdr, bam1_t **bams, int nbams, int *new_pos,
		char *ref_seq, int ref_len, int ref_start) {
    int i, kmer = KMER, ret = -1;
    dgraph_t *g;
    haps_t *haps = NULL, *cons = NULL;
    haps_t *ref = NULL, *ref_ = NULL;

    dump_input(hdr, bams, nbams, ref_seq, ref_len, ref_start);

    init_W128_score(-4,1);
    init_X128_score(-4,4);

    // Convert BAMS to seqs instead, so we can perform edits on them
    // while retaining the original data.
    int nhaps = nbams;
    haps = bam2haps(bams, nhaps);

    // Successive rounds allows for fixing more than 1 error in a read.
    fprintf(stderr, "Correcting\n");
    // 6- poor (3)
    // 7..10 good (3)
    // 11+ poor (3)
//    for (i = 0; i < 3; i++)
//	correct_errors(haps, nhaps, ERRK, 3);

    correct_errors(haps, nhaps, 27, 3, 0);
    correct_errors(haps, nhaps, 25, 1, 0);
    correct_errors(haps, nhaps, 20, 1, 0);
    correct_errors(haps, nhaps, 14, 1, 0);
//    for (i = 0; i < 10; i++)
//	if (correct_errors(haps, nhaps, 14, 2, 5) <= 0)
//	    break;
    //trim_adapters(haps, nhaps, "/nfs/srpipe_references/adapters/adapters.fasta", 14);
    // Also trim low quality terminal bases
    trim_adapters(haps, nhaps, "adapter.fa", ADAPTER_KMER, 10);

 bigger_kmer:

    // If our mapped sequences can be processed with a sensible kmer,
    // but the unmapped cannot, then maybe it's just the unmapped is
    // a broken homo-polymer.  In which case keep it unmapped and ignore
    // it.
    for (; kmer < MAX_KMER; kmer += 10) {
	fprintf(stderr, "Building graph with kmer=%d\n", kmer);
	g = graph_create(kmer);
	for (i = 0; i < nhaps; i++) {
	    if (add_seq(g, haps[i].seq, 0, 0) != 0) {
		// loop within a sequence.  If unmapped, we'll just skip using it
		if (bams[i]->core.flag & BAM_FUNMAP) {
		    del_seq(g, haps[i].seq, 0, 0);
		    continue;
		}
		fprintf(stderr, "Loop detected within mapped seq %s\n", haps[i].seq);
		break;
	    }
	}

	if (i == nhaps && loop_check(g, 0) == 0)
	    break;

	// FIXME: we get loops too often.  How and why?
	fprintf(stderr, "Loop detected, increasing kmer\n");

	graph_destroy(g);
	g = NULL;
    }
    if (kmer >= MAX_KMER) {
	fprintf(stderr, "No suitable kmer found\n");
	goto err;
    }
    fprintf(stderr, "Using kmer %d\n", kmer);

    //graph2dot(g, "f.dot", 0);

    // Merge bubbles excluding reference
    find_bubbles(g, 0);
    assert(loop_check(g, 0) == 0);
    //graph2dot(g, "_F.dot", 0);

    // Also try this after find_insertions() call
#if 0
    cons = compute_consensus(g); 
    if (add_seq(g, cons->seq, 0, 0) < 0 || loop_check(g, 0)) {
	fprintf(stderr, "Loop when adding consensus\n");
	graph_destroy(g);
	if ((kmer += 10) < MAX_KMER)
	    goto bigger_kmer;
	g = NULL;
	goto err;
    }

    // Merge in head & tail tips
    int merged;
    do {
	merged = 0;
	//puts("Merging\n");
	merged += merge_head_tips(g, cons->seq, strlen(cons->seq));
	graph2dot(g, "H.dot", 0);

	merged += merge_tail_tips(g);
	graph2dot(g, "T.dot", 0);
    } while (merged);
#endif

    int shift = bams[0]->core.pos+1;

    if (ref_seq) {
	ref_ = ref = malloc(sizeof(*ref));
	if (!ref)
	    goto err;
	fprintf(stderr, "Alloc %p\n", ref);
	ref->name = "Ref";
	ref->pos = ref_start;
	fprintf(stderr, "ref=%.*s\n", ref_len, ref_seq);

	ref->seq = malloc(ref_len + g->kmer-1 + 1);
	ref->seq[ref_len + g->kmer-1] = 0;
	memcpy(ref->seq + g->kmer-1, ref_seq, ref_len);
	memset(ref->seq, 'N', g->kmer-1);
	if (add_seq(g, ref->seq, 0, IS_REF) < 0 || loop_check(g, 0)) {
	    fprintf(stderr, "Loop when adding reference\n");
	    graph_destroy(g);
	    if ((kmer += 10) < MAX_KMER) {
		free_haps(ref_, 1);
		goto bigger_kmer;
	    }
	    g = NULL; goto err;
	}

	shift = ref_start + 1;
    }

    //graph2dot(g, "_g.dot", 0);

    if (!(cons = compute_consensus(g)))
	goto err;

    fprintf(stderr, "cons=%s\n", cons->seq);
    if (add_seq(g, cons->seq, 0, (ref?0:IS_REF)|IS_CON) < 0 || loop_check(g, 0)) {
	fprintf(stderr, "Loop when adding consensus\n");
	graph_destroy(g);
	if ((kmer += 10) < MAX_KMER)
	    goto bigger_kmer;
	g = NULL;
	goto err;
    }
    if (!ref)
	ref = cons;

    // Now also merge in reference bubbles
    find_bubbles(g, 1);
    assert(loop_check(g, 0) == 0);

    graph2dot(g, "_G.dot", 0);


    // Strings of bases inserted between reference coords
    // permit us to generate alignments starting or ending
    // inside the insertions.
    find_insertions(g);
    //graph2dot(g, "_I.dot", 0);

    //fprintf(stderr, "Pruning\n");
    //prune(g, argc > 3 ? atoi(argv[3]) : 2);

    for (i = 0; i < nhaps; i++) {
        if (seq2cigar_new(g, ref->seq, shift, bams[i], haps[i].seq, &new_pos[i]) < 0)
	    goto err;
    }

    trim_cigar_STR(ref->seq+g->kmer, shift, bams, nhaps, new_pos);

    ret = 0;
 err:
    graph_destroy(g);
    free_haps(haps, nhaps);
    free_haps(cons, 1);
    fprintf(stderr, "Free %p\n", ref_);
    free_haps(ref_, 1);

    fprintf(stderr, "Finished.\n");

    return ret;
}

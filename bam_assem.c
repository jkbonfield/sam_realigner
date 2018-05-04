// Copy of assem_bam3.c, but for use within bam_problem_regions instead of
// as a standalone tool.


// cc -g -I$HOME/work/samtools_master/htslib assem_bam.c hash_table.c pooled_alloc.c string_alloc.c -DKMER=70 -lm -L$HOME/work/samtools_master/htslib -lhts -lz -pthread

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
 *
 * - Track nodes on cons/ref path, but not in cons/ref.  Ie if we go from
 *   ref through non-ref node X,Y,Z and back to ref, then X,Y,Z are on-ref.
 *   This will greatly simplify alignment of cons vs ref later on.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <limits.h>
#include <ctype.h>

#include "hash_table.h"
#include "string_alloc.h"
#include "str_finder.h"
#include "bam_assem.h"
#include "ksw.h"

#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"

#define IS_REF 1
#define IS_CON 2

#define GOPEN 5
#define GEXT  1

#define GOPEN_REF 5
#define GEXT_REF  1

// Use most frequent 95% of words for error correcting.
#define CORRECT_PERC 0.1

// // But only if their combined confidence is high enough
// #define CORRECT_MIN_CONF 15

// And only when the confidence of the new kmer is MULT times
// larger than the confidence of the old kmer.
//#define CORRECT_MULT 1

// Appropriate when summing kmer-qual over all matching kmers.
// #define CORRECT_MIN_CONF 30
// #define CORRECT_MULT 5

//---------------------------------------------------------------------------

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

typedef unsigned char uc;

// Ambiguity codes include base-pad combinations.
// We set the scores so that given the choice of aligning
// A vs AC ambig or A* ambig, the AC scores higher.  Thus by
// elimination pads preferentially align to A* than AC.
int8_t X128[128][128];
int8_t X128_ref[128][128];
void init_X128_score(int8_t X128[128][128], int mis_base, int mis, int mat) {
    init_base_val();

    // This works well as 0 for seq vs seq alignments,
    // but may be best as eg -3 or -4 for cons vs ref in order
    // to get neighbouring snps labelled as compound indels.
    memset(X128, mis_base, 128*128*sizeof(X128[0][0]));
    int i, j;

    // Matches
    for (i = 0; i < 4; i++)
	X128[(uc)"ACGT"[i]][(uc)"ACGT"[i]] = mat;

    // Ambiguity codes; partial match
    for (i = 0; i < 3; i++)
	X128['A'][(uc)"MRW"[i]] = X128[(uc)"MRW"[i]]['A'] = mat/2;
    for (i = 0; i < 3; i++)
	X128['C'][(uc)"MSY"[i]] = X128[(uc)"MSY"[i]]['C'] = mat/2;
    for (i = 0; i < 3; i++)
	X128['G'][(uc)"RSK"[i]] = X128[(uc)"RSK"[i]]['G'] = mat/2;
    for (i = 0; i < 3; i++)
	X128['T'][(uc)"WYK"[i]] = X128[(uc)"WYK"[i]]['T'] = mat/2;
    // Ambiguity codes; mismatch (B = not A; D = not C; etc)
    X128['A']['B'] = X128['B']['A'] = mis/2;
    X128['C']['D'] = X128['D']['C'] = mis/2;
    X128['G']['H'] = X128['H']['G'] = mis/2;
    X128['T']['V'] = X128['V']['T'] = mis/2;

    // Ambiguity vs pads.
    // All to start with, to set mismatch (A vs G*) baseline.
    for (i = 0; i < 11; i++) {
	for (j = 0; j < 11; j++) {
	    X128[(uc)"AMRWCSYGKTN"[i]][tolower("AMRWCSYGKTN"[j])] = mis;
	    X128[tolower("AMRWCSYGKTN"[i])][(uc)"AMRWCSYGKTN"[j]] = mis;
	    // A* vs C* is a partial match (pad wise). How to score?
	    X128[tolower("AMRWCSYGKTN"[i])][tolower("AMRWCSYGKTN"[j])] = mis/2;
	}
    }
    // Now fix up the partial matches (A vs A*).
    for (i = 0; i < 4; i++) {
	X128[(uc)"ACGT"[i]][tolower("ACGT"[i])] = mat/2-1;
	X128[tolower("ACGT"[i])][(uc)"ACGT"[i]] = mat/2-1;
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
	X128[(uc)"AMRWCSYGKTN"[i]]['-'] = X128['-'][tolower("AMRWCSYGKTN"[i])] = mis*2;
	X128[tolower("AMRWCSYGKTN"[i])]['-'] = X128['-'][(uc)"AMRWCSYGKTN"[i]] = mis*2;
    }
}

// Convert a 6-way ACGT*N vector to X128 lookup value.
int vec2X(int V[6]) {
#if 0 // running as 1
    char b = "NTGKCYSBAWRDMHVN"[(!!V[0]<<3)+(!!V[1]<<2)+(!!V[2]<<1)+(!!V[3]<<0)];
    if (V[4]>0) b = tolower(b);
    return b;
#else
    // ambig if at least 1/4tr.
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

//---------------------------------------------------------------------------
// Interface to Heng Li's ksw code
typedef struct haps {
    int   pos;
    char *seq;
    char *name;
    uint8_t *qual; // not our copy; do not free
} haps_t;

//#include "hap2s.h"
//#include "haps.h"

#ifndef KMER
#  define KMER 40
#endif

#ifndef KMER_INC_FINAL
#  define KMER_INC_FINAL 10
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

#ifndef MIN_STR
#  define MIN_STR 0
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

    int above, below;
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

    //    node_hash_t *_node_hash;  // "seq" -> node
    HashTable *node_hash;  // "seq" -> node
    HashTable *edge_hash;  // id[2] -> edge

    string_alloc_t *spool;
} dgraph_t;

void validate_graph(dgraph_t *g);

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

    //kh_destroy(node_hash, g->_node_hash);

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
    //g->_node_hash = kh_init(node_hash);
    //g->node_hash = HashTableCreate(8, HASH_DYNAMIC_SIZE | HASH_POOL_ITEMS);
    g->node_hash = HashTableCreate(8, HASH_DYNAMIC_SIZE | HASH_POOL_ITEMS | HASH_NONVOLATILE_KEYS | HASH_FUNC_TCL);
    if (!g->node_hash) {
	graph_destroy(g);
	return NULL;
    }

    g->edge_hash = HashTableCreate(8, HASH_DYNAMIC_SIZE | HASH_POOL_ITEMS | HASH_FUNC_TCL);
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
	    node->bases[j][(uc)base_val[seq[j] & 0x7f]]++;
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
    if (n1->n_in >= MAX_EDGE)
	return NULL;
    n2->in[n2->n_in++] = edge;

    return edge;
}

// Move edge e with new source (n2->n3) to n1 (n1->n3)
void move_edge_in(dgraph_t *g, edge_t *e, node_t *n1, node_t *n2) {
    HashItem *hi = HashTableSearch(g->edge_hash, (char *)e->n, 2*sizeof(*e->n));
    assert(hi);
    HashTableDel(g->edge_hash, hi, 0);

    node_t *n3 = g->node[e->n[1]];
    int j;
    if (n1->n_out >= MAX_EDGE)
	return;
    n1->out[n1->n_out++] = e;
    for (j = 0; j < n3->n_in; j++)
	if (n3->in[j]->n[0] == n2->id)
	    n3->in[j]->n[0] = n1->id;

    // Correct edge_hash too.
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
    HashTableDel(g->edge_hash, hi, 0);

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

    if (1) {
	int i;
	for (i = 0; i < len1; i++)
	    n1->bases[i][(uc)base_val[seq1[i] & 0x7f]]++;
    }

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
	    if (n2->in[i]->n[0] != n1->id)
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
		memmove(&last->out[last_edge], &last->out[last_edge+1],
			(last->n_out - last_edge) * sizeof(*last->out));
		last->n_out--;
		int z,k;
		for (z = k = 0; z < to->n_in; z++) {
		    if (to->in[z]->n[0] != last->id)
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

    // pad strip
    for (i = j = 0; i < len; i++)
	if (seq[i] != '*')
	    seq[j++] = seq[i];
    len = j;

    if (!(ref & IS_REF)) {
	// Prune leading/trailing Ns.
	while (*seq == 'N' && len > 0)
	    seq++, len--;
	while (seq[len-1] == 'N' && len > 0)
	    len--;

	if (len == 0)
	    return 0;
    } else {
	// Ref: prune to no more than kmer-10 Ns
	char *s2 = seq;
	int l2 = len;
	while (*s2 == 'N' && l2 > g->kmer)
	    s2++, l2--;
//	if (l2 == 0)
//	    return 0;
//	if (len-l2 > g->kmer-10) {
//	    seq = s2 - (g->kmer-10);
//	    len = l2 - (g->kmer-10);
//	}
    }

    counter++;

    char *s = string_alloc(g->spool, len);
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
	    return -1; // loop
	}

	g->node[e->n[1]]->visited = counter;

	if (ref) {
	    node_t *n1 = g->node[e->n[0]];
	    node_t *n2 = g->node[e->n[1]];
	    // FIXME: detect possible loops here. Ref must be a clean pass through!
	    //n1->pos = i+g->kmer-1;
//	    if (ref & IS_REF)
//		n1->pos = i;
	    n1->ref |= ref;

	    //n2->pos = i+g->kmer;
//	    if (ref & IS_REF)
//		n2->pos = i+1;
	    n2->ref |= ref;
	}
    }

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
    return best_count != INT_MIN ? g->node[n->in[best_i]->n[0]] : NULL;
}

static void node_common_ancestor_match(dgraph_t *g, node_t *n1, node_t *n2, node_t *n_end,
				       int start, int (*vn)[5], node_t **path1, int np1, node_t **path2, int np2,
				       int *p, int *x1, int *x2, node_t **l1, node_t **l2) {
    int i, j;
    // March along both paths, so merge n2 into n1.
    //printf("M %2d %2d\n", n1->id, n2->id);

    // Merge coordinates
    if (n1->pos == INT_MIN) {
	n1->pos = n2->pos;
	n1->posa = n2->posa;
    }

    // Base frequencies
    //memcpy(n1->bases, vn[(*p)++], g->kmer*5*sizeof(int));
    for (j = 0; j < g->kmer; j++)
	memcpy(n1->bases[j], vn[*p+(g->kmer-j)-1], 5*sizeof(int));
    (*p)++;
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
	    if (n1->in[j] == n2->in[j] ||
		(*x2+1 < np2 && n2->in[i]->n[0] == path2[*x2+1]->id))
		break;
	if (j == n1->n_in && n2->in[i]->n[0] != start) { // not already a parent
	    if (n1->n_in >= MAX_EDGE)
		continue;
	    n1->in[n1->n_in++] = n2->in[i];
	    node_t *n3 = g->node[n2->in[i]->n[0]];
	    move_edge_out(g, n3, n2, n1);
	    //for (j = 0; j < n3->n_out; j++)
	    //    if (n3->out[j]->n[1] == n2->id)
	    //	n3->out[j]->n[1] = n1->id;
	} else if (j != n1->n_in && n2->in[i]->n[0] != start &&
		   ((*x2+1 >= np2) ||
		    (*x2+1 < np2 && path2[*x2+1]->id != n2->in[i]->n[0]))) {
	    // shared node, parent(n1) is also parent(n2) and
	    // not on path, so ensure we delete the edge linking
	    //to n2 // although we don't need to add it to n1 as
	    // it's already there.
	    node_t *n3 = g->node[n2->in[i]->n[0]];
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
	    if (n1->out[j] == n2->out[i] ||
		(*x2 > 0 && path2[*x2-1] && n2->out[i]->n[1] == path2[*x2-1]->id))
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
	    if (*x2 > 0 && path2[*x2-1] && n1->out[j]->n[1] == path2[*x2-1]->id)
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
	    if (t->in[k]->n[0] != n2->id)
		t->in[l++] = t->in[k];
	}
	t->n_in = l;
    }
    (*x1)++, (*x2)++;
    *l1 = n1;
    *l2 = n2;
}

static void node_common_ancestor_ins(dgraph_t *g, node_t *n1, node_t *n2, node_t *n_end, node_t *l1, node_t *l2,
				    int start, int (*vn)[5], node_t **path1, int np1, node_t **path2, int np2,
				    int op, int *p, int *x1, int *x2,
				    node_t **l1_p, node_t **l2_p, node_t **n2_p) {
    int i,j;

    // Insertion in path2, link into path1
    if (!n1)
	n1 = g->node[start]; // if inserting to end of alignment
    int first_ins = 1;
    while (op-- && *x2 < np2) {
	//printf("I  - %2d\n", n2->id);

	if (first_ins) {
	    // Replace l1 in from n1 with in from n2.
	    for (j = 0; j < l1->n_in; j++) {
		if (l1->in[j]->n[0] == n1->id) {
		    // FIXME: plus edge hash
		    l1->in[j]->n[0] = n2->id;
		    // FIXME: do we need to ensure this is in[0]?
		}
	    }

	    // Link n2 out to l1
	    for (j = 0; j < n2->n_out; j++)
		if (n2->out[j]->n[1] == l1->id)
		    break;
	    if (j == n2->n_out)
		// Move out n2->l2 to n2->l1
		move_edge_out(g, n2, l2, l1);

	    memcpy(&n2->bases[g->kmer-1], vn[*p], 5*sizeof(int));
	    first_ins = 0;
	}

	if (op == 0 || *x2 == np2-1) { // last
	    // Link n2 in to n1
	    // Link n1 out to n2
	    // Cull parent(n2) out (to n2)
	    node_t *p2 = g->node[n2->in[0]->n[0]];
	    n2->in[0]->n[0] = n1->id;
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
	}

	(*p)++;

	//n2->bases[g->kmer-1][4]++;

	path2[*x2] = NULL;
	l2 = n2;
	n2 = path2[++(*x2)];
    }

    *l2_p = l2;
    *l1_p = l2;
    *n2_p = n1;
}

// Recursively scan down from curr looking for visited to other_path.
// If found, prune the link to avoid this path from being reached.
int prune_extra_recurse(dgraph_t *g, node_t *last, node_t *curr, node_t *end, int other_path,
			int *vis, int *nvis) {
 tail_loop:
    if (curr->id == end->id)
	return 0;

    if (curr->visited & (1<<30))
	return 0;

    if ((curr->visited & ~(1<<30)) == (other_path & ~(1<<30))) {
	fprintf(stderr, "Found internal bubble. Prune link %d->%d\n", last->id, curr->id);
	int i, j;
	for (i = j = 0; i < last->n_out; i++) {
	    if (last->out[i]->n[1] != curr->id)
		last->out[j++] = last->out[i];
	}
	last->n_out = j;

	for (i = j = 0; i < curr->n_in; i++) {
	    if (curr->in[i]->n[0] != last->id)
		curr->in[j++] = curr->in[i];
	}
	curr->n_in = j;
	return 1;
    }

    vis[(*nvis)++] = curr->id;
    curr->visited |= (1<<30);

    if (curr->n_out == 1) {
	node_t *next = g->node[curr->out[0]->n[1]];
	last = curr;
	curr = next;
	goto tail_loop;
    }

    int i, r = 0;
    for (i = 0; i < curr->n_out; i++)
	r |= prune_extra_recurse(g, curr, g->node[curr->out[i]->n[1]], end, other_path, vis, nvis);

    return r;
}

// Hunt through path looking for branch points that stray from the path.
// When found, check if these could lead to other other path forming
// an second (inner) bubble.
int prune_extra_bubbles(dgraph_t *g, node_t **path, int np, node_t *n_end, int other_path,
			 int *vis, int *nvis) {
    int i, j, r = 0;
    int p_id = path[0]->visited;

    for (i = np-1; i >= 0; i--) {
	node_t *n = path[i];
	if (n->n_out == 1)
	    continue;
	for (j = 0; j < n->n_out; j++) {
	    node_t *n2 = g->node[n->out[j]->n[1]];
	    if ((n2->visited & ~(1<<30)) == (p_id & ~(1<<30)) || n2->id == n_end->id)
		continue;
	    //fprintf(stderr, "Check branch on path %d at %d -> %d\n", p_id, n->id, n2->id);
	    r |= prune_extra_recurse(g, n, n2, n_end, other_path, vis, nvis);
	}
    }

    return r;
}

void ksw_print_aln(FILE *fp, int len1, char *seq1, int len2, char *seq2, int ncigar, uint32_t *cigar) {
    int i1 = 0, i2 = 0, i =0;
    while (i1 < len1 || i2 < len2) {
	int op = cigar[i] & BAM_CIGAR_MASK;
	int oplen = cigar[i] >> BAM_CIGAR_SHIFT;
	switch(op) {
	case BAM_CMATCH:
	    fprintf(fp,
		    "vv %dM\t%.*s\nvv %dM\t%.*s\n", oplen, oplen, seq1+i1, oplen, oplen, seq2+i2);
	    i1+=oplen;
	    i2+=oplen;
	    break;

	case BAM_CINS:
	    fprintf(fp, "vv %dI\t%.*s\n", oplen, oplen, seq1+i1);
	    i1 += oplen;
	    fprintf(fp, "vv %dI\t", oplen);
	    while (oplen--)
		fputc('-', fp);
	    fputc('\n', fp);
	    break;

	case BAM_CDEL: {
	    fprintf(fp, "vv %dD\t", oplen);
	    int z = oplen;
	    while (z--)
		fputc('-', fp);
	    fputc('\n', fp);
	    fprintf(fp, "vv %dD\t%.*s\n", oplen, oplen, seq2+i2);
	    i2 += oplen;
	    break;
	}

	default:
	    abort();
	}

	i++;
    }
}

// Local maximum snp count in a given window size
int ksw_snp_count(int len1, char *seq1, int len2, char *seq2, int ncigar, uint32_t *cigar, int window) {
    int max_snp = 0;
    int snp = 0;

    char *mis = calloc(len2, 1);
    char *indel = calloc(len2, 1);
    if (!mis || !indel)
	return -1;

#define SNP_WITHIN 20

    int i1 = 0, i2 = 0, i =0, j;
    while (i1 < len1 || i2 < len2) {
	int op = cigar[i] & BAM_CIGAR_MASK;
	int oplen = cigar[i] >> BAM_CIGAR_SHIFT;
	switch(op) {
	case BAM_CMATCH:
	case BAM_CEQUAL:
	case BAM_CDIFF:
	    while (oplen--) {
		if (seq1[i1] != 'N' && seq2[i2] != 'N' &&
		    toupper(seq1[i1]) != toupper(seq2[i2]))
		    mis[i2]=1;
		i1++;
		i2++;
	    }
	    //i1+=oplen;
	    //i2+=oplen;
	    break;

	case BAM_CINS:
	    i1 += oplen;
	    for (j = i2>SNP_WITHIN?i2-SNP_WITHIN:0; j < len2 && j < i2+SNP_WITHIN; j++)
		indel[j]=1;
	    break;

	case BAM_CDEL: {
	    for (j = i2>SNP_WITHIN?i2-SNP_WITHIN:0; j < len2 && j < i2+oplen+SNP_WITHIN; j++)
		indel[j]=1;
	    i2 += oplen;
	    break;
	}

	default:
	    abort();
	}

	i++;
    }

    for (i2 = 0; i2 < len2 && i2 < window; i2++) {
	//if (mis[i2] && indel[i2])
	if (mis[i2])
	    snp++;
    }
    max_snp = snp;

    for (; i2 < len2; i2++) {
	//snp += (mis[i2] && indel[i2]) - (mis[i2-window] && indel[i2-window]);
	snp += mis[i2] - mis[i2-window];
	if (max_snp < snp)
	    max_snp = snp;
    }

    free(mis);
    free(indel);

    fprintf(stderr, "new max snp count %d\n", max_snp);
    return max_snp;
}

// Node n_end has parents p1 and p2 which meet up again at some common
// node n_start.  Find n_start.
int node_common_ancestor(dgraph_t *g, node_t *n_end, node_t *p1, node_t *p2, int *vis, int *nvis,
			 int use_ref, int shrink_align) {
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
	n = n->n_in ? g->node[n->in[0]->n[0]] : NULL;
    }

    // Backtrack up p2 until we find 'visited' set
    n = p2;
    while (n) {
	if (n->visited & (1<<31))
	    break;
	n = n->n_in ? g->node[n->in[0]->n[0]] : NULL;
    }
    int start = n->id;

    // Clear visited status
    n = p1;
    while (n) {
	n->visited &= ~(1<<31);
	n = n->n_in ? g->node[n->in[0]->n[0]] : NULL;
    }

    // Create (reversed) paths
    n = p1;
    while (n && n->id != start) {
	//fprintf(stderr, "1[%d]: %d\n", np1, n->id);
	path1[np1++] = n;
	n->visited = (1<<29);
	n = n->n_in ? g->node[n->in[0]->n[0]] : NULL;
    }

    n = p2;
    while (n && n->id != start) {
	//fprintf(stderr, "2[%d]: %d\n", np2, n->id);
	path2[np2++] = n;
	n->visited = (1<<29)-1;
	n = n->n_in ? g->node[n->in[0]->n[0]] : NULL;
    }

    // Recurse down path1 looking for anything that hits path2 and
    // vice versa, excluding the obvious final bubble merge.
    // In theory our bubble will be small and not contain self-bubbles,
    // but there can be occasions where the smallest bubble contains a
    // longer one.  Eg:
    //       |
    //      .A.       Bubble A to B, but another path from C to D
    //    ./   \.     exists that is longer and not found by the
    //   /       \    initial bubble finding algorithm.
    //  C--O-O-O--D
    //   \.     ./    Path 1 = (B),C,(A)
    //     \. ./      Path 2 = (B),D,(A)
    //       B
    //       |

    // This works and fixes a crash with CHM1_CHM13_2.bam -r 19:3168000-3169500
    // but we have this situation often without it being detrimental
    // (and indeed it's sometimes an improvement) plus it's slower.
    //
    // Alternatively spot we're going wrong and just give up.  It's once in a
    // blue moon that it matters and if we don't realign there it isn't the end
    // of the world.

    if (np1 && np2) {
	if (prune_extra_bubbles(g, path1, np1, n_end, path2[0]->visited, vis, nvis))
	    return -1;
	if (prune_extra_bubbles(g, path2, np2, n_end, path1[0]->visited, vis, nvis))
	    return -1;
    }

    if (!np1 && !np2) {
	free(path1);
	free(path2);
	return -1;
    }


//    // Report
//    n = p1;
//    printf(" => path 1: ");
//    while (n && n->id != start) {
//	printf(" %d", n->id);
//	n = n->n_in ? g->node[n->in[0]->n[0]] : NULL;
//    }
//    printf(" %d, len=%d\n", start, np1);
//
//    n = p2;
//    printf(" => path 2: ");
//    while (n && n->id != start) {
//	printf(" %d", n->id);
//	n = n->n_in ? g->node[n->in[0]->n[0]] : NULL;
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
	n = n->n_in ? g->node[n->in[0]->n[0]] : NULL;
    }
    for (i = g->kmer-2; l && i >= 0; i--)
	// with SHRINK_ALIGN we only end up using n->bases[g->kmer-1] as
	// these get discarded again.
	memcpy(v1[len1++], l->bases[i], 5*sizeof(int));

    n = p2;
    int len2 = 0;
    while (n && n->id != start) {
	memcpy(v2[len2++], n->bases[g->kmer-1], 5*sizeof(int));
	l = n;
	n = n->n_in ? g->node[n->in[0]->n[0]] : NULL;
    }
    for (i = g->kmer-2; l && i >= 0; i--)
	// with SHRINK_ALIGN we only end up using n->bases[g->kmer-1] as
	// these get discarded again.
	memcpy(v2[len2++], l->bases[i], 5*sizeof(int));

    uint32_t *cigar = NULL;
    int ncigar = 0;
    // Extra context (g->kmer-1 bases) to the alignment still helps the alignment itself,
    // despite potentially causing impossible to correct alignment regions if we get
    // gaps within the first kmer.
    char *vs1 = vec2seq(v1, len1);
    char *vs2 = vec2seq(v2, len2);
    if (shrink_align) {
	len1-=g->kmer-1;
	len2-=g->kmer-1;
    }
    fprintf(stderr, "use_ref=%d\n%.*s\n%.*s\n", use_ref, len1, vs1, len2, vs2);
    int s;
    fprintf(stderr, "X128[A][G]=%d,%d mat=%d\n", X128['A']['G'], X128['G']['A'], X128['A']['A']);
    if (use_ref)
	s=ksw_global_end(len1, (uint8_t *)vs1, len2, (uint8_t *)vs2,
			 128, (int8_t *)X128_ref, GOPEN_REF, GEXT_REF, 0,
			 &ncigar, &cigar,  1,1,1,1);
    else
	s=ksw_global_end(len1, (uint8_t *)vs1, len2, (uint8_t *)vs2,
			 128, (int8_t *)X128, GOPEN, GEXT, 0,
			 &ncigar, &cigar,  1,1,1,1);
    fprintf(stderr, "Score=%d\n", s);

    ksw_print_aln(stderr, len1, vs1, len2, vs2, ncigar, cigar);
    if (shrink_align) {
	len1+=g->kmer-1;
	len2+=g->kmer-1;
	cigar=realloc(cigar,++ncigar*sizeof(*cigar));
	cigar[ncigar-1]=BAM_CMATCH + ((g->kmer-1)<<BAM_CIGAR_SHIFT);
    }
//    {
//	fprintf(stderr, "v1 %.*s\nv2 %.*s\n", len1, vs1, len2, vs2);
//	fprintf(stderr, "vx %d;", ncigar);
//	int i;
//	for (i = 0; i < ncigar; i++)
//	    fprintf(stderr, " %d%c", (int)(cigar[i] >> BAM_CIGAR_SHIFT), BAM_CIGAR_STR[cigar[i] & BAM_CIGAR_MASK]);
//	fprintf(stderr, ";\n");
//    }
    free(vs1);
    free(vs2);

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
	if (n->in[i]->n[0] == (np2 ? path2[0]->id : start))
	    continue;
	else
	    n->in[j++] = n->in[i];
    }
    n->n_in = j;

    // Merge v1/v2 into vn.  FIXME v1/v2 into v2 is simpler
    {
	int x1 = 0, x2 = 0, p = 0;

	int cig_ind = 0;
	while (x1 < len1 && x2 < len2) {
	    if (cig_ind >= ncigar)
		abort();

	    int op = cigar[cig_ind] & BAM_CIGAR_MASK;
	    int oplen = cigar[cig_ind++] >> BAM_CIGAR_SHIFT;

	    switch (op) {
	    case BAM_CMATCH:
		while (oplen-- && x1 < len1 && x2 < len2) {
		    // match
		    memcpy(vn[p], v1[x1], sizeof(v1[x1]));
		    vn[p][0] += v2[x2][0];
		    vn[p][1] += v2[x2][1];
		    vn[p][2] += v2[x2][2];
		    vn[p][3] += v2[x2][3];
		    vn[p][4] += v2[x2][4];
		    x1++; x2++;
		    p++;
		}
		break;

	    case BAM_CINS:
		// Already in path1 only, nothing to do with path2;
		while (oplen-- && x1 < len1) {
		    memcpy(vn[p], v1[x1++], sizeof(v1[p]));
		    vn[p][4]++;
		    p++;
		}
		break;

	    case BAM_CDEL:
		// Insertion in path2, link into path1
		while (oplen--) {
		    memcpy(vn[p], v2[x2++], sizeof(v1[p]));
		    vn[p][4]++;
		    p++;
		}
		break;
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
//	    fprintf(stderr, "vn[%02d]={%d,%d,%d,%d,%d}\n",
//		   x1, vn[x1][0], vn[x1][1], vn[x1][2], vn[x1][3], vn[x1][4]);
//	}
//	char *vs1 = vec2seq(vn, p);
//	fprintf(stderr, "vn %.*s\n", p, vs1);
//	free(vs1);
    }


    // Merge path2 into path1.  These are in reverse order.
    // Ie. path1[0] and path2[0] both have *a* shared end (n_end);
    // path1[np1-1] and path2[np2-1] both have *a* shared parent (start).
    // => path 1:  5 4 3
    // => path 2:  14 13 12

    {
	int cig_ind = 0, oplen = 0, op = 0;
	int x1 = 0, x2 = 0, p = 0;

	node_t *l1 = n_end;
	node_t *l2 = n_end;
	while (x1 < np1 && x2 < np2) {
	    node_t *n1 = path1[x1];
	    node_t *n2 = path2[x2];

	    if (oplen == 0) {
		if (cig_ind >= ncigar)
		    abort();
		op = cigar[cig_ind] & BAM_CIGAR_MASK;
		oplen = cigar[cig_ind++] >> BAM_CIGAR_SHIFT;
		//fprintf(stderr, "op %c %d at %d\n", BAM_CIGAR_STR[op], oplen, p);
	    }

	    if (op == BAM_CMATCH) {
		//fprintf(stderr, "M %4d %4d\n", n1->id, n2->id);
		oplen--;
		node_common_ancestor_match(g, n1, n2, n_end, start, vn,
					   path1, np1, path2, np2,
					   &p, &x1, &x2, &l1, &l2);
	    } else if (op == BAM_CINS) {
		// Already in path1 only, nothing to do with path2

		// fixme: inc x1 here again?  Or is it 1 short from above?
		while (oplen-- && x1 < np1) {
		    //fprintf(stderr, "D %4d    -\n", n1->id);
		    //n1->bases[g->kmer-1][4]++;
		    //memcpy(n1->bases, vn[p++], 5*sizeof(int)); // FIXME?
		    memcpy(&n1->bases[g->kmer-1], vn[p++], 5*sizeof(int));
		    n1 = path1[++x1];
		}
		l1 = n1;
		oplen = 0;
	    } else {
		//fprintf(stderr, "I     - %4d\n", n2->id);
		node_common_ancestor_ins(g, n1, n2, n_end, l1, l2,
					 start, vn, path1, np1, path2, np2, oplen,
					 &p, &x1, &x2, &l1, &l2, &n2);
		oplen = 0;
	    }
	}

	if (!l1) {
	    free(path1);
	    free(path2);
	    return -1;
	}

	// gap at end
	// TODO: Consider just culling these nodes so cigar generation turns
	// into soft-clips.
	while (x1 < np1) {
	    // Already in path1 only, nothing to do with path2
	    node_t *n1 = path1[x1];
	    //fprintf(stderr, "d %2d  -\n", n1->id);
	    //n1->bases[g->kmer-1][4]++; // FIXME?
	    memcpy(&n1->bases[g->kmer-1], vn[p++], 5*sizeof(int));
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

	    fprintf(stderr, "i  - %2d\n", n2->id);

	    // tail:
	    // link n2 out to l1 in
	    for (i = 0; i < n2->n_out; i++)
		if (n2->out[i]->n[1] == l1->id)
		    break;
	    if (i == n2->n_out)
		move_edge_out(g, n2, l2, l1);
	    for (i = 0; i < l1->n_in; i++)
		if (l1->in[i]->n[0] == n1->id)
		    l1->in[i]->n[0] = n2->id; // FIXME: plus edge hash?

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

    free(cigar);
    free(v1);
    free(v2);
    free(vn);
    free(path1);
    free(path2);

//	static int n=0;
//	char buf[100];
//	sprintf(buf, "_bub%d.dot", n++);
//	graph2dot(g, buf, 0);

    return 0;
}

// Breadth first search starting from id.
// We keep a track of nodes and paths from every fork.
// As we progress one node at a time per path, we check 'visited'.
// If this is set, then we know we have a bubble and can tell
// which paths are involved.
int find_bubble_from2(dgraph_t *g, int id, int use_ref, int min_depth, int *vis, int *nvis) {
    int p_id = 1;
    int i;
    path_t *head = path_create(p_id++);
    int active = 1;
    int ret = 1;

    // Initialise with one path
    //printf("find_bubble_from2 node %d\n", id);
    path_add(head, g->node[id]);
//    printf("A: Node %d visisted by %d\n", id, head->id);
    g->node[id]->visited = head->id;
    vis[(*nvis)++]=id;

    // Iterate
    while (head && active) {
	path_t *p, *lp = NULL, *pnext = NULL;

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

		    if (g->node[n->out[i]->n[1]]->count < min_depth)
			continue; // skip pointless paths

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
//	    printf(" %d(%d,%d)", n->id, n->in[0]->n[0], p->pn->id);
//	}
//	printf("\n");

	// Step two: check & mark all active nodes for visited status
	for (lp = NULL, p = head; p; lp = p, p = pnext) {
	    pnext = p->next;
	    node_t *n = p->n, *pn = p->pn;

	    if (!n) continue;

	    if (n->visited) {
		int merge_id = n->in[0]->n[0];
//		printf("Bubble detected with node %d (%d, %d)\n",
//		       n->id, merge_id, pn->id);
		//graph2dot(g, "_before.dot", 0);

		// Note this could form a loop, but if so it gets broken for us.
#if 1
		if (node_common_ancestor(g, n, pn, g->node[merge_id], vis, nvis, use_ref, 1) < 0) {
		    ret = -1;
		    goto err;
		}
#else
		if (node_common_ancestor(g, n, pn, g->node[merge_id], vis, nvis, use_ref, 0) < 0) {
		    if (node_common_ancestor(g, n, pn, g->node[merge_id], vis, nvis, use_ref, 1) < 0) {
			ret = -1;
			goto err;
		    }
		}
#endif
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
			printf(" %d(%d,%d)", n->id, n->in[0]->n[0], p->pn->id);
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
	    if (n->in[0]->n[0] != pn->id) {
		int i;
		for (i = 0; i < n->n_in; i++)
		    if (n->in[i]->n[0] == pn->id)
			break;
		assert(i < n->n_in);
		struct edge *tmp = n->in[0];
		n->in[0] = n->in[i];
		n->in[i] = tmp;
	    }

	    n->visited = p->id;
	    vis[(*nvis)++]=n->id;
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

int find_bubbles(dgraph_t *g, int use_ref, int min_depth) {
    int i, found;
    // One for this loop and 1 each for the two prune_extra_bubbles loops
    int *v = malloc(g->nnodes*3 * sizeof(int)), nv = 0;

    for (i = 0; i < g->nnodes; i++)
	g->node[i]->visited = 0;

    do {
	found = 0;

	for (i = 0; i < g->nnodes; i++) {
	    int j;

	    if (g->node[i]->n_in == 0 && !g->node[i]->pruned) {
		//printf("Graph start at %d\n", i);

		nv = 0;
		int b = find_bubble_from2(g, i, use_ref, min_depth, v, &nv);
		if (b < 0)
		    return -1;

		// Only clear the nodes that we visited.
		for (j = 0; j < nv; j++)
		    g->node[v[j]]->visited = 0;

		if (b) {
		    found += b;
		    loop_check(g, 1);
		}
		//break; // only really need main start point?  Unknown...
	    }
	}
    } while (found);

    free(v);
    return 0;
}

void validate_graph(dgraph_t *g) {
    int i, j, k;

    for (i = 0; i < g->nnodes; i++) {
	node_t *n = g->node[i];
	if (n->pruned)
	    continue;

	int above = -1;
	for (j = 0; j < n->n_in; j++) {
	    //assert(n->in[j]->n[1] == n->id);
	    node_t *p = g->node[n->in[j]->n[0]];
	    for (k = 0; k < p->n_out; k++)
		if (p->out[k]->n[1] == n->id)
		    break;
	    assert(k < p->n_out);
	    if (above < p->above)
		above = p->above;
	}
	//assert(n->above == above+1 || (above <= 0 && n->above == 0));
	assert(n->above > above || (above <= 0 && n->above == 0));


	for (j = 0; j < n->n_out; j++) {
	    //assert(n->out[j]->n[0] == n->id);
	    node_t *p = g->node[n->out[j]->n[1]];
	    for (k = 0; k < p->n_in; k++)
		if (p->in[k]->n[0] == n->id)
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
#if 1
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


	fprintf(fp, "  n_%d [label=\"%d @ %d.%d x %d, ^%d", n->id, n->id, n->pos, n->ins, n->count, n->above);
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
	    fprintf(fp, "  n_%d -> n_%d [color=\"blue\"]\n", n->id, n->in[j]->n[0]);
	}
#endif
    }

#ifndef MIN
#  define MIN(a,b) ((a)<(b)?(a):(b))
#endif

#define MAX_GRAPH_KEY 10

#if 1
	// Hash key -> node map
	HashIter *iter = HashTableIterCreate();
	HashItem *hi;
	while (hi = HashTableIterNext(g->node_hash, iter)) {
	    node_t *n = hi->data.p;
	    fprintf(fp, " n_%.*s [label=\"%.*s\", color=\"grey\", group=%d]\n",
	    	    hi->key_len, hi->key,
		    MIN(MAX_GRAPH_KEY, hi->key_len), hi->key, n->id);
		    //1, hi->key + hi->key_len-1, n->id);
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


	//fprintf(fp, "  n_%d [label=\"%d @ %d.%d x %d, ^%d", n->id, n->id, n->pos, n->ins, n->count, n->above);
	fprintf(fp, "  n_%d [label=\"%d\\n", n->id, n->id);
	{
	    int x;
	    for (x=0; x<n->n_hi; x++)
		//fprintf(fp, "%s%.*s", x?",\\n":"",n->hi[x]->key_len, n->hi[x]->key);
		fprintf(fp, "%s%.*s", x?",\\n":"", 5, n->hi[x]->key);
	}
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
	    fprintf(fp, "  n_%d -> n_%d [color=\"blue\"]\n", n->id, n->in[j]->n[0]);
	}
#endif
    }

#ifndef MIN
#  define MIN(a,b) ((a)<(b)?(a):(b))
#endif

#define MAX_GRAPH_KEY 10

#if 0
	// Hash key -> node map
	HashIter *iter = HashTableIterCreate();
	HashItem *hi;
	while (hi = HashTableIterNext(g->node_hash, iter)) {
	    node_t *n = hi->data.p;
	    fprintf(fp, " n_%.*s [label=\"%.*s\", color=\"grey\", group=%d]\n",
	    	    hi->key_len, hi->key,
		    MIN(MAX_GRAPH_KEY, hi->key_len), hi->key, n->id);
		    //1, hi->key + hi->key_len-1, n->id);
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
#endif

int graph2dot_simple(dgraph_t *g, char *fn, int min_count) {
    FILE *fp = stdout;

    char *seq = malloc(g->nnodes+1);
    int *nidx = malloc(g->nnodes * sizeof(*nidx));
    int *start = malloc(g->nnodes * sizeof(*start));
    int nstart = 0;

    if (fn) {
	fp = fopen(fn, "w");
	if (!fp) {
	    perror(fn);
	    return -1;
	}
    }

    // Interesting nodes
    int i, j;
    for (i = 0; i < g->nnodes; i++) {
	node_t *n = g->node[i];
	if (n->pruned)
	    continue;

	if (n->n_in != 1 || n->n_out > 1 ||
	    (n->n_in == 1 && g->node[n->in[0]->n[0]]->n_out != 1)) {
	    start[nstart] = n->id;
	    nidx[n->id] = nstart++;
	}
    }

    fprintf(fp, "digraph g {\n");
    for (i = 0; i < nstart; i++) {
	node_t *n = g->node[start[i]];
	int seq_len = 0;

	if (n->count < min_count)
	    continue;

	for (;;) {
	    char base = 4;
	    int bcount = 0;
	    int k;
	    for (k = 0; k < 4; k++) {
		if (bcount < n->bases[g->kmer-1][k]) {
		    bcount = n->bases[g->kmer-1][k];
		    base = k;
		}
	    }
	    seq[seq_len++] = "ACGTN"[(uint8_t)base];

	    if (n->n_out != 1 || (n->n_out == 1 && g->node[n->out[0]->n[1]]->n_in != 1))
		break;

	    n = g->node[n->out[0]->n[1]];
	}
	seq[seq_len++] = 0;

	fprintf(fp, "  n_%d [label=\"", i);
	char *sp = seq;
	int len = strlen(seq);
	while (len > 0) {
	    fprintf(fp, "%.*s", len>20?20:len, sp);
	    len -= 20;
	    sp += 20;
	    if (len)
		fprintf(fp, "\\n");
	}
	fprintf(fp, "\"]\n");
	for (j = 0; j < n->n_out; j++) {
	    if (g->node[n->out[j]->n[1]]->count < min_count)
		continue;
	    fprintf(fp, "  n_%d -> n_%d [penwidth=%f]\n", i, nidx[n->out[j]->n[1]], log(n->out[j]->count+1)+2);
	}
    }

    fprintf(fp, "}\n");
    if (fn)
	return fclose(fp);

    free(seq);
    free(nidx);
    free(start);

    return 0;
}

typedef struct hseq {
    struct hseq *next;
    char *seq;  // ambig for alignments
    char *seq2; // real kmers, for matching to graph
    int len;
    int score;
} hseqs;

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

// Sequence bases used, barring clips
int mapped_bases(uint32_t *cig, int ncig) {
    int i, n = 0;
    for (i =0; i < ncig; i++) {
	if (bam_cigar_op(cig[i]) != BAM_CSOFT_CLIP &&
	    (bam_cigar_type(bam_cigar_op(cig[i])) & 1))
	    n += bam_cigar_oplen(cig[i]);
    }

    return n;
}

void bam_set_seqi(uint8_t *s, int i, char b) {
    // 0=top 4, 1=bot 4
    // 2=top 4, 2=bot 4 etc
    s[i>>1] = (s[i>>1] & (0xf0 >> ((~i&1)<<2))) | (b << ((~i&1)<<2));
}

void bam_set_seqi_base(bam1_t *b, int pos, char base) {
    static char L[256] = {0};
    static int init = 0;
    if (!init) {
	int i;
	for (i = 0; i < 16; i++)
	    L[(uc)"=ACMGRSVTWYHKDBN"[i]] = i;
	init = 1;
    }

    bam_set_seqi(bam_get_seq(b), pos, L[(uint8_t)base]);
}

// A tail tip which has gone off into unpositioned nodes.
// Try a trivial match back to the aligned path to see how it matches
// and whether we should adjust seq in order to get aligned coords again.
#define TIP_MIN_COUNT 3
#define TIP_MISMATCH 5
int fix_tail(dgraph_t *g, node_t *last, char *seq, int len) {
    int i;

    //fprintf(stderr, "Try tail fix on %.*s\n", len, seq);

    if (!last || last->pos < 0)
	return 0; // cannot fix.

    int score = 0, best_score = 0, best_i = 0;
    for (i = 0; i < len-(g->kmer-1); i++) {
	node_t *n = last, *p = NULL;
	int j;
	if (n->count < TIP_MIN_COUNT)
	    break;

	for (j = 0; j < n->n_out; j++) {
	    if (g->node[n->out[0]->n[1]]->pos >= 0) {
		p = g->node[n->out[0]->n[1]];
		break;
	    }
	}
	if (!p)
	    break;

	for (j = 0; j < p->n_hi; j++) {
	    if (p->hi[j]->key[g->kmer-1] == seq[i+g->kmer-1])
		break;
	}
	if (j == p->n_hi)
	    j = 0;

	//fprintf(stderr, "i=%d score=%d %.*s vs %.*s\n", i, score, p->hi[j]->key_len, p->hi[j]->key, g->kmer, seq+i);

	if (seq[i+g->kmer-1] == p->hi[j]->key[g->kmer-1]) {
	    score++;
	} else {
	    score-=TIP_MISMATCH;
	    seq[i+g->kmer-1] = p->hi[j]->key[g->kmer-1];
	}
	if (best_score < score) {
	    best_score = score;
	    best_i = i;
	}
	last = p;
    }

    if (best_i+g->kmer < len)
	seq[best_i+g->kmer]='x'; // force kmer to break and switch to soft-clip from here on
    //fprintf(stderr, "best_i=%d, best_score=%d\nAlign up to %s\n", best_i, best_score, seq);

    return best_i;
}

//#define NO_QUAL_FIX
// seq2cigar based on the newer find_bubbles and common_ancestor output.
int seq2cigar_new(dgraph_t *g, char *ref, int shift, bam1_t *b, char *seq, int *new_pos,
		  uint32_t **new_cig, uint32_t *new_ncig) {
    int i;
    node_t *n = NULL, *last = NULL;
    int cig_op = 999, cig_len = 0;
    int first_go = 1;
    int seq_start = 0;
#ifndef NO_QUAL_FIX
    char *orig_seq = seq;
#endif
    int len = strlen(seq);
    int ref_len = strlen(ref);

    char *sub = malloc(g->kmer + len + 1), *sub_k = sub + g->kmer;
    memcpy(sub_k, seq, len);
    sub_k[len] = 0;

    *new_pos = b->core.pos; // default to original

    uint32_t cig_a[MAX_CIGAR];
    int cig_ind = 0;

    //fprintf(stderr, "Orig seq = %.*s\n", len, seq);

    int last_ins = 0;
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

    int pos = 0, j;

    if (n->pos >= 0 /*&& (first_go || n->ins == 0)*/ /*|| (n->ref & IS_CON)*/) {
//	int np_dist = 0;
//	if (n->pos < 0) {
//	    // last base of 1st kmer is in an insertion.  Backtrack to see
//	    // if the sequence keeps matching the consensus (todo - we need
//	    // to check the actual sequence does!).
//	    node_t *np = n;
//	    while (np->n_in > 0 /*&& np->pos < 0*/ && (np->ref & IS_CON)) {
//		np = g->node[np->in[0]->n[0]];
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

	    // The first g->kmer j are just moving the node back so the first base
	    // is found and we can match the start of the kmer instead of the end.
	    // Beyond that (g->kmer <= j < g->kmer+i) we're faking up sequence in
	    // and attempt to resolve head-tips.  Do this only provided the alignment
	    // is sufficiently good to original seq.

	    int best_score = 0, score = 0, best_score_j = 1;
	    //for (j = 1; j < g->kmer; j++) {
	    for (j = 1; j - i <= g->kmer; j++) {
		int k;
		node_t *sn = NULL;
		if (seq[-(j-1)+i+g->kmer-1] == sub_k[g->kmer-j]) {
		    if (j<=g->kmer || n->count>=TIP_MIN_COUNT)
			score++;
		    if (best_score < score) {
			best_score = score;
			best_score_j = j;
		    }
		} else {
		    score-=TIP_MISMATCH;
		}
		//fprintf(stderr, "j=%d %d %s %s\n", j, score, seq-(j-1)+i+g->kmer-1, sub_k + g->kmer-j);
	    up_one:
		//printf("Search for %.*s", g->kmer, sub_k-j);
		for (k = 0; k < n->n_in; k++) { // added first only here
		    node_t *s2 = g->node[n->in[k]->n[0]];
		    if (s2->pos < 0)
			continue;
		    // and if j >= g->kmer then s2->count is high?
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

	    j = best_score_j;
	    //fprintf(stderr, "Best score %d at j %d, seq %s %s\n", best_score, j, seq-(j-1)+i+g->kmer-1, sub_k + g->kmer-j);

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
#ifndef NO_QUAL_FIX
	uint8_t *bam_seq = bam_get_seq(b);
	uint8_t *bam_qual = bam_get_qual(b);

	// Copy quality anyway, even if we later "fail".  The alternative
	// is to write to a new qual buffer here and memcpy it later,
	// as we do for cigars.
	//
	// Even when we fail to realign, modifying the qualities based on
	// het indels within STRs is a key way of reducing false positives
	// when calling existing alignments (albeit at the cost of a similar
	// growth in false negatives).
	//
	// Eg chr1:24M-30M CHM1_CHM13_2.bam goes from
	// FP/FN 142/167 to 85/217 (about 60 swing).
#endif
	for (; i <= len - g->kmer; i++) {
#ifndef NO_QUAL_FIX
	    if (orig_seq[i+g->kmer-1] == "=ACMGRSVTWYHKDBN"[bam_seqi(bam_seq, i+g->kmer-1)])
		bam_qual[i+g->kmer-1] += 5;
	    else
		bam_qual[i+g->kmer-1] = bam_qual[i+g->kmer-1]-10>0 ?bam_qual[i+g->kmer-1]-10 :0;
#endif

#define NO_SEQ_FIX
#ifndef NO_SEQ_FIX
	    bam_set_seqi_base(b, i+g->kmer-1, seq[i+g->kmer-1]);
#endif

	    //fprintf(stderr, "base %c vs %c vs %c\n", orig_seq[i+g->kmer-1], seq[i+g->kmer-1], "=ACMGRSVTWYHKDBN"[bam_seqi(bam_get_seq(b), i+g->kmer-1)]);
	    last = n;
	    if (i == len-g->kmer) {
		memcpy(s2, seq+i, g->kmer);
		if (pos+2 + g->kmer >= ref_len) {
		    memset(s2+g->kmer, 'N', g->kmer);
		    if (pos+2 < ref_len)
			memcpy(s2+g->kmer, ref+pos+2, ref_len - (pos+2));
		} else {
		    memcpy(s2+g->kmer, ref+pos+2, g->kmer);
		}
		s1 = s2;
	    }
	    if (!(n = find_node(g, s1++, g->kmer, 0)) || n->pruned)
		break;

	    int fixed_suffix = 0;
	fix_suffix:

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

		fprintf(stderr, "Unmapped kmer %.*s in seq %.*s\n", n->hi[0]->key_len, n->hi[0]->key, len - i, s1-1);
		int new_suffix = fixed_suffix ? 0 : fix_tail(g, last, s1-1, len-i);
		if (new_suffix) {
		    if (!(n = find_node(g, s1-1, g->kmer, 0)) || n->pruned)
			break;
		    fixed_suffix = 1;
		    goto fix_suffix;
		}
		break;

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

	if (i+g->kmer-1 < len) {
	    // trailing unaligned. FIXME: align it or soft-clip as appropriate.
	    ADD_CIGAR(BAM_CSOFT_CLIP, len-(i+g->kmer-1));
	}

	ADD_CIGAR(999, 0); // flush
    } else {
    unmapped:
	// Unmapped
	cig_ind = 0;
    }

    //int old_len = mapped_bases(bam_get_cigar(b), b->core.n_cigar);
    //int new_len = mapped_bases(cig_a, cig_ind);
    //if (new_len < 0.7*old_len) cig_ind = 0;

    *new_ncig = cig_ind;
    if (cig_ind) {
	*new_cig = malloc(cig_ind * sizeof(uint32_t));
	if (!*new_cig)
	    return -1;
	memcpy(*new_cig, cig_a, cig_ind * sizeof(uint32_t));
    } else {
	*new_cig = NULL;
    }

    free(sub);
    return 0;
}

// After a realignment, produce consensus from cigar strings to get indel
// locations and consensus.  Couple this with simple STR detection to
// adjust qualities in regions where heterozygous indels coincide with
// STRs.
int calc_BAQ(bam1_t **bams, int nbam, int ref_start, int ref_len, char *ref) {
    int i, j;
    int (*cons)[5] = calloc(ref_len, sizeof(*cons));
    int *indel = calloc(ref_len, sizeof(*indel));
    if (!cons || !indel)
	return -1;

    //             1 2   4       8
    //             A C   G       T
    int L[16] = {4,0,1,4,2,4,4,4,3,4,4,4,4,4,4,4}; // bam_seqi to 0,1,2,3

    // Compute a basic consensus and indel marker.
    for (i = 0; i < nbam; i++) {
	bam1_t *b = bams[i];
	uint32_t *cig = bam_get_cigar(b), cig_ind = 0;
	uint32_t ncig = b->core.n_cigar;
	uint64_t rp = b->core.pos, sp = 0;
	uint8_t *seq = bam_get_seq(b);

	for (j = 0; j < ncig; j++) {
	    uint32_t op = cig[cig_ind] & BAM_CIGAR_MASK;
	    uint32_t oplen = cig[cig_ind++] >> BAM_CIGAR_SHIFT;

	    switch (op) {
	    case BAM_CMATCH:
	    case BAM_CEQUAL:
	    case BAM_CDIFF:
		while (oplen--) {
		    if (rp >= ref_start && rp < ref_start + ref_len)
			cons[rp-ref_start][L[bam_seqi(seq, sp)]]++;
		    sp++;
		    rp++;
		}
		break;

	    case BAM_CINS:
		if (rp >= ref_start && rp < ref_start + ref_len)
		    indel[rp-ref_start]++;
		if (rp > ref_start && rp <= ref_start + ref_len)
		    indel[rp-ref_start-1]++;
		sp += oplen;
		break;

	    case BAM_CSOFT_CLIP:
		sp += oplen;
		break;

	    case BAM_CDEL:
		while (oplen--) {
		    if (rp >= ref_start && rp < ref_start + ref_len)
			indel[rp-ref_start]++;
		    sp++;
		    rp++;
		}
		break;

	    case BAM_CREF_SKIP:
		rp += oplen;
		break;
	    }
	}
    }

    // Now combine to get consensus string + bit status per pos.
    // 1 STR in cons
    // 2 STR in ref
    // 4 homozygous indel
    // 8 heterozygous indel
    char *seq = malloc(ref_len+1);
    uint8_t *str_status = calloc(ref_len, 1);
    if (!seq || !str_status)
	return -1;
    
    for (i = 0; i < ref_len; i++) {
	int c = 0, b = 4, j = 0;
	for (j = 0; j < 4; j++) {
	    if (c < cons[i][j]) {
		c = cons[i][j];
		b = j;
	    }
	}
	seq[i] = "ACGTN"[b];
	if (indel[i] >= .9*c)
	    str_status[i] |= 4; // hom indel present
	else if (indel[i] >= .1*c)
	    str_status[i] |= 8; // het indel present
    }
    seq[i] = 0;

    // Find STRs within seq.
    rep_ele *reps, *elt, *tmp;
    reps = find_STR(seq, ref_len, 0);
#define STR_MARGIN 5
    DL_FOREACH_SAFE(reps, elt, tmp) {
	int left  = elt->start-STR_MARGIN < 0 ? 0 : elt->start-STR_MARGIN;
	int right = elt->end+STR_MARGIN >= ref_len ? ref_len-1 : elt->end+STR_MARGIN;

	for (i = left; i <= right; i++)
	    if (str_status[i])
		break;

	if (i <= right) {
	    for (i = left; i <= right; i++)
		if (seq[i] != 'N') str_status[i] |= 1; // cons
	}

	DL_DELETE(reps, elt);
	free(elt);
    }

    // Find STRs within ref
    if (ref) {
	reps = find_STR(ref, ref_len, 0);
	DL_FOREACH_SAFE(reps, elt, tmp) {
	    int left  = elt->start-STR_MARGIN < 0 ? 0 : elt->start-STR_MARGIN;
	    int right = elt->end+STR_MARGIN >= ref_len ? ref_len-1 : elt->end+STR_MARGIN;

	    for (i = left; i <= right; i++)
		if (str_status[i])
		    break;

	    if (i <= right) {
		for (i = left; i <= right; i++)
		    str_status[i] |= 2;
	    }

	    DL_DELETE(reps, elt);
	    free(elt);
	}
    }

//    for (i = 0; i < ref_len; i++) {
//	fprintf(stderr, "%5d %c %x %d  %d %d %d %d\n", i, ref[i], str_status[i], indel[i], cons[i][0], cons[i][1], cons[i][2], cons[i][3]);
//    }

    uint8_t nq[256];
    for (i = 0; i < 256; i++) {
	nq[i] = i>7 ?i-5 : 2;
    }
    // old:    37/69, 30/137    68/218, 46/355
    
    // MARGIN 2
    // *= .7:  36/74, 33/138    62/220, 52/351
    // -5      37/72, 33/137    64/215, 52/349 // single subtract
    // -10     36/74, 33/138    61/224, 51/351
    // -7-7    35/75, 32/138    60/220, 53/351 // sub, and again if het-indel
    // -5-5    37/72, 33/137    62/215, 52/351
    // -5-5,+T 35/74, 30/139    59/225, 47/359 // with old STR trim baq too and calc_baq always
    // +T-5-5  37/70, 30/137    62/219, 47/355 // with old STR trim, and calq_baq on err only
    // +T-7-7  36/72, 30/138    61/220, 47/355 // "
    // -T-5-5  37/71, 32/135    62/214, 51/353 // no old STR trimming at all; STR helps FP indel
    // -T-7-7  35/74, 32/136    60/219, 52/354 // "
    //
    // MARGIN 0
    // -T-5-5   37/69, 32/135    64/214, 52/353
    //
    // MARGIN 5
//>>// -T-5-5   37/72, 32/135    60/215, 51/353
    // -T-5-55  36/74, 32/134    57/221, 51/353 // indel=nq[x] and het_indel=nq[nq[nq[x]]
    // -T-5-5+2 37/72, 32/135    61.215, 51/353 // ie T-5-3
    // -T-5-5-5 37/72, 32/135    60/215, 51/353 // indel, +het, +nospan; identical to T-5-5
    // -T-7-7   36/76, 32/136    56/226, 51/354
    
    nq[0] = 0; nq[1] = 1;

    // Now reprocess again, identifying bases in reads that neighbour hard STRs.
    for (i = 0; i < nbam; i++) {
	bam1_t *b = bams[i];
	uint32_t *cig = bam_get_cigar(b), cig_ind = 0;
	uint32_t ncig = b->core.n_cigar;
	uint64_t rp = b->core.pos, sp = 0;
	uint8_t *qual = bam_get_qual(b);

	for (j = 0; j < ncig; j++) {
	    uint32_t op = cig[cig_ind] & BAM_CIGAR_MASK;
	    uint32_t oplen = cig[cig_ind++] >> BAM_CIGAR_SHIFT;

	    switch (op) {
	    case BAM_CMATCH:
	    case BAM_CEQUAL:
	    case BAM_CDIFF:
		while (oplen--) {
		    if (rp >= ref_start && rp < ref_start + ref_len) {
			if ((str_status[rp-ref_start] & 3) && sp < b->core.l_qseq) {
			    qual[sp] = nq[qual[sp]];
			    if (str_status[rp-ref_start] & 8)
				qual[sp] = nq[qual[sp]];
			}
		    }
		    sp++;
		    rp++;
		}
		break;

	    case BAM_CINS:
		if (rp >= ref_start && rp < ref_start + ref_len) {
		    if (str_status[rp-ref_start] & 3) {
			while (oplen-- && sp < b->core.l_qseq) {
			    qual[sp] = nq[qual[sp]];
			    if (str_status[rp-ref_start] & 8)
				qual[sp] = nq[qual[sp]];
			    sp++;
			}
			oplen = 0;
		    }
		}
		sp += oplen;
		break;

	    case BAM_CSOFT_CLIP:
		sp += oplen;
		break;

	    case BAM_CDEL:
		while (oplen--) {
		    if (rp >= ref_start && rp < ref_start + ref_len) {
			if ((str_status[rp-ref_start] & 3) && sp < b->core.l_qseq) {
			    qual[sp] = nq[qual[sp]];
			    if (str_status[rp-ref_start] & 8)
				qual[sp] = nq[qual[sp]];
			}
		    }
		    sp++;
		    rp++;
		}
		break;

	    case BAM_CREF_SKIP:
		rp += oplen;
		break;
	    }
	}
    }

    free(cons);
    free(indel);
    free(seq);
    free(str_status);

    return 0;
}

int int64_compar(const void *vp1, const void *vp2) {
    return ((*(const uint64_t *)vp2)>>32) - ((*(const uint64_t *)vp1)>>32);
}

//-----------------------------------------------------------------------------
// Error correction via kmer hashing.
//
// The basic algorithm is:
//
// 0. Split sequence into a series of overlapping kmers.
// 1. Count kmers
// 2. Observe kmer frequency distribution to find a typical depth
//    for correct sequence.
// 3. Assume very rare kmers correspond to sequencing errors, vs frequent
//    ones which are good.  (NB: VERY frequent means repetitive sequence)
// 4. Produce a mapping of all single base edits from "good" kmers
//    (kmer neighberhood) and link them to their original kmer.
//    NB: neighbours can clash. If so exclude from neighbour mapping.
// 5. Loop through all kmers in sequences, and if bad (rare) and matching
//    a neighbour of a good kmer, then correct it.
//
// Note errk may be anything from 8 (array) to 30 (needs sparse
// array / hashing).  Hence we use hash tables here.
int correct_errors_fast(haps_t *h, int n, int errk, int min_count) {
    HashTable *kmer_hash = NULL, *neighbours = NULL;
    HashItem *hi;
    int i, counth = 0, countw = 0, tbases = 0;
    string_alloc_t *sp = string_pool_create(errk*10000);

    char **new_seq = calloc(n, sizeof(char *));
    if (!new_seq)
	return -1;

    HashTable *hash = kmer_hash;

    hash = HashTableCreate(8192, HASH_DYNAMIC_SIZE | HASH_POOL_ITEMS |
			   HASH_NONVOLATILE_KEYS | HASH_FUNC_TCL);
    kmer_hash = hash;

    //-----------
    // Hash words
    for (i = 0; i < n; i++) {
	char *seq = h[i].seq;
	int len = strlen(seq), j;

	tbases += len;
	for (j = 0; j < len-errk; j++) {
	    HashData hd;
	    HashItem *hi;
	    int nw;

	    hd.i = 0;
	    hi = HashTableAdd(hash, seq+j, errk, hd, &nw);
	    hi->data.i++;

	    counth++;
	    countw+=nw;
	}
    }

    //-----------
    // Compute kmer depth threshhold that corresponds to correcting ~90%
    // of the data set, and assume these represent correct kmers.
    int thresh = 0, count_good = 0;
    {
	// Discard 10% of words and error correct those to remaining 90%
	int F[256] = {0};
	HashIter *hiter = HashTableIterCreate();
	while ((hi = HashTableIterNext(hash, hiter))) {
	    //fprintf(stderr, "KMER %.*s freq %d\n", errk, hi->key, hi->data.i);
	    F[hi->data.i > 255 ? 255 : hi->data.i]++;
	}
	HashTableIterDestroy(hiter);

	// Compute mean of values > 2 and s.d.
	int cnt = 0, sum = 0, sum_sq = 0;
	for (i = 3; i < 256; i++) {
	    if (!F[i])
		continue;
	    //fprintf(stderr, "Dist[%d]=%d\n", i, F[i]);
	    cnt += F[i];
	    sum += i*F[i];
	    sum_sq += i*i*F[i];
	    if (F[i] >= min_count)
		count_good += F[i]; // approximation
	}
	if (!cnt) {
	    free(new_seq);
	    return 0;
	}
	double mean = (double)sum/cnt;
	double sd = sqrt(sum_sq/cnt - mean*mean);
	fprintf(stderr, "Mean %f sd %f => %d..%d\n", mean, sd,
		(int)(mean-sd*2), (int)(mean+sd*2+.999));

	thresh = mean-sd*2 > 3 ? mean-sd*2 : 3;
	if (thresh > 3) {
	    // 2nd pass on the bulk of the data
	    int i_end = mean+sd*2+.999;
	    cnt = sum = sum_sq = 0;
	    for (i = thresh; i < i_end; i++) {
		if (!F[i])
		    continue;
		cnt += F[i];
		sum += i*F[i];
		sum_sq += i*i*F[i];
	    }
	    if (!cnt) {
		free(new_seq);
		return 0;
	    }
	    mean = (double)sum/cnt;
	    sd = sqrt(sum_sq/cnt - mean*mean);

	    fprintf(stderr, "Mean %f sd %f => %d..%d\n", mean, sd,
		    (int)(mean-sd*2), (int)(mean+sd*2+.999));
	}

	thresh = mean-sd*2 > min_count ? mean-sd*2 : min_count;
    }


    //-----------
    // Find common words and produce neighbourhoods
    //
    // The slow bit...
    // Consider building neighbour at single point (so 4 neighbours per kmer
    // rather than 4 x kmer) and comparing all words vs neighbour.  Sliding
    // all seq kmers against this ought to work, but in practice does not.
    //
    // Instead we build the complete neighbourhood and compare every
    // kmer/2 within sequence vs all neighbours.

    fprintf(stderr, "%d unique words, %d est. good, %d total words, threshold %d\n",
	    countw, count_good, counth, thresh);
    neighbours = HashTableCreate(count_good*16, HASH_DYNAMIC_SIZE | HASH_POOL_ITEMS |
				 HASH_NONVOLATILE_KEYS | HASH_FUNC_TCL);

    int nn = 0;
    HashIter *hiter = HashTableIterCreate();
    while ((hi = HashTableIterNext(hash, hiter))) {
        if (hi->data.i >= thresh) {
	    int j;
	    //for (j = 0; j < errk; j++) {
	    for (j = 0; j < errk; j += errk-1) { // 1st and last base only
		int nw, k;
		HashData hd;
		hd.p = hi->key;
		int base = hi->key[j];
		for (k = 0; k < 5; k++) {
		    HashItem *hi2;
		    if ("ACGTN"[k] == base) continue;
		    char *N = string_alloc(sp, errk);
		    memcpy(N, hi->key, errk);
		    N[j] = "ACGTN"[k];
		    hi2 = HashTableAdd(neighbours, N, errk, hd, &nw);
		    //fprintf(stderr, "Add neighbour %.*s -> %.*s: had=%d\n", errk, (char *)N, errk, (char *)hi->key, !nw);
		    if (!nw) {
			// This deletion step works better if we have full neighbour
			// analysis via j++ loop above, instead of j += errk-1.
			//
			//fprintf(stderr, "Del neighbour %.*s\n", errk, (char *)N);
			hi2->data.p = NULL; // mark the clash
		    }
		    nn++;
		}
	    }
	}
    }
    HashTableIterDestroy(hiter);
    //fprintf(stderr, "%d neighbours hashed, est. %d\n", nn, count_good*8);


    //-----------
    // Now scan through all kmers in our sequences to see if they're rare but
    // with a close match to a known good kmer (neighbour).  Correct if so.
    int nc = 0;
    int nbases = 0;
    for (i = 0; i < n; i++) {
	char *seq = h[i].seq;
	int len = strlen(seq), j;
	char *s2 = NULL;

#undef EDGE_DIST
//#define EDGE_DIST 3
#define EDGE_DIST 0
	for (j = EDGE_DIST; j < len-errk-EDGE_DIST; j++) {
	    HashItem *hi, *hi2;

	    // Ditch common words
	    hi = HashTableSearch(hash, seq+j, errk);
	    //fprintf(stderr, "Word %.*s freq %d\n", errk, seq+j, hi?(int)hi->data.i:-1);
	    if (hi && hi->data.i >= min_count) {
		//fprintf(stderr, "Word %.*s not in neighbour hash %p %d, min %d\n", errk, seq+j,
		//	hi, hi?(int)hi->data.i:-1, min_count);
		continue;
	    }

	    // Not common, but also not a neighbour of a common kmer?
	    hi2 = HashTableSearch(neighbours, seq+j, errk);
	    if (!hi2 || !hi2->data.p)
		continue;

	    //fprintf(stderr, "Word %.*s matches %.*s\n", errk, seq+j, errk, (char *)hi2->data.p);

	    // Found a match, find diff and correct it.
	    int k;
	    for (k = 0; k < errk; k++)
		if (seq[j+k] != ((char *)hi2->data.p)[k])
		    break;
	    if (k == errk)
		continue;

	    // We modify a new copy here so we don't change the hash keys.
	    //
	    // However we should consider pointing seq to this temporary copy
	    // so we can fix neighbouring mismatching bases (and remove the
	    // j+=k at the loop end).
	    if (!s2) {
		new_seq[i] = s2 = string_alloc(sp, len+1);
		memcpy(s2, seq, len);
		s2[len] = '\0';
	    }

	    s2[j+k] = ((char *)hi2->data.p)[k];
	    nc++;
	    //fprintf(stderr, "Correct %.*s %d -> %.*s\n",
	    //        errk, hi->key, (int)hi->data.i, errk, (char *)hi2->data.p);

	    j += k; // skip to the next word not containing the corrected base
	}

	// debug
	if (s2) {
	    for (j = 0; j < len; j++)
		if (s2[j] != seq[j])
		    nbases++;
	}
    }
    fprintf(stderr, "Corrected %d bases (%5.2f%%)\n",
	    nbases, 100.0*nbases/tbases);

    // Finally replace the original sequences with the edited ones.
    // We delay this to here so our hash table keys don't change under
    // us, allowing us to avoid costly strdups in the hash table code.
    for (i = 0; i < n; i++) {
	if (new_seq[i])
	    strcpy(h[i].seq, new_seq[i]);
    }

    string_pool_destroy(sp);
    HashTableDestroy(hash, 0);
    HashTableDestroy(neighbours, 0);
    free(new_seq);

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
    static HashTable *hash = NULL;
    if (!hash) {
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
    }

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

    //HashTableDestroy(hash, 0);
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
		if (n2->in[k]->n[0] == n->id) {
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
	    //printf("Best prev = %d\n", n->in[best_i]->n[0]);
	    n = g->node[n->in[best_i]->n[0]];
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

// Use the n->above fields to compute the consensus from the path
// with longest/most use.
//
// Returns SNP cluster count, defined as the maximum number of
// SNPs within a given sliding window.
int compute_consensus_above(dgraph_t *g, char *ref, int window) {
    int i, j;
    node_t *n;
    char *seq = malloc(g->nnodes + g->kmer + 1);
    int *nnum = malloc((g->nnodes + g->kmer + 1)*sizeof(int));

    // Find trailing node with highest above count.
    int above = 0, end_n = 0, start_n = 0;
    for (i = 0; i < g->nnodes; i++) {
	n = g->node[i];
	if (n->pruned || n->n_out)
	    continue;
	if (above < n->above)
	    above = n->above, end_n = n->id;
    }

    // Recurse up
    n = g->node[end_n];
    j = g->nnodes + g->kmer + 1;
    seq[--j] = 0;
    nnum[j] = -1;

    while (n) {
	seq[--j] = vec2X(n->bases[g->kmer-1]);
	nnum[j] = n->id;

	for (above = i = 0; i < n->n_in; i++) {
	    node_t *p = g->node[n->in[i]->n[0]];
	    if (above <= p->above)
		above = p->above, start_n = n->in[i]->n[0];
	}
	n = n->n_in ? g->node[start_n] : NULL;
    }

    ref += g->kmer-1;
    //fprintf(stderr, "vc Cons from node %d to %d: %s\nRef %s\n", start_n, end_n, &seq[j], ref);

    // Compare vs ref.
    char *cons = seq+j;
    int *cid = nnum+j;
    int cons_len = strlen(cons);
    int ref_len = strlen(ref);

    uint32_t *cigar = NULL;
    int ncigar = 0;
    ksw_global_end(cons_len, (uint8_t *)cons, ref_len, (uint8_t *)ref,
		   //128, (int8_t *)X128, GOPEN_REF, GEXT_REF, 0,
		   128, (int8_t *)X128, 4,1, 0,
		   &ncigar, &cigar, 0,0,0,0);

    //ksw_print_aln(stderr, cons_len, cons, ref_len, ref, ncigar, cigar);

    int snp_count = ksw_snp_count(cons_len, cons, ref_len, ref, ncigar, cigar, window);

//    {
//	fprintf(stderr, "v1 %.*s\nv2 %.*s\n", cons_len, cons, ref_len, ref);
//	fprintf(stderr, "VX %d;", ncigar);
//	int i;
//	for (i = 0; i < ncigar; i++)
//	    fprintf(stderr, " %d%c", (int)(cigar[i] >> BAM_CIGAR_SHIFT), BAM_CIGAR_STR[cigar[i] & BAM_CIGAR_MASK]);
//	fprintf(stderr, ";\n");
//    }

    int x1 = 0, x2 = 0, k = 0;
    while (x1 < cons_len || x2 < ref_len) {
	if (k >= ncigar)
	    break;
	int op = cigar[k] & BAM_CIGAR_MASK;
	int oplen = cigar[k++] >> BAM_CIGAR_SHIFT;

	switch(op) {
	case BAM_CMATCH:
	    while (oplen--) {
		// Match
		//fprintf(stderr, "%c %c @ %d %d\n", cons[x1], ref[x2], cid[x1], x2);
		g->node[cid[x1]]->pos=x2;
		g->node[cid[x1]]->ins=0;
		x1++, x2++;
	    }
	    break;

	case BAM_CINS: {
	    // Insertion in consensus
	    int ival = 0;
	    while (oplen--) {
		g->node[cid[x1]]->pos=x2;
		g->node[cid[x1]]->ins=++ival;
		//fprintf(stderr, "%c - @ %d %d.%d\n", cons[x1], cid[x1], x2, ival);
		x1++;
	    }
	    break;
	}

	case BAM_CDEL:
	    // Deletion in consensus
	    while (oplen--) {
		//fprintf(stderr, "- %c () %d\n", ref[x2], x2);
		x2++;
	    }
	    break;

	default:
	    abort();
	}
    }

    free(cigar);
    free(seq);
    free(nnum);

    return snp_count;
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

#if 0
static void dump_input(bam_hdr_t *hdr, bam1_t **bams, int nbams, char *ref, int ref_len, int ref_start,
		       char *cons1, char *cons2, int cons_len) {
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
}
#endif

// Add 'number of nodes above' and 'number of nodes below' figures to
// every node.  If we have an incoming fork, number of nodes above is
// the MAX of left & right incoming forks.  Similar for nodes below
// and outgoing fork.
//
// Then once we have that, starting from any node in the graph, head
// all the way down taking the highest number going down, then all the
// way up. At this point we're now furthest up the graph.  Head down
// once more and we're now furthest point at bottom.  Ie down, up,
// down with the final "up, down" being the longest consensus path.

// Recurse down from n incrementing n->above as we go.
// On an out-fork follow both routes.
// On an in-fork, push the node onto the queue of nodes to-do (if
// it appears we have a higher count) and return.
void number_nodes_above_recurse(dgraph_t *g, node_t *n, int *queue, int *nqueue, int qn) {
    // Starting point, 1 more than previously.
    int i, count = 0, first = 1;

    // FIXME: try summing n->count instead of +1 per layer
    for (i = 0; i < n->n_in; i++)
	//if (count < g->node[n->in[i]->n[0]]->above + 1)
	//    count = g->node[n->in[i]->n[0]]->above + 1;
	if (count < g->node[n->in[i]->n[0]]->above + g->node[n->in[i]->n[0]]->count+1)
	    count = g->node[n->in[i]->n[0]]->above + g->node[n->in[i]->n[0]]->count+1;

    //fprintf(stderr, "node %d, count %d, in %d, out %d\n", n->id, count, n->n_in, n->n_out);

    while (n) {
	first--;
	int curr_above = n->above;
	if (n->above < count)
	    n->above = count;
	//count++;
	count += n->count+1;
	if (n->n_in > 1) {
	    int c2 = n->above;
	    for (i = 0; i < n->n_in; i++) {
		//if (c2 < g->node[n->in[i]->n[0]]->above + 1)
		//    c2 = g->node[n->in[i]->n[0]]->above + 1;
		if (c2 < g->node[n->in[i]->n[0]]->above + g->node[n->in[i]->n[0]]->count+1)
		    c2 = g->node[n->in[i]->n[0]]->above + g->node[n->in[i]->n[0]]->count+1;
	    }
	    if (c2 > curr_above) {
		n->above = c2;
		for (i = qn; i < *nqueue; i++)
		    if (queue[i] == n->id)
			break;
		if (i == *nqueue)
		    queue[(*nqueue)++] = n->id;
		assert(*nqueue < g->nnodes);
	    }
	    if (first < 0)
		return;
	}

	if (n->n_out <= 0)
	    break;

	if (n->n_out == 1) {
	    n = g->node[n->out[0]->n[1]];
	    continue;
	} else {
	    for (i = 0; i < n->n_out; i++)
		number_nodes_above_recurse(g, g->node[n->out[i]->n[1]], queue, nqueue, qn);
	    break;
	}
    }
}

int number_nodes_above(dgraph_t *g) {
    int i;
    int *queue = malloc(g->nnodes * sizeof(int));
    int nqueue = 0;

    // Build a queue of starting points
    for (i = 0; i < g->nnodes; i++) {
	node_t *n = g->node[i];
	if (n->pruned || n->n_in != 0)
	    continue;
	queue[nqueue++] = i;
	n->above = 0;
    }

    // Recurse down each node in queue, adding new ones as we see fit.
    for (i = 0; i < nqueue; i++) {
	node_t *n = g->node[queue[i]];
	number_nodes_above_recurse(g, n, queue, &nqueue, i);
    }

    free(queue);
    return 0;
}



static int default_kmer = KMER;
int bam_realign(bam_hdr_t *hdr, bam1_t **bams, int nbams, int *new_pos,
		char *ref_seq, int ref_len, int ref_start,
		char *cons1, char *cons2, int cons_len,
		int max_snp, int window) {
    int i, kmer = default_kmer, ret = -1;
    dgraph_t *g = NULL;
    haps_t *haps = NULL, *cons = NULL;
    haps_t *ref = NULL, *ref_ = NULL;
    int plus10 = 1;

    //dump_input(hdr, bams, nbams, ref_seq, ref_len, ref_start, cons1, cons2, cons_len);

    init_X128_score(X128,      0,-4,4); // for seqs against each other to form pileup
    init_X128_score(X128_ref, -3,-3,1); // for seq pileup against reference

    // Convert BAMS to seqs instead, so we can perform edits on them
    // while retaining the original data.
    int nhaps = nbams;
    haps = bam2haps(bams, nhaps);

    // Successive rounds allows for fixing more than 1 error in a read.
    fprintf(stderr, "Correcting\n");
    correct_errors_fast(haps, nhaps, 29, 3);
    correct_errors_fast(haps, nhaps, 27, 3);
    correct_errors_fast(haps, nhaps, 25, 2);
    correct_errors_fast(haps, nhaps, 20, 2);
    correct_errors_fast(haps, nhaps, 14, 2);

    //    correct_errors(haps, nhaps, 14, 2, 0);

    // tried:
    // 27, 25, 20, 14 (good; 82B cycles)
    // 25, 20, 14     (?  think close to top one. RECHECK parse error...)
    // 25, 20         (lots of unmapped / unaligned, but ok result)
    // 25, 14         (66B cycle; many unmapped / unaligned, but ...)

//    for (i = 0; i < 10; i++)
//	if (correct_errors(haps, nhaps, 14, 2, 5) <= 0)
//	    break;
    //trim_adapters(haps, nhaps, "/nfs/srpipe_references/adapters/adapters.fasta", 14);
    // Also trim low quality terminal bases
    trim_adapters(haps, nhaps, "adapter.fa", ADAPTER_KMER, 10);

 bigger_kmer:
    if (kmer > MAX_KMER)
	goto err;

    // If our mapped sequences can be processed with a sensible kmer,
    // but the unmapped cannot, then maybe it's just the unmapped is
    // a broken homo-polymer.  In which case keep it unmapped and ignore
    // it.
    for (; kmer < MAX_KMER; kmer += 10) {
	fprintf(stderr, "Building graph with kmer=%d, nseqs=%d\n", kmer, nhaps);
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

    graph2dot(g, "f.dot", 0);
    //graph2dot_simple(g, "x.dot", 5);

    // Merge bubbles excluding reference
    if (find_bubbles(g, 0, 10) < 0) {
	graph_destroy(g);
	if ((kmer += 10) < MAX_KMER)
	    goto bigger_kmer;
	g = NULL;
	goto err;
    }
    if (loop_check(g, 0)) {
	graph_destroy(g);
	if ((kmer += 10) < MAX_KMER)
	    goto bigger_kmer;
	g = NULL;
	goto err;
    }
    if (find_bubbles(g, 0, 2) < 0) {
	graph_destroy(g);
	if ((kmer += 10) < MAX_KMER)
	    goto bigger_kmer;

	g = NULL;
	goto err;
    }
    if (loop_check(g, 0)) {
	graph_destroy(g);
	if ((kmer += 10) < MAX_KMER)
	    goto bigger_kmer;
	g = NULL;
	goto err;
    }
    loop_check(g, 0);
    graph2dot(g, "_F.dot", 0);

    int shift = bams[0]->core.pos+1;

    if (ref_seq) {
	ref_ = ref = malloc(sizeof(*ref));
	if (!ref)
	    goto err;
	//fprintf(stderr, "Alloc %p\n", ref);
	ref->name = "Ref";
	ref->pos = ref_start;
	//fprintf(stderr, "ref=%.*s\n", ref_len, ref_seq);

	ref->seq = malloc(ref_len + g->kmer-1 + 1);
	ref->seq[ref_len + g->kmer-1] = 0;
	memcpy(ref->seq + g->kmer-1, ref_seq, ref_len);
	memset(ref->seq, 'N', g->kmer-1);
	if (add_seq(g, ref->seq, 0, IS_REF) < 0 || loop_check(g, 0)) {
	    fprintf(stderr, "Loop when adding reference\n");
	    graph_destroy(g);
	    if ((kmer += 10) < MAX_KMER) {
		free_haps(ref_, 1); ref_ = NULL;
		goto bigger_kmer;
	    }
	    fprintf(stderr, "No suitable kmer found\n");
	    g = NULL; goto err;
	}

	shift = ref_start + 1;
    }

    graph2dot(g, "_g.dot", 0);

    if (!(cons = compute_consensus(g)))
	goto err;

    fprintf(stderr, "cons=%s\n", cons->seq);
    if (add_seq(g, cons->seq, 0, (ref?0:IS_REF)|IS_CON) < 0 || loop_check(g, 0)) {
	fprintf(stderr, "Loop when adding consensus\n");
	graph_destroy(g);
	if ((kmer += 10) < MAX_KMER) {
	    free_haps(cons, 1); cons = NULL;
	    free_haps(ref_, 1); ref_ = NULL;
	    goto bigger_kmer;
	}
	g = NULL;
	goto err;
    }

    // Auto tune kmer + extra 10 to make alignments a bit better
    if (plus10 && kmer+KMER_INC_FINAL < MAX_KMER && KMER_INC_FINAL) {
	plus10--;
	graph_destroy(g);
	g = NULL;
	kmer += KMER_INC_FINAL;
	free_haps(cons, 1); cons = NULL;
	free_haps(ref_, 1); ref_ = NULL;
	goto bigger_kmer;
    }

    if (!ref)
	ref = cons;

    // Now also merge in reference bubbles
    if (find_bubbles(g, 1, 0) < 0) {
	graph_destroy(g);
	if ((kmer += 10) < MAX_KMER)
	    goto bigger_kmer;

	g = NULL;
	goto err;
    }
    assert(loop_check(g, 0) == 0);

    graph2dot(g, "_h.dot", 0);

    // Finally try merging in the consensus, both primary and secondary,
    // incase this makes further bubbles, which during collapse will merge
    // in some heads/tips.
    //
    // NB we have no phasing on cons1/2.  Possibly a real allele may
    // switch from cons1 to cons2 and back again, but it's a simple
    // approximation which resolves some issues.
    //
    // FIXME: cons should be actual consensus with indels, rather than
    // based on reference coordinates.  See gap5 pileup for this instead?
    //
    // As it stands, this doesn't improve things.  We get around 10% fewer
    // FN, but at a cost of 10% more FP.  This doesn't fix the ~40% more FN
    // we have vs original (unrealigned) data, so we need something more than
    // that anyway to fix this.  Find that first. (Suspect it's left vs right
    // alignment justification causing that and SNP vs Indel preferences.)

    // c2 only as c1 should be represented better by main assembly?
    // FIXME: make optional
    fprintf(stderr, "Cons len = %d\n", cons_len);
    if (cons1 || cons2) {
	if (cons1) add_seq(g, cons1, cons_len, 0);
	if (cons2) add_seq(g, cons2, cons_len, 0);
	if (loop_check(g, 0)) {
	    // adding consensus caused loop; give up on that idea.
	    cons1 = cons2 = NULL;
	    free_haps(cons, 1); cons = NULL;
	    free_haps(ref_, 1); ref_ = NULL;
	    graph_destroy(g);
	    if ((kmer += 10) < MAX_KMER)
		goto bigger_kmer;
	    g = NULL;
	    goto err;
	}
	if (find_bubbles(g, 1, 0) < 0) {
	    graph_destroy(g);
	    if ((kmer += 10) < MAX_KMER)
		goto bigger_kmer;
	    g = NULL;
	    goto err;
	}
    }

    // Map the graph to the reference.  We previously did this by
    // adding the reference as a sequence and collapsing bubbles
    // again, but this sometimes causes more issues.  Now we try
    // the alternative approach of producing the primary linear path
    // through the graph representing the bulk of the data and then
    // align that as a whole to the reference.
    number_nodes_above(g);
    if (compute_consensus_above(g, ref->seq, window) - max_snp >= 2) { // eg 2 in 15
	// Reject as too many extra clustered variants.
	// Likely we've made an error somewhere.

	fprintf(stderr, "Rejecting realignment as too many new SNPs\n");
	goto err;
    }
    graph2dot(g, "_G.dot", 0);

    // Compute new cigar strings
    uint32_t **new_cig = calloc(nhaps, sizeof(*new_cig));
    uint32_t *new_ncig = calloc(nhaps, sizeof(*new_ncig));
    if (!new_cig || !new_ncig)
	goto err;

    int orig_unmapped = 0, unmapped = 0;
    for (i = 0; i < nhaps; i++) {
	if (seq2cigar_new(g, ref->seq, shift, bams[i], haps[i].seq, &new_pos[i],
			  &new_cig[i], &new_ncig[i]) < 0)
	    goto err;
	orig_unmapped += (bams[i]->core.flag & BAM_FUNMAP) ? 1 : 0;
	unmapped += (new_ncig[i] == 0);
    }

//#define CHECK_DEPTH // seem detrimental overall.
#ifdef CHECK_DEPTH
    if ((double)unmapped/nhaps <= 0.3 || unmapped <= 1.3*orig_unmapped)
#endif
	{
	    for (i = 0; i < nhaps; i++) {
		bam1_t *b = bams[i];
		if (new_ncig[i]) {
		    b->core.flag &= ~BAM_FUNMAP;
		    replace_cigar(b, new_ncig[i], new_cig[i]);
		    free(new_cig[i]);
		} else {
		    new_pos[i] = b->core.pos; // ie. don't change it
		    b->core.flag |= BAM_FUNMAP;
		}
	    }
	}
#ifdef CHECK_DEPTH
    else goto err;
#endif
    free(new_cig);
    free(new_ncig);

    ret = 0;
 err:
    calc_BAQ(bams, nhaps, ref_start, ref_len, ref_seq);

    graph_destroy(g);
    free_haps(haps, nhaps);
    free_haps(cons, 1);
    free_haps(ref_, 1);

    fprintf(stderr, "Finished.\n");

    return ret;
}

//=============================================================================
#ifdef TEST_MAIN
//=============================================================================
//
// Basic load function to read an entire BAM file.
// This is just a test for small subsets, rather than streaming
// a large file.
bam1_t **load_bam(char *fn, int *nrecs, bam_hdr_t **hdr_p) {
    samFile *in = sam_open(fn, "r");
    bam_hdr_t *hdr;
    int nalloc = 128, nr = 0;
    bam1_t **bams = calloc(128, sizeof(*bams));

    if (!in)
	return NULL;

    if (!(*hdr_p = hdr = sam_hdr_read(in)))
	return NULL;

    bams[nr] = calloc(1, sizeof(bam1_t));
    while (sam_read1(in, hdr, bams[nr]) >= 0) {
	if (++nr >= nalloc) {
	    nalloc *= 2;
	    bams = realloc(bams, nalloc * sizeof(*bams));
	    memset(&bams[nalloc/2], 0, nalloc/2*sizeof(*bams));
	}
	bams[nr] = calloc(1, sizeof(bam1_t));
    }

    sam_close(in);

    *nrecs = nr;
    return bams;
}

#define MAX_LINE 1000000

haps_t *load_fasta(char *fn, int *nhaps) {
    haps_t *h = NULL;
    int nh = 0, nalloc = 0;
    char line[MAX_LINE];
    FILE *fp = fopen(fn, "r");

    if (!fp) {
	perror(fn);
	return NULL;
    }

    // Simple 2 line affair only
    while (fgets(line, MAX_LINE, fp)) {
	int l = strlen(line);
	if (line[l-1] == '\n')
	    line[l-1] = 0;

	if (line[0] == 0)
	    continue;

	if (line[0] != '>') {
	    fprintf(stderr, "Unknown format\n");
	    return NULL;
	}

	if (nh >= nalloc) {
	    nalloc = nh ? nh*2 : 64;
	    h = realloc(h, nalloc*sizeof(*h));
	    if (!h)
		return NULL;
	}

	h[nh].name = strdup(line+1);
	if (!fgets(line, MAX_LINE, fp)) {
	    if (nhaps) *nhaps = nh;
	    fclose(fp);
	    return h;
	}

	l = strlen(line);
	if (line[l-1] == '\n')
	    line[l-1] = 0;
	char *line_start = line;
	while (isspace(*line_start))
	    line_start++;
	h[nh].seq = strdup(line_start);
	h[nh].pos = 0;
	h[nh].qual = NULL;

	nh++;
    }

    fprintf(stderr, "Loaded %d seqs\n", nh);

    if (nhaps) *nhaps = nh;
    fclose(fp);
    return h;
}

int main(int argc, char **argv) {
    bam_hdr_t *hdr;
    int nbams, *newpos, start = 1, i;

    if (argc > 1 && strncmp(argv[1], "-k", 2) == 0) {
	default_kmer = atoi(argv[1]+2);
	argc--;
	argv++;
    }

    bam1_t **bams = load_bam(argv[1], &nbams, &hdr);

    if (!(newpos = calloc(nbams, sizeof(*newpos))))
	return 1;

    char *ref = NULL;
    if (argc > 2) {
	fprintf(stderr, "===> Adding reference\n");
	haps_t *h = load_fasta(argv[2], 0);
	ref = h->seq;
	start = atoi(h->name)-1;
	free(h);
    }

    if (bam_realign(hdr, bams, nbams, newpos, ref, ref?strlen(ref):0, start, NULL, NULL, 0, 0, 15) < 0) {
	fprintf(stderr, "Realign failed\n");
	return 1;
    }

    // FIXME/TODO: sort output
    samFile *fp;
    if (!(fp = sam_open("-", "w")))
	return 1;
    if (sam_hdr_write(fp, hdr) < 0)
	return 1;
    for (i = 0; i < nbams; i++) {
	bams[i]->core.pos = newpos[i];
	if (sam_write1(fp, hdr, bams[i]) < 0)
	    return 1;
    }
    if (sam_close(fp) < 0)
	return 1;


    free(ref);
    for (i = 0; i < nbams; i++)
	bam_destroy1(bams[i]);
    free(bams);
    bam_hdr_destroy(hdr);

    return 0;
}
#endif

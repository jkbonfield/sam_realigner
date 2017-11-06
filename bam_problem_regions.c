// Basic operation.
//
// - Process the BAM file in pileup fashion, column by column.
// - Identify regions of interest (het indels, mass soft-clipping, etc).
// - If within N bases on interesting region
//       Add reads to queue until sufficient outside region.
//   Else if leaving an interesting region
//   (was within N bases, but not now or changed reference)
//       Realign queue  [main algorithm goes here]
//   Else
//       Add to queue
//       Flush queue up to N bases way from current location.
//       (Or when we change reference ID.)
//
// We need to keep track of "start of interest" and "end of interest",
// merging interesting locations into a single region if they're
// within N bases.
//
// Also track min_interest and max_interest for reads that *overlap*
// start_interest and end_interest.  Set a flag on these reads
// ("realign me"). Additionally any read that doesn't itself overlap
// start_interest to end_interest but is contained entirely within
// min_interest to max_interest needs to be labelled as "realign me"
// too in order to maintain sort order.
//
// Then leaving region of interest is simply defined as a new read
// starting beyond max_interest (or switching reference).  If the
// number of reads in the queue is too large, flush anyway and
// consider the problem too hard to realign.

//#define DEBUG

#define REALIGNER_VERSION "0.1"

/*
 * Prunes quality based on snp calling score.
 * 
 * Bi-allelic consensus is computed (possibly as hom).
 * If call has high confidence, then bin to 2 qualities.
 * - high if base is one of the 1-2 consensus calls.
 * - low otherwise.
 *
 * Exceptions.
 * - Within 'D' bases of any read indel.
 * - Within 'D' bases of a soft-clip (possible clipped indel).
 * - Any sequence with mapping quality of <= M.
 * - Any column with high discrepancy, implying possibly tri-allelic.
 *
 */

/*
 * TODO: Het indels may mean we quantise the indel but not the
 * non-indel.  This gives a bias.  Need to treat both the same.
 *
 * Identify low complexity regions.  If we're preserving confidence
 * then it needs to be for all bases in that low-complexity section so
 * that local realignment doesn't change qualities.
 *
 * Consider using STR method for detecting start/end range of bases to
 * keep for SNPs as well as indels, particularly when the sequence is
 * soft-clipped during an STR.  (Possibility of misaligned bases
 * then.)  Or just a low complexity filter? Or concordant soft-clips?
 *
 * For indels, consider the soft-clip adjustment on reads ending in
 * STR adjacent to indels as these bases don't confirm the count.
 */

// Default params
// FIXME: ideally we should have MARGIN+STR size, so an indel at pos
// 1234 in an STR spanning 1234 to 1250 should consider the problem to
// be 1234-MARGIN to 1250+MARGIN instead of to 1234+MARGIN.
#define MARGIN 100       // margin around suspect regions
#define CON_MARGIN 100   // margin to add on to consensus when realigning

#define MIN_INDEL 2      // FIXME: parameterise
#define MAX_DEPTH 20000  // FIXME: parameterise

// Prevalence of low mapping quality, > PERC => store all
// Lower => larger files
#define LOW_MQUAL_PERC 0.5

// Amount of variable sized insertion
#define INS_LEN_PERC 0.1

// Percentage of seqs with large soft-clips
#define CLIP_PERC 0.1

// Amount of over-depth to consider this as as suspect region
#define OVER_DEPTH 3.0

// Percentage of reads spanning indel.
#define INDEL_OVERLAP_PERC 0.5

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <getopt.h>
#include <inttypes.h>
#include <unistd.h>

#include <htslib/sam.h>
#include <htslib/khash.h>
#include <htslib/faidx.h>

#include "tree.h"

#ifndef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))
#endif

#ifndef MAX
#define MAX(a,b) ((a)>(b)?(a):(b))
#endif

#ifndef ABS
#define ABS(a) ((a)>=0?(a):-(a))
#endif

KHASH_SET_INIT_INT(aux_exists)
typedef khash_t(aux_exists) *auxhash_t;

typedef struct {
    char *ref;
    char *region;
    int verbose;
    int margin, cons_margin;
    double low_mqual_perc;
    double clip_perc;
    double ins_len_perc;
    double over_depth;
    double indel_ov_perc;
    int clevel;
    int cons1, cons2;
} cram_realigner_params;

//-----------------------------------------------------------------------------
// Tree of bam objects, sorted by chromosome & position

typedef struct bam_sorted_item {
    RB_ENTRY(bam_sorted_item) link;
    bam1_t *b;
    uint64_t id;
    int end_pos;
} bam_sorted_item;

RB_HEAD(bam_sort, bam_sorted_item);

typedef struct bam_sort bam_sorted_list;

static int bam_item_cmp(bam_sorted_item *b1, bam_sorted_item *b2) {
    int d;
    if (b2->b->core.tid == -1)
	return -1;
    if ((d = b1->b->core.tid - b2->b->core.tid))
	return d;
    if ((d = b1->b->core.pos - b2->b->core.pos))
	return d;
    return b1->id - b2->id;
}

RB_PROTOTYPE(bam_sort, bam_sorted_item, link, bam_item_cmp);
RB_GENERATE(bam_sort, bam_sorted_item, link, bam_item_cmp);

bam_sorted_list *bam_sorted_list_new(void) {
    bam_sorted_list *bl = calloc(1, sizeof(bam_sorted_list));
    if (!bl)
        return NULL;
        
    RB_INIT(bl);

    return bl;
}

void bam_sorted_list_destroy(bam_sorted_list *bl) {
    bam_sorted_item *node, *next;
    
    if (!bl)
        return;

    for (node = RB_MIN(bam_sort, bl); node; node = next) {
	next = RB_NEXT(bam_sort, bl, node);
	RB_REMOVE(bam_sort, bl, node);
	free(node);
    }

    free(bl);
}

bam_sorted_item *bam_sorted_item_new(void) {
    return calloc(1, sizeof(bam_sorted_item));
}

/*
 * Inserts a bam structure into the bam_sorted_list, in positional
 * order. 
 *
 * Returns the newly created bam_sorted_list element on success.
 *         NULL on failure.
 */
bam_sorted_item *insert_bam_list2(bam_sorted_list *bl, bam_sorted_item *ele) {
    RB_INSERT(bam_sort, bl, ele);

    return ele;
}

bam_sorted_item *insert_bam_list_id(bam_sorted_list *bl, bam1_t *b, uint64_t id) {
    bam_sorted_item *ele = bam_sorted_item_new();
    if (!ele)
	return NULL;
    
    ele->b = b;
    ele->id = id;
    ele->end_pos = bam_endpos(ele->b);

    return insert_bam_list2(bl, ele);
}

static uint64_t global_id = 0;
bam_sorted_item *insert_bam_list(bam_sorted_list *bl, bam1_t *b) {
    return insert_bam_list_id(bl, b, global_id++);
}

/*
 * Removes an item from the bam_sorted_list.
 */
void remove_bam_list(bam_sorted_list *bl, bam_sorted_item *ele) {
    RB_REMOVE(bam_sort, bl, ele);

    free(ele);
}


//-----------------------------------------------------------------------------
typedef struct {
    // sam, bam AND hts! All the weasles in a single bag :)
    samFile *fp;
    bam_hdr_t *header;
    hts_itr_t *iter;
    bam_sorted_list *b_hist;
    bam1_t *b_unmap;
    int64_t count_in, count_out;
    cram_realigner_params *params;
} pileup_cd;

int flush_bam_list(pileup_cd *cd, bam_sorted_list *bl, int before_tid, int before, samFile *out, bam_hdr_t *header) {
    bam_sorted_item *bi, *next;

    // Sanity check
    int last_pos = 0;
    RB_FOREACH(bi, bam_sort, bl) {
	if (!(bi->b->core.flag & BAM_FUNMAP)) {
	    assert(bi->b->core.pos >= last_pos);
	    last_pos = bi->b->core.pos;
	}
    }

    for (bi = RB_MIN(bam_sort, bl); bi; bi = next) {
	next = RB_NEXT(bam_sort, bl, bi);
	if (bi->end_pos >= before ||
	    (bi->b->core.tid >= 0 && bi->b->core.tid >= before_tid))
	    break;

	cd->count_out++;
	if (sam_write1(out, header, bi->b) < 0)
	    return -1;
	bam_destroy1(bi->b);
	remove_bam_list(bl, bi);
    }

    return 0;
}

//-----------------------------------------------------------------------------
// Interface to the realigner code
int check_overlap(bam_sorted_list *bl, int start, int *start_ovl, int *end_ovl) {
    bam_sorted_item *bi, *next;
    int count = 0;
    RB_FOREACH(bi, bam_sort, bl) {
	if (bi->end_pos < start)
	    continue;
	    
	if (*start_ovl > bi->b->core.pos)
	    *start_ovl = bi->b->core.pos;
	if (*end_ovl < bi->end_pos)
	    *end_ovl = bi->end_pos;
    }

    return 0;
}

int realign_list(pileup_cd *cd, bam_hdr_t *hdr, bam_sorted_list *bl,
		 char *cons, char *cons2, int cons_len,
		 int start, int end, int start_ovl, int end_ovl,
		 faidx_t *fai) {
    bam_sorted_item *bi, *next;
    int count = 0, ba_sz = 0;
    bam1_t **ba = NULL;

    // FIXME start bi at >= start_ovl and then filter until core.pos > end.
    RB_FOREACH(bi, bam_sort, bl) {
	if (bi->b->core.pos < start_ovl) {
	    //fprintf(stderr, "   [%d..%d (of %d..%d)]\n", bi->b->core.pos, bi->end_pos, start, end);
	    continue;
	}

	if (bi->b->core.pos > end) {
	    //fprintf(stderr, "   [%d..%d (of %d..%d)] END\n", bi->b->core.pos, bi->end_pos, start, end);
	    break;
	}

	// Do we want to try realigning unmapped data too? Maybe...
	if (bi->b->core.flag & BAM_FUNMAP)
	    continue;

	if (count >= ba_sz) {
	    ba_sz = ba_sz ? ba_sz*2 : 256;
	    ba = realloc(ba, ba_sz * sizeof(*ba));
	    if (!ba)
		return -1;
	}
	ba[count++] = bi->b;
	//fprintf(stderr, "    %d..%d (of %d..%d)\n", bi->b->core.pos, bi->end_pos, start, end);
    }
    fprintf(stderr, "    Realign %d reads\n", count);

    if (count == 0)
	return 0;

    // Realign the selection
    int *new_pos = malloc(count * sizeof(*new_pos));
    if (!new_pos) {
	free(ba);
	return -1;
    }

    int i;
    for (i = 0; i < count; i++)
	new_pos[i] = ba[i]->core.pos;
    extern int bam_realign(bam_hdr_t *hdr, bam1_t **bams, int nbams, int *new_pos, char *ref, int ref_len, int ref_pos, char *cons1, char *cons2, int len);

    char region[1024];
    snprintf(region, 1023, "%s:%d-%d", hdr->target_name[ba[0]->core.tid], start_ovl+1, end_ovl);
    int seq_len;

    char *ref = NULL;
    if (fai)
	ref = fai_fetch(fai, region, &seq_len);
    if (!ref) {
	ref = cons;
	seq_len = cons_len;
    }

    if (bam_realign(hdr, ba, count, new_pos, ref, seq_len, start_ovl,
		    cd->params->cons1 ? cons : NULL,
		    cd->params->cons2 ? cons2 : NULL,
		    cons_len) < 0) {
	free(new_pos);
	free(ba);
	if (fai)
	    free(ref);
	return -1;
    }
    if (fai)
	free(ref);

    // Resort the list as realignment can move reads around.
    bi = RB_MIN(bam_sort, bl);
    i = 0;
    while (bi) {
	next = RB_NEXT(bam_sort, bl, bi);
	if (i < count && bi->b == ba[i]) {
	    // Could optimise this, but simpler to remove, update, insert.
	    RB_REMOVE(bam_sort, bl, bi);
	    bi->b->core.pos = new_pos[i];
	    RB_INSERT(bam_sort, bl, bi);
	    i++;
	}

	bi = next;
    }
    assert(i == count); // FIXME: CON_MARGIN 50 triggers this sometimes.  Why?
    
    free(new_pos);
    free(ba);

    return 0;
}

//-----------------------------------------------------------------------------
// Main pileup iterator

int pileup_callback(void *vp, bam1_t *b) {
    pileup_cd *cd = (pileup_cd *)vp;
    int ret = (cd->iter)
	? sam_itr_next(cd->fp, cd->iter, b)
	: sam_read1(cd->fp, cd->header, b);

    if (ret >= 0) {
	cd->count_in++;

	// Unmapped chromosome => end of pileup.  Record the read we've
	// already read and then feign EOF so we can handle these outside
	// of the pileup interface.
	if (b->core.tid == -1) {
	    cd->b_unmap = bam_dup1(b);
	    return -1;
	}

	int unmap = b->core.flag & BAM_FUNMAP;
	if (!unmap) {
	    // Mapped reads that have no matching ref location, eg
	    // they are entirely insertion, will not appear in any pileup
	    // column.  Therefore treat them as unmapped to avoid errors.
	    uint32_t *cig = bam_get_cigar(b);
	    int i, n = b->core.n_cigar;
	    for (i = 0; i < n; i++)
		if (bam_cigar_type(bam_cigar_op(cig[i])) & 2)
		    break;
	    if (i == n) {
#ifdef DEBUG
		printf("Note: %s mapped at #%d,%d but has no cigar ref op!\n",
		       bam_get_qname(b), b->core.tid, b->core.pos);
#endif
		unmap = 1;
	    }
	}

	if (unmap)
	    insert_bam_list(cd->b_hist, bam_dup1(b));
    }

    return ret;
}

// Turns an absolute reference position into a relative query position within the seq.
int ref2query_pos(bam1_t *b, int pos) {
     uint32_t *cig = bam_get_cigar(b);
     int i, n = b->core.n_cigar;
     int p = b->core.pos, q = 0;

     for (i = 0; i < n; i++) {
 	int op = bam_cigar_op(cig[i]);
 	int op_len = bam_cigar_oplen(cig[i]);
	if (p + ((bam_cigar_type(op) & 2) ? op_len : 0) < pos) {
	    if (bam_cigar_type(op) & 1) // query
		q += op_len;
	    if (bam_cigar_type(op) & 2) // ref
		p += op_len;
	    continue;
	}

	if (bam_cigar_type(op) & 1) // query
	    q += (pos - p); // consume partial op_len

	return q >= 0 ? q : 0;
     }

     return q;
}


// Converts a position on a query sequence to a position on the reference.
int bam_qpos2rpos(bam1_t *b, int qpos) {
    const uint32_t *cigar = bam_get_cigar(b);
    int k, rpos = b->core.pos, aqpos = 0; // accumulated qpos
    for (k = 0; k < b->core.n_cigar && aqpos < qpos; k++) {
	if (bam_cigar_type(bam_cigar_op(cigar[k]))&2) {
	    if (bam_cigar_oplen(cigar[k]) <= qpos - aqpos)
		rpos += bam_cigar_oplen(cigar[k]);
	    else
		rpos += qpos - aqpos;
	}
	if (bam_cigar_type(bam_cigar_op(cigar[k]))&1)
	    aqpos += bam_cigar_oplen(cigar[k]);
    }
    return rpos;
}

static char *tid_name(bam_hdr_t *h, int tid) {
    return h->target_name[tid];
}

// Cut down from Crumble's main loop
int transcode(cram_realigner_params *p, samFile *in, samFile *out,
	      bam_hdr_t *header, hts_itr_t *h_iter) {
    bam_plp_t p_iter;
    int tid, pos, last_tid = -2;
    int n_plp;
    const bam_pileup1_t *plp;
    pileup_cd cd = {0};
    bam_sorted_list *b_hist = bam_sorted_list_new();
    int counter = 0;

    cd.fp = in;
    cd.header = header;
    cd.iter = h_iter;
    cd.b_hist = b_hist;
    cd.count_in = 0;
    cd.count_out = 0;
    cd.params = p;

    p_iter = bam_plp_init(pileup_callback, &cd);
    bam_plp_set_maxcnt(p_iter, INT_MAX);
    int min_pos = INT_MAX, max_pos = 0;
    int min_pos2 = INT_MAX, max_pos2 = 0;

    int64_t total_depth = 0, total_col = 0;

    // Flip between OK and in-PROBLEM status, for tracking regions
    // to realign.
    enum status_t {S_OK, S_PROB};
    enum status_t status = S_OK;

    int start_reg = INT_MAX, end_reg = INT_MIN;
    int start_ovl = INT_MAX, end_ovl = INT_MIN;

    int right_most = 0;
    int flush_pos = 0;

    // Primary and secondary consensus computed on-the-fly.  Optionally
    // used in assembly graph generation.
    char *cons = malloc(1024), *cons2 = malloc(1024);
    assert(p->cons_margin < 1024);
    assert(p->margin < 1024);
    memset(cons,  'N', 1024);
    memset(cons2, 'N', 1024);
    int cons_sz = 1024;
    int start_cons = -1;

    if (!cons || !cons2)
	return -1;

    faidx_t *fai = p->ref ? fai_load(p->ref) : NULL;
    while ((plp = bam_plp_auto(p_iter, &tid, &pos, &n_plp))) {
	int i, preserve = 0, indel = 0;
	unsigned char base;

	if (tid != last_tid) {
	    // Ensure b_hist is only per chromosome
	    if (flush_bam_list(&cd, b_hist, tid, INT_MAX, out, header) < 0)
		return -1;
	    last_tid = tid;
	    start_reg = INT_MAX, end_reg = INT_MIN;
	    start_ovl = INT_MAX, end_ovl = INT_MIN;
	    total_depth = 0;
	    total_col = 0;

	    start_cons = -1;
	    memset(cons,  'N', cons_sz);
	    memset(cons2, 'N', cons_sz);

	    right_most = 0;
	    flush_pos = 0;
	}

	if (start_cons == -1)
	    start_cons = pos - p->cons_margin;

	// For avg depth; see depth renorm later on.
	total_depth += n_plp;
	total_col++;

	if (n_plp > MAX_DEPTH) {
	    if (p->verbose > 1)
		fprintf(stderr, "Excessive depth at tid %d, pos %d, depth %d\n", tid, pos, n_plp);
	    goto too_deep;
	}

	if (counter++ == 100000) {
	    if (p->verbose)
		fprintf(stderr, "Processing %s:%d\n", tid_name(header,tid), pos);
	    counter = 0;
	}

	if (h_iter) {
	    if (pos < h_iter->beg)
		continue;
	    if (pos >= h_iter->end)
		break;
	}

	// Counter number of inserted and deleted bases to pileup.
	// Homozygous cases aren't an issue, nor are very rare ones.
	// The mixed case may require realginment.
	//
	// TODO: consider counting separately for indels in middle of
	// read vs end of reads.  If they're in the same proportion
	// then it's probably correct, but indels near end of read
	// are usually misaligned, similarly lack of indels at end of
	// read when everything else has the indel.
	int has_ins = 0, has_del = 0, suspect = 0;
	int left_most = n_plp ? plp[0].b->core.pos : 0;
	int rm_tmp = right_most;
	int freq[256] = {0};
	for (i = 0; i < n_plp; i++) {
	    if (rm_tmp == 0 || plp[i].is_head) {
		int end_pos = bam_endpos(plp[i].b);
		if (right_most < end_pos)
		    right_most = end_pos;
	    }
	    // FIXME: check location along read.
	    // score more for end than middle.
	    if (plp[i].indel < 0 || plp[i].is_del)
		has_del++;
	    if (plp[i].indel > 0)
		has_ins++;

	    unsigned char base = bam_seqi(bam_get_seq(plp[i].b), plp[i].qpos);
	    base = seq_nt16_str[base]; // fixme, only helpful for debugging
	    freq[plp[i].is_del ? '*' : base] += bam_get_qual(plp[i].b)[plp[i].qpos];
	}
	if (has_ins > MIN_INDEL && has_ins < n_plp-MIN_INDEL)
	    suspect |= 16;
	if (has_del > MIN_INDEL && has_del < n_plp-MIN_INDEL)
	    suspect |= 16;

	// Build a consensus, used when we don't have a reference.
	unsigned char call = 'N', call2 = 'N';
	if (1 || !has_del) {
	    int mv = 0, mv2 = 0;
	    if      (mv  < freq['A']) call2 = call, mv2 = mv, mv  = freq[call  = 'A'];
	    else if (mv2 < freq['A'])                         mv2 = freq[call2 = 'A'];
	    if      (mv  < freq['C']) call2 = call, mv2 = mv, mv  = freq[call  = 'C'];
	    else if (mv2 < freq['C'])                         mv2 = freq[call2 = 'C'];
	    if      (mv  < freq['G']) call2 = call, mv2 = mv, mv  = freq[call  = 'G'];
	    else if (mv2 < freq['G'])                         mv2 = freq[call2 = 'G'];
	    if      (mv  < freq['T']) call2 = call, mv2 = mv, mv  = freq[call  = 'T'];
	    else if (mv2 < freq['T'])                         mv2 = freq[call2 = 'T'];
	    if (mv2 < freq['*'])
		mv2 = freq[call2 = '*'];
	    if (mv2 < 20 || mv2 < .2*mv)
		//call2 = 'N';
		call2 = call;
	}
	//fprintf(stderr, "%d %c\n", pos, call);
	while (cons_sz < pos - start_cons+1) {
	    cons_sz *= 2;
	    cons  = realloc(cons,  cons_sz);
	    cons2 = realloc(cons2, cons_sz);
	    if (!cons || !cons2)
		return -1;
	}
	cons [pos - start_cons] = call;
	cons2[pos - start_cons] = call2;

	// Check for unexpectedly deep regions.
	if (n_plp*(total_col+1) > p->over_depth * (total_depth+1))
	    suspect |= 1;

	// Depth renorm: keep average depth to within the last few Mb only.
	if (total_col > 1024*1024) {
	    total_col   >>= 1;
	    total_depth >>= 1;
	}

	// Look for excess of soft-clipping.
	// Also gather stats on the frequency of each indel size.
	// We expect only 2 alleles.  If more then we had better
	// double check the alignments.
	int indel_sz = 0;
	int indel_depth[101];
	indel_depth[0] = 0;
	int clipped = 0, n_overlap = 0;
	for (i = 0; i < n_plp; i++) {
	    int is_indel = (plp[i].indel || plp[i].is_del);

	    if ((plp[i].is_head && plp[i].qpos > 0) ||
		(plp[i].is_tail && plp[i].qpos+1 < plp[i].b->core.l_qseq))
		clipped++;

	    if (!plp[i].is_tail && !plp[i].is_head)
		n_overlap++;

	    if (!plp[i].is_head && !plp[i].is_tail && (plp[i].indel > 0 || (has_ins + has_del))) {
		while (indel_sz < plp[i].indel && indel_sz < 100)
		    indel_depth[++indel_sz] = 0;
		if (plp[i].indel >= 0)
		    indel_depth[MIN(plp[i].indel, 99)]++;
	    }
	}

	if ((clipped - 1.0) >= p->clip_perc * n_overlap && clipped > 1) {
	    if (p->verbose > 1)
		fprintf(stderr, "%s %d\tUnexpected high clip rate, %d of %d\n",
			tid_name(header,tid), pos, clipped, n_overlap);
	    suspect |= 2;
	}

	// Over the span of an indel, the ratio of top to total should
	// be consistently 1 or 2 things, similarly the depth.
	if (indel_sz) {
	    //printf("Indel at %d, max size %d\n", pos, indel_sz);
	    //int qv1 = 0, qv2 = 0;
	    int qd1 = 0, qd2 = 0;
	    int indel_overlap = 0;
	    for (i = 0; i <= indel_sz && i < 100; i++) {
		if (!indel_depth[i])
		    continue;
	    	//printf("%3d\t%2d\n", i, indel_depth[i]);
		indel_overlap += indel_depth[i];
		if (qd1 < indel_depth[i]) {
		    qd2 = qd1; //qv2 = qv1;
		    qd1 = indel_depth[i];
		    //qv1 = i;
		} else if (qd2 < indel_depth[i]) {
		    qd2 = indel_depth[i];
		    //qv2 = i;
		}
	    }
	    //printf("Top 2 = %d x %d,  %d x %d, out of %d, ov/n_plp=%f\n", qv1, qd1, qv2, qd2, indel_overlap, (double)indel_overlap / n_plp);
	    
	    if ((indel_overlap - qd1 - qd2) > p->ins_len_perc * (indel_overlap + .1)) {
		if (p->verbose > 1)
		    fprintf(stderr, "%s %d\tSuspect indel, depth %d / %d, common %d+%d\n",
			    tid_name(header,tid), pos, n_plp, indel_overlap, qd1, qd2);
		suspect |= 4;
	    }

	    if ((double)indel_overlap < p->indel_ov_perc * n_plp) {
		if (p->verbose > 1)
		    fprintf(stderr, "%s %d\tSuspect drop in indel overlap %d vs %d\n",
			    tid_name(header,tid), pos, indel_overlap, n_plp);
		suspect |= 8;
	    }
	}


	// If suspect, expand the regions out and compute overlap
	// extents.
	switch(status) {
	case S_OK:
	    if (suspect) {
		start_reg = pos - p->margin;
		end_reg = pos + p->margin;
		start_ovl = left_most;
		end_ovl = right_most;
		check_overlap(b_hist, start_reg, &start_ovl, &end_ovl);
		fprintf(stderr, "PROB start %d .. %d (%d .. %d)\n", start_reg, end_reg, start_ovl, end_ovl);
		status = S_PROB;
	    } else {
		// Periodic flush here of reads ending before  pos-MARGIN
		flush_pos = left_most - p->margin;
	    }
	    break;

	case S_PROB:
	    if (suspect) {
		end_reg = pos + p->margin;
		if (end_ovl < right_most)
		    end_ovl = right_most;
		//fprintf(stderr, "PROB extnd %d .. %d (%d .. %d)\n", start_reg, end_reg, start_ovl, end_ovl);
	    } else {

		int start_ovl2 = MAX(1, start_ovl - p->cons_margin);
		int end_ovl2 = MAX(1, end_ovl + p->cons_margin);

		if (left_most > end_reg && pos > end_ovl + p->cons_margin) {
		    fprintf(stderr, "PROB end %d %d .. %d (%d .. %d)\n", pos, start_reg, end_reg, start_ovl, end_ovl);
		    realign_list(&cd, header, b_hist,
				 cons  + start_ovl2 - start_cons,
				 cons2 + start_ovl2 - start_cons,
				 end_ovl2 - start_ovl2,
				 start_reg, end_reg,
				 start_ovl2, end_ovl2,
				 fai);
		    start_reg = end_reg = left_most;
		    start_ovl = end_ovl = 0;
		    status = S_OK;
		    flush_pos = left_most - p->margin;
		}
	    }

	    // FIXME: consider the case of building too large a
	    // backlog to realign.  Realign subsets and hope it's
	    // better than nothing?
	    break;
	}

    too_deep:

	// Migrate any finished sequence from the bl to b_hist lists.
	for (i = 0; i < n_plp; i++) {
	    bam1_t *b = plp[i].b;

	    if (!plp[i].is_tail)
		continue;

	    // Note: this may reorder seqs that start at the same coord,
	    // so we give it the read-id to preserve the order.
	    insert_bam_list(b_hist, bam_dup1(b));
	}

	// Flush history (preserving sort order).
	if (flush_bam_list(&cd, b_hist, INT_MAX, flush_pos, out, header) < 0)
	    return -1;
    }
    fai_destroy(fai);

    if (plp) {
	int i;
	for (i = 0; i < n_plp; i++) {
	    bam1_t *b = plp[i].b;
	    insert_bam_list(b_hist, bam_dup1(b));
	}
    }

    if (flush_bam_list(&cd, b_hist, INT_MAX, INT_MAX, out, header) < 0)
	return -1;

    // Handle trailing unmapped reads
    if (cd.b_unmap) {
	int next = 0;
	do {
	    cd.count_out++;
	    if (sam_write1(out, header, cd.b_unmap) < 0)
		return -1;
	    next = (cd.iter
		    ? sam_itr_next(cd.fp, cd.iter, cd.b_unmap)
		    : sam_read1(cd.fp, cd.header, cd.b_unmap)) >= 0;
	    if (next) cd.count_in++;
	} while (next);

	bam_destroy1(cd.b_unmap);
    }

    bam_plp_destroy(p_iter);
    bam_sorted_list_destroy(b_hist);

    if (cd.count_in != cd.count_out) {
	fprintf(stderr, "ERROR: lost a read?\n");
	fprintf(stderr, "Read  %"PRId64" reads\n",   cd.count_in);
	fprintf(stderr, "Wrote %"PRId64" reads\n\n", cd.count_out);
	return 1;
    }

    free(cons);
    free(cons2);

    return 0;
}

int parse_aux_list(auxhash_t *h, char *optarg) {
    if (!*h)
        *h = kh_init(aux_exists);

    while (strlen(optarg) >= 2) {
        int x = optarg[0]<<8 | optarg[1];
        int ret = 0;
        kh_put(aux_exists, *h, x, &ret);

        optarg += 2;
        if (*optarg == ',') // allow white-space too for easy `cat file`?
            optarg++;
        else if (*optarg != 0)
            break;
    }

    if (strlen(optarg) != 0) {
        fprintf(stderr, "main_samview: Error parsing option, "
                "auxiliary tags should be exactly two characters long.\n");
        return -1;
    }

    return 0;
}

void usage(FILE *fp) {
    fprintf(fp, "Realigner version %s\n\n", REALIGNER_VERSION);
    fprintf(fp, "Usage: realigner [options] in-file out-file\n");
    fprintf(fp, "\nOptions:\n"
                "-I fmt(,opt...)   Input format and format-options [auto].\n"
                "-O fmt(,opt...)   Output format and format-options [SAM].\n");
    fprintf(fp, "-u                Write uncompressed output\n");
    fprintf(fp, "-v                Increase verbosity\n");
    fprintf(fp, "-R filename       Reference FASTA file\n");
    fprintf(fp, "-r region         Sub region of input file in chr:start-end format. [all]\n");
    fprintf(fp, "-m int            Realignment margin surrounding suspect regions\n");
    fprintf(fp, "-c int            Consensus generation margin surrounding suspect regions\n");
    fprintf(fp, "-P float          Suspect if depth locally >= [%.1f] times deeper than expected\n", OVER_DEPTH);
    fprintf(fp, "-C float          Suspect if >= [%.2f] reads have soft-clipping\n", CLIP_PERC);
    fprintf(fp, "-M float          Suspect if >= [%.2f] reads have low mapping quality\n", LOW_MQUAL_PERC);
    fprintf(fp, "-Z float          Suspect if >= [%.2f] indel sizes do not fit bi-modal dist.\n", INS_LEN_PERC);
    fprintf(fp, "-V float          Suspect if <  [%.2f] reads span indel\n", INDEL_OVERLAP_PERC);
    fprintf(fp, "-X int            Whether to add primary (val&1) or secondary (val&2) consensus to graph\n");
    fprintf(fp, "\n");
}

int main(int argc, char **argv) {
    samFile *in, *out = NULL;
    htsFormat in_fmt = {0};
    htsFormat out_fmt = {0};
    bam_hdr_t *header;
    hts_itr_t *h_iter = NULL;
    int opt;

    cram_realigner_params params = {
	.verbose       = 0,                 // -v
	.ref           = NULL,		    // -R
	.region        = NULL,		    // -r
	.margin        = MARGIN,            // -m
	.cons_margin   = CON_MARGIN,        // -c
        .clip_perc     = CLIP_PERC,         // -C
        .low_mqual_perc= LOW_MQUAL_PERC,    // -M
        .ins_len_perc  = INS_LEN_PERC,      // -Z
        .over_depth    = OVER_DEPTH,        // -P
        .indel_ov_perc = INDEL_OVERLAP_PERC,// -V
	.clevel        = 6,                 // -u
	.cons1         = 0,                 // -X
	.cons2         = 0,                 // -X
    };

    while ((opt = getopt(argc, argv, "I:O:m:c:v:M:Z:P:V:r:R:uX:")) != -1) {
	switch (opt) {
	case 'u':
	    params.clevel = 0;
	    break;

	case 'I':
	    hts_parse_format(&in_fmt, optarg);
	    break;

	case 'O':
	    hts_parse_format(&out_fmt, optarg);
	    break;

	case 'r':
	    params.region = optarg;
	    break;

	case 'R':
	    params.ref = optarg;
	    break;

	case 'v':
	    params.verbose++;
	    break;

	case 'm':
	    params.margin = atoi(optarg);
	    break;
	case 'c':
	    params.cons_margin = atoi(optarg);
	    break;

        case 'C':
            params.clip_perc = atof(optarg);
            break;

        case 'M':
            params.low_mqual_perc = atof(optarg);
            break;

        case 'Z':
            params.ins_len_perc = atof(optarg);
            break;

        case 'P':
            params.over_depth = atof(optarg);
            break;

        case 'V':
            params.indel_ov_perc = atof(optarg);
            break;

	case 'X':
	    params.cons1 = atoi(optarg) & 1;
	    params.cons2 = atoi(optarg) & 2;
	    break;

	default: /* ? */
	    usage(stderr);
	    return 1;
	}
    }

    char *fnin = NULL;
    if ( optind>=argc ) {
        if ( !isatty(STDIN_FILENO) ) fnin = "-";  // reading from stdin
        else { usage(stdout); return 1; }
    }
    else fnin = argv[optind++];

    if (!(in = sam_open_format(fnin, "r", &in_fmt))) {
	perror(argv[optind]);
	return 1;
    }

    char mode[6] = "w";
    char *fnout = optind < argc ? argv[optind++] : "-";
    sam_open_mode(mode+1, fnout, NULL);
    if (params.clevel == 0)
	strcat(mode, "0");

    if (!(out = sam_open_format(fnout, mode, &out_fmt))) {
	perror("(stdout)");
	return 1;
    }

    if (!(header = sam_hdr_read(in))) {
	fprintf(stderr, "Failed to read file header\n");
	return 1;
    }
    if (out && sam_hdr_write(out, header) != 0) {
	fprintf(stderr, "Failed to write file header\n");
	return 1;
    }

    if (params.region) {
        hts_idx_t *idx = sam_index_load(in, fnin);
    	h_iter = idx ? sam_itr_querys(idx, header, params.region) : NULL;
    	if (!h_iter || !idx) {
    	    fprintf(stderr, "Failed to load index and/or parse iterator.\n");
    	    return 1;
    	}
    	hts_idx_destroy(idx);
    }

    if (transcode(&params, in, out, header, h_iter) != 0) {
	fprintf(stderr, "Error while reducing file\n");
	return 1;
    }

    bam_hdr_destroy(header);
    if (h_iter) hts_itr_destroy(h_iter);

    if (sam_close(in) != 0) {
	fprintf(stderr, "Error while closing input fd\n");
	return 1;
    }

    if (out && sam_close(out) != 0) {
	fprintf(stderr, "Error while closing output fd\n");
	return 1;
    }

    return 0;
}

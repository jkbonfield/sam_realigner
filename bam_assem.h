#ifndef _BAM_ASSEM_H_
#define _BAM_ASSEM_H_

#include <htslib/sam.h>

extern int bam_realign(bam_hdr_t *hdr, bam1_t **bams, int nbams, int *new_pos, char *ref, int ref_len, int ref_pos, char *cons1, char *cons2, int len, int max_snp, int window, int min_mqual);

#endif

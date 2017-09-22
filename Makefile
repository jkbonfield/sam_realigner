ALL=assem_bam_k3 assem_bam_k4 assem_bam_k13 assem_bam_k18 bam_problem_regions

CC=cc
CFLAGS=-g
INCLUDES=-I$(HOME)/work/samtools_master/htslib
DEFINES=
LDFLAGS=
#LIBS=-lm -L$(HOME)/work/samtools_master/htslib -lhts -lz -lbz2 -llzma -pthread
LIBS=-lm $(HOME)/work/samtools_master/htslib/libhts.a -lz -lbz2 -llzma -pthread

all: $(ALL)

.c.o:
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDES) -c $< -o $@

# OBJS1=assem.o hash_table.o pooled_alloc.o string_alloc.o align_ss.o align_sv.o
# assem: $(OBJS1)
# 	$(CC) -o $@ $(OBJS) $(LDFLAGS) $(LIBS)

OBJS=hash_table.o pooled_alloc.o string_alloc.o align_ss.o align_sv.o str_finder.o
assem_bam_k3: $(OBJS) assem_bam3_3.o
	$(CC) -o $@ $(OBJS) assem_bam3_3.o $(LDFLAGS) $(LIBS)

assem_bam_k4: $(OBJS) assem_bam3_4.o
	$(CC) -o $@ $(OBJS) assem_bam3_4.o $(LDFLAGS) $(LIBS)

assem_bam_k13: $(OBJS) assem_bam3_13.o
	$(CC) -o $@ $(OBJS) assem_bam3_13.o $(LDFLAGS) $(LIBS)

assem_bam_k18: $(OBJS) assem_bam3_18.o
	$(CC) -o $@ $(OBJS) assem_bam3_18.o $(LDFLAGS) $(LIBS)

assem_bam3_3.o: bam_assem.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDES) -DKMER=3 -DTEST_MAIN -c $< -o $@

assem_bam3_4.o: bam_assem.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDES) -DKMER=4 -DTEST_MAIN -c $< -o $@

assem_bam3_13.o: bam_assem.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDES) -DKMER=13 -DTEST_MAIN -c $< -o $@

assem_bam3_18.o: bam_assem.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDES) -DKMER=18 -DTEST_MAIN -c $< -o $@

clean:
	-rm assem_bam_k[0-9]* *.o

TEST_KMERS=test_kmer3 test_kmer4 test_kmer13 test_kmer18
.PHONY: test test3 test4 test13 test18 $(TEST_KMERS)
test check: test3 test4 test13 test18


test3: KMER=3
test3: assem_bam_k3 test_kmer3

test4: KMER=4
test4: assem_bam_k4 test_kmer4

test13: KMER=13
test13: assem_bam_k13 test_kmer13

test18: KMER=18
test18: assem_bam_k18 test_kmer18

$(TEST_KMERS):
	for d in eg$(KMER)?.sam; do \
	    echo "Testing $$d"; \
	    r=`echo $$d | sed 's/sam/ref/'`; \
	    if [ -e $$r ]; \
	    then ./assem_bam_k$(KMER) $$d $$r 2>/dev/null | sed -n '/^@HD/,$$p' > _$(KMER);\
	    else ./assem_bam_k$(KMER) $$d     2>/dev/null | sed -n '/^@HD/,$$p' > _$(KMER);\
	    fi; \
	    diff _$(KMER) `echo $$d | sed 's/sam/out/'` || exit 1; \
	    rm _$(KMER); \
	done

UPDATE_KMERS=update_test_kmer3 update_test_kmer4 update_test_kmer13 update_test_kmer18
.PHONY: update_tests update_test3 update_test4 update_test13 update_test18 $(UPDATE_KMERS)
update_tests: update_test3 update_test4 update_test13 update_test18

update_test3: KMER=3
update_test3: assem_bam_k3 update_test_kmer3

update_test4: KMER=4
update_test4: assem_bam_k4 update_test_kmer4

update_test13: KMER=13
update_test13: assem_bam_k13 update_test_kmer13

update_test18: KMER=18
update_test18: assem_bam_k18 update_test_kmer18

$(UPDATE_KMERS):
	for d in eg$(KMER)?.sam; do \
	    r=`echo $$d | sed 's/sam/ref/'`; \
	    o=`echo $$d | sed 's/sam/out/'`; \
	    if [ -e $$r ]; \
	    then ./assem_bam_k$(KMER) $$d $$r 2>/dev/null | sed -n '/^@HD/,$$p' > $$o; \
	    else ./assem_bam_k$(KMER) $$d     2>/dev/null | sed -n '/^@HD/,$$p' > $$o; \
	    fi; \
	done

# A test program for finding regions to realign.
# Could also be used to replace the GATK component?
bam_problem_regions: $(OBJS) bam_problem_regions.o bam_assem.o
	$(CC) -o $@ $(OBJS) bam_problem_regions.o bam_assem.o $(LDFLAGS) $(LIBS)

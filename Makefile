ALL=assem_bam

CC=cc
CFLAGS=-g
INCLUDES=-I$(HOME)/work/samtools_master/htslib
KMER=4
DEFINES=-DKMER=$(KMER)
LDFLAGS=
LIBS=-lm -L$(HOME)/work/samtools_master/htslib -lhts -lz -pthread

.c.o:
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDES) -c $< -o $@

# OBJS1=assem.o hash_table.o pooled_alloc.o string_alloc.o align_ss.o align_sv.o
# assem: $(OBJS1)
# 	$(CC) -o $@ $(OBJS) $(LDFLAGS) $(LIBS)

#OBJS2=assem_bam2.o hash_table.o pooled_alloc.o string_alloc.o align_ss.o align_sv.o
OBJS3=assem_bam3.o hash_table.o pooled_alloc.o string_alloc.o align_ss.o align_sv.o
assem_bam: $(OBJS3)
	$(CC) -o $@ $(OBJS3) $(LDFLAGS) $(LIBS)

clean:
	-rm assem_bam *.o
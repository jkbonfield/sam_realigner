Bam_realign
-----------

Note: bam_realign is still very much a WORK IN PROGRESS.  It is also
poorly named, as it can read SAM, BAM and CRAM due to being built on
top of htslib.  The plan is to see whether it is possible to improve
bcftools calling by realigning up-front.

The tool identifies regions of a BAM file that contain potentially
suspect alignments where it thinks a multiple sequence alignment would
yield a better result than discrete pairwise alignments (as produced
by bwa, bowtie, etc).  The main indicators of these are the presence
of heterozygous indels and concordant soft clipping.

Foreach region of suspect alignments, a local re-assembly is made of
the sequences in that region and then their consensus is realigned
back to the reference (if specified) to produce new CIGAR strings.

Some basic performance figures, on a CHM1 + CHM13 mix:
https://www.biorxiv.org/content/early/2017/11/22/223297
https://github.com/lh3/CHM-eval

NOTE: you will need to renormalise the output afterward. Eg with
    bcftools norm -m +both -f $HREF a.vcf > a.norm.vcf

Whole genome: CHM1_CHM13_2.bam, evaluation via compare_vcf.sh

Original

    Bcftools     All   >=Q30   / Filtered       Freebayes    All   >=Q30
    SNP   TP 3456155 / 3442048 / 3440987        SNP   TP 3335893 / 3319363
    SNP   FP   99681 /   73723 /   52589        SNP   FP  195468 /   48621
    SNP   FN   71240 /   85347 /   86408        SNP   FN  130780 /  147310

    InDel TP  439484 /  417136 /  436037        InDel TP  402276 /  392699
    InDel FP  181793 /  163720 /   20708        InDel FP   11194 /    6984
    InDel FN  116849 /  139197 /  120296        InDel FN  148600 /  158177

Realigned

    Bcftools     All   >=Q30   / Filtered       Freebayes    All   >=Q30
    SNP   TP 3435107 / 3423440 / 3422671        SNP   TP 3317394 / 3302687
    SNP   FP   72830 /   58245 /   44217        SNP   FP  112879 /   39771
    SNP   FN   91307 /  102974 /  103743        SNP   FN  150155 /  164862

    InDel TP  435745 /  423432 /  434033        InDel TP  405887 /  402291
    InDel FP   26572 /   23722 /   20462        InDel FP   13881 /   11627
    InDel FN  121122 /  133435 /  122834        InDel FN  147457 /  151053

All is all calls, Q>=30 is calls with QUAL field >= 30, and Filtered
is a set of custom filters for bcftools and this particular data set.

SNP:   TYPE='snp' && QUAL >= $qual && DP<90
Indel: TYPE='indel' && IDV >= 3 && IMF >= 0.03

The simple quality filtering on bcftools shows a slight shift from
favouring accuracy over recall, with a reduction in false positives
(FP) and an increase in false negatives (FN).  For SNPs this is fairly
balanced, but with Indels the FP reduction is huge, due to
realignment.  With the more aggressive filtering the realignment has
minimal difference.  One could argue that this is a problem of
bcftools setting of QUAL, not taking into acount IDV and IMF
parameters, and I wish I'd spotted this issue much earlier on as I'd
convinced myself it was doing a great job on fixing bcftools
overcalling!

For Freebayes, the changes to FP and FN are smaller, and perhaps not
beneficial.  It is expected that the indel FN/FP rate would not change
much as Freebayes already contains a (better?) built-in realigner.

The above was performed by a script using bcftools isec and counting.
It is very strict requiring both identical positions and variant, but
sometimes incorrectly labels FP/FN due to the multiple ways in which
the same variant can be described.  A more liberal approach is taken
by Heng Li's CHM-eval package, which more aggressively filters before
evaluation and has a more forgiving definition of equality which
requires SNPs at the same point but not necessarily the same variant
and Indels within 10bp (again not necessarily the same indel).

CHM-eval method, original

                                Bcftools        FreeBayes
    distEval    SNP     N+P     2719649102      2719649102
    distEval    SNP     TP      3489836         3499725   
    distEval    SNP     FN      41008           31119     
    distEval    SNP     TPc     3482311         3545352   
    distEval    SNP     FP      74585           210308    
    distEval    SNP     %FNR    1.16            0.88      
    distEval    SNP     %FDR    2.10            5.60      
    distEval    SNP     FPpM    27.424          77.329    
    distEval    INDEL   N+P     2719649102      2719649102
    distEval    INDEL   TP      348707          325374    
    distEval    INDEL   FN      28875           52208     
    distEval    INDEL   TPc     437500          292975    
    distEval    INDEL   FP      62120           4821      
    distEval    INDEL   %FNR    7.65            13.83     
    distEval    INDEL   %FDR    12.43           1.62      
    distEval    INDEL   FPpM    22.841          1.773     

CHM-eval method, realigned
         
    norm -m +both               Bcftools        Freebayes
    distEval    SNP     N+P     2719649102      2719649102
    distEval    SNP     TP      3473805         3486321   
    distEval    SNP     FN      57039           44523     
    distEval    SNP     TPc     3456888         3571504   
    distEval    SNP     FP      52072           166964    
    distEval    SNP     %FNR    1.62            1.26
    distEval    SNP     %FDR    1.48            4.47
    distEval    SNP     FPpM    19.147          61.392    
    distEval    INDEL   N+P     2719649102      2719649102
    distEval    INDEL   TP      339475          328953    
    distEval    INDEL   FN      38107           48629     
    distEval    INDEL   TPc     329481          322373    
    distEval    INDEL   FP      5995            9580      
    distEval    INDEL   %FNR    10.09           12.88
    distEval    INDEL   %FDR    1.79            2.89
    distEval    INDEL   FPpM    2.204           3.523     

Here we are trading an increase in false negatives for a decrease in
false positives.  With bcftools indel the FP decrease is huge, but
this comes down to the filtering not being ideal.

The main cause for false negative growth is the realigner sometimes
makes reads unmapped when it fails to find a new alignment for them.
Generally this is beneficial as such poorly aligned reads are more
likely to add to false positives than reducing false negatives, but
not universally so.  Generally it is the very rare reads which
erroneously (*assuming* standard allele frequences and a single
individual) confirm an indel which have become unmapped, showing this
is essentially equivalent to adding the IMF and IDV bcftools filtering
options.  Samtools flagstat respectively shows:

    Original   1026664938 + 0 mapped (99.19% : N/A)
    Realigned  1011455845 + 0 mapped (97.72% : N/A)

In summary, it does some good things, but needs more work to overcome
the loss of around 2% of sequence alignments.

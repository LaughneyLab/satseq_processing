#!/bin/bash

# This is for 4T1 and EO771.lmb, samples, other satseq 10x samples were processed similarly

export PATH=/workdir/mt269/cellranger-8.0.1/bin:$PATH


# Point to 10x reference sequence
refpath=/workdir/refdata-gex-mm10-2020-A

# Point to description of barcodes for each variant
# Note: for the paper we used the v2 file, but it had some minor errors
# that were corrected downstream. The errors are fixed in the v3 file,
# we recommend using that going forward
# featureFile=../data/twist_mSTING_full_feature_ref_v2.csv
featureFile=../data/twist_mSTING_full_feature_ref_v3.csv

## in this case we had directory structure like:
## 4T1/transcripts/*.fastq.gz:        fastqs from RNA for 8 4T1 samples
## 4T1/feature_barcodes/*.fastq.gz:   fastqs from satseq barcode reads for 4T1 samples
## E0771/transcripts/*.fastq.gz:      fastqs from RNA for 8 EO771.lmb samples
## EO771/feature_barcodes/*.fastq.gz: fastqs from satseq barcode reads for the EO771.lmb samples

# The fastq files in each directory are named such that they all start with the same prefix,
# except a number in the third field (separated by _) gives the sample number
# So, sorting by filename in each directory works to get the sample ID.
# (The important thing here being that you run cellranger using the matching set of
# feature_reads and transcript reads, from the same sample. Though, if you get it wrong, 
# you will get very poor feature mapping, which will be a red flag).

for d in 4T1 EO771; do
    cd $d
    numsamples=`ls feature_barcodes/ | awk -F '_' -v OFS="_" '{NF=3; print $0}' | sort | uniq | wc -l`
    for i in `seq 1 $numsamples`; do
        feat=`ls feature_barcodes/ | awk -F '_' -v OFS="_" '{NF=3; print $0}' | sort | uniq | awk -v i=$i 'NR==i'`
        tran=`ls transcripts/ | awk -F '_' -v OFS="_" '{NF=3; print $0}' | sort | uniq | awk -v i=$i 'NR==i'`
        sname=${d}_$i
        echo "$feat $tran $sname"
        echo -e "fastqs,sample,library_type
`pwd`,$tran,Gene Expression
`pwd`,$feat,Custom" > libraries_${sname}.csv
         cellranger count --id=$sname \
               --libraries=libraries_${sname}.csv \
               --create-bam=true \
               --localcores=20 \
               --localmem=128 \
               --transcriptome=$refpath \
               --feature-ref=$featureFile \
               --check-library-compatibility=false
    done
    cd ..
done


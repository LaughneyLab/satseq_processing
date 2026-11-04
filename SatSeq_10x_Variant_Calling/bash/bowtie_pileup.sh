#!/bin/bash

## Get the feature barcode associated with each cell. Uses cellranger bam output and a barcode reference

# adds two files to the cellranger 'outs' dir:
# -- barcode_align.bam: an alignment of R2 reads to a barcode template. Maps all reads that are assigned a
#        cell barcode by cellranger
#    Note: I have tried only mapping reads that have flag=4 (i.e., removing ones marked as 'PCR duplicate').
#        But, we seem to lose a significant fraction of mappings this way
# -- barcode_consensus_simple.fa.gz: A fasta file where sequence name is cell name, and sequence is the
#        'consensus' R2 sequence after all reads are mapped to barcode template.
#        The actual barcode should be found in this sequence in positions 20:56 (zero-based).
#    Note: I have tried calling consensus without --mode=simple, it produced fewer overall mappings. But
#        with the output bam file, samtools consensus can be re-run easily with different parameters, if desired.


bamfile=$1
if [[ -z $bamfile ]]; then
    bamfile=/workdir/satseq/Laughney-MD-16637_2024_06_11/4T1/4T1_1_with_bam/outs/possorted_genome_bam.bam
fi

outdir=`dirname $bamfile`
tmpdir=`mktemp -d /tmp/tmp_XXXXXX`
echo bamfile=$bamfile
echo tmpdir=$tmpdir

barcoderef=bref
if [[ ! -e bref.1.bt2 ]]; then
    barcoderef=$tmpdir/barcode_ref
     echo "creating barcode template $barcoderef.fa"
     ( echo ">barcode_template_CGCT"
       echo "AAGAACCTTGCGGGTAAATCNNNNNNNNNNNNNNNNCGCTNNNNNNNNNNNNNNNNCGCTTTAAGGCCGGTCCTAGCAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
       echo ">barcode_template_AGAA"
       echo "AAGAACCTTGCGGGTAAATCNNNNNNNNNNNNNNNNAGAANNNNNNNNNNNNNNNNCGCTTTAAGGCCGGTCCTAGCAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN") > $barcoderef.fa
    echo "indexing barcode_ref.fa"
    bowtie2-build $barcoderef.fa $barcoderef &>/dev/null
fi


echo "calling bowtie2"
time samtools view $bamfile '*' | awk '{for (i=12; i<= NF; i++) {
     if ($i ~ /^CB:Z:/) {
       cell=substr($i,6);
       print "@"$1":CB:"cell"\n"$10"\n+\n"$11
      }
   }}' | bowtie2 -p 30 --very-sensitive --end-to-end -x $barcoderef -U - -S $tmpdir/out.sam --np 0 --n-ceil 36


# after aligning to refseq, change refseq name to cell name so that pileup/consensus groups reads by cell
echo "making bam file"
( head -n 1 $tmpdir/out.sam
  zcat $outdir/raw_feature_bc_matrix/barcodes.tsv.gz | awk -v OFS="\t" '{print "@SQ","SN:"$0,"LN:79"}'
awk -v OFS="\t" 'NR > 3 && ($3 ~ /^barcode/)  {rseq=$3;cell=split($1,a,":"); cell=a[9]; $3=cell; print $0"\bt:Z:"rseq}' $tmpdir/out.sam ) > $tmpdir/out2.sam
samtools sort $tmpdir/out2.sam > $outdir/barcode_align_bothref.bam
samtools index  $outdir/barcode_align_bothref.bam

echo "Calling consensus"
samtools consensus $outdir/barcode_align_bothref.bam -l 1000 --mode simple -d 1 | bgzip > $outdir/barcode_consensus_simple_bothref.fa.gz
samtools coverage $outdir/barcode_align_bothref.bam | gzip > $outdir/barcode_coverage_bothref.txt.gz


rm -rf $tmpdir

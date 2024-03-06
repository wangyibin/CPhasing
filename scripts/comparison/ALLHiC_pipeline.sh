#!/bin/bash
#PBS -j oe
#PBS -q share
#PBS -V
#PBS -l nodes=1:ppn=10

if [[ ! -z $PBS_O_WORKDIR ]]; then
    cd $PBS_O_WORKDIR
fi

name=ploidy-
ploidy=
n50=REPLACE
fasta=${name}.$n50.fasta
bam=${name}.bam
cds=/public/home/agis_wangyb/project/0.ALLPhase/0.pore_c_dev/simulate/AT/TAIR10.cds
bed=/public/home/agis_wangyb/project/0.ALLPhase/0.pore_c_dev/simulate/AT/TAIR10.bed
nchrom=$((ploidy*5))
threads=10
restrict_site=GATC

source activate /public/home/agis_wangyb/miniconda3/envs/popCNV
ln -s ../$fasta .
ln -s ../$bam .


gmap_build -D . -db DB  $fasta 
gmap -d DB -D . -f 2 -n $ploidy -t $threads $cds > gmap.gff3
gmap2AlleleTableBED.pl $bed 

~/software/ALLHiC_components/Prune/ALLHiC_prune  -i Allele.ctg.table -b $bam
allhic extract prunning.bam $fasta --RE $restrict_site
allhic partition prunning.counts_${restrict_site}.txt prunning.pairs.txt $nchrom --minREs 25

for txt in *.counts_${restrict_site}.${nchrom}g*.txt; do 
    echo "allhic optimize $txt prunning.clm" 
done | parallel -j 10

ALLHiC_build $fasta
rm -rf groups.asm.fasta

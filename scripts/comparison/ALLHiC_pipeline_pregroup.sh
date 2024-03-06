
ploidy=2
n50=REPLACE
name=ploidy-
fasta=${name}.${n50}.fasta
bam=${name}.bam
nchrom=$((ploidy*5))
threads=10
cds=/public/home/agis_wangyb/project/0.ALLPhase/0.pore_c_dev/simulate/AT/TAIR10.cds
bed=/public/home/agis_wangyb/project/0.ALLPhase/0.pore_c_dev/simulate/AT/TAIR10.bed
restrict_site=GATC

ln -s ../$fasta .
ln -s ../$bam .

gmap_build -D . -db DB  $fasta 
gmap -d DB -D . -f 2 -n $ploidy -t $threads $cds > gmap.gff3
gmap2AlleleTableBED.pl $bed 

rm -rf DB 
rm -rf gmap.gff3

source activate /public/home/agis_wangyb/software/anaconda3/envs/popCNV
samtools sort -@ $threads $bam -o ${bam%%.bam}.sorted.bam
partition_gmap.py -r $fasta -g Allele.ctg.table -b ${bam%%.bam}.sorted.bam -t $threads 

cd wrk_dir

for i in {1..5}; do 
    echo "grep "^Chr$i" ../Allele.ctg.table > Chr$i/Allele.ctg.table"
done | parallel -j $threads

for i in {1..5}; do 
    echo "cd Chr$i; ~/software/ALLHiC_components/Prune/ALLHiC_prune -i Allele.ctg.table -b Chr${i}.bam; cd .."
done | parallel -j $threads 

for i in {1..5}; do 
    echo "cd Chr$i; allhic extract prunning.bam Chr$i.fa --RE $restrict_site; cd .."
done | parallel -j $threads

for i in {1..5}; do 
    echo "cd Chr$i; allhic partition prunning.counts_${restrict_site}.txt prunning.pairs.txt $ploidy --minREs 25; cd .."
done | parallel -j $threads

for i in {1..5}; do 
    cd Chr$i; 
    for txt in *.counts_${restrict_site}.${ploidy}g*.txt; do 
        echo "allhic optimize $txt prunning.clm"; done | parallel -j $threads;
    
    for tour in *.tour; do 
        new_tour=`echo $tour | cut -d "." -f 3,4 | sed "s/^$ploidy/Chr${i}g/"`
        mv $tour $new_tour ;
    done 
    
    cd ..

done

find . -name "*.tour" | xargs -I {} mv {} .

ALLHiC_build ../$fasta
rm groups.asm.fasta
cd ..

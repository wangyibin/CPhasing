
# Dependency

*   seqkit
*   python
*   minimap2 (version grater than 2.24)
*   python packages
*
    *   joblib
*
    *   scipy
*
    *   numpy
*
    *   matplotlib.pyplot

# Modules

## sliceReads

### example

```bash
$ python split_reads_v2.py -h
usage: split_reads_v2.py [-h] [-i FASTQ] [-w WINDOW] [-p NPART] [-t THREADS]

This is the script for split fastq file into window reads. For accelerating this process, we split fastq into N parts at first, then slice reads with window.

optional arguments:
  -h, --help            show this help message and exit
  -i FASTQ, --fastq FASTQ
                        <filepath> fastq file.
  -w WINDOW, --window WINDOW
                        <int> window size, default is 5000.
  -p NPART, --npart NPART
                        <int> number of part, default is 10
  -t THREADS, --threads THREADS
                        <int> number of threads, default is 10

```

### output result

*   split\_fastq

## correctPAF

The pipeline including two steps. First, the 5000bp window reads are compared to the draft assembly using minimap2.Because the short window reads are prone to occur alignment errors, we further correct errors based on the longest increasing subsequence (LIS) method.

### example

```bash
$ python correctPAF.py -h
usage: correctPAF.py [-h] [-f FASTA] [-i FASTQ] [-a MINAS] [-m MINMAPQ]
                          [-n NHAP] [-w WINDOW] [-t THREADS] [-o OUTPUT]

This is the script for maaping window UL-ONT reads into fasta, and identifying
chimeric contig through split alignment.

optional arguments:
  -h, --help            show this help message and exit
  -f FASTA, --fasta FASTA
                        <filepath> draft contig assembly fasta file.
  -i FASTQ, --fastq FASTQ
                        <filepath> window fastq file.
  -a MINAS, --minAS MINAS
                        <int> minimun Alignment Score, default is 2000.
  -m MINMAPQ, --minMapq MINMAPQ
                        <int> minimun Alignment Score, default is 10.
  -n NHAP, --nhap NHAP  <int> maximum supplement aligment records.
  -w WINDOW, --window WINDOW
                        <int> window size, default is 5000.
  -t THREADS, --threads THREADS
                        <int> number of threads used when running minimap2,
                        default is 10
  -o OUTPUT, --output OUTPUT
                        <str> output file prefix, default is output

```

### output

*   LIS
*
    *   \*.LIS.gtf
        All corrcted alignment result
*
    *   \*.mapq.LIS.gtf
        Corrected alignment result with mapping quality greater than thredshold.


```
    # It similar to '.gtf' format. including 9 colmuns.
    # seq_id	seq_length	type	start	end	Alignment-score	strand	alignment-type	attributes
    ctg12	400000	LIS/alignment	25132483	25132543	8000	+	P/S/E/M	reads_id "reads1"; window_id "reads1_1";

    ## alignment-type:
    ### P: Primary alignment
    ### S: Secondary alignment
    ### E: Alignments in contig Ends
    ### M: Alignments with mapping quality >= thredshold.
```
*   corrected paf
- - *.corrected.paf
All corrected alignment results. PAF format.
- - *.mapq.corrected.paf
Corrected alignment result with mapping quality greater than thredshold. PAF format.

```
    # PAF format, Cotaining 13 colmuns.
    # 1. window_id
    # 2. window_length
    # 3. alignment start site in window read
    # 4. alignment end site in window read
    # 5. string
    # 6. query name
    # 7. query length
    # 8. alignment start site in query
    # 9. alignment end site in query
    # 10. alignment length of window read
    # 11. alignment length of query
    # 12. mapping quality
    # 13. alignment score
```

## findChimeric

### STEP1. Identyfing splited alignment

#### Example

```
$ python find_splitAlignment.py -h
usage: find_splitAlignment.py [-h] [-l LIS] [-w WINDOW] [-m MINSA] [-e EDGE]
                              [-o OUTPUT]

This is the script for maaping window UL-ONT reads into fasta, and identifying
chimeric contig through split alignment.

optional arguments:
  -h, --help            show this help message and exit
  -l LIS, --lis LIS     <filepath> draft contig assembly fasta file.
  -w WINDOW, --window WINDOW
                        <int> window size, default is 5000.
  -m MINSA, --minSA MINSA
                        <int> Number of minimum split alignment in windows,
                        default is 3.
  -e EDGE, --edge EDGE  <int> Number of minimum split alignment in windows,
                        default is 20000.
  -o OUTPUT, --output OUTPUT
                        <str> output file prefix, default is output

```

#### output

*   \*.splitAlign.txt


```
    # 4 colmun file
    # contig_Name    start end number_Of_Split_Alignment split_alignment_reads_id 
    Chr4_hap1       74000   90000   4       Chr4-hap2_69239;aligned_3589_F_25_61855_34,Chr4-hap2_46259;
```
*   \*.mergedSplitAlign.txt
    Merging adjacent intervals


```
    # 4 colmuns file
    # Merging adjacent intervals
    # contig_Name    start end number_Of_Split_Alignment split_alignment_reads_id 
    Chr4_hap1       74000   90000   4       Chr4-hap2_69239;aligned_3589_F_25_61855_34,Chr4-hap2_46259;
```
### STEP2. Calculating ONT reads depth

This scripts calculating sequencing depth in a 5000bp window with 1000bp step. (default window size is 5000bp, step length is 1/5 window size)

#### Example
```
    $ python paf2depth.py -h
    usage: paf2depth.py [-h] [-p PAF] [-f FASTA] [-w WIN] [-o OUTPUT]

    This is the script for filter genome region.

    optional arguments:
      -h, --help            show this help message and exit
      -p PAF, --paf PAF     <filepath> the corrected UL-ONT reads/Hifi reads
                            mapping result, must be .paf format.
      -f FASTA, --fasta FASTA
                            <filepath> break-point of chimeric contigs.
      -w WIN, --win WIN     <int> window size when calculating depth.
      -o OUTPUT, --output OUTPUT
                            <str> output file prefix, default is output
```
#### output

*   \*.depth


```
    # 4 colmuns file
    # contig_name start end depth
    Chr1_hap1       0       5000    0.000
    Chr1_hap1       1000    6000    0.000
```
### STEP3. Normalization number of split alignment with sequencing depth.

> propotion of split alignment reads = (Number of split alignment in window) / (Sequencing depth in window)

When propotion of split alignment reads higher than threshold, we think that a chimeric error has occurred in this interval.

#### Example
```
    $ python normMergeBP.py -h
    usage: normMergeBP.py [-h] [-b BREAKPOINT] [-d DEPTH] [-min MINDEPTH]
                          [-c CUTOFF] [-w WIN] [-o OUTPUT]

    This is the script for filter genome region.

    optional arguments:
      -h, --help            show this help message and exit
      -b BREAKPOINT, --breakpoint BREAKPOINT
                            <filepath> the merged breakpoint file.
      -d DEPTH, --depth DEPTH
                            <filepath> depth file.
      -min MINDEPTH, --minDepth MINDEPTH
                            <int> minimum depth of windows.
      -c CUTOFF, --cutoff CUTOFF
                            <float> cutoff of identifying chimeric contigs, which
                            equal (count of splited alignment reads)/(avarage
                            depth in chimeric region).
      -w WIN, --win WIN     <int> window size when calculating depth.
      -o OUTPUT, --output OUTPUT
                            <str> output file prefix, default is output
```
#### output

*   \*.norm.result


```
    # 7 colmun file
    # contig_name   start   end type    number_of_split_alignment depth propotion
    Chr4_hap1       1050000 1060000 chimeric        3       5.31    0.56
```
*   \*.breakPos.txt


```
    # 2 colmun file
    # contig_name   break_point_of_chimeric
    Chr4_hap1       82000,143000,205000
```
### STEP4. correct chimeric contig

#### Example
```
    $ python correctFasta.py -h
    usage: correctFasta.py [-h] [-f FASTA] [-bp BREAKPOS] [-o OUTPUT]

    This is the script for chimeric contig correction.

    optional arguments:
      -h, --help            show this help message and exit
      -f FASTA, --fasta FASTA
                            <filepath> the raw fasta file.
      -bp BREAKPOS, --breakpos BREAKPOS
                            <filepath> chimeric position of contigs.
      -o OUTPUT, --output OUTPUT
                            <str> output file prefix, default is output
```
## highConfidence

### STEP1. filter region with low or too-high sequencing depth

#### Example

```bash
## need parameter: -d -f -w
## optional parameter: -w -m -b
$ python bed2depth_2.py -h
usage: bed2depth_2.py [-h] [-d DEPTHFILE] [-b BREAKPOINT] [-f FASTA] [-w WIN]
                      [-M MAX] [-o OUTPUT]

This is the script for filter genome region.

optional arguments:
  -h, --help            show this help message and exit
  -d DEPTHFILE, --depthFile DEPTHFILE
                        <filepath> UL-ONT reads/Hifi depth file, 4colum, Chr
                        start end depth.
  -b BREAKPOINT, --breakpoint BREAKPOINT
                        <filepath> break-point of chimeric contigs, optional.
  -f FASTA, --fasta FASTA
                        <filepath> fasta file.
  -w WIN, --win WIN     <int> window size when calculating depth, default=5000.
  -M MAX, --max MAX     <int> maximum depth, default=100.
  -o OUTPUT, --output OUTPUT
                        <str> output file prefix, default is output
```

#### output
```
    # 4colmun file
    # contig_name start end depth
    Chr1_hap1       0       5000    0.000
```
### STEP2. Identifing high confidence region.

#### Example
```
    usage: High_confidence_region.py [-h] [-d DEPTHFILE] [-sa SPLIAALIGN] [-l LIS] [-m MINCOUNT] [-M MINMAPQ] [-o OUTPUT]

    This is the script for filter genome region.

    optional arguments:
      -h, --help            show this help message and exit
      -d DEPTHFILE, --depthFile DEPTHFILE
                            <filepath> UL-ONT reads/Hifi depth file, 4colum, Chr start end depth.
      -sa SPLIAALIGN, --spliaAlign SPLIAALIGN
                            <filepath> splitAlign file, 4 col: <chr> <start> <end> <count of split-align> <split-aligned reads>.
      -l LIS, --lis LIS     <filepath> Corrected alignment result with mapping quality greater than thredshold.
      -m MINCOUNT, --minCount MINCOUNT
                            <int> minimum count of windows in LIS.
      -M MINMAPQ, --minMapq MINMAPQ
                            <float> minmum ratio of (mapq10 window)/(window counts).
      -o OUTPUT, --output OUTPUT
                            <str> output file prefix, default is output
```
#### output

```bash
# 3 colmun file
# contig_name   start   end
chr4_hap1   1000    506642
```

## Useful scripts

### extract\_hicLink\_fromBam.py

Extract hic links located in high confidence region from hic mapping bam file.
```
    Usage: python extract_hicLink_fromBam.py <high_confidence_region|bed_file> <bam_file> <working_directory> <prefix_of_output> <threads>
```

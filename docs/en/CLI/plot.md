---
icon: lucide/image
tags:
    - Plot
--- 

## Examples
### Plot chromosome-level heatmap from contig-level mapping result
- Convert pairs to cool file  

=== "pairs.pqs" 
    ```shell
    cphasing-rs contigsizes contigs.fa > contigs.contigsizes
    cphasing pairs2cool sample.pairs.pqs contigs.contigsizes sample.10k.cool -bs 10k
    ```

=== "pairs.gz" 
    ```shell
    cphasing-rs contigsizes contigs.fa > contigs.contigsizes
    cphasing pairs2cool sample.pairs.gz contigs.contigsizes sample.10k.cool -bs 10k
    ```

- Adjust contig-level matrix to chromosome-level by agp file

    ```shell
    cphasing plot -a groups.agp -m sample.10k.cool -o groups.500k.wg.png   
    ```

!!! note
    This function will generate two intermedia `sample.10000.chrom.cool` and `sample.500k.chrom.cool`.


<center>
    <img src="https://pic.superbed.cc/item/6791ff5efa9f77b4dce26304.png" width=50%>
    <br>
    <font color="red">groups.500k.wg.png</font>
</center>


### Directly plot heatmap from `.cool` file
```shell
cphasing plot -m sample.100k.chrom.cool -o groups.100k.wg.png 
```
!!!note 
    `--binsize` can be specified when you want to plot a larger resolution, e.g., 500k, `plot` will automatically generate a 1m matrix from the input matrix. Still, the binsize of the output matrix should be in integer multiples than the binsize of the input matrix.

    ```shell
    cphasing plot -m sample.100k.chrom.cool -o groups.500k.wg.png -bs 500k
    ```

<center>
    <img src="https://pic.superbed.cc/item/679201c8fa9f77b4dce27b64.png" width=40% >
    <img src="https://pic.superbed.cc/item/6791ff5efa9f77b4dce26304.png" width=40% >
    <br>
    <font color="red">groups.100k.wg.png</font>
    &emsp; &emsp; &emsp; &emsp; &emsp; &emsp; &emsp; 
    <font color="red">groups.500k.wg.png</font>
    
</center>



### Gallery of different parameters

![different_parameters](https://pic.superbed.cc/item/67921867fa9f77b4dce3785c.png)

#### Add line to highlight the haplotype border 
```shell
cphasing plot chrom.500k.cool -o chrom.500k.png --whitered --add-hap-border
```
<center>
<img src="https://pic.superbed.cc/item/67b72bfc01c078e31052c7a0.png" width=50% >
</center>

#### Plot specified chromosome
![chrom_plot](https://pic.superbed.cc/item/67921640fa9f77b4dce3654f.png)
#### Plot all chromosome split in one picture
![per_chrom](https://pic.superbed.cc/item/67921812fa9f77b4dce37540.png)


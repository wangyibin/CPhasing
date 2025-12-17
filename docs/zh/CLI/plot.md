---
icon: lucide/image
tags:
    - Plot
--- 

## 示例
### 从contig水平比对结果绘制热图
- 将`pairs`文件转成`.cool`矩阵文件
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
- 根据AGP文件绘制热图
    根据AGP文件将contig水平的`.cool`矩阵转成染色体水平，并绘制热图

    ```shell
    cphasing plot -a groups.agp -m sample.10k.cool -o groups.500k.wg.png   
    ```

!!! note
    这步会产生两个文件一个是`sample.10k.chrom.cool` and `sample.500k.chrom.cool`，后续在agp文件不变的情况下，可以使用这两个文件进行热图的调整.


<center>
    <img src="https://pic.superbed.cc/item/6791ff5efa9f77b4dce26304.png" width=30%>
    <br>
    <font color="black">groups.500k.wg.png</font>
</center>


### 直接使用`.cool` 文件绘制热图
```shell
cphasing plot -m sample.100k.chrom.cool -o groups.100k.wg.png 
```
!!!note 

    用户可以指定`--binsize`参数来调整热图的分辨率，如果输入矩阵的分辨率跟设定的分辨率不一致，程序会自动生成新的分辨率绘制热图，但需要设定的分辨率和输入矩阵的分辨率呈整数倍的关系。例如，输入100k.cool，想绘制500k分辨率的热图:

    ```shell
    cphasing plot -m sample.100k.chrom.cool -o groups.500k.wg.png -bs 500k
    ```

<center>
    <img src="https://pic.superbed.cc/item/679201c8fa9f77b4dce27b64.png" width=30% >
    <img src="https://pic.superbed.cc/item/6791ff5efa9f77b4dce26304.png" width=30% >
    <br>
    <font color="black">groups.100k.wg.png</font>
    &emsp; &emsp; &emsp; &emsp; 
    <font color="black">groups.500k.wg.png</font>
    
</center>




### 不同参数绘制的热图

![different_parameters](https://pic.superbed.cc/item/67921867fa9f77b4dce3785c.png)

#### 同源染色体组加边框
```shell
cphasing plot chrom.500k.cool -o chrom.500k.png --whitered --add-hap-border
```
<center>
<img src="https://pic.superbed.cc/item/67b72bfc01c078e31052c7a0.png" width=50% >
</center>

#### 指定染色体绘制
![chrom_plot](https://pic.superbed.cc/item/67921640fa9f77b4dce3654f.png)
#### 将每条染色体分开绘制
![per_chrom](https://pic.superbed.cc/item/67921812fa9f77b4dce37540.png)


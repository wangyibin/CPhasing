---
icon: lucide/circle-question-mark
---

### 第一轮分组的结果不好:   
在我们的两轮聚类算法中，第一轮聚类依赖于同源染色体之间比对错误；如果用户输入低水平Switch error的contigs或输入高精度的Pore-c数据，h-*trans*将不足以支撑将来自同源染色体的contig聚到一起，这容易导致结果不理想。用户可以为`hyperpartition`或`pipeline`设置`-q1 0`以增加h-*trans*错误率。但是，当您在porec.gz或pairs.gz中输入大量的互作数据时，此参数可能会引发内存不足的错误。 


### 挂载上的染色体总大小远低于预估基因组大小
如果存在以下两种情况，可以通过调整 `cphasing pipeline`的模式至敏感（`--preset sensitive`）或者超敏感（`--preset very-sensitive`）
1. 输入的数据量低。2. 输入的基因组较为复杂，存在大量的纯合或者近乎纯合的区域。 需要注意的是，以上两种模式会让部分较碎的contig聚类或者排序错误。同时如果属于第二种情况，容易发生贪婪的聚类，即两条高度同源的染色体组被分到一组里面。


### 如何在组装非整倍体基因组时设置`-n`参数:  
非整倍体基因组，如现代栽培的甘蔗，包含数目不相等的同源染色体。我们建议`-n`参数可以设置为零（`-n 0:0`），让程序自动判别分组数。
此外，我们也允许用户输入一个包含两列的文件：第一列是第一轮分区的索引（1-base），第二列是每个同源染色体的染色体数目。然后指定`- n10:second.number.tsv`。在`cphasing pipeline`或`cphasing hyperpartition`中使用。

```text title="second.number.tsv"
1    13
2    12
3    12
4    11
5    10
6    12
7    12
8    10
9    12
10    12
```

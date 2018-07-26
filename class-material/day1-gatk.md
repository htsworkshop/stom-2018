---
layout: page
---

---

### III-1. Variant calling and joint genotuping using GATK HaplotypeCaller

In this part, we will follow the best practice of GATK
HaplotypeCaller. We used [GATK version
4.0.5.1](https://github.com/broadinstitute/gatk/releases/tag/4.0.5.1)
released on Jun 12, 2018.

---

#### a. Preambles for joint calling

- We will first create a directory to store output files.
> <pre>
$ cd ~/
$ mkdir --p stom2018/s3.1 
$ cd stom2018/s3.1/ </pre>

- We will use the PUR trio deeply sequenced repeatedly from two
  different centers.
> <pre>
$ ls /ws_data/topmed/crams/ </pre>

- Because they are CRAM files, it is important to set environment
  variable `REF_PATH`. 
> <pre>
$ export REF_PATH=~/md5/%2s/%2s/%s 
$ env | grep REF_PATH </pre>

- To list the sample IDs, we can combine our UNIX skill sets as
  follows:
> <pre>
$ ls /ws_data/topmed/crams | grep cram$ | cut -d . -f 1 
NWD230091
NWD231092
NWD259170
NWD315403
NWD495157
NWD684137 </pre>
Make sure that you observe the same output

- You may use `xargs` to apply a series of command to each of these
  sample with a single command. For example,
> <pre>
$ ls /ws_data/topmed/crams/ | grep cram$ | cut -d . -f 1 | xargs -I {} ls -l /ws_data/topmed/crams/{}.recab.chr20.cram 
-rw-rw-r-- 1 user1 user1 921367429  7월 25 19:59 /ws_data/topmed/crams/NWD230091.recab.chr20.cram
-rw-rw-r-- 1 user1 user1 497682425  7월 25 19:56 /ws_data/topmed/crams/NWD231092.recab.chr20.cram
-rw-rw-r-- 1 user1 user1 524656359  7월 25 19:56 /ws_data/topmed/crams/NWD259170.recab.chr20.cram
-rw-rw-r-- 1 user1 user1 561420513  7월 25 19:57 /ws_data/topmed/crams/NWD315403.recab.chr20.cram
-rw-rw-r-- 1 user1 user1 556571604  7월 25 19:57 /ws_data/topmed/crams/NWD495157.recab.chr20.cram
-rw-rw-r-- 1 user1 user1 543732572  7월 25 19:57 /ws_data/topmed/crams/NWD684137.recab.chr20.cram</pre>
  * Q. Can you explain what the command did?
  
- Another example of using `xargs`:
> <pre>
$ ls /ws_data/topmed/crams/ | grep cram$ | cut -d . -f 1 | xargs -I {} echo "My name is {}"
My name is NWD230091
My name is NWD231092
My name is NWD259170
My name is NWD315403
My name is NWD495157
My name is NWD684137 </pre>
  * Q. Can you explain how this command works?

---

#### b. Creating per-sample GVCF files
	
- GATK HaplotypeCaller best practice recommends to create GVCF files
for each sample for joint variant calling. Due to time constrains, we
will focus the 200kb region starting from 10Mb of chromosome 20.

- Using `xargs`, we can run across all six samples with a single
command.
> <pre>
$ ls /ws_data/topmed/crams/ | grep cram$ | cut -d . -f 1 | \
xargs -I {} gatk HaplotypeCaller -R /ws_data/ref/hs38DH.fa \
-I /ws_data/topmed/crams/{}.recab.chr20.cram \
-O {}.g.vcf.gz \
-ERC GVCF \
-L chr20:10000000-10200000 </pre>
Note that the command will take a while to finish.

- When everything finishes, take a look at what was generated in the
  current directory.
> <pre>
$ ls -l </pre>

- Take a look at what a GVCF look like..
> <pre>
$ zcat NWD230091.g.vcf.gz | grep -v ^## | less -S </pre>
  * Q. What is the difference between GVCF files from typical VCF files?

---

#### c. Consolidating the variant sites across multiple samples.

- The next step for joint calling is to merge the variant sites across
multiple samples. The new version of GATK HaplotypeCaller uses
_GenomicsDB_ to store collections of variant sites.

- Run the following command to merge the GVCF files to create a GenomicsDB
> <pre>$ gatk GenomicsDBImport \
-V NWD230091.g.vcf.gz \
-V NWD231092.g.vcf.gz \
-V NWD259170.g.vcf.gz \
-V NWD315403.g.vcf.gz \
-V NWD495157.g.vcf.gz \
-V NWD684137.g.vcf.gz \
--genomicsdb-workspace-path PURtrio \
--intervals chr20:10000000-10200000 </pre>

- Upon completion, take a peek at the GenomicsDB directory.
> <pre>$ ls PURtrio/ 
$ du -sh PURtrio/
</pre>
  * Do you expect the size of GenomicsDB will be large or small? How
    much, if you call whole genomes across 10,000 samples?

---

#### d. Joint genoptyping of the variant sites using GenomicsDB

- The joint genotyping can be performed using the GenomicsDB, without
requiring any other parameters, as follows:
> <pre>$ gatk GenotypeGVCFs \
-R /ws_data/ref/hs38DH.fa \
-V gendb://PURtrio \
-O PURtrio.joint.vcf.gz \
-L chr20:10000000-10200000 </pre>

- Take a look at the output VCF file.
> <pre>$ zcat PURtrio.joint.vcf.gz | grep -v ^## | less -S </pre>
Use left and right arrows to navigate the columns across samples.
  * Q. Can you understand what each column means? 
  * Q. How about INFO fields? 
  * Q. How about FORMAT fields?

- To extract GT fields only across the 6 samples, we will use some
  hacky solution...
> <pre>$ zcat PURtrio.joint.vcf.gz | grep -v ^# | cut -f 10- | perl \
  -lane 'print join("\t",substr($F[0],0,3),substr($F[1],0,3),substr($F[2],0,3),substr($F[3],0,3),substr($F[4],0,3),substr($F[5],0,3),)' \
  | sort | uniq -c | sort -n </pre>
  * Q. Considering these 6 samples are actually duplicate of a single
    trio, can you figure out which pairs of samples are identical? 
  * Q. Can you also figure out which one is offspring and which ones
    are parents?
  * Q. Does the genotypes mostly look Mendelian concordant or not?
  
---

Feel free to ask questions to your instructor(s) if you are stuck. 
Otherwise, let's 
go to next step 
[III-2 Variant calling and joint genotuping using vt and cramore](../class-material/day1-vt-cramore.html), or 
go back to [Day 1 Overview](../day1).

---

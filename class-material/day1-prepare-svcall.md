---
layout: page
---

---

### IV-1.  Preparing input files for Lumpy

In this part, we will briefly learn how to prepare for SV calling using 
Lumpy, a breakpoint detection tool for structural variations.

---

#### a. Prepare BAM files 

- Lumpy uses information of split reads and discordant read pairs, stored as separate BAM files
Recent versions of Lumpy also can handle CRAM files, but it still requires split reads and 
discordant read pairs in BAM format, so we will use BAM files from TOPMed whole genome sequenced
samples. The BAM versions of a 1-Mbp region (chr22:36,000,000-37,000,000) are generated at

> <pre>
/BiO/data/topmed/bams</pre>
- These BAMs are generated from original CRAM files using samtools:
> <pre>
samtools view -h -b NWD230091.recab.cram chr22:36000000-37000000  > NWD230091.bam
samtools index NWD230091.bam</pre>
- There are six (6) BAM files in the directory. 
 
LUMPY Express expects BWA-MEM aligned BAM files as input.
It automatically parses sample, library, and read group information using the @RG
tags in the BAM header.
Each BAM file is expected to contain exactly one sample.

The minimum input is a coordinate-sorted BAM file (-B), from which LUMPY Express
extracts splitters and discordants using SAMBLASTER before running LUMPY.
Optionally, users may supply coordinate-sorted splitter (-S) and discordant (-D)
BAM files which will bypass SAMBLASTER extraction for faster analysis.

If you have not yet aligned your data, you can generate these additional data 
with [SpeedSeq](https://github.com/cc2qe/speedseq), which
performs BWA-MEM alignment, marks duplicates and extracts split and discordant
read-pairs. For example, 
```
speedseq align -R "@RG\tID:id\tSM:sample\tLB:lib" \
    human_g1k_v37.fasta \
    sample.1.fq \
    sample.2.fq
```

- Otherwise, split-read and discordant read pair data may be generated from BAM/CRAM files 
already aligned with BWA-MEM.

```
# Extract the discordant paired-end alignments.
samtools view -b -F 1294 sample.bam > sample.discordants.unsorted.bam

# Extract the split-read alignments
samtools view -h sample.bam \
    | scripts/extractSplitReads_BwaMem -i stdin \
    | samtools view -Sb - \
    > sample.splitters.unsorted.bam
```



---

Feel free to ask questions to your instructor(s) if you are stuck. 
Otherwise, let's go to next step 
[IV-2. SV discovery using Lumpy](../class-material/day1-lumpy)
, or go back to [Day 1 Overview](../day1).

---

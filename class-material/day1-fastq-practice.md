---
layout: page
---

---

### II-2. Processing raw sequence reads

In this part, we will briefly learn how to understand the format of
FASTQ files, and how to align them using the `bwa` alignment software
tool.

---

#### a. Browse raw sequence reads

- The first task today is to read the example FASTQ files extracted
from the 1000 Genomes Project. Type the following command
> <pre>
zless ~/data/1000g/fastqs/HG00406.R1.fastq.gz </pre>
  * Q. What did appear in the screen? 
  * Q. Can you tell what is the base call and quality of the very
    first base of the first read?
  * Q. What happens if you press [Tab] in the middle of typing the
    file path?
  * Q. Why was `zless` used instead of `less`?
	
- Next, get a basic summary statistics on the FASTQ file using the
  following command
> <pre>
zcat ~/data/1000g/fastqs/HG00406.R1.fastq.gz | wc -l </pre>
  * Q. What did appear in the screen?
  * Q. What do you think the output means? 
  * Q. How many reads are contained in the FASTQ file?
  * Q. How many reads are contained in the other pairs, contained in
       `~/data/1000g/fastqs/HG00406.R2.fastq.gz`?
  
- Next, let's check whether the two paired FASTQ files have matching
  read names in each line, using the following command.
> <pre>
paste <(zcat ~/data/1000g/fastqs/HG00406.R1.fastq.gz) <(zcat ~/data/1000g/fastqs/HG00406.R2.fastq.gz) | less -S </pre>
   * Q. What do you see in the output?
   * Q. Can you explain what the command is doing?
   * Q. What is the difference between `less` and `less -S`?
   * Q. How can you verify that the FASTA files have matching read names down to the end of the file?	 
     * If you can't figure out the answer, feel free to __[open a clue](day1-clue-ii-a-1)__


---

#### b. Understand and access the reference genome

- Next, we will align the raw sequence reads in FASTQ format to the
  human reference genome using `bwa mem` software tool. To do so, we
   will use GRCh38 reference genome. Let's see if we have 

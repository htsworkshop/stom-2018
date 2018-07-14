---
layout: page
title: Day 1 
permalink: /day1/
---

# Day 1. Analysis of High-throughput Sequence Reads

<hr>
<br>

Sequencing human genomes is becoming more and more afforable.

Most bioinformatics happens on Unix/Linux platforms but why and how do we use Unix?

Increasingly, the raw output of biological research exists as _in silico_ data, usually in the form of large text files. Unix is particularly suited to working with such files and has many powerful (and flexible) commands that can process your data for you. 

The real strength of learning Unix is that most of these commands can be combined in an almost unlimited fashion. So if you can learn just five Unix commands, you will be able to do a lot more than just five things. Our objective here is to learn a subset of Unix and to become a productive Unix user without knowing or using every program and feature.

<br>
<hr>
<br>

### Schedule:

| Session | Time           | Topics                   | 
| :-----: |:--------------:| :----------------------- | 
| I       | 9:20-10:50 AM  | **Design and analysis of sequencing studies in population scale** | 
|         | 10:50-11:00 AM | Coffee Break             | 
| II      | 11:00-12:30 PM | **Alignment and quality control of sequenced genomes**       | 
|         | 12:30-2:00 PM  | Lunch                    | 
| III     | 2:00-3:30 PM   | **Calling short variants (SNPs and Indels) from sequence reads**    | 
|         | 3:30-3:40 PM   | Coffee Break             | 
| IV      | 3:40-5:10 PM   | **Calling structural variants (large deletions and CNVs) from sequence reads**   | 

<br>
<hr>
<br>

### Instructors:
Hyun Min Kang (HMK)
Goo Jun (GJ)

<br>
<hr>
<hr>
<br>

## Topics:

### 1)   Design and analysis of sequencing studies in population scale

- This session consists of lectures only, and there will be no
  hands-on practice after the lecture.
  
<br>
<hr>
<hr>
<br>

### 2) Alignment and quality control of sequenced genomes

- The hands-on session will begin right after the completion of the
  lecture part.
  
---
  
#### a. Make connection to server

- _IMPORTANT NOTE:_ Make sure that you have received the login information to access the KOBIC
  server.
- _IMPORTANT NOTE:_ For security reasons, we will denote the IP address of the
  server as `IP.AD.DR.ESS` and the login ID as `kobic00`. Please
  replace these information with your own login ID and IP address.
- If you have Windows laptop, open __MobaXTerm__ and start a local
  terminal. 
- For Mac OS X users, open the __Terminal__ utility.
- To connect to the server, type `ssh kobic00@IP.AD.DR.ESS [Enter]` in your shell.
```
$ ssh kobic00@IP.AD.DR.ESS
kobic@@IP.AD.DR.ESS's password: [Need to type password here]
Last login: Sun Jul 15 05:59:10 2018 from IP.AD.DR.ESS
[kobic00@localhost ref]$ 
```

---

#### b. Make a symbolic link to the data directory

- Depending on users, the directory structures are slightly different.
  Notably, there are two replicated data directories, `/BiO/data`, and
  `/Bio2/data`. To avoid performance slowdown, it is important to
  use designated directories for each user.
- To ensure that each user access their designated data directory,
  we will create a symbolic link under your home directory to point
  the designated data directory.
- If the last two digits of your username is between `01` and `25`,
  type the following commands. The second command is for verification
  purpose. Make sure that `BiO` is spelled correctly with lowercase
  `i`.
```
[kobic00@localhost ~]$ pwd
/BiO/home/kobic00
[kobic00@localhost ~]$ ln -s /BiO/data/ data
[kobic00@localhost ~]$ ls -l
total 0
lrwxrwxrwx 1 kobic00 kobic 11 Jul 15 08:08 data -> /BiO/data/
```
- If the last two digits of your username is between `26` and `50`,
  type the following commands. The second command is for verification
  purpose. Make sure that `BiO` is spelled correctly with lowercase `i`
```
[kobic00@localhost ~]$ pwd
/BiO2/home/kobic00
[kobic00@localhost ~]$ ln -s /BiO2/data/ data
[kobic00@localhost ~]$ ls -l
total 0
lrwxrwxrwx 1 kobic00 kobic 11 Jul 15 08:08 data -> /BiO2/data/
```

---

#### c. Browse raw sequence reads

- The first task today is to read the example FASTQ files extracted
from the 1000 Genomes Project. Type the following command
```
[kobic00@localhost ~]$ zless ~/data/1000g/fastqs/HG00406.R1.fastq.gz
```
  * Q. What did appear in the screen? 
  * Q. Can you tell what is the base call and quality of the very
    first base of the first read?
  * Q. What happens if you press [Tab] in the middle of typing the
    file path?
  * Q. Why was `zless` used instead of `less`?
	
- Next, get a basic summary statistics on the FASTQ file using the
  following command
```
[kobic00@localhost ~]$ zcat ~/data/1000g/fastqs/HG00406.R1.fastq.gz | wc -l
```
  * Q. What did appear in the screen?
  * Q. What do you think the output means? 
  * Q. How many reads are contained in the FASTQ file?
  * Q. How many reads are contained in the other pairs, contained in
       `~/data/1000g/fastqs/HG00406.R2.fastq.gz`?
  
- Next, let's check whether the two paired FASTQ files have matching
  read names in each line, using the following command.
```
[kobic00@localhost ~]$ paste <(zcat ~/data/1000g/fastqs/HG00406.R1.fastq.gz) \
  <(zcat ~/data/1000g/fastqs/HG00406.R2.fastq.gz) | less -S
```
   * Q. What do you see in the output?
   * Q. Can you explain what the command did?
   * Q. 

—- Lunch Break [1:30 hr] —

#### III) Calling short variants (SNPs and Indels) from sequence reads
[1:30hr] HMK


—- Coffee Break [10 mins] —

#### IV) Calling structural variants (large deletions and CNVs) from
sequence reads GJ

- Project organization   
- Compiling software


—- End/Wrap-Up —

<br>

### Reference material
[Example data: bootcamp_01_unix.tar.gz](../class-material/bootcamp_01_unix.tar.gz) [13M]  
[Slides-1.1](../class-material/slides_day1-1_unix-motivation.pdf)  
[Slides-1.2](../class-material/slides_day1-2_unix-intro.pdf)  
[Slides-1.3](../class-material/slides_day1-3_unix-work.pdf)  

<!--- files dont exist yet... 
[Slides-1.4](../class-material/slides_day1-4_unix-compiling.pdf)  
-->

[A Quick Guide to Organizing Computational Biology Projects (Noble 2009)](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424)  

[Unix Referene Commands and Glossary](../class-material/unix-reference.html)  

[Introduction to Bash Shell Scripting](https://en.wikibooks.org/wiki/Bash_Shell_Scripting)  

[Make Install tutorial](http://www.ee.surrey.ac.uk/Teaching/Unix/unix7.html)  


  


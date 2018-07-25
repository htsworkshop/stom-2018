---
layout: page
---

--- 

### II-1. Setting up the access to the compute server.

In this part, we will briefly learn how to access the Won/Park lab compute
servers and to initially set up environments.

---

#### a. Make a connection to the server

- _IMPORTANT NOTE:_ Make sure that you have the login information to
	access the Won/Park lab server. If you do not have an account,
	please let the instructor know.
- _IMPORTANT NOTE:_ For security reasons, we will denote the IP address of the
  server as `IP.AD.DR.ESS` and the login ID as `username`. Please
  replace this information with your own login ID and IP address (or
  hostname) of the server.
- The first thing you need to do is open a local terminal to access a
  shell prompt that you can type in.
  * If you have Windows laptop, open `MobaXTerm` and start a local
  terminal. Follow
  __[this link](https://mobaxterm.mobatek.net/download-home-edition.html)__ to
  download the software.
  * For Mac OS X users, just open the `Terminal` utility. If you do
  not know where it is located, open `Finder`, and click
  `Applications` on the left panel, and then open the `Utilities` folder,
  then you will be able to open `Terminal.app`. Keeping the
  `Terminal` application in the Dock would probably good for the next
  few days.
- To connect to the server, type `ssh -p 9001 username@IP.AD.DR.ESS [Enter]` in
  your shell.
> <pre>
$ ssh -p 9001 username@IP.AD.DR.ESS
username@@IP.AD.DR.ESS's password: [Need to type password here]
Last login: Sun Jul 15 05:59:10 2018 from IP.AD.DR.ESS
username@HOST:~$  </pre>
  * Q. Were you able to successfully connect?
  * Q. What are the same and different from the example output shown
    above?

---

#### b. Check access to the data directory

- Depending on users, the directory structures are slightly different.
  Notably, there are two replicated data directories, `/ws_data`. Make
  sure that you have access to the directory.

> <pre>
$ pwd
/home/username
$ ls -l /ws_data
total 24
drwxr-xr-x 6 user1 user1 4096  7월 16 16:55 1000g
drwxr-xr-x 3 user1 user1 4096  7월 17 14:44 10x
drwxrwxr-x 4 user1 user1 4096  7월 25 20:04 bin
drwxrwxr-x 2 user1 user1 4096  7월 16 22:43 bravo
drwxrwxr-x 2 user1 user1 4096  7월 25 19:23 ref
drwxrwxr-x 5 user1 user1 4096  7월 25 19:36 topmed

---

Is everything OK? Otherwise, do not hesitate to ask questions to your
instructors, as this step is extremely important to be able to continue all
the next steps. 

If everything is okay, then let's go to next step [II-2 Processing raw sequence reads](../class-material/day1-fastq-practice.html)
, or go back to [Day 1 Overview](../day1).

---

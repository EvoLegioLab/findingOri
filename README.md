![findingOri logo](https://github.com/EvoLegioLab/findingOri/blob/master/logo/findingOriLogo.png)

Ori and Ter prediction for prokaryotic circular chromosomes.

findingOri("accession",,) # NCBI accession number
findingOri(,"file path" ,) # Local fasta file path, UNIX standard
findingOri("accession","v") # verbose, with plots
findingOri(,"file path" ,"v") # verbose, with plots

### **New in 1.0.0 (October 2017)**

> <sub>* First release.

### **Getting Started**

FindingOri consists of a set of R files, with a wrapper that calls R packages and other source files. To run it you need to first install Prodigal under bash [ [Installing Prodigal](https://www.github.com/hyattpd/Prodigal/wiki/installation)]. Called R packages are: seqinr, ape, stringr and pracma. You can then install from source as follows:

$ source("wrapper.R")

### **Features**

> <sub>* Predicts Ori and Ter of small circular chromosomes from FASTA DNA sequences.

> <sub>* The user need to ensure that the target chromosome is the only sequence.

> <sub>* Runs quickly: findingOri analyzes the *E. coli K-12* genome in 60 seconds on a modern MacBook Pro.

> <sub>* Runs unsupervised.

### **More Information**

> <sub>* [Website](https://github.com/EvoLegioLab/findingOri)

#### **Contributors**

> <sub>* Author: Torbjörn Larsson

#### **License**

[GPL](https://github.com/EvoLegioLab/findingOri/blob/master/LICENSE)


![findingOri logo](https://github.com/EvoLegioLab/findingOri/blob/master/logo/rsz_findingorilogo.png)

Ori and Ter prediction for prokaryotic circular chromosomes.

```bash
findingOri("accession",,) # NCBI accession number
findingOri(,"file path" ,) # Local fasta file path, UNIX standard
findingOri("accession","v") # verbose, with plots
findingOri(,"file path" ,"v") # verbose, with plots
```

### New in 1.0.0 (October 2017)

* First release.

### Getting Started

FindingOri consists of a set of R files, with a wrapper that calls R packages and other source files. To run it you need to first install Prodigal under bash [ [Installing Prodigal](https://www.github.com/hyattpd/Prodigal/wiki/installation)]. Called R packages are: seqinr, ape, stringr and pracma. You can then install from source as follows:

```bash
$ source("wrapper.R")
```

### Features

* **Predicts Ori and Ter of small circular chromosomes** from FASTA DNA sequences.

* The user need to ensure that the target chromosome is **the only sequence**.

* **Runs quickly:** findingOri analyzes the *E. coli K-12* genome in 60 seconds on a modern MacBook Pro.

* **Runs unsupervised.**

### More Information

* [Website](https://github.com/EvoLegioLab/findingOri)

#### Contributors

* Author: Torbjörn Larsson

#### License

[GPL](https://github.com/EvoLegioLab/findingOri/blob/master/LICENSE)


#**idability.py: build and evaluate hitting-set-based codes**

**Authors**: 

* Eric A. Franzosa (<mailto:franzosa@hsph.harvard.edu>)
* Lauren McIver
* Curtis Huttenhower

**Contents**:

[TOC]

##**Description**

`idability.py` is a standalone python program which processes tabular data to build and evaluate hitting-set-based codes. Table rows represent features, and table columns represent "samples" or "subjects". The entry (*F*, *S*) represents subject *S*'s value for feature *F*. In their simplest form, values are 1s and 0s, representing the presence or absence (respectively) of a feature within a given subject. In this simple example:

```
   S1 S2 S3
F1  0  1  1
F2  1  0  1
F3  1  1  0
F4  1  1  1
```

Subject *S1* has features *F2*, *F3*, and *F4* but is missing feature *F1*. 

A *code* is a set of features that uniquely identify a given subject, in the sense that the features collectively distinguish that subject from the background population. Put another way, the subject's code features are never observed together in another subject. Features can be prioritized for inclusion in a code in a number of ways. In the simplest method, subject *S*'s code is built from the rarest features in the population (more specifically, the code is grown iteratively by adding the rarest feature among non-*S* subjects that are not already excluded by a feature in the code). Under this formulation, the code is of minimal size, and its construction is analogous to the classical [Hitting-Set Problem](http://en.wikipedia.org/wiki/Set_cover_problem#Hitting_set_formulation). In the simple example above, {*F2*, *F3*} is a minimal code for *S1*.

`idability.py` has been additionally optimized to build personalized codes from microbiome data (as derived from metagenomic sequencing experiments). In this context, minimal codes are suboptimal as they tend not to be robust to temporal variation, and so an alternative feature prioritization scheme is employed. Please refer to the [Advanced Configuration](#Advanced Configuration) and [Metagenomic Codes Demo](#Metagenomic Codes Demo) sections below for further details.

##**Citation**

If you use `idability.py` or the datasets provided here in a publication, please cite:

Franzosa EA, Katherine H, Meadow JF, Gevers D, Lemon KP, Bohannan BJM, Huttenhower C. [Identifying personal microbiomes using metagenomic codes.](http://www.pnas.org/content/early/2015/05/08/1423854112) Proceedings of the National Academy of Sciences (2015): 201423854.

##**Prerequisites**

``idability.py`` requires python 2.7+.

##**Installation**

Clone the repository via the following command:

```
$ hg clone https://bitbucket.org/biobakery/idability
```

or [directly download the repository](https://bitbucket.org/biobakery/idability/get/default.tar.gz) as an archive and expand its contents.

##**Minimal Codes Demo**

The repository contains a simple demo to illustrate basic use of the software (similar in spirit, but larger, than the minimal example above). To run the demo, execute:

```
$ ./idability.py demo1.pcl
```

This produces a codes file called `demo1.codes.txt`. The first few lines of this file look like this:

```
#SAMPLE CODE
S01     F15     F02     F13
S02     F14     F10     F05
S03     F08     F13
...
```

This indicates, for example, that the set of features {*F02*, *F13*, *F15*} were a unique code for subject *S01*. You can verify this by inspecting the input file `demo1.pcl`. If `idability.py` is passed a set of codes in addition to a table, it will apply the codes to the table and report which codes matched which subjects:

```
$ ./idability.py demo1.pcl --codes demo1.codes.txt
```

This produces a file called `demo1.demo1.hits.txt`. (The general form of the "hits" file is `INPUT_TABLE.CODES_FILE.hits.txt`; this can be configured using the program's `--output` flag.) The first few lines of the hits file look like:

```
# 1|TP: 12
# 2|TP+FP: 0
# 3|FN+FP: 0
# 4|FN: 0
# 5|NA: 3
S01     matches S01
S02     matches S02
S03     matches S03
...
```

The first five lines (beginning with `#`) represent a confusion matrix: they summarize which codes matched the correct subjects (true positives, TPs), which codes failed to match the correct subject (false negatives, FNs), and which codes spuriously matched the wrong subjects (false positives, FPs). Combinations of these categories are allowed. The "NA" category indicates the presence of a "null code": three of the subjects from `demo1.pcl` did not have a unique code as their features were a subset of some other subject's features.

The demo contains a second file, `demo2.pcl`, which represents a perturbation of the original table. Running the following command:

```
$ ./idability.py demo2.pcl --codes demo1.codes.txt
```

Produces a file called `demo2.demo1.hits.txt`. The first few lines of this file look like:

```
# 1|TP: 7
# 2|TP+FP: 3
# 3|FN+FP: 1
# 4|FN: 1
# 5|NA: 3
S01     matches S01
S02     no_matches
S03     matches S03     S05     S06
...
```

Note that the results now contain instances of FNs and FPs in addition to TPs. For example, *S02* has been perturbed and its original code no longer matches the features detected in *S02* in the second dataset (a FN). *S03*'s original code continues to match *S03*, but now two other subjects (*S05* and *S06*) also match this code (FPs).

##**Metagenomic Codes Demo**

`idability.py` was originally developed to explore individual-specific adaptation of human microbiomes. Specifically, we were interested to know if body sites within an individual contain a collection of microbial taxa or genes that uniquely distinguish that individual from the population. This is equivalent to the hitting-set based code construction process described above. However, our investigation revealed that minimal metagenomic codes were unstable over time. Hence, we adapted the classical greedy approach to minimal hitting-set construction to instead prioritize identification of stable metagenomic codes. The repository contains a demo based on microbial marker genes sampled from individuals involved in the [Human Microbiome Project](http://www.hmpdacc.org/) as surveyed by the [MetaPhlAn](http://huttenhower.sph.harvard.edu/metaphlan) software package. (A *marker gene* is a gene that is consistently found in isolate genomes from a given clade [here, bacterial and archaeal species] and not found outside that clade.)

To begin, unzip the two data files, which contain marker measurements for a set of 50 individuals' gut microbiomes (as represented from stool samples) sampled ~6 months apart.

```
$ gunzip markers-stool-visit1.pcl.gz 
$ gunzip markers-stool-visit2.pcl.gz 
```

Try running the default code construction process used above on the visit1 file, and then applying the visit1 codes to the visit2 table:

```
$ ./idability.py markers-stool-visit1.pcl
$ ./idability.py markers-stool-visit2.pcl --codes markers-stool-visit1.codes.txt
```

This yields:

```
# 1|TP: 4
# 2|TP+FP: 1
# 3|FN+FP: 22
# 4|FN: 23
# 5|NA: 0
...
```

The results are less than stellar due to the prioritization of minimal (unstable) codes. Repeat this process by running the program in "meta_mode", which fine-tunes the code construction process to identify sets of features that are more likely to be stable over time:

```
$ ./idability.py markers-stool-visit1.pcl --meta_mode rpkm
$ ./idability.py markers-stool-visit2.pcl --codes markers-stool-visit1.codes.txt --meta_mode rpkm
```

This yields:

```
# 1|TP: 43
# 2|TP+FP: 0
# 3|FN+FP: 1
# 4|FN: 6
# 5|NA: 0
...
```

The results are much better: the majority of individuals' visit2 samples still match their visit1 codes, and spurious matches are rare.

##**Advanced Configuration**

The `--meta_mode` parameter used in the Microbial Community Demo above encapsulates a variety of advanced parameter settings to the `idability.py` program. These serve to increase code robustness and specificity in metagenomics applications. Settings include:

* A novel feature prioritization scheme (*abundance gap sorting*)
* More flexible definitions of feature *presence* and *absence*
* A minimum code size threshold
* Methods to exclude redundant features

These settings can be individually fine-tuned for user-specific applications. Consult the program's help menu to learn more about them:

```
$ ./idability.py -h
```

##**Preprocessed Metagenomics Datasets**

The following datasets were used in the publication and are available for download here. There is one gzipped tarball for each type of metagenomic features (OTUs, species, marker genes [markers], and kilobase windows [kbwindows]). Each tarball contains pairs of files for each body site: one file is a table reflecting subjects' feature measurements at their first visit (visit1), and a second file contains values from the follow-up visit (visit2). For example, `otus-tables.tar.gz` contains:

* `otus-anterior_nares-visit1.pcl`
* `otus-anterior_nares-visit2.pcl`

(Along with similar pairings for 17 other body sites.) The table formats are as described above (tab-delimitted with feature rows and subject columns). All data are derived from raw sequencing reads and analysis conducted during the [Human Microbiome Project (HMP)](http://www.hmpdacc.org/). Additional details of the features and datasets are available in the publication. Subjects are identified by HMP-issued `RANDSID` identifiers (column headers). Tables have been pre-processed to ensure that each subject appears in both tables with the same identifier. The OTU and species tables contain values in relative abundance units, while the MetaPhlAn markers and kilobase windows tables contain values in reads mapped per kilobase of genomic feature per million total reads (RPKM).

* [Download OTU paired tables](https://bitbucket.org/biobakery/idability/downloads/otus-tables.tar.gz)
* [Download MetaPhlAn species paired tables](https://bitbucket.org/biobakery/idability/downloads/species-tables.tar.gz)
* [Download MetaPhlAn marker genes paired tables](https://bitbucket.org/biobakery/idability/downloads/markers-tables.tar.gz)
* [Download Kilobase Windows paired tables](https://bitbucket.org/biobakery/idability/downloads/kbwindows-tables.tar.gz)
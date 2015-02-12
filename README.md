[TOC]

#**idability.py: software for building and evaluating hitting-set-based codes**#

AUTHORS: Eric A. Franzosa (franzosa@hsph.harvard.edu), Curtis Huttenhower (chuttenh@hsph.harvard.edu)

##**Description**##

idability.py processes tabular data to build and evaluate hitting-set-based codes. Table rows represent features, and table columns represent "samples" or "subjects". The entry (F, S) represents subject S's value for feature F. In their simplest form, values are 1s and 0s, representing the presence or absence (respectively) of a feature within a given subject. In this simple example:

```
#!cmd
   S1 S2
F1  0  1
F2  1  0
```

Subject *S1* has feature *F2* only, and subject *S2* has feature *F1* only. A "code" is a set of features that uniquely identify a given subject. Features can be prioritized for inclusion in a code in a number of ways. In the simplest method, subject S's code is built from the rarest features in the population (more specifically, the code is grown iteratively by adding the rarest feature among non-S subjects that are not already excluded by a feature in the code). Under this formulation, the code is of minimal size, and its construction is analogous to the classical [**Hitting-Set Problem**](http://en.wikipedia.org/wiki/Set_cover_problem#Hitting_set_formulation). In the simple example above, *F2* is a code for *S1* and *F1* is a code for *S2*.

##**Basic Demo**##

The repository contains a demo to illustrate using the software (similar in spirit, but more advanced, than the minimal example above). To run the demo, execute:

```
#!cmd
./idability demo1.dat
```

This produces a codes file called ``demo1.codes.txt``. The first few lines of this file look like:

```
#!cmd
#SAMPLE CODE
S01     F15     F02     F13
S02     F14     F10     F05
S03     F08     F13
...
```

This indicates, for example, that the set of features {*F02*, *F13*, *F14*} were a unique code for subject *S01*. You can verify this by inspecting the input file, ``demo1.dat``. If ``idability.py`` is based a set of codes in addition to a table, it will apply the codes to the table and report which codes were hit:

```
#!cmd
./idability demo1.dat --codes demo1.codes.txt
```

This produces a file called ``demo1.demo1.hits.txt``. The general form of the "hits" file is ``INPUT_TABLE.CODES_FILE.hits.txt" (this can be configured using the program's ``-o, --output`` flag). The first few lines of the hits file look like:

```
#!cmd
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

The first five lines (beginning with ``#``) represent a confusion matrix: they summarize which codes matched the correct subjects (true positives, TPs), which codes failed to match the correct subject (false negatives, FNs), and which codes spurious matched the wrong subjects (false positives, FPs). Combinations of these categories are allowed. The "NA" category indicates the presence of a "null code" -- in this case, three of the subjects from demo1.dat did not have a unique code (their features were a subset of some other subject's features).

The demo contains a second file, ``demo2.dat``, which represents a perturbation of the original table. Running the following command:

```
#!cmd
./idability demo2.dat --codes demo1.codes.txt
```

Produces a file called ``demo2.demo1.hits.txt``. The first few lines of this file look like:

```
#!cmd
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

Note that the results now contain instances of FNs and FPs in addition to TPs. For example, *S02* has been perturbed and its original code no longer matches the features detected in S02 in the second dataset (a FN). *S03*'s original code continue to match *S03*, but now two other subjects (*S05* and *S06*) also match this code (FPs).


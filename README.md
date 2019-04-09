# Tax4Fun2

Welcome to the homepage of Tax4Fun2 (under development).
Older versions are also available under https://sourceforge.net/projects/tax4fun2/

**Tax4Fun2 requirements**

BLAST+ (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
R packages seqinr, ape
diamond v0.9.24 is needed for functional annotation (binaries for Windows and Linux are downloaded as part of the reference data), Mac users need to compile diamond. Please see the wiki for instructions.

**Installation instructions**

1) Please download the latest release (version 1.1)
https://github.com/bwemheu/Tax4Fun2/releases/download/V1.1/Tax4fun2_1.1.tar.gz

2) In RStudio
Select 'Install packages ..' under Tools

3) Switch install from _CRAN_ to _Package archive file_

4) Click on _Browse_ and select the Tax4Fun2 file

5) Click _Install_


**Build the default reference database and download the example data**

In order to provide a straight-forward solution, we implemented a function in Tax4Fun2 v1.1 which will download and build the reference database. A second function will download the example data. The example script with these two commands and all other functions can be found here:
https://github.com/bwemheu/Tax4Fun2/releases/download/V1.1/Tax4Fun2_1.1_example.R

Have fun!

# Autotrim v0.6.1

## Authors
Waldvogel A-M, Schell T

## Description
__Trimming (overrepresented *k*-mers from) multiple fastq files.__

After each trimming a test for overrepresented *k*-mers will be executed, if not switched off. Any overrepresented *k*-mers found will be added to a global list and the trimming will be repeated on the original input files together with the overrepesented *k*-mers until no overrepresentation is detected.  
Autotrim uses Trimmomatic for trimming, FastQC for overrepresentaion screening and MultiQC to generate summary reports.  
The dependencies will be automatically detected if `java`, `trimmomatic`, `fastqc` or `multiqc` is in your `$PATH`. Execution of MultiQC is optionally and will be skipped if not found in your `$PATH` and `-mqcp` is not specified.  
Autotrim distinguishes automatically between paired and single end data.

## Dependencies
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

__Optional__
- [MultiQC](http://multiqc.info/)

## Usage
### Prerequisites
- From FastQC version 0.11.6 the *k*-mer plot is disabled by default but is needed in Autotrim. To enable it, edit `FastQC/Configuration/limits.txt` by changing the `1` in the end of the line starting with `kmer` to `0`.

__Optional__
- If you don't want to specify `-tp`, `trimmomatic` should be in your `$PATH`. To do so, create an executable file named `trimmomatic` with the content `java -jar /path/to/trimmomatic.jar $*` in a directory that is your `$PATH`. Be sure to write the actual path to the `jar` file.
### Input
There are two different ways to specify the input files:

1. Multiple folders containing one or two fastq files are in the same directory. The fastq file(s) need to end with `.f(ast)q(.gz)`. For paired data the two fastq files have to end with `_1.f(ast)q(.gz)` and `_2.f(ast)q(.gz)`. Folders containing more than two fastq files or paired data without matching the criteria are skipped.
2. File of file names containing a path to a single end fastq or two tab seperated paths for paired end data per line. The files for paired end data don't have to be in the same folder nor any special file names.

__While choosing a root direcory will skip the subdirectory when it contains more than two fastq files (e.g. output), the file of file names option will overwirte existing files without asking!__

### Output
Fastq output files are created in the same directory as the corrensponding input file.
The `.f(ast)q(.gz)` file extension will be removed and `_autotrim.fq` for single end data and `_autotrim.paired.fq` and `_autotrim.unpaired.fq` will be added to the end of the file name respectively.  
FastQC reports will be created for single end input data and the paired output fastq files for paired end input (not for the unpaired output).  
Standard out `*.trimmomatic_log` and standard error `*.trimmomatic_err` from each Trimmomatic run will be saved in the same folder as the input.  
The MultiQC report will be placed either in the root directory if `-d` is used or in the directory of the log and overrepresented *k*-mers files if `-fofn` is used.

```
autotrim.pl [-d <root_dir> | -fofn <file_of_file_names> -log <dir_to_place_the_log-file>]
            -to <trimmomatic_options> -trim <trimmomatic_trimmer>

Options: [default]
	-d STR		Root directory for automatic input of multiple data sets. The file containing
			overrepresented k-mers (kmer.fa) and the autotrim log file (autotrim.log) will be
			saved in this directory.
	-fofn STR	File of file names containing tab-seperated paths to one data set per line.
	-log STR	Specify the path to save the file containing overrepresented k-mers (kmer.fa) and
			the log file of autotrim (autotrim.log) if using -fofn.
	-to STR		A file containing Trimmomatic options that should be used. All options need to be
			in the first line of the file.
			Create a trimlog for every data set writing "-trimlog" without a path, it will be
			saved in the same folder as the single or "_1" input file.
	-tt INT		Trimmomatic threads. Specify either within -to or with -tt. [1]
	-trim STR	A file containing Trimmomatic trimmer in the first line in the particular order
			they should be executed.
			If trimming overrepresented k-mers, the "ILLUMINACLIP" will be inserted after the
			last specified "ILLUMINACLIP".
	-k INT		K-mer length to screen for overrepresentaion with FastQC between 2 and 10. [7].
	-nok		No overrepresentation screening of k-mers. Trim each set once and create a FastQC
			report.
			Recommended for RNA-seq. [off]
	-v		Verbose. Print executed commands of Trimmomatic, FastQC and MultiQC to STDOUT and
			log file. [off]
	-tp STR		Trimmomatic path. The whole path to the Trimmomatic jar file. Specify if
			"trimmomatic" is not in your $PATH.
	-fqcp STR	FastQC path. The whole path to the FastQC executable. Specify if "fastqc" is not
			in your $PATH.
	-mqcp STR	MultiQC path. The whole path to the MultiQC executable. Specify if "multiqc" is
			not in your $PATH and you want MultiQC automatically executed.
	-rn		Rename files according the folder they are placed in. [off]
	-version	Print version number and exit.
	-h or -help	Print this help and exit.
```

## Citation
If you use this tool please cite:

- Waldvogel A-M, Wieser A, Schell T, Patel S, Schmidt H et al. (2018). The genomic footprint of climate adaptation in *Chironomus riparius*. *Molecular Ecology*, 27(6):1439–1456. [https://doi.org/10.1111/mec.14543](https://doi.org/10.1111/mec.14543)

Additional to the dependencies:

- Trimmomatic  
Bolger AM, Lohse M, Usadel B (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics, 30(15):2114–2120. [https://doi.org/10.1093/bioinformatics/btu170](https://doi.org/10.1093/bioinformatics/btu170)
- FastQC  
Andrews S (2010). FastQC: a quality control tool for high throughput sequence data. [http://www.bioinformatics.babraham.ac.uk/projects/fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc).
- MultiQC  
Ewels P, Magnusson M, Lundin S, Käller M (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*, 32(19):3047–3048. [https://doi.org/10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354)

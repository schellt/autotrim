Autotrim v0.5

Authors:
	Waldvogel AM, Schell T

Citation:
	If you use this tool please cite:
	Waldvogel AM, Wieser A, Schell T, Patel S, Schmidt H, Hankeln T, Feldmeyer B,
	Pfenninger M. (in press, status: 2018-02-15). The genomic footprint of
	climate adaptation in Chironomus riparius. Molecular Ecology.

Description:
	Trimming (overrepresented k-mers from) multiple fastq files.

	After each trimming a test for overrepresented k-mers will be executed, if not switched
	off. Any overrepresented k-mers found will be added to a global list and the trimming will
	be repeated on the original input files together with the overrepesented k-mers until no
	overrepresentation is detected. Autotrim uses Trimmomatic for trimming and FastQC for
	overrepresentaion screening. The dependencies will be automatically detected if
	"trimmomatic" and "fastqc" is in your $PATH. Autotrim distinguishes automatically between
	paired and single end data.
	
	Trimmomatic:	http://www.usadellab.org/cms/?page=trimmomatic
	FastQC:		https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

	INPUT
	There are two different ways to specify the input files:
		1) Multiple folders containing one or two fastq files are in the same directory.
		   The fastq files need to end with .fastq(.gz) or .fq(.gz).
		   For paired data the two fastq files have to end with "_1.f(ast)q(.gz)" or
		   "_2.f(ast)q(.gz)".
		   Folders containing more than two fastq files or paired data without ending in
		   "_1.f(ast)q(.gz)" and "_2.f(ast)q(.gz)" are skipped.
		2) File of file names containing a path to a single end fastq or two tab seperated
		   paths for paired end data per line.
		   The files for paired end data don't have to be in the same folder nor any
		   special file names.
	=> While choosing a root direcory will skip the subdirectory when containg more than two
	fastq files (e.g. output), the file of file names option will overwirte existing files
	without asking!

	OUTPUT
	Fastq output files are created in the same directory as the corrensponding input file.
	The f(ast)q(.gz) file extension will be removed and "_autotrim.fq" for single end data and
	"_autotrim.paired.fq" respective "_autotrim.unpaired.fq" will be added to the end of the
	file name. FastQC reports will be created for single end data and the paired output fastq
	files for paired end input (not for the unpaired output)
	Standard out "trimmomatic_log" and standard error "trimmomatic_err" from each Trimmomatic
	run will be saved in the same folder as the input.

Usage:
	autotrim_0.5.pl [-d <root_dir> | -fofn <file_of_file_names>
	-log <dir_to_place_the_log-file> -kfa <file_to_save_overrepresented_kmers>]
	-to <trimmomatic_options> -trim <trimmomatic_trimmer>
	
Options: [default]
	-d STR		Root directory for automatic input of multiple data sets. The file
			containing overrepresented k-mers (kmer.fa) and the autotrim log file
			(autotrim.log) will be saved in this directory.
	-fofn STR	File of file names containing paths to one data set per line.
	-log STR	Specify the path to save the file containing overrepresented k-mers
			(kmer.fa) and the log file of autotrim (autotrim.log) if using -fofn.
	-to STR		A file containing Trimmomatic options that should be used. All options
			need to be in the first line of the file.
			Create a trimlog for every data set writing "-trimlog" without a path,
			it will be saved in the same folder as the single or "_1" input file.
	-tt INT		Trimmomatic threads. Specify either with -to or with -tt. [1]
	-trim STR	A file containing Trimmomatic trimmer in the first line in the
			particular order they should be executed.
			If trimming overrepresented k-mers, the "ILLUMINACLIP" will be inserted
			after the last specified "ILLUMINACLIP".
	-k INT		Kmer length to screen for overrepresentaion with FastQC between 2 and 10.
			[7].
	-nok		No overrepresentation screening of kmers. Trim each set once and create
			a FastQC report.
			Recommended for RNA-seq. [off]
	-v		Verbose. Print executed commands of Trimmomatic and FastQC to STDOUT and
			log file. [off]
	-tp STR		Trimmomatic path. The whole path to the Trimmomatic jar file. Specify if
			"trimmomatic" is not in your $PATH.
	-fqcp STR	FastQC path. The whole path to the FastQC executable. Specify if "fastqc"
			is not in your $PATH.
	-rn		Rename files according the folder they are placed in. [off]
	-h or -help	Print this help and exit.

Known issues:
	- With -d empty folders count as jobs in the progessmeter
	
	- No error when Trimmomatic wasn't run successful but error from FastQC which can't find
	  the output files
	  
	- File of file names (-fofn) recognizes only tab separated but not space separated lines
	  as paired end
	  
	- No option for maximum number of iterations
	
	- No check where an overrepresented k-mer is located in the read
	
	- When no Trimmomatic options are specified with -to perl shows the
	  "Use of uninitialized..." error

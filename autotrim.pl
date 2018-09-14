#!/usr/bin/perl

use strict;
use warnings;
use Cwd 'abs_path';
use IPC::Cmd qw[can_run run];

sub print_help{
	print STDOUT "\n";
	print STDOUT "Autotrim v0.6\n";
	print STDOUT "\n";
	print STDOUT "Authors:\n";
	print STDOUT "\tWaldvogel A-M, Schell T\n";
	print STDOUT "\n";
	print STDOUT "Citation:\n";
	print STDOUT "\tIf you use this tool please cite:\n\tWaldvogel A‐M, Wieser A, Schell T, et al. The genomic footprint of climate adaptation in\n\tChironomus riparius. Mol Ecol. 2018;27:1439–1456. https://doi.org/10.1111/mec.14543\n";
	print STDOUT "\n";
	print STDOUT "Description:\n";
	print STDOUT "\tTrimming (overrepresented k-mers from) multiple fastq files.\n";
	print STDOUT "\n";
	print STDOUT "\tAfter each trimming a test for overrepresented k-mers will be executed, if not switched off.\n";
	print STDOUT "\tAny overrepresented k-mers found will be added to a global list and the trimming will be repeated\n\ton the original input files together with the overrepesented k-mers until no overrepresentation\n\tis detected.\n";
	print STDOUT "\tAutotrim uses Trimmomatic for trimming and FastQC for overrepresentaion screening and MultiQC to\n\tgenerate summary statistics.\n";
	print STDOUT "\tThe dependencies will be automatically detected if \"java\", \"trimmomatic\", \"fastqc\" or \"multiqc\"\n\tis in your \$PATH. Execution of MultiQC is optionally and will be skipped if not found in your\n\t\$PATH and -mqcp is not specified.\n";
	print STDOUT "\tAutotrim distinguishes automatically between paired and single end data.\n";
	print STDOUT "\n";
	print STDOUT "\tINPUT\n";
	print STDOUT "\tThere are two different ways to specify the input files:\n";
	print STDOUT "\t\t1) Multiple folders containing one or two fastq files are in the same directory.\n";
	print STDOUT "\t\t   The fastq files need to end with .fastq(.gz) or .fq(.gz).\n";
	print STDOUT "\t\t   For paired data the two fastq files have to end with \"_1.f(ast)q(.gz)\" or\n\t\t   \"_2.f(ast)q(.gz)\".\n";
	print STDOUT "\t\t   Folders containing more than two fastq files or paired data without ending in\n\t\t   \"_1.f(ast)q(.gz)\" and \"_2.f(ast)q(.gz)\" are skipped.\n";
	print STDOUT "\t\t2) File of file names containing a path to a single end fastq or two tab seperated paths\n\t\t   for paired end data per line.\n";
	print STDOUT "\t\t   The files for paired end data don't have to be in the same folder nor any special file\n\t\t   names.\n";
	print STDOUT "\t=> While choosing a root direcory will skip the subdirectory when containg more than two fastq\n\t   files (e.g. output), the file of file names option will overwirte existing files without\n\t   asking!\n";
	print STDOUT "\n";
	print STDOUT "\tOUTPUT\n";
	print STDOUT "\tFastq output files are created in the same directory as the corrensponding input file.\n";
	print STDOUT "\tThe f(ast)q(.gz) file extension will be removed and \"_autotrim.fq\" for single end data and\n\t\"_autotrim.paired.fq\" or \"_autotrim.unpaired.fq\" will be added to the end of the file name\n\trespectively.\n";
	print STDOUT "\tFastQC reports will be created for single end data and the paired output fastq files for paired\n\tend input (not for the unpaired output)\n";
	print STDOUT "\tStandard out \"trimmomatic_log\" and standard error \"trimmomatic_err\" from each Trimmomatic run\n\twill be saved in the same folder as the input.\n";
	print STDOUT "\tThe MultiQC report will be placed either in the root directory if -d used or in the directory of\n\tthe log and overrepresented k-mers files if -fofn is used\n";
	print STDOUT "\n";
	print STDOUT "Usage:\n";
	print STDOUT "\tautotrim_0.6.pl [-d <root_dir> | -fofn <file_of_file_names> -log <dir_to_place_the_log-file>]\n\t-to <trimmomatic_options> -trim <trimmomatic_trimmer>\n";
	print STDOUT "\n";
	print STDOUT "Options: [default]\n";
	print STDOUT "\t-d STR\t\tRoot directory for automatic input of multiple data sets. The file containing\n\t\t\toverrepresented k-mers (kmer.fa) and the autotrim log file (autotrim.log) will be\n\t\t\tsaved in this directory.\n";
	print STDOUT "\t-fofn STR\tFile of file names containing tab-seperated paths to one data set per line.\n";
	print STDOUT "\t-log STR\tSpecify the path to save the file containing overrepresented k-mers (kmer.fa) and\n\t\t\tthe log file of autotrim (autotrim.log) if using -fofn.\n";
	print STDOUT "\t-to STR\t\tA file containing Trimmomatic options that should be used. All options need to be\n\t\t\tin the first line of the file.\n";
	print STDOUT "\t\t\tCreate a trimlog for every data set writing \"-trimlog\" without a path, it will be\n\t\t\tsaved in the same folder as the single or \"_1\" input file.\n";
	print STDOUT "\t-tt INT\t\tTrimmomatic threads. Specify either within -to or with -tt. [1]\n";
	print STDOUT "\t-trim STR\tA file containing Trimmomatic trimmer in the first line in the particular order\n\t\t\tthey should be executed.\n";
	print STDOUT "\t\t\tIf trimming overrepresented k-mers, the \"ILLUMINACLIP\" will be inserted after the\n\t\t\tlast specified \"ILLUMINACLIP\".\n";
	print STDOUT "\t-k INT\t\tK-mer length to screen for overrepresentaion with FastQC between 2 and 10. [7].\n";
	print STDOUT "\t-nok\t\tNo overrepresentation screening of kmers. Trim each set once and create a FastQC\n\t\t\treport.\n";
	print STDOUT "\t\t\tRecommended for RNA-seq. [off]\n";
	print STDOUT "\t-v\t\tVerbose. Print executed commands of Trimmomatic, FastQC and MultiQC to STDOUT and\n\t\t\tlog file. [off]\n";
	print STDOUT "\t-tp STR\t\tTrimmomatic path. The whole path to the Trimmomatic jar file. Specify if\n\t\t\t\"trimmomatic\" is not in your \$PATH.\n";
	print STDOUT "\t-fqcp STR\tFastQC path. The whole path to the FastQC executable. Specify if \"fastqc\" is not\n\t\t\tin your \$PATH.\n";
	print STDOUT "\t-mqcp STR\tMultiQC path. The whole path to the MultiQC executable. Specify if \"multiqc\" is\n\t\t\tnot in your \$PATH and you want MultiQC automatically executed.\n";
	print STDOUT "\t-rn\t\tRename files according the folder they are placed in. [off]\n";
	print STDOUT "\t-h or -help\tPrint this help and exit.\n";
	exit;
}

my $root_dir = "";
my $trimmomatic_path = "";
my $fastqc_path = "";
my $multiqc_path = "";
my $fastqck = "";
my $fofn = "";
my $global_trim_opts = "";
my $global_trim_opts_file = "";
my $global_trimmer = "";
my $global_trimmer_file = "";
my $trimmomatic_threads = 1;
my $log_path = "";
my $kfa = "";
my $log = "";
my $verbose = 0;
my $nok = 0;
my $rename_switch = 0;

if(scalar(@ARGV) == 0){
	print_help;
}

my $input_error = 0;

for (my $i = 0; $i < scalar(@ARGV);$i++){
	if ($ARGV[$i] eq "-d"){
		$root_dir = abs_path($ARGV[$i+1]) . "/";
	}
	if ($ARGV[$i] eq "-tp"){
		$trimmomatic_path = abs_path($ARGV[$i+1]);
		$trimmomatic_path = "java -jar " . $trimmomatic_path;
	}
	if ($ARGV[$i] eq "-fqcp"){
		$fastqc_path = abs_path($ARGV[$i+1]);
	}
	if ($ARGV[$i] eq "-mqcp"){
		$multiqc_path = abs_path($ARGV[$i+1]);
	}
	if ($ARGV[$i] eq "-fofn"){
		$fofn = abs_path($ARGV[$i+1]);
	}
	if ($ARGV[$i] eq "-tt"){
		$trimmomatic_threads = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "-to"){
		$global_trim_opts_file = abs_path($ARGV[$i+1]);
	}
	if ($ARGV[$i] eq "-trim"){
		$global_trimmer_file = abs_path($ARGV[$i+1]);
	}
	if ($ARGV[$i] eq "-k"){
		$fastqck = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "-log"){
		if(not -d "$ARGV[$i+1]"){
			print STDERR "[autotrim] ERROR -log $ARGV[$i+1] is not a directory!\n";
			$input_error = 1;
		}
		$kfa = abs_path($ARGV[$i+1]) . "/kmer.fa";
		$log = abs_path($ARGV[$i+1]) . "/autotrim.log";
		$log_path = abs_path($ARGV[$i+1]) . "/";
	}
	if ($ARGV[$i] eq "-v"){
		$verbose = 1;
	}
	if ($ARGV[$i] eq "-nok"){
		$nok = 1;
	}
	if ($ARGV[$i] eq "-h" or $ARGV[$i] eq "-help"){
		print_help;
	}
	if ($ARGV[$i] eq "-rn"){
		$rename_switch = 1;
	}
}


if(not -d "$root_dir" and $root_dir ne ""){
	print STDERR "[autotrim] ERROR $root_dir is not a directory!\n";
	$input_error = 1;
}
if($root_dir ne "" and $fofn ne ""){
	print STDERR "[autotrim] ERROR Specify either a directory or a file of file names and kmer.fa path!\n";
	$input_error = 1;
}
if($root_dir ne "" and $kfa ne ""){
	print STDERR "[autotrim] ERROR Specify either a directory or a file of file names and kmer.fa path!\n";
	$input_error = 1;
}
if($root_dir eq "" and $fofn eq ""){
	print STDERR "[autotrim] ERROR No input specified!\n";
	$input_error = 1;
}
if($fofn ne "" and $kfa eq ""){
	print STDERR "[autotrim] ERROR No path for path for log file specified!\n";
	$input_error = 1;
}
if($global_trimmer_file eq ""){
	print STDERR "[autotrim] ERROR No trimmer specified!\n";
	$input_error = 1;
}
if($trimmomatic_threads !~ m/^\d+$/ or $trimmomatic_threads < 1){
	print STDERR "[autotrim] ERROR Trimmomatic threads is no integer >= 1!\n";
	$input_error = 1;
}
if($fastqck ne ""){
	if($fastqck !~ m/^\d+$/ or $fastqck < 2 or $fastqck > 10){
		print STDERR "[autotrim] ERROR FastQC k-mer length is no integer between 2 and 10!\n";
		$input_error = 1;
	}
}
if($global_trim_opts_file ne ""){
	my $opts_test = `head -n 1 $global_trim_opts_file`;
	chomp $opts_test;
	if ($trimmomatic_threads > 1){
		if ($opts_test =~ m/^-threads / or $opts_test =~ m/ -threads /){
			print STDERR "[autotrim] ERROR Trimmomatic threads specified twice! Use either -tt or -to\n";
			$input_error = 1;
		}
	}
}

if ($input_error == 1){
	print STDERR "[autotrim] ERROR Input error detected!\n";
	exit;
}


if($fastqck ne ""){
	$fastqck = "-k " . $fastqck . " ";
}
if($kfa eq ""){
	$kfa = $root_dir . "kmer.fa";
	$log = $root_dir . "autotrim.log";
}

can_run("java") or die "[autotrim] ERROR \"java\" is not in your path!\n";

if ($trimmomatic_path eq ""){
	$trimmomatic_path = can_run("trimmomatic") or die "[autotrim] ERROR \"trimmomatic\" is not in your path!\n[autotrim] Please specify the Trimmomatic path with -tp\n";
}
if ($fastqc_path eq ""){
	$fastqc_path = can_run("fastqc") or die "[autotrim] ERROR \"fastqc\" is not in your path!\n[autotrim] Please specify the FastQC path with -fqcp\n";
}

open (LOG, '>', $log);

print STDOUT "[autotrim] " . $0 . " " . join(" ",@ARGV) . "\n";
print LOG "[autotrim] " . $0 . " " . join(" ",@ARGV) . "\n";

my $trimmomatic_version = "";
if (`$trimmomatic_path -version 2> /dev/null` eq ""){
	$trimmomatic_version = "< 0.36";
}
else{
	$trimmomatic_version = `$trimmomatic_path -version`;
}
chomp $trimmomatic_version;

print STDOUT "[autotrim] Using Trimmomatic " . $trimmomatic_version . "\n";
print LOG "[autotrim] Using Trimmomatic " . $trimmomatic_version . "\n";

my $fastqc_version = `$fastqc_path -v`;
chomp $fastqc_version;

print STDOUT "[autotrim] Using " . $fastqc_version . "\n";
print LOG "[autotrim] Using " . $fastqc_version . "\n";

if ($multiqc_path eq ""){
	$multiqc_path = can_run("multiqc");
	if(not defined($multiqc_path)){
		$multiqc_path = "";
		print "[autotrim] MultiQC is not in your path!\n[autotrim] Please specify the MultiQC path with -mqcp if you want to use it\n";
	}
}
if ($multiqc_path ne ""){
	my $multiqc_version = `$multiqc_path --version`;
	chomp $multiqc_version;
	$multiqc_version=(split(/, /,$multiqc_version))[1];
	print STDOUT "[autotrim] Using MultiQC " . $multiqc_version . "\n";
	print LOG "[autotrim] Using MultiQC " . $multiqc_version . "\n";
}

if($global_trimmer_file ne ""){
	$global_trimmer = `head -n 1 $global_trimmer_file`;
	chomp $global_trimmer;
}

my @fulldirlist;
my %paths;
my %kmers;

if($root_dir ne ""){

	opendir (DIR, $root_dir) or die "[autotrim] ERROR Could not open root-directory " . $root_dir . "\n";        #opens the root dir
	
	while (my $subdir = readdir(DIR)) {					#reads the root dir
		next if ($subdir =~ m/^\./);					#ignores "." and ".." (dir up)
		next unless (-d "$root_dir$subdir");				#returns folders and ignores files in the root dir
		push (@fulldirlist,$root_dir . $subdir . "/");			#saves the full path of all found sub dirs in an array
	}
	
	for (my $i = 0; $i < scalar(@fulldirlist); $i++){
		opendir (SUBDIR, $fulldirlist[$i]) or die "[autotrim] ERROR Could not open sub-directory " . $fulldirlist[$i] . "!" . "\n";	#opens every found sub dir
		my @files = ();
		while (my $file = readdir(SUBDIR)) {				#reads the sub dir
			next unless (-f "$fulldirlist[$i]/$file");		#returns only files from the sub dir
			if ($file =~ m/\.f(ast)?q(\.gz)?$/){				#ignores files that don't match the pattern (end with .fq)
				push (@files, $fulldirlist[$i] . $file);
			}
		}
		if (scalar(@files) > 2){
			print STDERR "[autotrim] WARNING Skipping directory $fulldirlist[$i]\n";
			print LOG "[autotrim] WARNING Skipping directory $fulldirlist[$i]\n";
		}
		else{
			my $file_err = 0;
			if(scalar(@files) == 2){
				foreach(@files){
					if ($_ !~ m/_[1-2]\.f(ast)?q(\.gz)?$/){
						print STDERR "[autotrim] WARNING Unexpected filename $_\n";
						print LOG "[autotrim] WARNING Unexpected filename $_\n";
						$file_err = 1;
					}
				}
			}
			if($file_err == 0){
				if(scalar(@files) == 2 and $files[0] =~ m/_2\.f(ast)?q(\.gz)?$/){	#in rare cases the first file is not the _1 file - e.g. only the _1 is gzipped and _2 not
					my $f1 = $files[1];
					my $f2 = $files[0];
					@files = ($f1,$f2);
				}
				if(scalar(@files) != 0){
					$paths{$fulldirlist[$i]} = join(";",@files);
				}
			}
			else{
				print STDERR "[autotrim] WARNING Skipping directory $fulldirlist[$i]\n";
				print LOG "[autotrim] WARNING Skipping directory $fulldirlist[$i]\n";
			}
		}
	}

}

if($fofn ne ""){
	open (FOFN, $fofn) or die "[autotrim] ERROR Could not open file " . $fofn . "\n";
	
	my $fofn_index = 0;
	
	while (my $line = <FOFN>){
		chomp $line;
		$fofn_index++;
		my @files = split(/\t/,$line);
		if (scalar(@files) > 2){
			print STDERR "[autotrim] WARNING Skipping fofn line $fofn_index because of too many entries\n";
			print LOG "[autotrim] WARNING Skipping fofn line $fofn_index because of too many entries\n";
		}
		else{
			my $file_err = 0;
			foreach(@files){
				if(not -f $_){
					print STDERR "[autotrim] WARNING $_ does not exist. Skipping fofn line $fofn_index\n";
					print LOG "[autotrim] WARNING $_ does not exist. Skipping fofn line $fofn_index\n";
					$file_err = 1;
				}
			}
			if($file_err == 0){
				$paths{"fofn" . $fofn_index} = join(";",@files);
			}
		}
	}
}
	
my @path_keys = keys(%paths);

if(scalar(@path_keys) == 0){
	print STDERR "[autotrim] No jobs detected. Exiting.\n";
	print LOG "[autotrim] No jobs detected. Exiting.\n";
	exit;
}

my @ktrimmer;
my $status = 0;
my $frac = 0;

for (my $i = 0; $i < scalar(@path_keys); $i++){
	my $scalar_path_keys = scalar(@path_keys);
	$frac = sprintf "%.1f", ($status / scalar(@path_keys)) * 100;
	print STDOUT "[autotrim] " . $status . " / " . $scalar_path_keys . " [" . $frac . "%] jobs trimmed\r";
	my @in_files = split(/;/,$paths{$path_keys[$i]});
	if(scalar(@in_files) == 1){
		my $in = $in_files[0];
		my $out = $in;
		$out =~ s/\.f(ast)?q(\.gz)?$//;
		if($rename_switch == 1){
			my @rename = split(/\//,$out);
			$rename[-1] = $rename[-2];
			$out = join("/",@rename);
		}
		my @out = ($out . "_autotrim.fq");
		if($global_trim_opts_file ne ""){
			$global_trim_opts = `head -n 1 $global_trim_opts_file`;
			chomp $global_trim_opts;
			if($trimmomatic_threads > 1){
				$global_trim_opts = $global_trim_opts . " -threads " . $trimmomatic_threads;
			}
			$global_trim_opts = $global_trim_opts . " ";
		}
		else{
			if($trimmomatic_threads > 1){
				$global_trim_opts = "-threads " . $trimmomatic_threads . " ";
			}
		}
		$global_trim_opts = insert_trimlog($global_trim_opts,$in,$rename_switch);
		if($verbose == 1){
			print STDOUT "\e[K";
			print STDOUT "$trimmomatic_path SE $global_trim_opts$in $out[0] $global_trimmer > $out\_autotrim.trimmomatic_log 2> $out\_autotrim.trimmomatic_err\n";	#start first trimmomatic run
			print STDOUT "[autotrim] " . $status . " / " . $scalar_path_keys . " [" . $frac . "%] jobs trimmed\r";
			print LOG "$trimmomatic_path SE $global_trim_opts$in $out[0] $global_trimmer > $out\_autotrim.trimmomatic_log 2> $out\_autotrim.trimmomatic_err\n";
		}
		system("$trimmomatic_path SE $global_trim_opts$in $out[0] $global_trimmer > $out\_autotrim.trimmomatic_log 2> $out\_autotrim.trimmomatic_err");
		my $trimmomatic_error = `grep "Exception in thread" $out\_autotrim.trimmomatic_err`;
		chomp $trimmomatic_error;
		if($trimmomatic_error ne ""){
			die "[autotrim] ERROR from Trimmomatic: " . $trimmomatic_error . "\n";
		}
		$trimmomatic_error = `grep "Unknown option" $out\_autotrim.trimmomatic_err`;
		chomp $trimmomatic_error;
		if($trimmomatic_error ne ""){
			die "[autotrim] ERROR from Trimmomatic: " . $trimmomatic_error . "\n";
		}
		my ($match1,$match2) = run_fastqc(\@out,\%kmers,$status,$scalar_path_keys,$frac,$nok);	#call fastqc subroutine
		if($nok == 0){
			while ($match1 == 1 or $match2 == 1){			#if overrepresented kmers are left after trimming the trimming has to be repeated
				trimmomatic_loop($trimmomatic_path,\@in_files,$global_trim_opts,$global_trimmer,$verbose,$status,$scalar_path_keys,$frac,$rename_switch);   #call trimmomatic subroutine
				($match1,$match2) = run_fastqc(\@out,\%kmers,$status,$scalar_path_keys,$frac,$nok);       #call fastqc subroutine
			}
		}
		else{
			print STDOUT "\e[K";
			print STDOUT "[autotrim] " . ((split /\//,$in_files[0])[-2] . " trimmed successfully\n");
			print LOG "[autotrim] " . ((split /\//,$in_files[0])[-2] . " trimmed successfully\n");
			my $result = `tail -n 2 $out\_autotrim.trimmomatic_err | head -n 1`;
			chomp $result;
			print STDOUT "[Trimmomatic] " . $result . "\n";
			print LOG "[Trimmomatic] " . $result . "\n";
		}
	}
	if(scalar(@in_files) == 2){
		my $in1 = $in_files[0];
		my $in2 = $in_files[1];
		my $out1 = $in1;
		$out1 =~ s/\.f(ast)?q(\.gz)?$//;
		my $out1p;
		my $out1u;
		if($rename_switch == 1){
			my @rename1 = split(/\//,$out1);
			$rename1[-1] = $rename1[-2];
			$out1 = join("/",@rename1);
			$out1p = $out1 . "_1_autotrim.paired.fq";
			$out1u = $out1 . "_1_autotrim.unpaired.fq";
			$out1 = $out1 . "_1";
		}
		else{
			$out1p = $out1 . "_autotrim.paired.fq";
			$out1u = $out1 . "_autotrim.unpaired.fq";
		}
		my $out2 = $in2;
		$out2 =~ s/\.f(ast)?q(\.gz)?$//;
		my $out2p;
		my $out2u;
		if($rename_switch == 1){
			my @rename2 = split(/\//,$out2);
			$rename2[-1] = $rename2[-2];
			$out2 = join("/",@rename2);
			$out2p = $out2 . "_2_autotrim.paired.fq";
			$out2u = $out2 . "_2_autotrim.unpaired.fq";
		}
		else{
			$out2p = $out2 . "_autotrim.paired.fq";
			$out2u = $out2 . "_autotrim.unpaired.fq";
		}
		my @out = ($out1p,$out2p);
		if($global_trim_opts_file ne ""){
			$global_trim_opts = `head -n 1 $global_trim_opts_file`;
			chomp $global_trim_opts;
			if($trimmomatic_threads > 1){
				$global_trim_opts = $global_trim_opts . " -threads " . $trimmomatic_threads;
			}
			$global_trim_opts = $global_trim_opts . " ";
		}
		else{
			if($trimmomatic_threads > 1){
				$global_trim_opts = "-threads " . $trimmomatic_threads . " ";
			}
		}
		$global_trim_opts = insert_trimlog($global_trim_opts,$in1,$rename_switch);
		if($verbose == 1){
			print STDOUT "\e[K";
			print STDOUT "$trimmomatic_path PE $global_trim_opts$in1 $in2 $out1p $out1u $out2p $out2u $global_trimmer > $out1\_autotrim.trimmomatic_log 2> $out1\_autotrim.trimmomatic_err\n";	#start first trimmomatic run
			print STDOUT "[autotrim] " . $status . " / " . scalar(@path_keys) . " [" . $frac . "%] jobs trimmed\r";
			print LOG "$trimmomatic_path PE $global_trim_opts$in1 $in2 $out1p $out1u $out2p $out2u $global_trimmer > $out1\_autotrim.trimmomatic_log 2> $out1\_autotrim.trimmomatic_err\n";
		}
		system("$trimmomatic_path PE $global_trim_opts$in1 $in2 $out1p $out1u $out2p $out2u $global_trimmer > $out1\_autotrim.trimmomatic_log 2> $out1\_autotrim.trimmomatic_err");
                my $trimmomatic_error = `grep "Exception in thread" $out1\_autotrim.trimmomatic_err`;
		chomp $trimmomatic_error;
                if($trimmomatic_error ne ""){
                        die "[autotrim] ERROR from Trimmomatic: " . $trimmomatic_error . "\n";
                }
		$trimmomatic_error = `grep "Unknown option" $out1\_autotrim.trimmomatic_err`;
		chomp $trimmomatic_error;
		if($trimmomatic_error ne ""){
			die "[autotrim] ERROR from Trimmomatic: " . $trimmomatic_error . "\n";
		}
		my ($match1,$match2) = run_fastqc(\@out,\%kmers,$status,$scalar_path_keys,$frac,$nok);
		if($nok == 0){
			while ($match1 == 1 or $match2 == 1){
				trimmomatic_loop($trimmomatic_path,\@in_files,$global_trim_opts,$global_trimmer,$verbose,$status,$scalar_path_keys,$frac,$rename_switch);   #call trimmomatic subroutine
				($match1,$match2) = run_fastqc(\@out,\%kmers,$status,$scalar_path_keys,$frac,$nok);       #call fastqc subroutine
			}
		}
		else{
			print STDOUT "\e[K";
			print STDOUT "[autotrim] " . ((split /\//,$in_files[0])[-2] . " trimmed successfully\n");
			print LOG "[autotrim] " . ((split /\//,$in_files[0])[-2] . " trimmed successfully\n");
			my $result = `tail -n 2 $out1\_autotrim.trimmomatic_err | head -n 1`;
			chomp $result;
			print STDOUT "[Trimmomatic] " . $result . "\n";
			print LOG "[Trimmomatic] " . $result . "\n";
		}
}
	$status++;
}
$frac = sprintf "%.1f", ($status / scalar(@path_keys)) * 100;
print STDOUT "[autotrim] " . $status . " / " . scalar(@path_keys) . " [" . $frac . "%] jobs trimmed\n";
if($nok == 0){
	print STDOUT "[autotrim] " . scalar(keys(%kmers)) . " different overrepresented k-mers trimmed\n";
	print LOG "[autotrim] " . scalar(keys(%kmers)) . " different overrepresented k-mers trimmed\n";
}

my %multiqc_paths;
if ($multiqc_path ne ""){
	print STDOUT "[autotrim] Executing MultiQC\n";
	print LOG "[autotrim] Executing MultiQC\n";
	for(my $i = 0; $i < scalar(@path_keys); $i++){
		if($fofn ne ""){
			my @files = split(/;/,$paths{$path_keys[$i]});
			foreach(@files){
				my $mqc_path = $_;
				my @path = split(/\//,$_);
				splice(@path,-1);
				$mqc_path = join("/",@path) . "/";
				$multiqc_paths{$mqc_path} = "";
			}
                }
		else{
			$multiqc_paths{$root_dir} = ""; 
		}
	}
	my $mqc_outdir;
	if($root_dir ne ""){
		$mqc_outdir = $root_dir;
	}
	else{
		$mqc_outdir = $log_path;
	}
	if($verbose == 1){
		#print STDOUT "\e[K";
		print STDOUT $multiqc_path . " -o " . $mqc_outdir . " " . join(" ",keys(%multiqc_paths)) . "\n";
		print LOG $multiqc_path . " -o " . $mqc_outdir . " " . join(" ",keys(%multiqc_paths)) . "\n";
		system($multiqc_path . " -o " . $mqc_outdir . " " . join(" ",keys(%multiqc_paths)) . " 2>> " . $log);
	}
	else{
		system($multiqc_path . " -q  -o " . $mqc_outdir . " " . join(" ",keys(%multiqc_paths)) . " 2>> " . $log);
	}
}


close LOG;

sub insert_trimlog{
	my ($global_trim_opts,$in,$rename_switch) = @_;
	my $out = $in;
	$out =~ s/\.f(ast)?q(\.gz)?$//;
	if($rename_switch == 1){
		my @rename = split(/\//,$out);
		$rename[-1] = $rename[-2];
		$out = join("/",@rename);
	}
	if($global_trim_opts =~ m/^-trimlog$/){
		$global_trim_opts =~ s/^-trimlog$/-trimlog $out\_autotrim.trimlog/;
	}
	if($global_trim_opts =~ m/^-trimlog /){
		$global_trim_opts =~ s/^-trimlog /-trimlog $out\_autotrim.trimlog /;
	}
	if($global_trim_opts =~ m/ -trimlog$/){
		$global_trim_opts =~ s/ -trimlog$/ -trimlog $out\_autotrim.trimlog/;
	}
	if($global_trim_opts =~ m/ -trimlog /){
		$global_trim_opts =~ s/ -trimlog / -trimlog $out\_autotrim.trimlog /;
	}
	return $global_trim_opts;
}

sub run_fastqc{								#fastqc subroutine
	my ($out_ref,$kmers_ref,$status,$scalar_path_keys,$frac,$nok) = @_;	#give arguments to sub
	my @in_files = @{$out_ref};
	my $out1 = $in_files[0];
	$out1 =~ s/_autotrim.*\.f(ast)?q$//;
	if($verbose == 1){
		print STDOUT "\e[K";
		print STDOUT "$fastqc_path $fastqck-q -t " . scalar(@in_files) . " " . join(" ",@in_files) . "\n";			#fastqc creates output with _fastqc extension (_fastqc.zip and _fastqc.html)
		print LOG "$fastqc_path $fastqck-q -t " . scalar(@in_files) . " " . join(" ",@in_files) . "\n";
		print STDOUT "[autotrim] " . $status . " / " . $scalar_path_keys . " [" . $frac . "%] jobs trimmed\r";
	}
	system("$fastqc_path $fastqck-q -t " . scalar(@in_files) . " " . join(" ",@in_files));
	if($nok == 0){
		#read _fastqc.zip output
		my @fastqc_data_R1;
		my @fastqc_data_R2;
		my $r1pass = 0;
		my $r2pass = 0;
		if(scalar(@in_files == 1)){
			$r2pass = 1;
			my $fastqc_report = $in_files[0];
			$fastqc_report =~ s/\.f(ast)?q$/_fastqc.zip/;
			my $zip_root = (split /\//,$fastqc_report)[-1];
			$zip_root =~ s/\.zip$//;
			@fastqc_data_R1 = `unzip -p $fastqc_report $zip_root/fastqc_data.txt`;         #unzip the necessary fastqc-output to find kmers R1
		}
		else{
			my $fastqc_report_R1 = $in_files[0];
			$fastqc_report_R1 =~ s/\.f(ast)?q$/_fastqc.zip/;
			my $zip_root_R1 = (split /\//,$fastqc_report_R1)[-1];
			$zip_root_R1 =~ s/\.zip$//;
			@fastqc_data_R1 = `unzip -p $fastqc_report_R1 $zip_root_R1/fastqc_data.txt`;         #unzip the necessary fastqc-output to find kmers R1
			my $fastqc_report_R2 = $in_files[1];
			$fastqc_report_R2 =~ s/\.f(ast)?q$/_fastqc.zip/;
			my $zip_root_R2 = (split /\//,$fastqc_report_R2)[-1];
			$zip_root_R2 =~ s/\.zip$//;
			@fastqc_data_R2 = `unzip -p $fastqc_report_R2 $zip_root_R2/fastqc_data.txt`;         #unzip the necessary fastqc-output to find kmers R2
		}
		
		chomp @fastqc_data_R1;
		chomp @fastqc_data_R2;
	
		my $match1 = 0;
		my $match2 = 0;
		#R1
		for(my $j = 0; $j < scalar(@fastqc_data_R1); $j++){
			if ($match1 == 1){
				if($fastqc_data_R1[$j] =~ m/^[AGTC]/i){
					${$kmers_ref}{(split /\t/,$fastqc_data_R1[$j])[0]}=1;
				}
			}
			if($fastqc_data_R1[$j] =~ m/^>>Kmer Content/ and $fastqc_data_R1[$j] !~ m/\tpass$/){		#if trimming left overrepresented kmers
				#print $fastqc_data_R1[$j] . "\n";
				$match1 = 1;
			}
			if($fastqc_data_R1[$j] =~ m/^>>Kmer Content/ and $fastqc_data_R1[$j] =~ m/\tpass$/){		#if trimming was successful
				$r1pass = 1;
			}
		}
#		#R2
		for(my $j = 0; $j < scalar(@fastqc_data_R2); $j++){
			if ($match2 == 1){
				if($fastqc_data_R2[$j] =~ m/^[AGTC]/i){
					${$kmers_ref}{(split /\t/,$fastqc_data_R2[$j])[0]}=1;
				}
			}
			if($fastqc_data_R2[$j] =~ m/^>>Kmer Content/ and $fastqc_data_R2[$j] !~ m/\tpass$/){		#if trimming left overrepresented kmers
				$match2 = 1;
			}
			if($fastqc_data_R2[$j] =~ m/^>>Kmer Content/ and $fastqc_data_R2[$j] =~ m/\tpass$/){		#if trimming was successful
				$r2pass = 1;
			}
		}
	
		if($r1pass == 1 and $r2pass == 1){
			print STDOUT "\e[K";
			print STDOUT "[autotrim] " . ((split /\//,$in_files[0])[-2] . " trimmed successfully\n");
			print LOG "[autotrim] " . ((split /\//,$in_files[0])[-2] . " trimmed successfully\n");
			my $result = `tail -n 2 $out1\_autotrim.trimmomatic_err | head -n 1`;
			chomp $result;
			print STDOUT "[Trimmomatic] " . $result . "\n";
			print LOG "[Trimmomatic] " . $result . "\n";
		}
	
		if($match1 == 1 or $match2 == 1){
			open (KMERS, '>', $kfa);
			my @kmers_keys = keys(%kmers);
			for(my $k = 0; $k < scalar(@kmers_keys); $k++){					#save kmer-sequences in %kmers
				print KMERS ">kmer_$k\n";
				print KMERS $kmers_keys[$k] . "\n";					#write overrepresented kmer-sequences in fasta
			}
			close KMERS;
		}
		return ($match1,$match2);
	}
}

sub trimmomatic_loop{										#trimmomatic run2+ subroutine
	my ($trimmomatic_path,$in_files_ref,$global_trim_opts,$global_trimmer,$verbose,$status,$scalar_path_keys,$frac,$rename_switch) = @_;
	my @ktrimmer = split(/ /,$global_trimmer);
	for (my $i = scalar(@ktrimmer)-1; $i > -1; $i--){
		if($ktrimmer[$i] =~ m/^ILLUMINACLIP:/){
			splice(@ktrimmer,$i+1,0,"ILLUMINACLIP:$kfa:0:1:1:1:true");
			last;
		}
		if($i == 0 and $ktrimmer[$i] !~ m/^ILLUMINACLIP:/){
			splice(@ktrimmer,0,0,"ILLUMINACLIP:$kfa:0:1:1:1:true");
		}
	}
	my $k_trimmer = join(" ",@ktrimmer);

	my @in_files = @{$in_files_ref};
	if(scalar(@in_files) == 1){
		my $in = $in_files[0];
		my $out = $in;
		$out =~ s/\.f(ast)?q(\.gz)?$//;
		if($rename_switch == 1){
			my @rename = split(/\//,$out);
			$rename[-1] = $rename[-2];
			$out = join("/",@rename);
                }
		my $outfq = $out . "_autotrim.fq";
		if($verbose == 1){
			print STDOUT "\e[K";
			print STDOUT "$trimmomatic_path SE $global_trim_opts$in $outfq $k_trimmer > $out\_autotrim.trimmomatic_log 2> $out\_autotrim.trimmomatic_err\n";
			print LOG "$trimmomatic_path SE $global_trim_opts$in $outfq $k_trimmer > $out\_autotrim.trimmomatic_log 2> $out\_autotrim.trimmomatic_err\n";
			print STDOUT "[autotrim] " . $status . " / " . $scalar_path_keys . " [" . $frac . "%] jobs trimmed\r";
		}
		system("$trimmomatic_path SE $global_trim_opts$in $outfq $k_trimmer > $out\_autotrim.trimmomatic_log 2> $out\_autotrim.trimmomatic_err");
	}
	else{
		my $in1 = $in_files[0];
		my $in2 = $in_files[1];
		my $out1 = $in1;
		$out1 =~ s/\.f(ast)?q(\.gz)?$//;
		my $out1p;
		my $out1u;
		if($rename_switch == 1){
			my @rename1 = split(/\//,$out1);
			$rename1[-1] = $rename1[-2];
			$out1 = join("/",@rename1);
			$out1p = $out1 . "_1_autotrim.paired.fq";
			$out1u = $out1 . "_1_autotrim.unpaired.fq";
			$out1 = $out1 . "_1";
                }
		else{
			$out1p = $out1 . "_autotrim.paired.fq";
			$out1u = $out1 . "_autotrim.unpaired.fq";
		}
		my $out2 = $in2;
		$out2 =~ s/\.f(ast)?q(\.gz)?$//;
		my $out2p;
		my $out2u;
		if($rename_switch == 1){
			my @rename2 = split(/\//,$out2);
			$rename2[-1] = $rename2[-2];
			$out2 = join("/",@rename2);
			$out2p = $out2 . "_2_autotrim.paired.fq";
			$out2u = $out2 . "_2_autotrim.unpaired.fq";
                }
		else{
			$out2p = $out2 . "_autotrim.paired.fq";
			$out2u = $out2 . "_autotrim.unpaired.fq";
		}
		if($verbose == 1){
			print STDOUT "\e[K";
			print STDOUT "$trimmomatic_path PE $global_trim_opts$in1 $in2 $out1p $out1u $out2p $out2u $k_trimmer > $out1\_autotrim.trimmomatic_log 2> $out1\_autotrim.trimmomatic_err\n";
			print LOG "$trimmomatic_path PE $global_trim_opts$in1 $in2 $out1p $out1u $out2p $out2u $k_trimmer > $out1\_autotrim.trimmomatic_log 2> $out1\_autotrim.trimmomatic_err\n";
			print STDOUT "[autotrim] " . $status . " / " . $scalar_path_keys . " [" . $frac . "%] jobs trimmed\r";
		}
		system("$trimmomatic_path PE $global_trim_opts$in1 $in2 $out1p $out1u $out2p $out2u $k_trimmer > $out1\_autotrim.trimmomatic_log 2> $out1\_autotrim.trimmomatic_err");
	}
}

exit;

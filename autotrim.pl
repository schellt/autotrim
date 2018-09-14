#!/usr/bin/perl

use strict;
use warnings;

opendir (DIR, $ARGV[0]) or die "Could not open root-directory!" . "\n";	#opens the root dir

my @fulldirlist;
my @dirlist;
my @pathlist;
my %filelist;
my %kmers;
#my $kmers_ref = \%kmers;

while (my $subdir = readdir(DIR)) {					#reads the root dir
	next if ($subdir =~ m/^\./);					#ignores "." and ".." (dir up)
	next unless (-d "$ARGV[0]/$subdir");				#returns folders and ignores files in the root dir
	push (@dirlist,$subdir);	
	push (@fulldirlist,$ARGV[0] . $subdir . "/");			#saves the full path of all found sub dirs in an array
}

for (my $i = 0; $i < scalar(@fulldirlist); $i++){
	opendir (SUBDIR, $fulldirlist[$i]) or die "Could not open sub-directory " . $fulldirlist[$i] . "!" . "\n";	#opens every found sub dir
	my @arr = ();
	while (my $file = readdir(SUBDIR)) {				#reads the sub dir
		next unless (-f "$fulldirlist[$i]/$file");		#returns only files from the sub dir
		if ($file =~ m/_[1-2]\./){				#ignores files that don't match the pattern (end with .fq)
			push (@pathlist,$fulldirlist[$i] . $file);	#saves the full path of all found files in an array
			$file =~ s/\.fq$//;				#removes the file extension
			push (@arr,$file);				#saves all filenames without extension from a sub dir into an array
		}
	$filelist{$dirlist[$i]} = join("\t",@arr);			#saves all files from a sub dir into a hash
	}
}

#print "dirlist:\n";
#print join("\n",@dirlist) . "\n";
#print "\n";
print "folders with containing files:\n";
for (my $i = 0; $i < scalar(@dirlist); $i++){
	print $dirlist[$i] . "\t" . $filelist{$dirlist[$i]} . "\n";
}
#print "\n";
#print "fulldirlist:\n";
#print join("\n",@fulldirlist) . "\n";
#print "\n";

print "pathlist:\n";
print join("\n",@pathlist) . "\n";
print "\n";



for (my $i = 0; $i < scalar(@pathlist); $i = $i+2) {		#take corresponding forward and reverse reads for trimming from list
	#print $pathlist[$i] . "\t" . $pathlist[$i+1] . "\n";
	my $R1paired = $ARGV[0] . (split /\//,$pathlist[$i])[-2] . "/" . (split /\//,$pathlist[$i])[-2] . "_R1.paired.fastq_true";
	my $R1unpaired = $ARGV[0] . (split /\//,$pathlist[$i])[-2] . "/" . (split /\//,$pathlist[$i])[-2] . "_R1.unpaired.fastq_true";
	my $R2paired = $ARGV[0] . (split /\//,$pathlist[$i])[-2] . "/" . (split /\//,$pathlist[$i+1])[-2] . "_R2.paired.fastq_true";
	my $R2unpaired = $ARGV[0] . (split /\//,$pathlist[$i])[-2] . "/" . (split /\//,$pathlist[$i+1])[-2] . "_R2.unpaired.fastq_true";
#start first trimmomatic run
	#print ("trimmomatic.sh PE -phred33 $pathlist[$i] $pathlist[$i+1] $R1paired $R1unpaired $R2paired $R2unpaired ILLUMINACLIP:/home/amoppold/Ind_Reseq/adapter_2.0.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 TOPHRED33\n");
	system ("trimmomatic PE -threads 30 -phred33 $pathlist[$i] $pathlist[$i+1] $R1paired $R1unpaired $R2paired $R2unpaired ILLUMINACLIP:/home/tilman/results/2015-02-19_adapter/adapter_combined.fa:2:30:10:8:true SLIDINGWINDOW:4:20 MINLEN:50 TOPHRED33");

#start fastqc-trimming-loop
	my ($match1,$match2) = run_fastqc($R1paired,$R2paired,$ARGV[0],$pathlist[$i],\%kmers);		#call fastqc subroutine
	while ($match1 == 1 or $match2 == 1){									#if overrepresented kmers are left after trimming $match=1 and the trimming has to be repeated
		trimmomatic_loop($pathlist[$i],$pathlist[$i+1],$R1paired,$R1unpaired,$R2paired,$R2unpaired,$ARGV[0]);	#call trimmomatic subroutine
		($match1,$match2) = run_fastqc($R1paired,$R2paired,$ARGV[0],$pathlist[$i],\%kmers);
	}
	if($match1 == 0 and $match2 == 0){
		run_fastqc($R1paired,$R2paired,$ARGV[0],$pathlist[$i],\%kmers);
	}
	#print join ("\n",keys(%kmers)) . "\n";
	#print $fastqc_report_R1 . "\n";
	#print $zip_root_R1 . "\n";
	#print join("\n", @fastqc_data_R1) . "\n";

}

sub run_fastqc{								#fastqc subroutine
	my ($R1paired,$R2paired,$root_path,$path,$kmers_ref) = @_;	#give arguments to sub
	#print ("fastqc $R1paired $R2paired\n");			#fastqc creates output with _fastqc extension (_fastqc.zip and _fastqc.html)
	system ("fastqc -t 2 $R1paired $R2paired");
#read _fastqc.zip output
	my $fastqc_report_R1 = $R1paired . "_fastqc.zip";		#input R1.zip
	my $fastqc_report_R2 = $R2paired . "_fastqc.zip";		#input R2.zip

	my $zip_root_R1 = (split /\//,$fastqc_report_R1)[-1];
	$zip_root_R1 =~ s/\.zip$//;
	my @fastqc_data_R1 = `unzip -p $fastqc_report_R1 $zip_root_R1/fastqc_data.txt`;		#unzip the necessary fastqc-output to find kmers R1
	
	my $zip_root_R2 = (split /\//,$fastqc_report_R2)[-1];
	$zip_root_R2 =~ s/\.zip$//;
	my @fastqc_data_R2 = `unzip -p $fastqc_report_R2 $zip_root_R2/fastqc_data.txt`;		#unzip the necessary fastqc-output to find kmers R2
	
	chomp @fastqc_data_R1;
	my $match1 = 0;
	my $match2 = 0;
	#R1
	for(my $j = 0; $j < scalar(@fastqc_data_R1); $j++){
		if ($match1 == 1){
			#print $fastqc_data_R1[$j] . "\n";
			if($fastqc_data_R1[$j] =~ m/^[AGTC]/i){
				${$kmers_ref}{(split /\t/,$fastqc_data_R1[$j])[0]}=1;
			}
		}
		if($fastqc_data_R1[$j] =~ m/^>>Kmer Content/ and $fastqc_data_R1[$j] !~ m/\tpass$/){		#if trimming left overrepresented kmers
			#print $fastqc_data_R1[$j] . "\n";
			$match1 = 1;
		}
		if($fastqc_data_R1[$j] =~ m/^>>Kmer Content/ and $fastqc_data_R1[$j] =~ m/\tpass$/){		#if trimming was successful
			print ((split /\//,$path)[-2] . " trimmed successfully \n");
		}
	}
	#R2
	for(my $j = 0; $j < scalar(@fastqc_data_R2); $j++){
		if ($match2 == 1){
			#print $fastqc_data_R1[$j] . "\n";
			if($fastqc_data_R2[$j] =~ m/^[AGTC]/i){
				${$kmers_ref}{(split /\t/,$fastqc_data_R2[$j])[0]}=1;
			}
		}
		if($fastqc_data_R2[$j] =~ m/^>>Kmer Content/ and $fastqc_data_R2[$j] !~ m/\tpass$/){		#if trimming left overrepresented kmers
			#print $fastqc_data_R1[$j] . "\n";
			$match2 = 1;
		}
		if($fastqc_data_R2[$j] =~ m/^>>Kmer Content/ and $fastqc_data_R2[$j] =~ m/\tpass$/){		#if trimming was successful
			print ((split /\//,$path)[-2] . " trimmed successfully \n");
		}
	}
	#foreach(keys(%{$kmers_ref})){
		#print $_ . "\n";
	#}

	if($match1 == 1 or $match2 == 1){
		open (KMERS, '>', $root_path . "/kmers.fa");
		my @kmers_keys = keys(%kmers);
		for(my $k = 0; $k < scalar(@kmers_keys); $k++){					#save kmer-sequences in %kmers
			print KMERS ">kmer_$k\n";
			print KMERS $kmers_keys[$k] . "\n";					#write overrepresented kmer-sequences in fasta
		}
		close KMERS;
	}
	return ($match1,$match2);

}

sub trimmomatic_loop{										#trimmomatic run2 subroutine
	my ($pathlistR1,$pathlistR2,$R1paired,$R1unpaired,$R2paired,$R2unpaired,$path) = @_;
	#print ("trimmomatic.sh PE -phred33 $pathlistR1 $pathlistR2 $R1paired $R1unpaired $R2paired $R2unpaired ILLUMINACLIP:/home/amoppold/Ind_Reseq/adapter_2.0.fa:2:30:10:8:true ILLUMINACLIP:$path\kmers.fa:0:1:1:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 TOPHRED33\n");
	system ("trimmomatic PE -threads 30 -phred33 $pathlistR1 $pathlistR2 $R1paired $R1unpaired $R2paired $R2unpaired ILLUMINACLIP:/home/tilman/results/2015-02-19_adapter/adapter_combined.fa:2:30:10:8:true ILLUMINACLIP:$path\kmers.fa:0:1:1:1:true SLIDINGWINDOW:4:20 MINLEN:50 TOPHRED33");
}

exit;

#README/usage
# -The pipeline, parameters and resources are based on appropriate best practices according to authors' recommendations (e.g. Broad Institute for GATK).
#	-There may be deviations from best practices where practical.
# -When specifying fastq directories, the entry should include a preceding tag containing the run-specific details (e.g. 180213_K00267_0143_AHMYNGBBXX::/scratch/c1154021/fastq). These are necessary for red group creation.
#	-The tag should be separated from the directory with two colons ('::').
#	-Multiple directories can be given, separated by a semi colon (';). A meta-data tag will be needed for each dir.
#	-Subdirectories within the directories are also searched for fastq files.

#!/usr/bin/perl -w
# use strict;
# use warnings;

use Getopt::Std;
use Term::ANSIColor;

our (%options, %tissues);

#Get the switches from the command line. ':' after an option specifies an argument is needed, absence of ':' specifies a boolean flag.
#In this case, I've used capital letters to represent required parameters.
getopts("a:b:c:d:fg:h:I:J:j:k:m:n:O:o:p:q:r:s:t:u:v:", \%options);

$options = parameter_input(\%options);

%options= %{$options};

if ($options{"k"}){$tissues = tissue_pairing(\%options);}

#FASTQ file inputs
#Specify the input directory containing fastq files. Search all subdirectories and find all FQ files within the parent dir.
#From the FQ files, establish sample names and group together paired-end reads and files from different lanes if necessary.
#       e.g. if sampleX was sequenced on lane 6 and 7, and has xxx_L6_xxx_1.fq, xxx_L6_xxx_2.fq, xxx_L7_xxx_1.fq and xxx_L7_xxx_2.fq,
#       then sample names for downstream files should be sampleX_L6_1, sampleX_L6_2, etc. If files from several lanes are to be merged for a sample,
#       then sample names after merging should resemble sampleX_merged_1, sampleX_merged_2, etc. After reads 1 & 2 are merged (usually at BAM creation),
#       then sample names should resemble sampleX_merged, or sampleX_L3, etc.
#our $pattern = "\.f(ast)?q(\.gz)?\@?\$";
our $pattern = "f(ast)?q(.gz)?\$";
our @fastq_dirs = split(/,/, $options{"I"});

our %fq_files;

#Find fastq files and return array of filepaths.
foreach my $fastq_dir (sort {"\L$a" cmp  "\L$b"} @fastq_dirs){
	my ($sub_dirs, $fq_files) = read_dir($fastq_dir, $pattern, \%fq_files);
	%fq_files = %{$fq_files};print keys %fq_files;
}

our $outdir= $options{"O"};

our %output_dirs = (
		mapping => "mapping",
		samtools_processing => "SAMtools_processing",
		samtools_mpileup => "SAMtools_mpileup",
		indel_realignment => "indel_realignment",
		BQSR => "BQSR",
		variants => "GATK_haplotypecaller_variants",
		variant_annotation => "variant_annotation",
		coverage => "SAMtools_mpileup",
		metrics => "bam_metrics",
		VQSR => "VQSR",
		Log => "Log",
		FastQC => "FastQC",
		spreadsheets => "spreadsheets",
		somatic_variants => "somatic_variants"
);

for my $dir (sort {"\L$a" cmp "\L$b"} keys %output_dirs){$output_dirs{$dir} = $outdir."/".$output_dirs{$dir};}

create_outdirs(\%output_dirs);

our $outfiles = establish_out_names(\%fq_files, \%output_dirs, \%options, $tissues);

our $log_file = captains_log(\%options);

our $commands = generate_commands(\%options, $outfiles, \%output_dirs, $tissues);

execution($commands, $outfiles, $log_file);

########################################################################################
#Subroutine parameter_input
#Obtains the necessary parameter choices for the pipeline, along with filepaths of supporting files 
#(reference genome fasta file, database VCFs, etc).

#Requires: %options hash containing command line parameters.
#Returns: %options hash.
########################################################################################
sub parameter_input{
	my $options = shift;
	
	my %options = %{$options};
	
	&usage unless (exists $options{"I"} && exists $options{"J"});
	
	#$options{"a"} = "/scratch/c1154021/bioinformatics_software/Trimmomatic-0.32/adapters/all_trimmomatic_adapter_sequences.fa" unless exists $options{"a"};
	$options{"b"} = "/wgp1/wgp/resources/human/targets/SureSelect_All_Exon_50mb_with_annotation_hg19.bed" unless exists $options{"b"};
	$options{"c"} = 1 unless exists $options{"c"};
	$options{"d"} = "/scratch/c1154021/db_resources/dbsnp_150_GRCh37p13.vcf" unless exists $options{"d"};
	$options{"e"} = "/scratch/c1154021/bioinformatics_software/snpEff_4_2/snpEff.jar";
	$options{"g"} = "/scratch/c1154021/bioinformatics_software/gatk-4.0.3.0/gatk-package-4.0.3.0-local.jar" unless exists $options{"g"};
	#$options{"h"} = "/scratch/c1154021/db_resources/VCF/hapmap_3.3.b37.sites.vcf" unless exists $options{"h"};
	$options{"i"} = "/scratch/c1154021/db_resources/af-only-gnomad.raw.sites.b37.vcf" unless exists $options{"i"};
	$options{"j"} = "/share/apps/picard-tools/MergeSamFiles.jar" unless exists $options{"j"};
	$options{"m"} = "/scratch/c1154021/db_resources/Mills_and_1000G_gold_standard.indels.b37_UCSC-reference.vcf" unless exists $options{"m"};
	#$options{"n"} = "/home/marc/Scripts/seq_stat_new.pl" unless exists $options{"n"};
	#$options{"o"} = "/scratch/c1154021/db_resources/VCF/1000G_omni2.5.b37.sites.vcf " unless exists $options{"o"};
	$options{"O"} = "/scratch/c1154021/Projects/$options{J}" unless exists $options{"O"};
	#$options{"p"} = "/scratch/c1154021/db_resources/VCF/1000G_phase1_release_v3.2001123.snps_indels_sv.sites.vcf" unless exists $options{"p"};
	#$options{"q"} = "/scratch/c1154021/Projects/common/db_resources/1000G_phase1.indels.b37_UCSC-reference.vcf" unless exists $options{"q"};
	$options{"r"} = "/scratch/c1154021/reference_genome/human_g1k_v37.fasta" unless exists $options{"r"};
	$options{"s"} = "ILLUMINA" unless exists $options{"s"};
	#$options{"t"} = "/scratch/c1154021/bioinformatics_software/Trimmomatic-0.32/trimmomatic-0.32.jar" unless exists $options{"t"};
	#$options{"u"} = "/scratch/c1154021/bioinformatics_software/muTect-1.1.4/muTect-1.1.4.jar" unless exists $options{"u"};
	
	#Options without defaults:
	#b,f,I,J,v.
	
	$options{"e2"} = $options{"e"};
	$options{"e2"} =~ s/\.jar$/\.config/i;
	
	print color ("yellow"), "Parameters:\n";
	
	for my $option (sort {"\L$a" cmp "\L$b"} keys %options){
		print color ("magenta"), "\t-$option\t$options{$option}\n";
	}
	
	print color 'reset';
	
	#other options used by steps:
	#bwa aln: -t, 
	#samtolls rmdup -s/-S
	#BaseRecalibrator -l INFO [-cov cyclecovariate is not needed]
	#HaplotypeCaller -rf BadCigar
	
	return \%options;
}	
########################################################################################
#End of subroutine create_outdirs.
########################################################################################

########################################################################################
#Subroutine usage
#Prints a message warning of usage if users do not supply necessary command line parameters.
#This subroutine then kills the script's execution.
########################################################################################
sub usage{
	die "Usage: perl exome_pipline.pl \\\n\t -I [input fastq dir] \\\n\t -o [output dir] \\\n\t -J [job name] \\\n\t -r [ref genome ] \\\n\t -t [trimmomatic jar] \\\n\t -g [GATK jar] \\\n\t -a [sequencing adapter file] \\\n\t -c [# cores to use] \\\n\t -d [dbSNP VCF] \\\n\t -O [1000G omni2.5 VCF] \\\n\t -p [1000G phase 1 VCF] \\\n\t -h [hapmap VCF] \\\n\t -m [Mills & Devine VCF] \\\n\t -b [BED file for regional calling] \\\n\t -P [sequencing platform] \\\n\t -k [normal:tumour, normal:-, -:tumour...]\n";
}
########################################################################################
#End of subroutine create_outdirs.
########################################################################################

########################################################################################
#Subroutine read_dir
#Locates files wth names containing a pattern within a specified directory.
#Subdirectories are also searched. 
#Requires: the parent directory, the pattern.
#Returns: an array (reference) of subdirectories, an array (reference) of matched files.

#Caveats
#The subroutine expects files to be contained within the /data/results/xx/inputs directory on wotan (or other server).
#Expected subdirectories are patient-ID/sequencer-lane/run-date. Regexes have been designed to this specification;
#changing the filesystem may therefore require updating the regexes accordingly.
########################################################################################
sub read_dir{
	#get the directory, patterns and files. Files will only exist if the subroutine is calling itself and not when first called by the main program.
	my ($dir, $pattern, $fastq_files) = @_;
	
	print "Reading $dir\n";
	
	#get run meta data
	my $run = $dir;
	$run =~ s/\:\:.+$//;
	my ($run_date, $sequencer, $run_no, $flowcell) = split (/_/, $run);
	
	$dir =~ s/^.+\:\://;print "run-dir = $run\t$dir\n";
	
	my %fastq_files = %{$fastq_files};
	
	#open the directory.
	opendir DIR, "$dir" or die "cannot open $dir!\n";
	
	#Read in all files from the directory (excluding the "." and ".." files which are present), and store the full filepaths in an array.
	my @filepaths = map{$dir."/".$_} grep {"$dir/$_" && !/^\.{1,2}$/} readdir DIR or die "crash on readdir\n";
	
	closedir DIR;

	#foreach file and sub directory within the current directory.
	foreach my $filepath (sort {"\L$a" cmp "\L$b"} @filepaths){		
		#If the current file matches a certain criteria ('-f' specifies a file).
		if (-f "$filepath" && $filepath =~ /$pattern/){
			my $id = $filepath;
			my $lane;
			
			$id =~ s/^.+\///g;
			
			if ($id =~ /\_L\d+\_/){
				$lane = $&;
				$lane =~ s/\_//g;
				$lane =~ s/L//;
			}
			
			else{$lane = "Unknown";}
			
			$id =~ s/\_.+$//;
			
			push @{$fastq_files{$id}{$run_date}{$sequencer}{$flowcell}{$lane}{"fastq"}}, $filepath;
			print "$id\t$lane\t$run\t$filepath\n";
			
		}
		
		#if the directory contains sub directories, store them in the array and call this sub routine on them.
		#the result is each sub directory within a directory will be searched before moving onto the next directory;
		#e.g. a/1/i and a/1/ii will be searched before a/2/i, etc.
		if (-d $filepath){print "nest run-dir: $run\t$dir\n";
			my ($returned_array, $fastq_files) = read_dir($run."::".$filepath, $pattern, \%fastq_files);
			
			%fastq_files = %{$fastq_files};
			
			my @returned_array = @$returned_array;
			
			push @filepaths, @returned_array;
		}
	}
	
	return (\@filepaths, \%fastq_files);
}
########################################################################################
#End of subroutine read_dir.
########################################################################################

sub tissue_pairing{
	my $options = shift;
	
	my %options = %{$options};
	
	my $pairings = $options{"k"};
		
	my @pairings = split (/,/, $pairings);
	
	my %tissues;
		
	foreach my $pairing (sort {"\L$a" cmp "\L$b"} @pairings){
		my @tissues = split (/:/, $pairing);
		
		my ($germline, $tumour) = @tissues;
		
		#have a 'paired tumour' for every 'germline' sample necessary.
		#When looping over subject ID, generate a paired-tumour entry for the germline ID and necessary MuTect output files.
		#Then generate MuTect commands and annotation steps as neccessary.
		$tissues{$tumour} = $germline;
	}
	
	return \%tissues;
}

########################################################################################
#Subroutine establish_out_names.
#For each patient/subject, take metadata from the fastq files (patient ID, lane, run, etc) and generate a hash
#containing filenames used in downstream analyses. As some patients have multiple samples
#(e.g. sequenced on different lanes or in different runs), the BAM files for these can be merged.
#The earliest stage this can occur is after mapping as mapped reads should include read group data,
#allowing for tracing of reads to original files.

#Requires: the %fq_files hash
#Returns: A hash containing downstream output file names.
########################################################################################
sub establish_out_names{
	my ($fq_files, $output_dirs, $options, $tissues) = @_;
	
	my %fq_files = %{$fq_files};
	my %output_dirs = %{$output_dirs};
	my %options = %{$options};
	my %tissues = %{$tissues};
	
	#check the names of the returned files.
	for my $id (sort {"\L$a" cmp "\L$b"} keys %fq_files){
		my $id_merged = $id."_merged";
		
		push @{$fastq_files{$id}{$run_date}{$sequencer}{$flowcell}{$lane}{"fastq"}}, $filepath;
		
		for my $run_date (sort {"\L$a" cmp "\L$b"} keys %{$fq_files{$id}}){
			for my $sequencer (sort {"\L$a" cmp "\L$b"} keys %{$fq_files{$id}{$run_date}}){
				for my $flowcell (sort {"\L$a" cmp "\L$b"} keys %{$fq_files{$id}{$run_date}{$sequencer}}){
					for my $lane (sort {"\L$a" cmp "\L$b"} keys %{$fq_files{$id}{$run_date}{$sequencer}{$flowcell}}){
						my $sample = $id."_".$run_date."_".$sequencer."_".$flowcell."_".$lane;
						
						#The number of FQ files for this sample. Skip if != 2.
						if (scalar @{$fq_files{$id}{$run_date}{$sequencer}{$flowcell}{$lane}{"fastq"}} != 2){
							print color ("red"), "Unexpected number of FQ files for sample $sample!\nSkipping sample.\n";
							print color 'reset';
							
							next;
						}
						
						%RG_header = (
							$sample => {
								ID => $sample,		#ID should be unique for every sequencing run & flowcell (the source of a group of reads). Present in RG header and as a tag in reads.
								PU => $sample,		#PU (Platform Unit) contains sequencing run, flowcell and sample/library info. The most precise object in an RG header.
											#Takes precedence over ID for BQSR, but is not essential (especially if you make ID as precise as PU by including sample info).
								SM => $id,		#The subject (e.g. patient_X, mouse_1, cell_A). Used as the sample column header in VCFs.
								Lane => $lane,		#Lane of the flowcell
								PL => $options{"s"},	#Sequencing platform (either ILLUMINA, SOLID, LS454, HELICOS or PACBIO).
								LB => $sample,		#If you have a library ID, you can use it. Otherwise use the sample ID (since this should be unique enough to identify the library).
											#Used by MarkDuplicates (and rmdup?)
								DT => $run_date		#Date of sequencing run.
							}
						);
				
						#The read group headline, implementing sample meta-data.
						$outfiles{$id}{$sample}{"RG_headline"} = "\\'\@RG\\\\tID:$RG_header{$sample}{\"ID\"}\\\\tSM:$RG_header{$sample}{\"SM\"}\\\\tLB:$RG_header{$sample}{\"LB\"}\\\\tLane:$RG_header{$sample}{\"Lane\"}\\\\tDT:$RG_header{$sample}{\"DT\"}\\\\tPL:$RG_header{$sample}{\"PL\"}\\'";
						
						#Filenames for files generated throughout the pipeline are created and held within a hash.
						#This is advantageous when handling multiple samples as the pipeline commands can be constructed with references to the hash,
						#and filenames do not have to be explicitly specified.
						$outfiles{$id}{$sample}{"R1_fq"} = ${$fq_files{$id}{$run_date}{$sequencer}{$flowcell}{$lane}{"fastq"}}[0];
						$outfiles{$id}{$sample}{"R2_fq"} = ${$fq_files{$id}{$run_date}{$sequencer}{$flowcell}{$lane}{"fastq"}}[1];
						$outfiles{$id}{$sample}{"R1_aln_sai"} = "$output_dirs{\"mapping\"}/$sample\_R1_aln.sai";
						$outfiles{$id}{$sample}{"R2_aln_sai"} = "$output_dirs{\"mapping\"}/$sample\_R2_aln.sai";
						$outfiles{$id}{$sample}{"sam"} = "$output_dirs{\"mapping\"}/$sample\.sam";
						$outfiles{$id}{$sample}{"view_bam"} = "$output_dirs{\"mapping\"}/$sample.bam";
						
						$outfiles{$id}{"merged"}{"starting_bam"} = $fq_files{$id}{$lane}{$run}{"fastq"};
						$outfiles{$id}{"merged"}{"header_sam"} = $outfiles{$id}{"merged"}{"starting_bam"};
						$outfiles{$id}{"merged"}{"header_sam"} =~ s/.bam/_header.sam/;
						$outfiles{$id}{"merged"}{"reheader_bam"} = $outfiles{$id}{"merged"}{"header_sam"};
						$outfiles{$id}{"merged"}{"reheader_bam"} =~ s/_header.sam/_reheadered.bam/;
						$outfiles{$id}{"merged"}{"reorder_bam"} = $outfiles{$id}{"merged"}{"reheader_bam"};
						$outfiles{$id}{"merged"}{"reorder_bam"} =~ s/_reheadered.bam/_reordered.bam/;
						
						$outfiles{$id}{"merged"}{"cat_rg_headers"} = "$output_dirs{\"mapping\"}/$id_merged\_RG_headers.txt";
						$outfiles{$id}{"merged"}{"merged_bam"} = "$output_dirs{\"mapping\"}/$id_merged.bam";
						$outfiles{$id}{"merged"}{"sorted_bam"} = "$output_dirs{\"samtools_processing\"}/$id_merged\_sorted";
						$outfiles{$id}{"merged"}{"fixmate_bam"} = "$output_dirs{\"samtools_processing\"}/$id_merged\_fixmate.bam";
						$outfiles{$id}{"merged"}{"post_fixmate_sorted_bam"} = "$output_dirs{\"samtools_processing\"}/$id_merged\_post_fixmate_sorted";
						$outfiles{$id}{"merged"}{"nodups_bam"} = "$output_dirs{\"samtools_processing\"}/$id_merged\_nodups.bam";
						$outfiles{$id}{"merged"}{"indels_intervals"} = "$output_dirs{\"indel_realignment\"}/$id_merged\_indels.intervals";
						$outfiles{$id}{"merged"}{"realigned_bam"} = "$output_dirs{\"indel_realignment\"}/$id_merged\_realigned.bam";
						$outfiles{$id}{"merged"}{"mpileup_bam"} = "$output_dirs{\"samtools_mpileup\"}/$id_merged\_mpileup_bam";
						$outfiles{$id}{"merged"}{"recalibrated_grp"} = "$output_dirs{\"BQSR\"}/$id_merged\_recalibrated.grp";
						$outfiles{$id}{"merged"}{"post_recalibrated_bam"} = "$output_dirs{\"BQSR\"}/$id_merged\_post_recalibrated.bam";
						$outfiles{$id}{"merged"}{"BQSR_plots"} = "$output_dirs{\"BQSR\"}/$id_merged\_BQSR_plots.pdf";
						$outfiles{$id}{"merged"}{"recalibrated_bam"} = "$output_dirs{\"BQSR\"}/$id_merged\_recalibrated.bam";
						$outfiles{$id}{"merged"}{"genotype_calls"} = "$output_dirs{\"variants\"}/$id_merged\HC_genotyping_variants.gvcf";
						$outfiles{$id}{"merged"}{"variant_calls"} = "$output_dirs{\"variants\"}/$id_merged\_haplotypecaller_variants.vcf";
						$outfiles{$id}{"merged"}{"SNPs_only"} = "$output_dirs{\"variants\"}/$id_merged\_SNPs_only.vcf";
						$outfiles{$id}{"merged"}{"indels_only"} = "$output_dirs{\"variants\"}/$id_merged\_indels_only.vcf";
						$outfiles{$id}{"merged"}{"SNPs_recal"} = "$output_dirs{\"VQSR\"}/$id_merged\_SNPs_filtered.csv";
						$outfiles{$id}{"merged"}{"indels_recal"} = "$output_dirs{\"VQSR\"}/$id_merged\_indels_filtered.csv";
						$outfiles{$id}{"merged"}{"SNPs_tranches"} = "$output_dirs{\"VQSR\"}/$id_merged\_SNPs.tranches";
						$outfiles{$id}{"merged"}{"indels_tranches"} = "$output_dirs{\"VQSR\"}/$id_merged\_indels.tranches";
						$outfiles{$id}{"merged"}{"SNPs_recalibrated"} = "$output_dirs{\"VQSR\"}/$id_merged\_SNPs_recalibrated.vcf";
						$outfiles{$id}{"merged"}{"indels_recalibrated"} = "$output_dirs{\"VQSR\"}/$id_merged\_indels_recalibrated.vcf";
						$outfiles{$id}{"merged"}{"variants_recalibrated"} = "$output_dirs{\"VQSR\"}/$id_merged\_variants_recalibrated.vcf";
						$outfiles{$id}{"merged"}{"variants_annotated"} = "$output_dirs{\"variant_annotation\"}/$id_merged\_variants_annotated.vcf";
						
						$outfiles{$id}{"merged"}{"snpeff_annotation"} = "$output_dirs{\"spreadsheets\"}/$id_merged\_snpeff_annotation.vcf";
						$outfiles{$id}{"merged"}{"final_properties"} = "$output_dirs{\"spreadsheets\"}/$id_merged\_final.properties";
						$outfiles{$id}{"merged"}{"final_txt"} = "$output_dirs{\"spreadsheets\"}/$id_merged\_final.txt";
						$outfiles{$id}{"merged"}{"final_xlsx"} = "$output_dirs{\"spreadsheets\"}/$id_merged\_final.xlsx";
						$outfiles{$id}{"merged"}{"headers_file"} = "$output_dirs{\"spreadsheets\"}/$id_merged\_headers_list.txt";
						
						if (exists $tissues{$id}){
							$outfiles{$id}{"merged"}{"somatic_variant_calls"} = "$output_dirs{\"somatic_variants\"}/$id_merged\_somatic_variants.vcf";
							$outfiles{$id}{"merged"}{"coverage_wig"} = "$output_dirs{\"somatic_variants\"}/$id_merged\_mutect_coverage.wig";
							$outfiles{$id}{"merged"}{"somatic_snpeff_annotation"} = "$output_dirs{\"spreadsheets\"}/$id_merged\_somatic_snpeff_annotation.vcf";
							$outfiles{$id}{"merged"}{"somatic_variants_annotated"} = "$output_dirs{\"variant_annotation\"}/$id_merged\_somatic_variants_annotated.vcf";
						}
					}							
				}
			}
		}
	}
	
	return \%outfiles;
}
########################################################################################
#End of subroutine establish_out_names.
########################################################################################

########################################################################################
#Subroutine create_outdirs.
#Create directories for storing output files.

#Requires: the parent output directory, and the array of output directory names.
########################################################################################
sub create_outdirs{
	my ($output_dirs) = @_;

	my %output_dirs = %{$output_dirs};

	for my $output_type (sort {"\L$a" cmp "\L$b"} keys %output_dirs){
		my $output_dir = $output_dirs{$output_type};
		
		unless (-d "$output_dir"){
			print color ("red"), "Making $output_dir\n";
			print color 'reset';
			
			mkdir "$output_dir";

			#print a warning message and exit if the directory could not be created.
			unless (-d "$out_dir/$output_dir"){
				die "Directory $out_dir/$output_dir could not be created. Exiting.\n";
			}
		}
	}
}
########################################################################################
#End of subroutine create_outdirs.
########################################################################################

########################################################################################
#Subroutine generate_commands.
#Generates an array of commands stored in a hash for each patient/subject.

#Requires: the %options hash, the %outfiles hash.
########################################################################################
sub generate_commands{
	my ($options, $outfiles, $output_dirs, $tissues) = @_;
	
	my %options = %{$options};
	my %outfiles = %{$outfiles};
	my %output_dirs = %{$output_dirs};
	my %tissues = %{$tissues};
	
	my %commands;
	
	#Create commands for the initial mapping of each individual sample.
	for my $id (sort {"\L$a" cmp "\L$b"} keys %outfiles){
		#Bam count is used to determined whether there are multiple BAMs for merging or not.
		my $bam_count = 0;
		
		#Generate mappng commands, and commands for removing intermediate files (i.e. SAM and SAI files).
		for my $sample (sort {"\L$a" cmp "\L$b"} keys %{$outfiles{$id}}){
			next if $sample eq "merged";
			
			$bam_count++;
			
			my $ref_sampe = $options{"r"};
			#$ref_sampe =~ s/\.fasta$//i;
			
			#=====================================================================#
			push @{$commands{$id}},
				"bwa aln -t $options{c} \"$options{r}\" \"$outfiles{$id}{$sample}{R1_fq}\" > \"$outfiles{$id}{$sample}{R1_aln_sai}\"",
				"bwa aln -t $options{c} \"$options{r}\" \"$outfiles{$id}{$sample}{R2_fq}\" > \"$outfiles{$id}{$sample}{R2_aln_sai}\"",
				"bwa sampe -r \"$outfiles{$id}{$sample}{RG_headline}\" \"$ref_sampe\" \"$outfiles{$id}{$sample}{R1_aln_sai}\" \"$outfiles{$id}{$sample}{R2_aln_sai}\" \"$outfiles{$id}{$sample}{R1_fq}\" \"$outfiles{$id}{$sample}{R2_fq}\" > \"$outfiles{$id}{$sample}{sam}\"",
				"samtools view -bSh \"$outfiles{$id}{$sample}{sam}\" > \"$outfiles{$id}{$sample}{view_bam}\"";
				"rm \"$outfiles{$id}{$sample}{R1_aln_sai}\"",
				"rm \"$outfiles{$id}{$sample}{R2_aln_sai}\"",
				"rm \"$outfiles{$id}{$sample}{sam}\"";
			#=====================================================================#
				
			#Perform a fastqc analysis on the mapped bam file if wanted.
			if ($options{"f"}){
				unless (-d "$output_dirs{FastQC}/$id"){`mkdir "$output_dirs{FastQC}/$id"`;}
				#=====================================================================#
				push @{$commands{$id}}, "fastqc -f bam_mapped \"$outfiles{$id}{$sample}{view_bam}\" -o \"$output_dirs{FastQC}/$id\"";
				#=====================================================================#
			}
		}
		
		#If there are multiple bam files per subject/patient (i.e. read from more than 1 sequencer lane or run), merge the bam files into a single file.
		#In order to retain knowledge of read origin, we will need the RG header tags from each source bam file. These will be concatenated into a header file and included in the merged bam.
		#Using samtools merge, we can merge the bams, include the headers and annotate each read in the merged file with unique original ID (therefore reads are traceable to their origin).
		#Caveat: GATK requires each read's RG tag to relate to one of the RG header's ID section (e.g. patient-lane-run), so the ID section must be different for each RG header. The SM section
		#for all RG headers must be the same. In terms of this perl script, ID is equvalent to $sample and SM is equivalent to $id. Definitions are important!
		if ($bam_count > 1){
			my $merge_command = "samtools merge -fr -h \"$outfiles{$id}{merged}{cat_rg_headers}\" \"$outfiles{$id}{merged}{merged_bam}\"";
			my $cat_command = "cat";
			# my $merge_command = "java -jar $options{j} O=\"$outfiles{$id}{merged}{merged_bam}\" USE_THREADING=true SORT_ORDER=queryname VALIDATION_STRINGENCY=LENIENT TMP_DIR=~/tmp";
			
			for my $sample (sort {"\L$a" cmp "\L$b"} keys %{$outfiles{$id}}){
				next if $sample eq "merged";
				
				my $header_file = $outfiles{$id}{$sample}{"view_bam"};
				$header_file =~ s/\.bam/\_RG_header.txt/;
				
				$merge_command .= " \"$outfiles{$id}{$sample}{view_bam}\"";
				$cat_command .= " \"$header_file\"";
				
				#=====================================================================#
				push @{$commands{$id}}, "samtools view -H \"$outfiles{$id}{$sample}{view_bam}\" | awk '/\@RG/' > \"$header_file\"";
				#=====================================================================#
			}
			
			$cat_command .= " > \"$outfiles{$id}{merged}{cat_rg_headers}\"";
			
			#=====================================================================#
			push @{$commands{$id}}, "$cat_command";
			push @{$commands{$id}}, "$merge_command";
			#=====================================================================#
			
		}
		
		else{
			for my $sample (sort {"\L$a" cmp "\L$b"} keys %{$outfiles{$id}}){
				next if $sample eq "merged";
				
				$outfiles{$id}{"merged"}{"merged_bam"} = $outfiles{$id}{$sample}{"view_bam"};
			}
		}
		
		my $sample = "merged";
		
		#=====================================================================#
		push @{$commands{$id}},
			"samtools sort -@ $options{c} -n \"$outfiles{$id}{$sample}{merged_bam}\" \"$outfiles{$id}{$sample}{sorted_bam}\"",
			"samtools fixmate \"$outfiles{$id}{$sample}{sorted_bam}.bam\" \"$outfiles{$id}{$sample}{fixmate_bam}\"",
			##"samtools view -H \"$outfiles{$id}{$sample}{starting_bam}\" | sed 's/16569/16571/' > \"$outfiles{$id}{$sample}{header_sam}\"",
			##"java -jar /Users/marc/Documents/bioinformatics_software/picard-tools-2.0.1/picard.jar ReplaceSamHeader I=\"$outfiles{$id}{$sample}{starting_bam}\" HEADER=\"$outfiles{$id}{$sample}{header_sam}\" OUTPUT=\"$outfiles{$id}{$sample}{reheader_bam}\"",
			##"java -jar /Users/marc/Documents/bioinformatics_software/picard-tools-2.0.1/picard.jar ReorderSam INPUT=\"$outfiles{$id}{$sample}{reheader_bam}\" OUTPUT=\"$outfiles{$id}{$sample}{reorder_bam}\" REFERENCE=\"$options{r}\"",
			"samtools sort -@ $options{c} \"$outfiles{$id}{$sample}{fixmate_bam}\" \"$outfiles{$id}{$sample}{post_fixmate_sorted_bam}\"",
			"samtools rmdup \"$outfiles{$id}{$sample}{post_fixmate_sorted_bam}\.bam\" \"$outfiles{$id}{$sample}{nodups_bam}\"",
			"samtools index \"$outfiles{$id}{$sample}{nodups_bam}\"",
			"rm \"$outfiles{$id}{$sample}{merged_bam}\"",
			"rm \"$outfiles{$id}{$sample}{sorted_bam}.bam\"",
			"rm \"$outfiles{$id}{$sample}{fixmate_bam}\"",
			"rm \"$outfiles{$id}{$sample}{post_fixmate_sorted_bam}.bam\"",
			
			#Human indel realignment and BQSR
			#"java -jar \"$options{g}\" RealignerTargetCreator -I \"$outfiles{$id}{$sample}{post_fixmate_sorted_bam}.bam\" -R \"$options{r}\" -o \"$outfiles{$id}{$sample}{indels_intervals}\" -nt $options{c} -filterNoBases --known \"$options{q}\" --known \"$options{m}\"",
			#"java -jar \"$options{g}\" IndelRealigner -I \"$outfiles{$id}{$sample}{post_fixmate_sorted_bam}.bam\" -R \"$options{r}\" -targetIntervals \"$outfiles{$id}{$sample}{indels_intervals}\" -o \"$outfiles{$id}{$sample}{realigned_bam}\" -filterNoBases -known \"$options{q}\" -known \"$options{m}\"",
			"java -jar \"$options{g}\" BaseRecalibrator -I \"$outfiles{$id}{$sample}{nodups_bam}\" -R \"$options{r}\" -O \"$outfiles{$id}{$sample}{recalibrated_grp}\" --known-sites \"$options{d}\" -RF GoodCigarReadFilter",
			"java -jar \"$options{g}\" ApplyBQSR -I \"$outfiles{$id}{$sample}{nodups_bam}\" -R \"$options{r}\" -O \"$outfiles{$id}{$sample}{post_recalibrated_grp}\" -bqsr \"$outfiles{$id}{$sample}{recalibrated_grp}\" -RF GoodCigarReadFilter --emit-original-quals",
			
			#Mouse indel realignment and BQSR
			# "java -jar \"$options{g}\" RealignerTargetCreator -I \"$outfiles{$id}{$sample}{nodups_bam}\" -R \"$options{r}\" -o \"$outfiles{$id}{$sample}{indels_intervals}\" -nt $options{c} -filterNoBases --known \"$options{l}\"",
			# "java -jar \"$options{g}\" IndelRealigner -I \"$outfiles{$id}{$sample}{nodups_bam}\" -R \"$options{r}\" -targetIntervals \"$outfiles{$id}{$sample}{indels_intervals}\" -o \"$outfiles{$id}{$sample}{realigned_bam}\" -filterNoBases -known \"$options{l}\"",
			# "java -jar \"$options{g}\" BaseRecalibrator -I \"$outfiles{$id}{$sample}{realigned_bam}\" -R \"$options{r}\" -l INFO -o \"$outfiles{$id}{$sample}{recalibrated_grp}\" -knownSites \"$options{k}\" -knownSites \"$options{l}\" -nct $options{c} -filterNoBases -rf BadCigar",
			# "java -jar \"$options{g}\" BaseRecalibrator -I \"$outfiles{$id}{$sample}{realigned_bam}\" -R \"$options{r}\" -l INFO -o \"$outfiles{$id}{$sample}{post_recalibrated_grp}\" -knownSites \"$options{k}\" -knownSites \"$options{l}\" -nct $options{c} -filterNoBases -BQSR \"$outfiles{$id}{$sample}{recalibrated_grp}\" -rf BadCigar",

			"java -jar \"$options{g}\" AnalyzeCovariates -before \"$outfiles{$id}{$sample}{recalibrated_grp}\" -after \"$outfiles{$id}{$sample}{post_recalibrated_grp}\" -plots \"$outfiles{$id}{$sample}{BQSR_plots}\"",
			"java -jar \"$options{g}\" PrintReads -I \"$outfiles{$id}{$sample}{post_recalibrated_bam}\" -R \"$options{r}\" -O \"$outfiles{$id}{$sample}{recalibrated_bam}\" -RF GoodCigarReadFilter";
			"java -jar \"$options{g}\" HaplotypeCaller -R \"$options{r}\" -D \"$options{d}\" -I \"$outfiles{$id}{$sample}{recalibrated_bam}\" -o \"$outfiles{$id}{$sample}{genotype_calls}\" --annotation Coverage --annotation AlleleBalance --annotation ReadPosRankSumTest --annotation MappingQualityRankSumTest --annotation RMSMappingQuality --annotation QualByDepth -ERC GVCF -nct $options{c}",
			"rm \"$outfiles{$id}{$sample}{recalibrated_grp}\"",
			"rm \"$outfiles{$id}{$sample}{post_recalibrated_grp}\"",
			"rm \"$outfiles{$id}{$sample}{nodups_bam}\"",
			"rm \"$outfiles{$id}{$sample}{realigned_bam}\"";
		#=====================================================================#
			
		#If wanted, perform on-target metric (Kevin) or in-depth (Marc) analysis here.
		#Note that other scripts called here are subject to change and there may be other issues.
		#see exomeMetrics_BWA_Illumina.pl for on-target analysis and Marc for in-depth analysis scripts.
		if ($options{"n"} =~ /seq_stat_new.pl/){
			#Call seq stat new with appropriate parameters.
		}
		
		else{
			my $html_out = $outfiles{$id}{$sample}{"recalibrated_bam"};
			$html_out =~ s/recalibrated\.bam/on_target_metrics\.html/;
			
			#=====================================================================#
			#Call exomeMetrics_BWA_Illumina.pl with appropriate parameters
			# push @{$commands{$id}}, "exomeMetrics_BWA_Illumina.pl -p $options{j} -i $id -b \"$outfiles{$id}{$sample}{recalibrated_bam}\" -f \"$options{b}\" -r \"$options{r}\" -o $html_out";
			#=====================================================================#
		};
		
		#=====================================================================#
		#push @{$commands{$id}},	
			#"java -jar \"$options{g}\" HaplotypeCaller -R \"$options{r}\" -D \"$options{d}\" -I \"$outfiles{$id}{$sample}{recalibrated_bam}\" -o \"$outfiles{$id}{$sample}{genotype_calls}\" --annotation Coverage --annotation AlleleBalance --annotation ReadPosRankSumTest --annotation MappingQualityRankSumTest --annotation RMSMappingQuality --annotation QualByDepth -rf BadCigar -ERC GVCF -variant_index_type LINEAR -variant_index_parameter 128000 -nct $options{c}",
			# "java -jar \"$options{g}\" HaplotypeCaller -R \"$options{r}\" -D \"$options{d}\" -I \"$outfiles{$id}{$sample}{recalibrated_bam}\" -o \"$outfiles{$id}{$sample}{variant_calls}\" --annotation Coverage --annotation AlleleBalance --annotation ReadPosRankSumTest --annotation MappingQualityRankSumTest --annotation RMSMappingQuality --annotation QualByDepth -rf BadCigar -stand_emit_conf 10 -stand_call_conf 30",
			#"java -jar \"$options{g}\" GenotypeGVCFs -R \"$options{r}\" -V \"$outfiles{$id}{$sample}{genotype_calls}\" -D \"$options{d}\" -nda -o \"$outfiles{$id}{$sample}{variant_calls}\"",
			#"java -jar \"$options{g}\" ValidateVariants -R \"$options{r}\" -V \"$outfiles{$id}{$sample}{variant_calls}\" -D \"$options{d}\"";
			#"java -jar \"$options{g}\" SelectVariants -R \"$options{r}\" -V \"$outfiles{$id}{$sample}{variant_calls}\" -selectType SNP -o \"$outfiles{$id}{$sample}{SNPs_only}\"",
			#"java -jar \"$options{g}\" SelectVariants -R \"$options{r}\" -V \"$outfiles{$id}{$sample}{variant_calls}\" -selectType INDEL -o \"$outfiles{$id}{$sample}{indels_only}\"",
			#"java -jar \"$options{g}\" VariantFiltration -R \"$options{r}\" -V \"$outfiles{$id}{$sample}{SNPs_only}\" --filterExpression \"QD \< 2.0 || FS \> 60.0 || MQ \< 35.0 || MQRankSum \< -12.5 || ReadPosRankSum \< -8.0\" --filterName \"HardFilter\" -o \"$outfiles{$id}{$sample}{SNPs_recalibrated}\"",
			#"java -jar \"$options{g}\" VariantFiltration -R \"$options{r}\" -V \"$outfiles{$id}{$sample}{indels_only}\" --filterExpression \"QD \< 2.0 || FS \> 200.0 || ReadPosRankSum \< -20.0\" --filterName \"HardFilter\" -o \"$outfiles{$id}{$sample}{indels_recalibrated}\"",
			#"java -jar \"$options{g}\" CombineVariants -R \"$options{r}\" -V:SNP \"$outfiles{$id}{$sample}{SNPs_recalibrated}\" -V:Indel \"$outfiles{$id}{$sample}{indels_recalibrated}\" -genotypeMergeOptions UNSORTED -o \"$outfiles{$id}{$sample}{variants_recalibrated}\"";
			
			# "rm \"$outfiles{$id}{$sample}{variant_calls}\"",
			# "rm \"$outfiles{$id}{$sample}{SNPs_only}\"",
			# "rm \"$outfiles{$id}{$sample}{indels_only}\"",
			# "rm \"$outfiles{$id}{$sample}{SNPs_recalibrated}\"",
			# "rm \"$outfiles{$id}{$sample}{indels_recalibrated}\"",
			
# # 			"java -jar \"$options{e}\" eff -c \"$options{e2}\" -i vcf -o gatk GRCh37.75 \"$outfiles{$id}{$sample}{variants_recalibrated}\" > \"$outfiles{$id}{$sample}{snpeff_annotation}\"",
			# "java -jar \"$options{g}\" VariantAnnotator -R \"$options{r}\" -A SnpEff -V \"$outfiles{$id}{$sample}{variants_recalibrated}\" -snpEffFile \"$outfiles{$id}{$sample}{snpeff_annotation}\" -o \"$outfiles{$id}{$sample}{variants_annotated}\"";
			
			#For multi-allelic variants, snpEff will annotate each variant allele but GATK will only include the most significant effect in the final VCF.
			#If you therefore want to check annotation for all variant alleles, you'll have to look in the VCF file output directly from snpEff.
			#"rm \"$outfiles{$id}{$sample}{snpeff_annotation}\"";
			
 			# "run_ANNOVAR.pl -i \"$outfiles{$id}{$sample}{variants_recalibrated}\" -a knownGene -d /data/db/annovar/humandb -b hg19 -g 1000g2014sep_all -s snp138 -S ljb26_sift -c cosmic70 -o \"$outfiles{$id}{$sample}{variants_annotated}\"",
			# "annotate_vcf_from_annovar_output.pl -i \"$outfiles{$id}{$sample}{variants_recalibrated}\" -d \"$output_dirs{variant_annotation}\" -h /data/db/human/HGNC/HGNC.txt > \"$outfiles{$id}{$sample}{variants_annotated}\"",
			
# #  			#Convert vcf into property file:
			# "multiple_sample_vcf_to_properties.pl -i \"$outfiles{$id}{$sample}{variants_annotated}\" -s gatk_haplotype_caller > \"$outfiles{$id}{$sample}{final_properties}\"",

# #  			#Convert properties to spreadsheet (tab delimited text file).
			# #Note that header.list defines which properties to include in table and in which order.
			# "trio_properties_to_txt.pl -i \"$outfiles{$id}{$sample}{final_properties}\" -h \"$outfiles{$id}{$sample}{headers_file}\" > \"$outfiles{$id}{$sample}{final_txt}\"",

# #   			#Create excel spreadsheet
			# "textToExcel.pl -i \"$outfiles{$id}{$sample}{final_txt}\" -o \"$outfiles{$id}{$sample}{final_xlsx}\"";
			
			
			
			#=====================================================================#
			#Work in progress#
			#Add annotation using Ensembl's Variant Effect Predictor (VEP).
			#Tested on a local machine (Shark) using a downloaded version (release 91) and a cache of GRCh37. Some versions (i.e. ensembl 92) don't work well with GRCh37.
			#	There are several options during install, make sure the ensembl VEP (e.g. 91) and genome versions (e.g. 91_GRCh37) are the same.
			#	It may also ask if you want to convert (format of) the genome/cache - this decreases run time so saying yes is recommended.
			#	The installation process may take several hours (large downloads/conversions, etc).
			#It is possible to connect to Ensembl via an API during run time to add various annotations, but this increases run time (run time of a simple offline command: ~2 minutes).
			#Install commands (copied directly from website, https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html):
			#	git clone https://github.com/Ensembl/ensembl-vep.git
			#	cd ensembl-vep
			#	git pull
			#	git checkout release/91		#Set version to 91
			#	perl INSTALL.pl -a acf -y GRCh37 --CACHE_VERSION 91 --NO_HTSLIB -s homo_sapiens -v 91	#Options added to hopefully install without prompts (untested).
			#=====================================================================#
			#/Users/marc/Documents/bioinformatics_software/ensembl-vep/vep --cache -i /Users/marc/Documents/work/Projects/E052/somatic_variants/E052-H-001_merged_somatic_variants.vcf -o /Users/marc/Documents/work/Projects/E052/variant_annotation/E052-H-001_merged_somatic_variants_annotated_offline-test.vcf --offline --force_overwrite
			#=====================================================================#
			
		#=====================================================================#
			
		if (exists $tissues{$id}){
			my $germline_pair = $tissues{$id};
			
			my $germline_bam = $outfiles{$germline_pair}{$sample}{"recalibrated_bam"};
			
			#=====================================================================#
			#These MuTetct2 commands are based on documentation for GATK 4.0.1.2, and may not work with earlier/later versions.
			push @{$commands{$id}},
				#for somatic samples (paired with germline samples).
				"java -jar \"$options{g}\" Mutect2 -R \"$options{r}\" -I \"$outfiles{$id}{$sample}{post_recalibrated_bam}\" -tumor $id -I \"$outfiles{$germline_pair}{$sample}{post_recalibrated_bam}\" -normal $germline_pair --germline-resource $options{i} -O \"$outfiles{$id}{$sample}{somatic_variant_calls}\"";
				
				#for germline samples (stand alone).
				"java -jar \"$options{g}\" Mutect2 -R \"$options{r}\" -I \"$outfiles{$germline_pair}{$sample}{post_recalibrated_bam}\" -normal $germline_pair --germline-resource $options{i} -O \"$outfiles{$id}{$sample}{somatic_variant_calls}\"";
				
				#snpEff is no longer supported by GATK tools
				#"java -jar \"$options{e}\" eff -c \"$options{e2}\" -i vcf -o gatk GRCh37.75 \"$outfiles{$id}{$sample}{somatic_variant_calls}\" > \"$outfiles{$id}{$sample}{somatic_snpeff_annotation}\"",
				#"java -jar \"$options{g}\" VariantAnnotator -R \"$options{r}\" -A SnpEff -V \"$outfiles{$id}{$sample}{somatic_variant_calls}\" -snpEffFile \"$outfiles{$id}{$sample}{somatic_snpeff_annotation}\" -o \"$outfiles{$id}{$sample}{somatic_variants_annotated}\"",
				
				#default to Kevin's annovar scripts - however I need to find them and finish these commands (2018-05-11).
				#"/wgp1/wgp/apps/nestor-bin/run_ANNOVAR.pl -i \"$outfiles{$id}{$sample}{somatic_variant_calls}\" -a knownGene -d /scratch/share_ngs/hg19/annovar/humandb -b hg19 -g 1000g2014sep_all -s snp138 -S ljb26_sift -c cosmic70 -o $sample\"",
				#"/wgp1/wgp/apps/nestor-bin/annotate_vcf_from_annovar_output.pl -i \"$outfiles{$id}{$sample}{somatic_variant_calls}\" -d . -h /wgp1/wgp/resources/human/HGNC/HGNC.txt
			
			#=====================================================================#
		}
			
	}
	
	return \%commands;
}
	
########################################################################
#Subroutine delorean
#Gets the local time and converts it into a date (format: yyyy-mm-dd) and time (hh:mm:ss).
#These are then returned,
########################################################################
sub delorean{
	my @flux_capacitor = localtime(time);
	
	#Set the variables for the time and date. Note: year is a count of years after 1900 (e.g. 2015 = 115), and month is 0-based integer (0 = Jan, 11 = Dec).
	my ($sec, $min, $hour, $date, $month, $year) = @flux_capacitor;
	
	$month += 1;
	$year += 1900;
	
	$month = sprintf ("%02d", $month);
	$date = sprintf ("%02d", $date);
	$day = sprintf ("%02d", $day);
	$sec = sprintf ("%02d", $sec);
	$min = sprintf ("%02d", $min);
	$hour = sprintf ("%02d", $hour);
	
	my $return_date = "$year-$month-$date";
	
	my $return_time = "$hour:$min:$sec";
	
	return ($return_date, $return_time);
}
########################################################################
#End of subroutine delorean.
########################################################################

########################################################################
#Subroutine captains_log
#Create a log file to record commands performed by the current execution of the script.
########################################################################
sub captains_log{
	my $options = shift;
	
	my %options = %{$options};
	
	my $out_dir = $options{"O"}."/Log/";
	my $job = $options{"J"};
	my $name = `whoami`;
	chomp $name;
	$name =~ s/\?//;
	
	mkdir "$out_dir" unless (-d "$out_dir");
	
	#Get the current date/time.
	my ($date, $time) = &delorean;
	
	my $log_file = $out_dir."VarScape-log_(".$date."_".$time."_".$job."_".$name.").txt";
	
	open LOG, ">$log_file" or warn "WARNING! Log file ($log_file) could not be created!\n";
	print LOG "exome_pipeline_raven.pl log.\nProject: $job\tUser: $name\tDate:\t$date\tTime: $time\n\n";
	print LOG "The following commands were executed at the stated times.\n\n";
	close LOG or warn "WARNING! Log file ($log_file) could not be closed!\n";
	
	return $log_file;
}
########################################################################
#End of subroutine captains_log.
########################################################################

########################################################################
#Subroutine fastqc
#Call fastqc and performed QC on mapped BAM files (prior to merging).
########################################################################
sub fastqc{
	my ($commands, $outfiles, $out_dirs, $id, $sample) = @_;
	
	my %commands = %{$commands};
	my %outfiles = %{$outfiles};
	my %out_dirs = %{$out_dirs};
	print "pushing fq command\n";
	unless (-d "$output_dirs{FastQC}/$id"){`mkdir "$output_dirs{FastQC}/$id"`;}
	push @{$commands{$id}}, "fastqc -f bam_mapped \"$outfiles{$id}{$sample}{view_bam}\" -o \"$output_dirs{FastQC}/$id\"";
	
	return \%commands;
}
########################################################################
#End of subroutine fastqc
########################################################################

########################################################################
#Subroutine execution
#Execute commands stored in a hash of arrays.
########################################################################
sub execution{
	my ($commands, $outfiles, $log_file) = @_;
	
	my %commands = %{$commands};
	my %outfiles = %{$outfiles};
	
	for my $id (sort {"\L$a" cmp "\L$b"} keys %commands){
		my $command_no = 0;
		
		my $no_commands = scalar @{$commands{$id}};
		
		print "cd $options{O}\n";
		`cd $options{"O"}`;
		
		#create new PBS and set nodes and reporting options.
		`newpbs -name $id\_variant_calling_pipeline -nodes 1 -cpus 12 -walltime 36 -debug`;
		`cat $id\_variant_calling_pipeline.pbs | sed 's/\#PBS -l nodes.*/\#PBS -l select=1:ncpus=12/' >> $id\_variant_calling_pipeline_new.pbs`;
		`rm $id\_variant_calling_pipeline.pbs`;
		`mv $id\_variant_calling_pipeline_new.pbs $id\_variant_calling_pipeline.pbs`;
		
		`echo module load java >> $id\_variant_calling_pipeline.pbs`;
		`echo module load perl >> $id\_variant_calling_pipeline.pbs`;
		`echo module load BWA/0.7.15 >> $id\_variant_calling_pipeline.pbs`;
		`echo module load samtools/0.1.19 >> $id\_variant_calling_pipeline.pbs`;
		`echo module load fastqc/0.11.2 >> $id\_variant_calling_pipeline.pbs`;
	
		
		open HEADERS, ">$outfiles{$id}{merged}{headers_file}" or die "headers file ($outfiles{$id}{$sample}{headers_file}) could not be opened.\n";
		print HEADERS "##HEADERS\n#Comment lines are ignored.\n##List headers they will appear in final table.\n#Comment out columns not required.\n###########################\n";
		print HEADERS "chromosome\nposition\nref.allele\nalt.allele\nfilter\nvariant.quality\nnumber.with.snp\n\ngene.id\ngene.description\ngenomic.function\ngenomic.annotation\ndbSNP.id\ndbSNP.build\n1000Genomes.frequency\nsift.outcome\n\n$id\_$sample.alt.allele.depth\n$id\_$sample.alt.allele.frequency\n$id\_$sample.ref.allele.depth\n$id\_$sample.ref.allele.frequency\n$id\_$sample.total.depth\n";
		close HEADERS;
		
		my ($date, $time) = &delorean;
		
		print color ("red"), "---------------------------------------------------------------\nStarting analysis of subject $id.\nDate: $date\tTime: $time\tID: $id\tNo. commands: $no_commands.\n---------------------------------------------------------------\n";
		print color 'reset';
		
		open LOG, ">>$log_file" or warn "Log file could not be opened during execution!\n";
		print LOG "---------------------------------------------------------------\nStarting analysis of subject $id.\nDate: $date\tTime: $time\tID: $id.\n---------------------------------------------------------------\n";
		close LOG;
			
		foreach my $command (@{$commands{$id}}){
			$command_no++;
			
			my ($date, $time) = &delorean;
			
			print color ("cyan"), "Date: $date\tTime: $time\tID: $id.\tCommand $command_no/$no_commands:\n";
			print color ("green"), "$command\n\n";
			print color 'reset';
			
			open LOG, ">>$log_file" or warn "Log file could not be opened during execution!\n";
			print LOG "Date: $date\tTime: $time\tID: $id.\tCommand:\n$command\n\n";
			close LOG;
			
			#print commands to PBS script.
			`echo "$command" >> $id\_variant_calling_pipeline.pbs`;
		}
		
		open LOG, ">>$log_file" or warn "Log file could not be opened during execution!\n";
		print LOG "---------------------------------------------------------------\nAnalysis of subject $id is finished.\n---------------------------------------------------------------\n\n";
		close LOG;
		
		print color ("red"), "---------------------------------------------------------------\nAnalysis of subject $id is finished.\n---------------------------------------------------------------\n\n";
		print color 'reset';
	}
}
########################################################################
#End of subroutine execution.
########################################################################
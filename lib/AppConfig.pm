#!/usr/bin/perl
###############################################################################
#
#    AppConfig.pl
#    
#    Makes more useful names for the fields in the config file
#    The app_* scripts should include this file first
#
#    Copyright (C) 2011 Michael Imelfort
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################
package AppConfig;
require Exporter;
use File::Basename;
use File::Temp qw(tempfile);
use Carp;

our @ISA = qw(Exporter);
our @EXPORT=qw(
    $VERSION
    IGNORE_FAILURE
    WARN_ON_FAILURE
    DIE_ON_FAILURE
    $global_analysis_config_filename
    $global_qiime_mapping_filename
    %FNB 
    %FNA 
    $FNB_HEADER 
    $FNA_HEADER 
    $FNA_LINE_FINISHER 
    $FNA_FOOTER 
    $FNA_SILVA_FOOTER
    $APP_ROOT 
    $APP_RAW 
    $APP_BYJOB 
    $APP_BYRUN 
    %global_samp_ID_list 
    %global_raw_counts 
    %global_chimer_counts 
    %global_acacia_counts 
    $global_acacia_config
    $default_trim_length
    $global_barcode_length 
    $QA_dir
    $proc_dir 
    $res_dir
    $image_dir
    $global_acacia_output_dir 
    $global_working_dir
    $global_QA_dir
    $global_processing_dir
    $global_results_dir
    $global_image_dir
    $global_TB_processing_dir
    $global_SB_processing_dir
    $global_TB_results_dir
    $global_SB_results_dir
    $ACACIA_out_file
    $tn_log_file
    $tn_R_log_file
    $sn_log_file
    $tn_dist_file
    $sn_dist_file
    $nn_prefix    
    $nn_fasta_file
    $nn_otus_file
    $nn_otu_table_file
    $nn_expanded_otu_table_file
    $nn_tree_file
    $tn_prefix    
    $tn_otu_table_file
    $tn_expanded_otu_table_file 
    $tn_tree_file
    $sn_prefix    
    $sn_fasta_file
    $sn_otus_file
    $sn_otu_table_file
    $sn_expanded_otu_table_file 
    $sn_tree_file
    $global_mapping_file
    %acacia_config_hash
    $QIIME_split_out
    $QIIME_split_out_qual
    $QIIME_TAX_tax_file
    $QIIME_TAX_blast_file
    $QIIME_TAX_aligned_blast_file
    $QIIME_map_file
    $QIIME_imputed_file 
    $SILVA_TAX_tax_file
    $SILVA_TAX_blast_file
    $SILVA_TAX_aligned_blast_file
    $SILVA_imputed_file
    $MERGED_TAX_tax_file
    $MERGED_TAX_blast_file
    $MERGED_TAX_aligned_blast_file
    $MERGED_imputed_file 
    $global_R_log_file
    checkFileExists
    checkAndRunCommand
    runExternalCommand
    logExternalCommand
    getWorkingDirs 
    makeOutputDirs 
    makeResultsDirs
    makeImageDirs 
    denoise 
    getReadCounts 
    parseConfigQA 
    updateConfigQA
    );

# Version of this APP (update for each release candidate/tag)
our $VERSION = '3.0.3';

# Failure modes when executing a command
use constant {
    IGNORE_FAILURE => 0,
    WARN_ON_FAILURE => 1,
    DIE_ON_FAILURE => 2
};

our $global_analysis_config_filename = "config.txt";
our $global_qiime_mapping_filename = "qiime_mapping.txt";


#
# A file is created in PyroDB which can be used to split the sff file and 
# make all the relavant job dirs etc... XXX.pdbm
# This file is basically a qiime mapping file, BUT it has the same name as
# the sff. (or fasta and qual). The format is given below:
#
# SampleID	BarcodeSequence	LinkerPrimerSequence	Description
# <JID.SID> MID             acgggcggtgtgtRc         <PDB sample name>
# ...
#
our $FNB_HEADER = "#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription";
our %FNB = ();
$FNB{'SampleID'} = 0;
$FNB{'BarcodeSequence'} = 1;
$FNB{'LinkerPrimerSequence'} = 2;
$FNB{'Description'} = 3;

#
# Once the sff has been munged, each job will be placed into a folder in the by_jobID dir
# and given an app config file. The format is given below:
#
# #SampleID	BarcodeSequence	LinkerPrimerSequence	Description	        RAW	CHIME	ACC	USE
# <SID>     MID             acgggcggtgtgtRc         <PDB sample name>   XX  XX      XX  XX
# ...
#
our $FNA_HEADER = "#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription\tRAW\tCHIME\tACC\tUSE";
our $FNA_LINE_FINISHER = "\tXX\tXX\tXX\t1\n";
our $FNA_FOOTER = "@@\
NORMALISE=\
DB=0\
MUL_RARE_M=\
MUL_RARE_X=\
MUL_RARE_S=\
MUL_RARE_N=";

our $FNA_SILVA_FOOTER = "@@\
NORMALISE=\
DB=SILVA\
MUL_RARE_M=\
MUL_RARE_X=\
MUL_RARE_S=\
MUL_RARE_N=";

our %FNA = ();
$FNA{'SampleID'} = 0;
$FNA{'BarcodeSequence'} = 1;
$FNA{'LinkerPrimerSequence'} = 2;
$FNA{'Description'} = 3;
$FNA{'RAW'} = 4;
$FNA{'CHIME'} = 5;
$FNA{'ACC'} = 6;
$FNA{'USE'} = 7;

#
# APP_ROOT should be set by the module system. there are a number of dirs we need to get from there
#
our $APP_ROOT = `echo \$APP_ROOT`;
chomp $APP_ROOT;
our $APP_RAW = $APP_ROOT."/raw";
our $APP_BYJOB = $APP_ROOT."/by_jobid";
our $APP_BYRUN = $APP_ROOT."/by_run";

#
# We make a number of directories during the process. Store their names here
#
our $QA_dir = "QA_CD_HIT_OTU";
our $proc_dir = "processing";
our $res_dir = "results";
our $image_dir = "images";

our $global_acacia_output_dir = "UNSET";
our $global_working_dir = "UNSET";
our $global_mapping_file = "UNSET";



#
# The is the default trim length during QA.
#
our $default_trim_length = 250;

our $global_barcode_length = "variable_length";

#
# Some programs make special output files. Store these filenames here
#
our $QIIME_GG_TAX_ROOT = "/srv/whitlam/bio/db/gg/from_www.secongenome.com/current";

our $QIIME_map_file = "qiime_mapping.txt";
our $QIIME_split_out = "seqs.fna";
our $QIIME_split_out_qual = "seqs_filtered.qual";
our $QIIME_TAX_tax_file = "$QIIME_GG_TAX_ROOT/gg_12_10_taxonomy.txt";
our $QIIME_TAX_blast_file = "$QIIME_GG_TAX_ROOT/gg_12_10_otus/rep_set/99_otus.fasta";
our $QIIME_TAX_aligned_blast_file = "$QIIME_GG_TAX_ROOT/gg_12_10_otus/rep_set/99_otus.fasta";
our $QIIME_imputed_file = "$QIIME_GG_TAX_ROOT/gg_12_10_aligned.fasta";

#SILVA DBs
our $SILVA_TAX_tax_file = "/srv/whitlam/bio/db/Silva/QIIME_files/taxonomy_mapping/Silva_taxa_mapping_104set_97_otus.txt";
our $SILVA_TAX_blast_file = "/srv/whitlam/bio/db/Silva/QIIME_files/rep_set/silva_104_rep_set.fasta";
our $SILVA_TAX_aligned_blast_file = "/srv/whitlam/bio/db/Silva/QIIME_files/rep_set/silva_104_rep_set_aligned.fasta";
our $SILVA_imputed_file = "/srv/whitlam/bio/db/Silva/QIIME_files/core_aligned_set/core_Silva_aligned.fasta";

# Merged DBs (GG Bacteria + Archaea with Silva Eukaryotes)
our $MERGED_TAX_tax_file = "/srv/whitlam/bio/db/merged_gg_silva/merged_gg_silva_taxo.txt";
our $MERGED_TAX_blast_file = "/srv/whitlam/bio/db/merged_gg_silva/merged_gg_silva.fna";



our $CHIME_good_file = "good.fasta";
our $CHIME_bad_file = "ch.fasta";
our $ACACIA_out_file = "acacia_out_all_tags.seqOut";

#
# Some global variables we use to store read counts and sample IDs
# 
our %global_samp_ID_list = ();          # list of sample IDs
our %global_raw_counts = ();            # raw reads in $QA_dir/$QIIME_split_out
our %global_chimer_counts = ();
our %global_acacia_counts = ();

#
# File and directory names made by app_make_results (or QIIME)
#
# values are set in getWorkingDirs below
#
#### DIRS made by app_do_QA.pl
our $global_QA_dir = "UNSET";
our $global_processing_dir = "UNSET";
our $global_results_dir = "UNSET";
our $global_image_dir = "UNSET";

#### PROCESSING DIRS
our $global_TB_processing_dir = "UNSET";
our $global_SB_processing_dir = "UNSET";

#### DIRS made by app_make_results.pl

#### RESULTS
our $global_TB_results_dir = "UNSET";
our $global_SB_results_dir = "UNSET";

our $tn_log_file = "UNSET";
our $tn_R_log_file = "UNSET";
our $sn_log_file = "UNSET";
our $tn_dist_file = "UNSET";
our $sn_dist_file = "UNSET";

our $nn_prefix = "UNSET";
our $nn_fasta_file = $nn_prefix."UNSET";
our $nn_otus_file = $nn_prefix."UNSET";
our $nn_otu_table_file = "UNSET";
our $nn_expanded_otu_table_file = "UNSET";
our $nn_tree_file = "UNSET";

our $tn_prefix = "UNSET";
our $tn_otu_table_file = "UNSET";
our $tn_expanded_otu_table_file = "UNSET";
our $tn_tree_file = "UNSET";

our $sn_prefix = "UNSET";
our $sn_fasta_file = $sn_prefix."UNSET";
our $sn_otu_table_file = "UNSET";
our $sn_expanded_otu_table_file = "UNSET";
our $sn_otus_file = $sn_prefix."UNSET";
our $sn_tree_file = "UNSET";

#### made by app_make_images.pl

our $global_R_log_file = "UNSET";

######################################################################
# SHARED SUBS
######################################################################
sub checkFileExists {
    my $file = shift;
    unless(-e $file) {
        die "ERROR!\n\nCannot find:\n$file\n";
    }
}
sub logExternalCommand {
    print "\t", shift, "\n\n";
}

sub runExternalCommand {
    my $cmd = shift;
    logExternalCommand($cmd);
    system($cmd);
}

sub isCommandInPath
{
    #-----
    # Is this command in the path?
    #
    my ($cmd, $failure_type) = @_;
    if (system("which $cmd |> /dev/null")) {
        handleCommandFailure($cmd, $failure_type);
    }
}

sub checkAndRunCommand
{
    #-----
    # Run external commands more sanelier
    #
    my ($cmd, $params, $failure_type) = @_;
    
    isCommandInPath($cmd, $failure_type);
    
    # join the parameters to the command
    my $param_str = join " ", map {formatParams($_)} @{$params};
    
    my $cmd_str = $cmd . " " . $param_str;
    
    logExternalCommand($cmd_str);

    if (system($cmd_str)) {
        handleCommandFailure($cmd_str, $failure_type)
    }
}

sub formatParams {

    #---------
    # Handles and formats the different ways of passing parameters to 
    # checkAndRunCommand
    #
    my $ref = shift;
    
    if (ref($ref) eq "ARRAY") {
        return join(" ", @{$ref});
    } elsif (ref($ref) eq "HASH") {
        return join(" ", map { $_ . " " . $ref->{$_}} keys %{$ref});
    }
    croak 'The elements of the $params argument in checkAndRunCommand can ' .
        'only contain references to arrays or hashes\n';
}


sub handleCommandFailure {
    #-----
    # What to do when all goes bad!
    #
    my ($cmd, $failure_type) = @_;
    if (defined($failure_type)) {
        if ($failure_type == DIE_ON_FAILURE) {
            croak "**ERROR: $0 : " . $! . "\n";
        } elsif ($failure_type == WARN_ON_FAILURE) {
            carp "**WARNING: $0 : " . $! . "\n";
        }
    }
}

sub getWorkingDirs
{
    #-----
    # Set a number of output directories
    #
    my ($base_directory, $results_output_dir, $analysis_type) = @_;
    
    # get the acacia denoised directory
    $global_acacia_output_dir = $acacia_config_hash{'OUTPUT_DIR'};
    
    if ($analysis_type eq "QIIME") {
        $QA_dir = "QA_QIIME"
    }

    # Get the working dir
    # working dir is the dir of the config file
    # get the present dir
    my $pwd = `pwd`;
    chomp $pwd;
    $global_working_dir = dirname("$pwd/$base_directory");
    
    # set the mapping file
    $global_mapping_file = "$global_working_dir/$QA_dir/$QIIME_map_file";
    # now we set these guys
    $global_QA_dir = "$global_working_dir/$QA_dir";
    $global_processing_dir = "$results_output_dir/$proc_dir";
    $global_results_dir = "$results_output_dir/$res_dir";
    $global_image_dir = "$results_output_dir/$image_dir";
    $global_TB_processing_dir = "$global_processing_dir/table_based";
    $global_SB_processing_dir = "$global_processing_dir/sequence_based";
    $global_TB_results_dir = "$global_results_dir/table_based";
    $global_SB_results_dir = "$global_results_dir/sequence_based";
    $tn_log_file = "$global_TB_results_dir/otu_table_normalisation.stats";
    $tn_R_log_file = "$global_TB_processing_dir/otu_table_normalisation_R_commands.log";
    $sn_log_file = "$global_SB_results_dir/sequence_normalisation.log";
    $tn_dist_file = "$global_TB_processing_dir/otu_table_normalisation_dist.txt";
    $sn_dist_file = "$global_SB_processing_dir/sequence_normalisation_dist.txt";
    $nn_prefix = "non_normalised";
    $nn_fasta_file = $nn_prefix.".fa";
    $nn_otus_file = $nn_prefix."_otus.txt";
    $nn_otu_table_file = "$global_TB_results_dir/$nn_prefix"."_otu_table.txt";
    $nn_expanded_otu_table_file = "$global_TB_results_dir/$nn_prefix"."_otu_table_expanded.tsv";
    $nn_tree_file = $nn_prefix."_tree.tre";
    $tn_prefix = "table_normalised";
    $tn_otu_table_file = "$global_TB_results_dir/$tn_prefix"."_otu_table.txt";
    $tn_expanded_otu_table_file = "$global_TB_results_dir/$tn_prefix"."_otu_table_expanded.tsv";
    $tn_tree_file = "$tn_prefix"."_tree.tre";
    $sn_prefix = "sequence_normalised";
    $sn_fasta_file = $sn_prefix.".fa";
    $sn_otu_table_file = "$global_SB_results_dir/$sn_prefix"."_otu_table.txt";
    $sn_expanded_otu_table_file = "$global_SB_results_dir/$sn_prefix"."_otu_table_expanded.tsv";
    $sn_otus_file = $sn_prefix."_otus.txt";
    $sn_tree_file = "$sn_prefix"."_tree.tre";    
}

sub makeOutputDirs
{
    #-----
    # Directories must be made before we can put files there
    #
    my ($job_dir) = @_;
    `mkdir -p $global_working_dir$job_dir/$QA_dir`;
    `mkdir -p $global_working_dir$job_dir/$QA_dir/$global_acacia_output_dir`;
}

sub makeResultsDirs
{
    #-----
    # Make directories needed suring app_make_results.pl
    #
    my ($do_sb) = @_;
    print "$global_TB_results_dir\n";
    `mkdir -p $global_TB_results_dir`;
    `mkdir -p $global_TB_processing_dir`;
    # we only need to do this as an extra
    if(1 == $do_sb)
    {
        `mkdir -p $global_SB_processing_dir`;
        `mkdir -p $global_SB_results_dir`;
    }
}

sub makeImageDirs
{
    #-----
    # Make directories needed by app_make_images.pl
    #
    `mkdir -p $global_image_dir`;
    $global_R_log_file = $global_image_dir."/R.log"; 
}

sub getReadCounts
{
    #-----
    # get the read counts for raw, chimera removed and acacia filtered sequences
    #
    # this guy is called from within the QA dir so the files are local
    #
    # the three files to parse are:
    # $QIIME_split_out
    # $CHIME_good_file
    # $global_acacia_output_dir/$ACACIA_out_file
    #
    my $raw_only = shift;
    open my $tmp_fh, "<", $QIIME_split_out or die $!;
    while(<$tmp_fh>)
    {
        next if(!($_ =~ /^>/));
        my @f1 = split /_/, $_;
        my @f2 = split />/, $f1[0];
        my $fl = $f2[1];
        foreach my $uid (keys %global_samp_ID_list)
        {
            # QIIME replaces _ with . in Sample IDs, need to make the correction
            # here or keys in the hash won't be found in the output sequence
            # file.
            my $qiimeified_uid = qiimeify_uid($uid);
            if($fl eq $qiimeified_uid)
            {
                # this guy begins with the exact MID
                $global_raw_counts{$uid}++;
                if ($raw_only) {
                    $global_chimer_counts{$uid}++;
                    $global_acacia_counts{$uid}++;
                }
                last;
            }
        }
    }
    if ($raw_only) {
        return;
    }
    open $tmp_fh, "<", $CHIME_good_file or die $!;
    while(<$tmp_fh>)
    {
        next if(!($_ =~ /^>/));
        my @f1 = split /_/, $_;
        my @f2 = split />/, $f1[0];
        my $fl =$f2[1];  
        foreach my $uid (keys %global_samp_ID_list)
        {
            my $qiimeified_uid = qiimeify_uid($uid);
            if($fl eq $qiimeified_uid)
            {
                # this guy begins with the exact MID
                $global_chimer_counts{$uid}++;
                last;
            }
        }
    }    

    open $tmp_fh, "<", "$global_acacia_output_dir/$ACACIA_out_file" or die "Could not open file: $global_acacia_output_dir/$ACACIA_out_file $!";
    while(<$tmp_fh>)
    {
        next if(!($_ =~ /^>/));
        my @f1 = split /_/, $_;
        my @f2 = split />/, $f1[0];
        my $fl =$f2[1];  
        foreach my $uid (keys %global_samp_ID_list)
        {
            my $qiimeified_uid = qiimeify_uid($uid);
            if($fl eq $qiimeified_uid)
            {
                # this guy begins with the exact MID
                $global_acacia_counts{$uid}++;
                last;
            }
        }
    }
} 

sub parseConfigQA
{
    #-----
    # parse the app config file and produce a qiime mappings file
    #
    my ($config_prefix, $return_params_hash) = @_;
    open my $conf_fh, "<", $config_prefix or die $!;
    open my $mapping, ">", "qiime_mapping.txt" or die $!;
    print $mapping "$FNB_HEADER\n";
    my $used_sample_count = 0;
    # Parse the sample information.
    while(<$conf_fh>)
    {
        next if($_ =~ /^#/);
        last if($_ =~ /^@/);
        chomp $_;
        my @fields = split /\t/, $_;
        if ($fields[0] =~ /[_\-]/) {
            die "ERROR: Sample ID ". $fields[0] . " contains an underscore or hyphen ".
            "(_) and this can cause APP to fail (as QIIME replaces the offending character).\n" .
            "Replace the underscore/hyphen (with a .) in the config file and rerun
            APP_do_QA.pl.\n\n";
        }
        # save the MID for later
        $global_samp_ID_list{$fields[$FNA{'SampleID'}]} = $fields[$FNA{'USE'}];
        $global_raw_counts{$fields[$FNA{'SampleID'}]} = 0;
        $global_chimer_counts{$fields[$FNA{'SampleID'}]} = 0;
        $global_acacia_counts{$fields[$FNA{'SampleID'}]} = 0;

        if(! ("0" eq $fields[$FNA{'USE'}]))
        {
            $used_sample_count++;
            print $mapping "$fields[$FNA{'SampleID'}]\t$fields[$FNA{'BarcodeSequence'}]\t$fields[$FNA{'LinkerPrimerSequence'}]\t$fields[$FNA{'Description'}]\n";
        }
    }
    if ($return_params_hash) {
        while(<$conf_fh>) {
            next if($_ =~ /^#/);
            chomp $_;
            if ($_ =~ /^(.*?)=(.*)$/) {
                $return_params_hash->{$1} = $2;
            }
        }    
    }
    close $conf_fh;
    close $mapping;
    return $used_sample_count;
}

sub updateConfigQA
{
    #-----
    # parse the app config  and update to include read counts
    # this guy is called from within the QA dir so we need to do a ../ on the file names
    #
    my ($config_prefix) = @_;
    open my $conf_fh, "<", "../$config_prefix" or die $!;
    open my $conf_fh_tmp, ">", "../$config_prefix.tmp" or die $!;
    use Data::Dumper;
    while(<$conf_fh>)
    {
        if($_ =~ /^#/) { print $conf_fh_tmp $_; next; }
        if($_ =~ /^@/) { print $conf_fh_tmp $_; last; }
        chomp $_;
        my @fields = split /\t/, $_;
        #print Dumper(@fields);
        #print Dumper(\%global_raw_counts);
        print $conf_fh_tmp join("\t", ($fields[$FNA{'SampleID'}],
                                       $fields[$FNA{'BarcodeSequence'}],
                                       $fields[$FNA{'LinkerPrimerSequence'}],
                                       $fields[$FNA{'Description'}],
                                       $global_raw_counts{$fields[$FNA{'SampleID'}]},
                                       $global_chimer_counts{$fields[$FNA{'SampleID'}]},
                                       $global_acacia_counts{$fields[$FNA{'SampleID'}]},
                                       $fields[$FNA{'USE'}])), "\n";
    }
    
    # just print out the rest of the file
    while(<$conf_fh>)
    {
        print $conf_fh_tmp $_;
    }

    close $conf_fh;
    close $conf_fh_tmp;
    
    my $mv_string  = "mv ../$config_prefix.tmp ../$config_prefix";
    `$mv_string`;
}

sub qiimeify_uid {
    my $qiimeified_uid = shift;
    $qiimeified_uid =~ s/_/./g;
    return $qiimeified_uid;
}

1;

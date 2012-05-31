#!/usr/bin/perl
###############################################################################
#
#    app_combine.pl
#    
#    Combine data from multiple jobs into one folder for processing
#
#    Copyright (C) 2011 Michael Imelfort and Paul Dennis
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

#pragmas
use strict;
use warnings;

#core Perl modules
use Getopt::Long;

#CPAN modules
use File::Basename;

#locally-written modules
#load the pretty names for the fields
use AppConfig;

BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# get input params and print copyright
printAtStart();

my $options = checkParams();

######################################################################
# CODE HERE
######################################################################
#### GLOBALS
print "Checking if all the config checks out...\t\t";
my @global_job_sample_list = split /,/, $options->{'jobs'};
# for storing jobs separately
my %global_job_list = ();

# check to see if we havre at least two samples to combine
if($#global_job_sample_list < 1)
{
    print "**ERROR: you MUST supply at least two job folders\n";
    die;
}

# base output folder
my $global_combined_base = "ac";
if(exists $options->{'output'})
{
    $global_combined_base = $options->{'output'};
}
$global_combined_base =~ s/\/$//;

# acacia config file
if (defined($options->{'acacia_conf'})) {
    updateAcaciaConfigHash($options->{'acacia_conf'})
}

if (defined($options->{'length'})) {
    $acacia_config_hash{TRIM_TO_LENGTH} = $options->{'length'};
}

print "All good!\n";

# get job and sample names
foreach my $js (@global_job_sample_list)
{
    my @split_name = split /\./, $js;
    if(!exists $global_job_list{$split_name[0]})
    {
        $global_job_list{$split_name[0]} = 1;
    }
    if(0 == $#split_name)
    {
        # whole job
        open my $conf_fh, "<", "$split_name[0]/app_$split_name[0].config" or die $!;
        while(<$conf_fh>)
        {
            next if($_ =~ /^#/);
            last if($_ =~ /^@/);
            my @fields = split /\t/, $_;
            $global_samp_ID_list{$fields[0]} = 1;
        }
        close $conf_fh;
    }
    else
    {
        # single sample
        $global_samp_ID_list{$split_name[1]} = 1;
    }
}

# let the user know we're on the case
print "Combining reads from samples:\n";
my $mod_counter = 0;
my $col_width = 8;
my $sep = "\n ";
foreach my $samp (keys %global_samp_ID_list)
{
    print $sep.$samp;
    $mod_counter++;
    if($mod_counter == 8)
    {
        $mod_counter = 0;
        $sep = "\n ";
    }
    else
    {
        $sep = "\t";
    }
}
print "\n\nIn new combined job folder: $global_combined_base\n----------------------------------------------------------------\n";

#### Extract the reads from each sample

# make the output folders
`mkdir -p $global_combined_base`;
getWorkingDirs($global_combined_base);
makeOutputDirs("/$global_combined_base");

# make some files to write to
open my $global_fna_out_fh, ">", "$global_combined_base/$QA_dir/$QIIME_split_out" or die $!;
open my $global_qual_out_fh, ">", "$global_combined_base/$QA_dir/$QIIME_split_out_qual" or die $!;
open my $global_qimme_mapping_file_fh, ">", "$global_combined_base/$QA_dir/$QIIME_map_file" or die $!;
print $global_qimme_mapping_file_fh "$FNB_HEADER\n";

# go through all the individual fna files
# this is roughly equivalent to running split libraries
foreach my $job (keys %global_job_list)
{
    recombine_fnas($job);
}

close $global_fna_out_fh;
close $global_qual_out_fh;
close $global_qimme_mapping_file_fh;

print "----------------------------------------------------------------\n";

# make the config environment
print "Making the config file\n";
open my $new_conf_fh, ">", "$global_combined_base/app_".$global_combined_base.".config" or die "$!: $global_combined_base\n" ;
print $new_conf_fh "$FNA_HEADER\n";

foreach my $job (keys %global_job_list)
{
    open my $conf_fh, "<", "$job/app_$job.config" or die $!;
    while(<$conf_fh>)
    {
        next if($_ =~ /^#/);
        last if($_ =~ /^@/);
        my @fields = split /\t/, $_;
        if(exists $global_samp_ID_list{$fields[0]})
        {
            print $new_conf_fh $_;
        }
    }
    close $conf_fh;
}

if(exists $options->{'SILVA'})
{
    print $new_conf_fh "$FNA_SILVA_FOOTER\n";
}
else
{
    print $new_conf_fh "$FNA_FOOTER\n";
}
close $new_conf_fh;

#### start the $QA_dir pipeline!
my $c_dir = `pwd`;
chomp $c_dir;
$global_working_dir = $c_dir."/$global_combined_base";
chdir "$global_working_dir/$QA_dir";
removeChimeras();
denoise();

#### Fix the config file
print "Fixing read counts...\n";
getReadCounts();
updateConfigQA("app_$global_combined_base.config");

print "QA complete!\n";

######################################################################
# CUSTOM SUBS
######################################################################
sub recombine_fnas
{
    #-----
    # given a list of jobs and samples, recombine to form one complete unified job dir
    # first we make a lookup of readIDs to job IDs
    
    my ($job_ID) = @_;
    my $split_fasta = "$job_ID/$QA_dir/$QIIME_split_out";
    my $split_qual = "$job_ID/$QA_dir/$QIIME_split_out_qual";
    my $mapping_file = "$job_ID/$QA_dir/$QIIME_map_file";
    print "Processing: $split_fasta...\t\t";

    my %seq_id_hash = ();
    
    # make a mapping file
    my $conf_file = "$job_ID/app_".$job_ID.".config";
    open my $c_fh, "<", $conf_file or die $!;
    # make this dir if it doesn't exist
    `mkdir -p $job_ID/$QA_dir`;
    open my $map_fh, ">", $mapping_file  or die $!;
    print $map_fh $FNB_HEADER."\n";
    while(<$c_fh>)
    {
        next if($_ =~ /^#/);
        last if($_ =~ /^@/);
        chomp $_;
        my @fields = split /\t/, $_;

        if(exists($global_samp_ID_list{$fields[$FNA{'SampleID'}]}))
        {
            print $map_fh "$fields[$FNA{'SampleID'}]\t$fields[$FNA{'BarcodeSequence'}]\t$fields[$FNA{'LinkerPrimerSequence'}]\t$fields[$FNA{'Description'}]\n";
            print $global_qimme_mapping_file_fh "$fields[$FNA{'SampleID'}]\t$fields[$FNA{'BarcodeSequence'}]\t$fields[$FNA{'LinkerPrimerSequence'}]\t$fields[$FNA{'Description'}]\n";
        }
    }
    close $c_fh;
    close $map_fh;
    
    # first we see if a raw file exists
    # make one if it doesn't
    #if(! -e $split_fasta)
    #{
        my $c_dir = `pwd`;
        chomp $c_dir;
        chdir("$job_ID/$QA_dir");
        splitLibraries($job_ID);
        chdir($c_dir);
    #}

    # open the raw file
    open my $fna_fh, "<", $split_fasta  or die $!;
    open my $qual_fh, "<", $split_qual  or die $!;
    while(<$fna_fh>)
    {
        # for the sequence
        my $header = $_;
        my $seq = <$fna_fh>;
        
        # for the qual
        my $qual_header = <$qual_fh>;
        my $qual_val = <$qual_fh>;
        
        # split on ">"
        my @tmp_1 = split />/, $_;
        # then split on spaces
        my @tmp_2 = split / /, $tmp_1[1];
        # then split on underscore
        my @tmp_3 = split /_/, $tmp_2[0];
        
        # check to see if this is one of the sequences we wanted
        my $samp_ID = $tmp_3[0];
        if(exists $global_samp_ID_list{$samp_ID})
        {
            print $global_fna_out_fh $header;
            print $global_fna_out_fh $seq;
            print $global_qual_out_fh $qual_header;
            print $global_qual_out_fh $qual_val;            
        }
    }
    close $fna_fh;
    close $qual_fh;
    
    print "Done\n";
}



######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    my @standard_options = ( "help|h+", "jobs|j:s", "output|o:s",  "length|l:i", "acacia_conf:s", "SILVA+");
    my %options;

    # Add any other command line options, and the code to handle them
    # 
    GetOptions( \%options, @standard_options );

    # if no arguments supplied print the usage and exit
    #
    exec("pod2usage $0") if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    # Compulsosy items
    if(!exists $options{'jobs'} ) { print "**ERROR: you MUST supply at least two job folders\n"; exec("pod2usage $0"); }
    if(!exists $options{'length'} ) { print "**ERROR: you MUST supply a length to trim reads to\n"; exec("pod2usage $0"); }

    return \%options;
}

sub printAtStart {
print<<"EOF";
---------------------------------------------------------------- 
 $0
 Copyright (C) 2011 Michael Imelfort and Paul Dennis
    
 This program comes with ABSOLUTELY NO WARRANTY;
 This is free software, and you are welcome to redistribute it
 under certain conditions: See the source for more details.
---------------------------------------------------------------- 
EOF
}

__DATA__

=head1 NAME

    app_combine.pl

=head1 COPYRIGHT

   copyright (C) 2011 Michael Imelfort and Paul Dennis

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 DESCRIPTION

   Does filtering, denoising and de-replication of 454 pyrotag datasets

=head1 SYNOPSIS

    app_combine.pl -jobs|j JOB_ID[.SAMPLE_ID],JOB_ID[.SAMPLE_ID][,JOB_ID[.SAMPLE_ID][,...]] [-acacia_conf CONFIG_FILE] -length|l TRIM_LENGTH [-output|o FOLDER] [-SILVA] [-help|h]

      -jobs -j                      APP job folders and possibly sample numbers (must be at aleast two)
      -length -l                    Trim all reads to this length (bp).
      [-acacia_conf CONFIG_FILE]    Alternate acacia config file (Full path!)
      [-output|o FOLDER]            Folder to write the results to [default: ac]
      [-SILVA]                      Use the SILVA database
      [-help -h]                    Displays basic usage information
      
      NOTES:
      
      This script ignores the "USE" flags in the app_XX.config files
      
      Each job and each sample are given a unique ID by PyroDB. Use this script if you 
      want to combine data from a number of samples from in number of different jobs. 
      
      If you specify a different acacia config file, then you must use
      the following values, or this script will break!
          
      FASTA_LOCATION=good.fasta
      OUTPUT_DIR=denoised_acacia
      OUTPUT_PREFIX=acacia_out_
      SPLIT_ON_MID=FALSE

      Examples:
      
      app_combine.pl -j 56.3,56.8,57.2 -l 230   # combines samples 3 and 8 from Job 56 and sample 2 from job 57 and trims them to 230 bp
      app_combine.pl -j 56.3,56.8,57,59.108 -l 350   # combines samples 3 and 8 from Job 56, all samples from job 57 and sample 108 from job 59 and trims them to 350 bp.
      
=cut


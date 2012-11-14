#!/usr/bin/perl
###############################################################################
#
#    app_create_analysis.pl
#    
#    Creates the analysis directory from the fna, qual and config files.
#
#    Copyright (C) 2011 Michael Imelfort and Paul Dennis,
#                  2012 Adam Skarshewski
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
use Carp;

#CPAN modules
use File::Basename;

#locally-written modules
#load the pretty names for the fields
use AppConfig;
use AppCommon;

BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

######################################################################
# CODE HERE
######################################################################

# get input params and print copyright
my $options = checkParams();

printAtStart();

# Setup the analysis folder
my $params_hash = {};

if (! -e $options->{'c'}) {
    die "Unable to find config file: " . $options->{'c'} . "\n";
}

my $job_name = basename($options->{'c'});

if ($job_name =~ /app_(.*).config$/) {
    $job_name = $1;
} else {
    croak "The app config file needs to be of the form app_<prefix>.config, ".
        "where <prefix> can be chosen by the user.\n" ;
}
    
my @argv_to_app_params = (['c', 'original_config_file'],
                          ['l', 'trim_length'],
                          ['o', 'output_folder'],
                          ['acacia_conf', 'acacia_conf_file'],
                          ['i', 'identity'],
                          ['e', 'e_value'],
                          ['b', 'custom_blast_database'],
                          ['t', 'custom_taxonomy_file'],
                          ['I', 'custom_imputed_file'],
                          ['a', 'trim_length'],
                          ['n', 'number_of_threads']);

foreach my $arg_pair (@argv_to_app_params) {
    my $arg = $arg_pair->[0];
    my $app_param = $arg_pair->[1];
    if (defined($options->{$arg})) {
        $params_hash->{$app_param} = $options->{$arg};
    }
}

use Data::Dumper;

my $output_dir = setup_analysis_folder($job_name, $params_hash);
print "\nAnalysis folder created. Use the following command to start the analysis:\n";
print "\tapp_run_analysis.pl -d $output_dir\n\n";


########################
# SUBS
########################

################################################################################
# Subroutine: setup_analysis_folder($job_name, $params_hash)
# 
# 
#
################################################################################

sub setup_analysis_folder {
    my ($job_name, $params_hash) = @_;
    
    my $output_dir =
        create_unique_analysis_dir($params_hash->{output_folder});
    
    create_data_file_links($output_dir, $job_name);
    
    my ($sample_array, $config_array) =
        parse_app_config_file($params_hash->{original_config_file});
    
    my @samples_to_use;
    open(my $fh, ">$output_dir/sample_exclusion.txt");
    print {$fh} "#NAME\tREAD_COUNTS\tUSE\n";
    foreach my $sample_count (@{$sample_array}) {
        if (exists($sample_count->[4]) && ($sample_count->[4] == 0)) {
            print {$fh} $sample_count->[0] . "\t??\t0\n";
        } else {
            print {$fh} $sample_count->[0] . "\t??\t1\n";
            push @samples_to_use, $sample_count->[0];
        }

    }
    close($fh);

    foreach my $param (keys %{$params_hash}) {
        if ($param eq 'original_config_file') {
            next;
        }
        push @{$config_array}, [uc($param), uc($params_hash->{$param})]
    }
    
    create_analysis_config_file("$output_dir/$global_analysis_config_filename",
                                $config_array, \@samples_to_use);
    
    `cp app_$job_name.config $output_dir/app.config`;
    
    return $output_dir;
    
}


######################################################################
# TEMPLATE SUBS
######################################################################

sub checkParams {
    my @standard_options = ( "help|h+", "c:s",  "l:i", "acacia_conf:s",
                            "o:s", "i:i", "e:f", "b:s", "t:s", "I:s",
                            "a:s", "n:s");
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
    if(!exists $options{'c'} ) { print "**ERROR: you MUST give a config file\n"; exec("pod2usage $0"); }
    #if(!exists $options{''} ) { print "**ERROR: \n"; exec("pod2usage $0"); }

    return \%options;
}

sub printAtStart {
print<<"EOF";
---------------------------------------------------------------- 
 $0
 Version $VERSION
 Copyright (C) 2011 Michael Imelfort and Paul Dennis
               2012 Adam Skarshewski
    
 This program comes with ABSOLUTELY NO WARRANTY;
 This is free software, and you are welcome to redistribute it
 under certain conditions: See the source for more details.
---------------------------------------------------------------- 
EOF
}

__DATA__

=head1 NAME

    app_create_analysis.pl

=head1 COPYRIGHT

   copyright (C) 2011 Michael Imelfort and Paul Dennis,
                 2012 Adam Skarshewski

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

    app_create_analysis.pl -c|config CONFIG_FILE -l|length TRIM_LENGTH [-help|h]

      -c CONFIG_FILE               app config file to be processed
      [-l TRIM_LENGTH]             Trim all reads to this length (default: 250bp)
      [-o OUTPUT_FOLDER]           Output folder name (default: app_analysis_<date>)
      [-acacia_conf CONFIG_FILE]   alternate acacia config file (Full path!)
      [-i identity VALUE]          Set blast identity for OTU clustering (pick_otus.py) [default: 97%]
      [-e EVALUE]                  Set e-value for blast (assign_taxonomy.py) [default 0.001]      
      [-b FILE]                    Path to a custom blast database / bwa database
      [-t FILE]                    Path to a custom taxonomy for otus
      [-I FILE]                    Path to a custom imputed file
      [-a FILE]                    assign_taxonomy method [default blast, alternative bwasw (for BWA)]
      [-n NUM_THREADS]             Use this many threads where possible [default 5]
      
      [-help -h]                   Displays basic usage information
         
    NOTE:
      
    If you specify a different acacia config file, then you must use
    the following values, or this script will break!
      
    FASTA_LOCATION=good.fasta
    OUTPUT_DIR=denoised_acacia
    OUTPUT_PREFIX=acacia_out_
    SPLIT_ON_MID=FALSE
=cut
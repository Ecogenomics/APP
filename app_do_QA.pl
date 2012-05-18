#!/usr/bin/perl
###############################################################################
#
#    app_do_QA.pl
#    
#    Uses qiime scripts + acacia to do mid splitting and denoising
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
use Carp;

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

if (defined($options->{'acacia_conf'})) {
    updateAcaciaConfigHash($options->{'acacia_conf'})
}

# get the Job_ID we're working on
my $job_ID = basename($options->{'config'});
if ($job_ID =~ /app_(.*).config$/) {
    $job_ID = $1;
} else {
    croak "The app config file needs to be of the form app_<prefix>.config, ".
        "where <prefix> can be chosen by the user.\n" ;
}

# get the working directories
getWorkingDirs($options->{'config'});

# make the output directories
makeOutputDirs("");

my $QA_params = {};
# parse the config file
if (! parseConfigQA($options->{'config'}, $QA_params)) {
    croak "No samples were selected for analysis, correct the USE column in " .
        "config file\n";    
};


my $split_library_params;
if ($QA_params->{ION_TORRENT}) {    
    $split_library_params = {-l => 150,
                             -s => 15};
    $acacia_config_hash{TRIM_TO_LENGTH} = 150;
}

print "All good!\n";

#### start the $QA_dir pipeline!
chdir "$global_working_dir/$QA_dir";
splitLibraries($job_ID, $split_library_params);
removeChimeras();
denoise();

#### Fix the config file
print "Fixing read counts...\n";
getReadCounts();
updateConfigQA($options->{'config'});

######################################################################
# CUSTOM SUBS
######################################################################

## SEE ./lib/AppConfig.pm

######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    my @standard_options = ( "help|h+", "config|c:s", "acacia_conf:s" );
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
    if(!exists $options{'config'} ) { print "**ERROR: you MUST give a config file\n"; exec("pod2usage $0"); }
    #if(!exists $options{''} ) { print "**ERROR: \n"; exec("pod2usage $0"); }

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

    app_do_QA.pl

=head1 COPYRIGHT

   copyright (C) 2011 Michael Imelfort and Paul Dennis, 2012 Adam Skarshewski

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

    app_do_QA.pl -c|config CONFIG_FILE [-help|h]

      -c CONFIG_FILE               app config file to be processed
      [-acacia_conf CONFIG_FILE]   alternate acacia config file (Full path!)
      [-help -h]                   Displays basic usage information
         
         
    NOTE:
      
    If you specify a different acacia config file, then you must use
    the following values, or this script will break!
      
    FASTA_LOCATION=good.fasta
    OUTPUT_DIR=denoised_acacia
    OUTPUT_PREFIX=acacia_out_
    SPLIT_ON_MID=FALSE
=cut


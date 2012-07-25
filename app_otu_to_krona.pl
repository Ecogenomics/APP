#!/usr/bin/env perl
###############################################################################
#
#    app_otu_to_krona.pl
#    
#    Make a pretty Krona graph
#
#    Copyright (C) 2012 Connor Skennerton
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
use Pod::Usage;
use IO::File;
#CPAN modules

#locally-written modules

BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# get input params and print copyright
printAtStart();
my $global_options = checkParams();
my $in;
my $output_file_name = &overrideDefault("text.krona.html", 'output');
if (defined $global_options->{'input'}) {
    $in = &openRead($global_options->{'input'});
} else {
    $in = \*STDIN;
}

my $number_of_samples = 0;
# the number of columns that is defined by pyrotagger
my $base_columns = 14;

# array to sample names
my @sample_names;
# an array of corresponding tmp files
my @tmp_files;
# with the qiime output this line is waste
# pyrotagger uses it as the header
my $first_line = <$in>;
if($global_options->{"pyrotagger"}) {
    $first_line =~ s/\s+$//;
    
    # figure out how many samples there are
    my @n = split(/\t/, $first_line);
    my $number_of_columns = scalar @n;
    $number_of_samples = $number_of_columns - $base_columns;
    @sample_names = splice(@n, 1, $number_of_samples);
} else {
    # the second line in qiime is the header
    my $second_line = <$in>;
    chomp $second_line;
    $second_line =~ s/\s+$//;
    @sample_names = split(/\t/, $second_line);
    # remove first and last columns
    shift @sample_names; pop @sample_names;
    # @names = map {$_.".krona.tmp.otus"} @z;
}


my @file_names = map{$_.".krona.tmp.otus"} @sample_names;
@tmp_files = map {IO::File->new($_, 'w')} @file_names;
while(<$in>) {
    chomp;
    my @otus;
    my @c = split(/\t/, $_);
    if($global_options->{"pyrotagger"}) {
        
        # if pyrotagger doesn't find a hit it outputs nothing
        # therefore the elements of the array 
        next if($#c <( $number_of_samples + $base_columns));

        @otus = @c[1 .. $number_of_samples];
        # pyrotagger outputs a blank rather than 0, need to fix that
        $_ ||= '0' for @otus;
        my $tax_string = join("\t", splice(@c, $number_of_samples + $base_columns));
        my @full_name = splice(@c,-3,1);
        $tax_string .= "\t".$full_name[0];

        foreach my $i (0 .. $#otus) {
            $tmp_files[$i]->print($otus[$i],"\t", $tax_string, "\n");
        }
    } else {
        # deal with the greengenes or silva formatted tax string at the end of a line
        my @l = split(/;?\w_{1,2}|;/,$c[-1]);
        
        # Add the OTU number on the end
        push @l, "otu".$c[0];
        
        foreach my $i (1 .. $#c - 1) {
            $tmp_files[$i-1]->print(
                                    $c[$i],"\t",join("\t", @l), "\n");
        }
    }
}
foreach (@tmp_files) {
    close;
}
close $in;

if (! defined $global_options->{'temp'}) {
    my $files = join (" ", @file_names);
    `ktImportText -o $output_file_name $files`
}

if (! defined $global_options->{'keep'}) {
    unlink @file_names;
}

exit;
######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help|h+", "keep|k+", "temp|t+", "input|i:s", "output|o:s", "pyrotagger|p+" );
    my %options;

    # Add any other command line options, and the code to handle them
    # 
    GetOptions( \%options, @standard_options );

    # if no arguments supplied print the usage and exit
    #
    #pod2usage() if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    pod2usage() if $options{'help'};

    # Compulsosy items
    #if(!exists $options{''} ) { print "**ERROR: $0 : \n"; exec("pod2usage $0"); }

    return \%options;
}

sub printAtStart {
print<<"EOF";
---------------------------------------------------------------- 
 $0
 Copyright (C) 2012 Connor Skennerton
    
 This program comes with ABSOLUTELY NO WARRANTY;
 This is free software, and you are welcome to redistribute it
 under certain conditions: See the source for more details.
---------------------------------------------------------------- 
EOF
}

sub openWrite
{
    #-----
    # Open a file for writing
    #
    my ($fn) = @_;
    open my $fh, ">", $fn or die "**ERROR: could not open file: $fn for writing $!\n";
    return $fh;
}

sub openRead
{   
    #-----
    # Open a file for reading
    #
    my ($fn) = @_;
    open my $fh, "<", $fn or die "**ERROR: could not open file: $fn for reading $!\n";
    return $fh;
}

sub overrideDefault
{
    #-----
    # Set and override default values for parameters
    #
    my ($default_value, $option_name) = @_;
    if(exists $global_options->{$option_name}) 
    {
        return $global_options->{$option_name};
    }
    return $default_value;
}

__DATA__

=head1 NAME

    app_otu_to_krona.pl

=head1 COPYRIGHT

   copyright (C) 2012 Connor Skennerton

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

   Take an OTU table outputted by APP and use krona to create an interactive
   pie chart layout thingy...

=head1 SYNOPSIS

    app_otu_to_krona.pl [-help|h] [-input|i] [-keep|k] [-temp|t] [-output|o]

      [-help -h]                   Displays basic usage information
      [-input -i]                  Input file name [default: STDIN]
      [-keep -k]                   Keep all temorary files created
      [-temp -t]                   Only produce the temorary files in Krona format without
                                   running krona
      [-output -o]                 Name of the output html file produced by krona [default: text.krona.html]
      [-pyrotagger -p]             Input OTU table is in pyrotagger format
         
=cut

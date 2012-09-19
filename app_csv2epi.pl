#!/usr/bin/perl
###############################################################################
#
#    app_csv2epi.pl
#
#    Convert a csv file to another csv file which is formatted for use in the robot
#    Input file should look like this!
#    A1,9.4,pyroL803Fmix
#    A2,8.1,pyroL803Fmix
#    Bleg...
#    
#    NOTE:
#    All units are in ng and uL
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

#locally-written modules
use AppPrimers;
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
#### Globals / Defaults
# set the upper limit for the pooled tube
my $global_pool_tube_capacity = 2000; # uL
my $global_dilutant_tube_capacity = 2000; # uL

# set the amount of DNA we need for each sample
my $global_sample_DNA_required_total = 4; #ng

# we need to monitor the maximum amount of fluid in each original tube
my $global_sample_tube_volume = 13; # uL
my $global_sample_tube_capacity = 300; # uL

# where is all this stuff going?
my $global_pool_tube_position = "A1";
my $global_total_pool_volume = 0;
# set the minium volume we should take from any one well
my $global_pool_minimum_volume = 1; # uL

# where is the dilutant stored
my $global_dilutant_tube_position = "A2";
my $global_total_dilutant_volume = 0; # uL
# we won't bother diluting if 
my $global_dilutant_minimum_volume = 2; # uL

# prefix for output files
my $global_out_prefix = "";

#### OVERRIDE DEFAULTS
if (exists $options->{'sample_DNA_required_total'}) { $global_sample_DNA_required_total = $options->{'sample_DNA_required_total'}; }
if (exists $options->{'sample_tube_volume'}) { $global_sample_tube_volume = $options->{'sample_tube_volume'}; }
if (exists $options->{'sample_tube_capacity'}) { $global_sample_tube_capacity = $options->{'sample_tube_capacity'}; }

if (exists $options->{'pool_tube_position'}) { $global_pool_tube_position = $options->{'pool_tube_position'}; }
if (exists $options->{'pool_minimum_volume'}) { $global_pool_minimum_volume = $options->{'pool_minimum_volume'}; }
if (exists $options->{'pool_tube_capacity'}) { $global_pool_tube_capacity = $options->{'pool_tube_capacity'}; }

if (exists $options->{'dilutant_tube_position'}) { $global_dilutant_tube_position = $options->{'dilutant_tube_position'}; }
if (exists $options->{'dilutant_minimum_volume'}) { $global_dilutant_minimum_volume = $options->{'dilutant_minimum_volume'}; }
if (exists $options->{'dilutant_tube_capacity'}) { $global_dilutant_tube_capacity = $options->{'dilutant_tube_capacity'}; }
if (exists $options->{'prefix'}) { $global_out_prefix= $options->{'prefix'}; }

#### PRINT INFORMATION FOR THE USER
print "sample_DNA_required_total: $global_sample_DNA_required_total ng\n";
print "sample_tube_volume: $global_sample_tube_volume uL\n";
print "sample_tube_capacity: $global_sample_tube_capacity uL\n";
print "pool_tube_position: $global_pool_tube_position\n";
print "pool_minimum_volume: $global_pool_minimum_volume uL\n";
print "pool_tube_capacity: $global_pool_tube_capacity uL\n";
print "dilutant_tube_position: $global_dilutant_tube_position\n";
print "dilutant_minimum_volume: $global_dilutant_minimum_volume uL\n";
print "dilutant_tube_capacity: $global_dilutant_tube_capacity uL\n";

# but which primers have we seen?
my %global_seen_primers_hash = ();

# concentration per well location
my %global_well_conc_hash = ();

# primer per well location
my %global_well_primer_hash = ();

# sometimes we need to dilute the stronger samples
my %well_dilute_hash = ();

# need a list of valid wells
my %global_well2id = ();
my %global_id2well = ();
populateVars();

my $global_csv_header = "Rack,Source,Rack,Destination,Volume,Tool\r\n";

# open files
open my $global_input_fh, "<", $options->{'in'} or die $!;
open my $global_pool_fh, ">", $global_out_prefix."_pool.csv" or die $!;
open my $global_dilute_fh, ">", $global_out_prefix."_dilute.csv" or die $!;

#### PARSE THE INPUT FILE
while(<$global_input_fh>)
{
    chomp $_;

    next if($_ eq "");
    
    # remove whitespace
    $_ =~ s/ //g;
    my @line_fields = split /,/, $_;
    
    # check we know this well
    if(!exists $global_well2id{$line_fields[0]})
    {
        die "**ERROR: Unkown well: \"$line_fields[0]\"\n";
    }
    
    # check we know this primer:
    if(!exists $APP_prim_len_hash{$line_fields[2]})
    {
        die "**ERROR: Unkown primer: \"$line_fields[2]\"\n";
    }
    
    # for now, just take the specified concentrations and store the primer
    $global_well_conc_hash{$line_fields[0]} = $line_fields[1];
    $global_well_primer_hash{$line_fields[0]} = $APP_prim_len_hash{$line_fields[2]};
    $global_seen_primers_hash{$line_fields[2]} = $APP_prim_len_hash{$line_fields[2]};
}
# close the input file
close $global_input_fh;

#### MUNGE STUFF UP
# multiple length primers? Find the shortest one.
my $min_primer_len = 1000000000;
my @primers = values %global_seen_primers_hash;
if(0 == $#primers)
{
    # only one primer!
    $min_primer_len = $primers[0];
}
else
{
    foreach my $primer (@primers)
    {
        if($primer < $min_primer_len) { $min_primer_len = $primer; }
    }
}
print "------------------------------\n";
print "File contains: " . ($#primers + 1) . " different primers with a min length of: $min_primer_len\n";
print "------------------------------\n";
print "SANITY CHECK\n";

# now normalise for the primer lengths
foreach my $key (keys %global_well_primer_hash)
{
    $global_well_primer_hash{$key} = $global_well_primer_hash{$key} / $min_primer_len;
}

#### SANITY CHECK
# Run through once and check to see if the parameters chosen match up nicely
my $over_warnings = "";
my $under_warnings = "";

foreach my $well (keys %global_well_conc_hash)
{
    # calculate how much volume to take from the well
    my $volume = ($global_sample_DNA_required_total / $global_well_conc_hash{$well}) * $global_well_primer_hash{$well};
    
    # is this less than the robot likes to take at a minumum?
    if($volume < $global_pool_minimum_volume)
    {
        $under_warnings .= "Warning: Sample in well $well with concentration $global_well_conc_hash{$well} has volume ".roundVolume($volume)." \n";
        
        # calculate the amount of silution needed
        my $dil_amnt = (($global_well_conc_hash{$well} * $global_sample_tube_volume) / $global_sample_DNA_required_total) - $global_sample_tube_volume;
        
        # no point in diluting petty amounts
        if($dil_amnt >= $global_dilutant_minimum_volume)
        {
            $under_warnings .= "\tWill dilute with ".roundVolume($dil_amnt)."\n";
            
            # we are starting with finite sized tubes. So we need to make sure we don't overflow.
            if(($dil_amnt + $global_sample_tube_volume) > $global_sample_tube_capacity)
            {
                $under_warnings .= "\tOVERFLOW!!!\n";
                die $under_warnings;
            }
            # store what we calculated
            $well_dilute_hash{$well} = $dil_amnt;
            $global_total_dilutant_volume += $dil_amnt;
        }
        else
        {
            $under_warnings .= "Dilution amount: ".roundVolume($dil_amnt)." is less than $global_dilutant_minimum_volume so won't dilute\n";
        }
        
        # set the concentration to match the amount needed...
        $global_well_conc_hash{$well} = $global_sample_DNA_required_total;
    }
    # do we require more than is available in the tube?
    # (we will fix this later)
    elsif($volume > $global_sample_tube_volume)
    {
        $over_warnings .= "Warning: Sample in well $well with concentration $global_well_conc_hash{$well} requires ".roundVolume($volume)." \n";
    }
    
    $global_total_pool_volume += $volume;
}

if($global_total_pool_volume > $global_pool_tube_capacity)
{
    die "**ERROR: Total volume (".roundVolume($global_total_pool_volume).") exceeds global limit of $global_pool_tube_capacity\n";
}

if($global_total_dilutant_volume > $global_dilutant_tube_capacity)
{
    die "**ERROR: Total dilutant volume ($global_total_dilutant_volume) exceeds global limit of $global_dilutant_tube_capacity\nTry setting a higher value for sample_DNA_required_total";
}

if("" ne $over_warnings)
{
    print "\n\tWARNING!\n\n\tTaking amounts which are more than contained in each tube: $global_pool_minimum_volume\n\n";
    print $over_warnings;
    print "\n";
}
if("" ne $under_warnings)
{
    print "\n\tWARNING!\n\n\tTaking amounts which are lower than min specified: $global_pool_minimum_volume\n\n";
    print $under_warnings;
    print "\n";
}

#### PRINT THE OUTPUT FILE
print $global_pool_fh $global_csv_header;
print $global_dilute_fh $global_csv_header;

$global_total_pool_volume = 0;
$global_total_dilutant_volume = 0;
foreach my $well_num (sort {$a <=> $b} keys %global_id2well)
{
    my $well = $global_id2well{$well_num};
    if(exists $global_well_conc_hash{$well})
    {
        my $volume = ($global_sample_DNA_required_total / $global_well_conc_hash{$well}) * $global_well_primer_hash{$well};
        if ($volume >  $global_sample_tube_volume) { $volume =  $global_sample_tube_volume; }
        if(exists $well_dilute_hash{$well})
        {
            # do a dilution first (if needed)
            print $global_dilute_fh printLine($global_dilutant_tube_position,$well,$well_dilute_hash{$well}); 
            $global_total_dilutant_volume += $well_dilute_hash{$well};
        }
        $global_total_pool_volume += $volume;        
        print $global_pool_fh printLine($well, $global_pool_tube_position,$volume);
    }
}

# close files
close $global_pool_fh;
close $global_dilute_fh;

#### TELL THE USER WHAT HAPPENED
print "------------------------------\n";
print "Pooling will produce: ".roundVolume($global_total_pool_volume)." uL\n";
print "Diluting will use: ".roundVolume($global_total_dilutant_volume)." uL\n";
print "------------------------------\n";


######################################################################
# CUSTOM SUBS
######################################################################
sub printLine
{
    #-----
    # print one line to output
    #
    my ($position_from, $position_to, $volume) = @_;
    return "1,$position_from,1,$position_to,".roundVolume($volume).",".chooseTool($volume)."\r\n";
}

sub roundVolume
{
    #-----
    # Round a volume
    # 
    my ($volume) = @_;
    return sprintf("%.2f", $volume);
}

sub chooseTool
{
    #-----
    # Choose the appropriate tool based on volume
    #
    my ($volume) = @_;
    if($volume <= 50)
    {
        return 1;
    }
    elsif($volume <= 300)
    {
        return 2;
    }
    else
    {
        die "OMG!!!! Wrong pipette for volume $volume\n";
    }
}


sub populateVars
{
    $global_id2well{1} = "A1";
    $global_id2well{2} = "B1";
    $global_id2well{3} = "C1";
    $global_id2well{4} = "D1";
    $global_id2well{5} = "E1";
    $global_id2well{6} = "F1";
    $global_id2well{7} = "G1";
    $global_id2well{8} = "H1";
    $global_id2well{9} = "A2";
    $global_id2well{10} = "B2";
    $global_id2well{11} = "C2";
    $global_id2well{12} = "D2";
    $global_id2well{13} = "E2";
    $global_id2well{14} = "F2";
    $global_id2well{15} = "G2";
    $global_id2well{16} = "H2";
    $global_id2well{17} = "A3";
    $global_id2well{18} = "B3";
    $global_id2well{19} = "C3";
    $global_id2well{20} = "D3";
    $global_id2well{21} = "E3";
    $global_id2well{22} = "F3";
    $global_id2well{23} = "G3";
    $global_id2well{24} = "H3";
    $global_id2well{25} = "A4";
    $global_id2well{26} = "B4";
    $global_id2well{27} = "C4";
    $global_id2well{28} = "D4";
    $global_id2well{29} = "E4";
    $global_id2well{30} = "F4";
    $global_id2well{31} = "G4";
    $global_id2well{32} = "H4";
    $global_id2well{33} = "A5";
    $global_id2well{34} = "B5";
    $global_id2well{35} = "C5";
    $global_id2well{36} = "D5";
    $global_id2well{37} = "E5";
    $global_id2well{38} = "F5";
    $global_id2well{39} = "G5";
    $global_id2well{40} = "H5";
    $global_id2well{41} = "A6";
    $global_id2well{42} = "B6";
    $global_id2well{43} = "C6";
    $global_id2well{44} = "D6";
    $global_id2well{45} = "E6";
    $global_id2well{46} = "F6";
    $global_id2well{47} = "G6";
    $global_id2well{48} = "H6";
    $global_id2well{49} = "A7";
    $global_id2well{50} = "B7";
    $global_id2well{51} = "C7";
    $global_id2well{52} = "D7";
    $global_id2well{53} = "E7";
    $global_id2well{54} = "F7";
    $global_id2well{55} = "G7";
    $global_id2well{56} = "H7";
    $global_id2well{57} = "A8";
    $global_id2well{58} = "B8";
    $global_id2well{59} = "C8";
    $global_id2well{60} = "D8";
    $global_id2well{61} = "E8";
    $global_id2well{62} = "F8";
    $global_id2well{63} = "G8";
    $global_id2well{64} = "H8";
    $global_id2well{65} = "A9";
    $global_id2well{66} = "B9";
    $global_id2well{67} = "C9";
    $global_id2well{68} = "D9";
    $global_id2well{69} = "E9";
    $global_id2well{70} = "F9";
    $global_id2well{71} = "G9";
    $global_id2well{72} = "H9";
    $global_id2well{73} = "A10";
    $global_id2well{74} = "B10";
    $global_id2well{75} = "C10";
    $global_id2well{76} = "D10";
    $global_id2well{77} = "E10";
    $global_id2well{78} = "F10";
    $global_id2well{79} = "G10";
    $global_id2well{80} = "H10";
    $global_id2well{81} = "A11";
    $global_id2well{82} = "B11";
    $global_id2well{83} = "C11";
    $global_id2well{84} = "D11";
    $global_id2well{85} = "E11";
    $global_id2well{86} = "F11";
    $global_id2well{87} = "G11";
    $global_id2well{88} = "H11";
    $global_id2well{89} = "A12";
    $global_id2well{90} = "B12";
    $global_id2well{91} = "C12";
    $global_id2well{92} = "D12";
    $global_id2well{93} = "E12";
    $global_id2well{94} = "F12";
    $global_id2well{95} = "G12";
    $global_id2well{96} = "H12";

    $global_well2id{"A1"} = 1;
    $global_well2id{"B1"} = 2;
    $global_well2id{"C1"} = 3;
    $global_well2id{"D1"} = 4;
    $global_well2id{"E1"} = 5;
    $global_well2id{"F1"} = 6;
    $global_well2id{"G1"} = 7;
    $global_well2id{"H1"} = 8;
    $global_well2id{"A2"} = 9;
    $global_well2id{"B2"} = 10;
    $global_well2id{"C2"} = 11;
    $global_well2id{"D2"} = 12;
    $global_well2id{"E2"} = 13;
    $global_well2id{"F2"} = 14;
    $global_well2id{"G2"} = 15;
    $global_well2id{"H2"} = 16;
    $global_well2id{"A3"} = 17;
    $global_well2id{"B3"} = 18;
    $global_well2id{"C3"} = 19;
    $global_well2id{"D3"} = 20;
    $global_well2id{"E3"} = 21;
    $global_well2id{"F3"} = 22;
    $global_well2id{"G3"} = 23;
    $global_well2id{"H3"} = 24;
    $global_well2id{"A4"} = 25;
    $global_well2id{"B4"} = 26;
    $global_well2id{"C4"} = 27;
    $global_well2id{"D4"} = 28;
    $global_well2id{"E4"} = 29;
    $global_well2id{"F4"} = 30;
    $global_well2id{"G4"} = 31;
    $global_well2id{"H4"} = 32;
    $global_well2id{"A5"} = 33;
    $global_well2id{"B5"} = 34;
    $global_well2id{"C5"} = 35;
    $global_well2id{"D5"} = 36;
    $global_well2id{"E5"} = 37;
    $global_well2id{"F5"} = 38;
    $global_well2id{"G5"} = 39;
    $global_well2id{"H5"} = 40;
    $global_well2id{"A6"} = 41;
    $global_well2id{"B6"} = 42;
    $global_well2id{"C6"} = 43;
    $global_well2id{"D6"} = 44;
    $global_well2id{"E6"} = 45;
    $global_well2id{"F6"} = 46;
    $global_well2id{"G6"} = 47;
    $global_well2id{"H6"} = 48;
    $global_well2id{"A7"} = 49;
    $global_well2id{"B7"} = 50;
    $global_well2id{"C7"} = 51;
    $global_well2id{"D7"} = 52;
    $global_well2id{"E7"} = 53;
    $global_well2id{"F7"} = 54;
    $global_well2id{"G7"} = 55;
    $global_well2id{"H7"} = 56;
    $global_well2id{"A8"} = 57;
    $global_well2id{"B8"} = 58;
    $global_well2id{"C8"} = 59;
    $global_well2id{"D8"} = 60;
    $global_well2id{"E8"} = 61;
    $global_well2id{"F8"} = 62;
    $global_well2id{"G8"} = 63;
    $global_well2id{"H8"} = 64;
    $global_well2id{"A9"} = 65;
    $global_well2id{"B9"} = 66;
    $global_well2id{"C9"} = 67;
    $global_well2id{"D9"} = 68;
    $global_well2id{"E9"} = 69;
    $global_well2id{"F9"} = 70;
    $global_well2id{"G9"} = 71;
    $global_well2id{"H9"} = 72;
    $global_well2id{"A10"} = 73;
    $global_well2id{"B10"} = 74;
    $global_well2id{"C10"} = 75;
    $global_well2id{"D10"} = 76;
    $global_well2id{"E10"} = 77;
    $global_well2id{"F10"} = 78;
    $global_well2id{"G10"} = 79;
    $global_well2id{"H10"} = 80;
    $global_well2id{"A11"} = 81;
    $global_well2id{"B11"} = 82;
    $global_well2id{"C11"} = 83;
    $global_well2id{"D11"} = 84;
    $global_well2id{"E11"} = 85;
    $global_well2id{"F11"} = 86;
    $global_well2id{"G11"} = 87;
    $global_well2id{"H11"} = 88;
    $global_well2id{"A12"} = 89;
    $global_well2id{"B12"} = 90;
    $global_well2id{"C12"} = 91;
    $global_well2id{"D12"} = 92;
    $global_well2id{"E12"} = 93;
    $global_well2id{"F12"} = 94;
    $global_well2id{"G12"} = 95;
    $global_well2id{"H12"} = 96;
}

######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    my @standard_options = ( "prefix|p:s", "help|h+", "in|i:s", "sample_tube_volume:i", "sample_tube_capacity:i", "sample_DNA_required_total:i", "pool_tube_position:s", "pool_minimum_volume:i", "pool_tube_capacity:i", "dilutant_tube_position:s", "dilutant_minimum_volume:i", "dilutant_tube_capacity:i");
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
    if(!exists $options{'in'} ) { print "**ERROR: You need to supply a csv file to parse\n"; exec("pod2usage $0"); }
    #if(!exists $options{''} ) { print "**ERROR: \n"; exec("pod2usage $0"); }

    return \%options;
}

sub printAtStart {
print<<"EOF";
---------------------------------------------------------------- 
 $0
 Version $VERSION
 Copyright (C) 2011 Michael Imelfort and Paul Dennis
    
 This program comes with ABSOLUTELY NO WARRANTY;
 This is free software, and you are welcome to redistribute it
 under certain conditions: See the source for more details.
---------------------------------------------------------------- 
EOF
}

__DATA__

=head1 NAME

    app_csv2epi.pl

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

   Convert a csv file to another csv file which is formatted for use in the robot

=head1 SYNOPSIS

    app_csv2epi.pl -in|i CSV_FILE [-help|h]

    Convert a csv file to another csv file which is formatted for use in the robot
    
      -in -i CSV_FILE                               File to parse -- DO NOT INCLUDE HEADER IN FILE
      -prefix -p NAME                               Prefix to attach to output files (default: NONE)
      [-help -h]                                    Displays basic usage information
      
    Source rack options:
      
      [-sample_tube_volume AMOUNT (uL) ]            The amount of material in the sample tube (default: 13 uL)
      [-sample_tube_capacity AMOUNT (uL) ]          The capacity of the sample tube (default: 300 uL)
      [-sample_DNA_required_total AMOUNT (ng) ]     The amount of DNA required in the pool for EACH SAMPLE (default: 13 uL)
     
    Pool options:
      
      [-pool_tube_position POSITION ]               The position in the rack of the pooling tube (default: A1)
      [-pool_minimum_volume AMOUNT (uL) ]           The minimum amount of material we take from each well during pooling  (default 1 uL)
      [-pool_tube_capacity AMOUNT (uL) ]            The capacity of the pooling tube (default: 2000 uL)
      
    Dilution options:
      
      [-dilutant_tube_position POSITION ]           The position in the rack of the dilution tube (default: A2)
      [-dilutant_minimum_volume AMOUNT (uL) ]       The minimum amount we add to dilute (default 2 uL)
      [-dilutant_tube_capacity AMOUNT (uL) ]        The capacity of the dilution tube (default: 2000 uL)
         
=cut


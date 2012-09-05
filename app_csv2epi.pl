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
    $global_id2well{2} = "A2";
    $global_id2well{3} = "A3";
    $global_id2well{4} = "A4";
    $global_id2well{5} = "A5";
    $global_id2well{6} = "A6";
    $global_id2well{7} = "A7";
    $global_id2well{8} = "A8";
    $global_id2well{9} = "A9";
    $global_id2well{10} = "A10";
    $global_id2well{11} = "A11";
    $global_id2well{12} = "A12";
    $global_id2well{13} = "B1";
    $global_id2well{14} = "B2";
    $global_id2well{15} = "B3";
    $global_id2well{16} = "B4";
    $global_id2well{17} = "B5";
    $global_id2well{18} = "B6";
    $global_id2well{19} = "B7";
    $global_id2well{20} = "B8";
    $global_id2well{21} = "B9";
    $global_id2well{22} = "B10";
    $global_id2well{23} = "B11";
    $global_id2well{24} = "B12";
    $global_id2well{25} = "C1";
    $global_id2well{26} = "C2";
    $global_id2well{27} = "C3";
    $global_id2well{28} = "C4";
    $global_id2well{29} = "C5";
    $global_id2well{30} = "C6";
    $global_id2well{31} = "C7";
    $global_id2well{32} = "C8";
    $global_id2well{33} = "C9";
    $global_id2well{34} = "C10";
    $global_id2well{35} = "C11";
    $global_id2well{36} = "C12";
    $global_id2well{37} = "D1";
    $global_id2well{38} = "D2";
    $global_id2well{39} = "D3";
    $global_id2well{40} = "D4";
    $global_id2well{41} = "D5";
    $global_id2well{42} = "D6";
    $global_id2well{43} = "D7";
    $global_id2well{44} = "D8";
    $global_id2well{45} = "D9";
    $global_id2well{46} = "D10";
    $global_id2well{47} = "D11";
    $global_id2well{48} = "D12";
    $global_id2well{49} = "E1";
    $global_id2well{50} = "E2";
    $global_id2well{51} = "E3";
    $global_id2well{52} = "E4";
    $global_id2well{53} = "E5";
    $global_id2well{54} = "E6";
    $global_id2well{55} = "E7";
    $global_id2well{56} = "E8";
    $global_id2well{57} = "E9";
    $global_id2well{58} = "E10";
    $global_id2well{59} = "E11";
    $global_id2well{60} = "E12";
    $global_id2well{61} = "F1";
    $global_id2well{62} = "F2";
    $global_id2well{63} = "F3";
    $global_id2well{64} = "F4";
    $global_id2well{65} = "F5";
    $global_id2well{66} = "F6";
    $global_id2well{67} = "F7";
    $global_id2well{68} = "F8";
    $global_id2well{69} = "F9";
    $global_id2well{70} = "F10";
    $global_id2well{71} = "F11";
    $global_id2well{72} = "F12";
    $global_id2well{73} = "G1";
    $global_id2well{74} = "G2";
    $global_id2well{75} = "G3";
    $global_id2well{76} = "G4";
    $global_id2well{77} = "G5";
    $global_id2well{78} = "G6";
    $global_id2well{79} = "G7";
    $global_id2well{80} = "G8";
    $global_id2well{81} = "G9";
    $global_id2well{82} = "G10";
    $global_id2well{83} = "G11";
    $global_id2well{84} = "G12";
    $global_id2well{85} = "H1";
    $global_id2well{86} = "H2";
    $global_id2well{87} = "H3";
    $global_id2well{88} = "H4";
    $global_id2well{89} = "H5";
    $global_id2well{90} = "H6";
    $global_id2well{91} = "H7";
    $global_id2well{92} = "H8";
    $global_id2well{93} = "H9";
    $global_id2well{94} = "H10";
    $global_id2well{95} = "H11";
    $global_id2well{96} = "H12";


    $global_well2id{"A1"} = 1;
    $global_well2id{"A2"} = 2;
    $global_well2id{"A3"} = 3;
    $global_well2id{"A4"} = 4;
    $global_well2id{"A5"} = 5;
    $global_well2id{"A6"} = 6;
    $global_well2id{"A7"} = 7;
    $global_well2id{"A8"} = 8;
    $global_well2id{"A9"} = 9;
    $global_well2id{"A10"} = 10;
    $global_well2id{"A11"} = 11;
    $global_well2id{"A12"} = 12;
    $global_well2id{"B1"} = 13;
    $global_well2id{"B2"} = 14;
    $global_well2id{"B3"} = 15;
    $global_well2id{"B4"} = 16;
    $global_well2id{"B5"} = 17;
    $global_well2id{"B6"} = 18;
    $global_well2id{"B7"} = 19;
    $global_well2id{"B8"} = 20;
    $global_well2id{"B9"} = 21;
    $global_well2id{"B10"} = 22;
    $global_well2id{"B11"} = 23;
    $global_well2id{"B12"} = 24;
    $global_well2id{"C1"} = 25;
    $global_well2id{"C2"} = 26;
    $global_well2id{"C3"} = 27;
    $global_well2id{"C4"} = 28;
    $global_well2id{"C5"} = 29;
    $global_well2id{"C6"} = 30;
    $global_well2id{"C7"} = 31;
    $global_well2id{"C8"} = 32;
    $global_well2id{"C9"} = 33;
    $global_well2id{"C10"} = 34;
    $global_well2id{"C11"} = 35;
    $global_well2id{"C12"} = 36;
    $global_well2id{"D1"} = 37;
    $global_well2id{"D2"} = 38;
    $global_well2id{"D3"} = 39;
    $global_well2id{"D4"} = 40;
    $global_well2id{"D5"} = 41;
    $global_well2id{"D6"} = 42;
    $global_well2id{"D7"} = 43;
    $global_well2id{"D8"} = 44;
    $global_well2id{"D9"} = 45;
    $global_well2id{"D10"} = 46;
    $global_well2id{"D11"} = 47;
    $global_well2id{"D12"} = 48;
    $global_well2id{"E1"} = 49;
    $global_well2id{"E2"} = 50;
    $global_well2id{"E3"} = 51;
    $global_well2id{"E4"} = 52;
    $global_well2id{"E5"} = 53;
    $global_well2id{"E6"} = 54;
    $global_well2id{"E7"} = 55;
    $global_well2id{"E8"} = 56;
    $global_well2id{"E9"} = 57;
    $global_well2id{"E10"} = 58;
    $global_well2id{"E11"} = 59;
    $global_well2id{"E12"} = 60;
    $global_well2id{"F1"} = 61;
    $global_well2id{"F2"} = 62;
    $global_well2id{"F3"} = 63;
    $global_well2id{"F4"} = 64;
    $global_well2id{"F5"} = 65;
    $global_well2id{"F6"} = 66;
    $global_well2id{"F7"} = 67;
    $global_well2id{"F8"} = 68;
    $global_well2id{"F9"} = 69;
    $global_well2id{"F10"} = 70;
    $global_well2id{"F11"} = 71;
    $global_well2id{"F12"} = 72;
    $global_well2id{"G1"} = 73;
    $global_well2id{"G2"} = 74;
    $global_well2id{"G3"} = 75;
    $global_well2id{"G4"} = 76;
    $global_well2id{"G5"} = 77;
    $global_well2id{"G6"} = 78;
    $global_well2id{"G7"} = 79;
    $global_well2id{"G8"} = 80;
    $global_well2id{"G9"} = 81;
    $global_well2id{"G10"} = 82;
    $global_well2id{"G11"} = 83;
    $global_well2id{"G12"} = 84;
    $global_well2id{"H1"} = 85;
    $global_well2id{"H2"} = 86;
    $global_well2id{"H3"} = 87;
    $global_well2id{"H4"} = 88;
    $global_well2id{"H5"} = 89;
    $global_well2id{"H6"} = 90;
    $global_well2id{"H7"} = 91;
    $global_well2id{"H8"} = 92;
    $global_well2id{"H9"} = 93;
    $global_well2id{"H10"} = 94;
    $global_well2id{"H11"} = 95;
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


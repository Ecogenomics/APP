package AppCommon;
require Exporter;
use File::Basename;
use File::Temp qw(tempfile);
use Carp;

use AppConfig;

our @ISA = qw(Exporter);
our @EXPORT=qw(
    %known_config_options
    create_unique_analysis_dir
    create_data_file_links
    parse_app_config_file
    create_analysis_config_file
    create_qiime_mapping_file
    convert_hash_to_array
    get_read_counts_from_cd_hit_otu
    get_read_counts_from_OTU_table
);

our %known_config_options = (
    NORMALISE => ['SEQ'],
    MUL_RARE_M => [],
    MUL_RARE_X => [],
    MUL_RARE_S => [],
    MUL_RARE_N => [],
    PIPELINE => []
);

################################################################################
# Subroutine: create_unique_analysis_dir($output_dir)
# 
################################################################################

sub create_unique_analysis_dir {
    
    my $output_dir = shift;
    
    if (! defined($output_dir)) {
        my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst)
            = localtime(time);
            
        $output_dir = sprintf("app_analysis_%04i%02i%02i", $year+1900, $mon+1, $mday);
        
        if (-e $output_dir) {
            my $folder_suffix = 1;
            while (-e $output_dir . "_" . $folder_suffix) {
                $folder_suffix++;
            }
            $output_dir .= "_" . $folder_suffix;
        }
    }

    mkdir($output_dir) or die "Unable to create output folder: $output_dir - $!";
    return $output_dir;
}

################################################################################
# Subroutine: create_data_file_links($output_dir, $job_ID)
#
################################################################################

sub create_data_file_links {
    my ($output_dir, $job_ID) = @_;
    
    my @fasta_files;
    foreach my $suffix ((".fa", ".fna", ".fasta")) {
        if (-e "$job_ID$suffix") {
            push @fasta_files, "$job_ID$suffix";
        }
    }
    if (scalar @fasta_files) {
        if (scalar @fasta_files > 1) {
            croak "ERROR: Too many possible FASTA files to choose, remove or rename " .
            "to stop ambiguity. Offending files: " . join(", ", @fasta_files) . "\n";
        }
    } else {
        croak "ERROR: Unable to find fasta file: $job_ID.fna\n";
    }
    
    if (! -e "$job_ID.qual") {
         croak "ERROR: Unable to find qual file: $job_ID.qual\n";
    }
    
    symlink "../$job_ID.qual", "$output_dir/raw_sequences.qual";
    symlink "../$fasta_files[0]", "$output_dir/raw_sequences.fasta";
}

################################################################################
# Subroutine: parse_app_config_file($config_file)
# 
################################################################################

sub parse_app_config_file {
    my $config_file = shift;
    open(my $fh, $config_file);
    
    my $in_config_section = 0;
    my $config_array = [];
    my $sample_array = [];
    while (my $line = <$fh>) {
        chomp $line;
        
        # Skip the line if starts with a hash(#)
        if ($line =~ /^(\s)*#/) {next};
        
        # Are we in the config section yet?
        if ($line =~ /@@/) {
            $in_config_section = 1;
            next;
        }
        
        # If we are in the config section...
        if ($in_config_section) {
            my @splitline = split /=/, $line;
            
            # Check if the option is known/allowed.
            if (exists($known_config_options{$splitline[0]})) {
                my @possible_choices = @{$known_config_options{$splitline[0]}};
                my %options = map { $_ => 1 } @possible_choices;
                
                # Check if the choices for this option are limited.
                if ((! scalar @possible_choices) || exists ($options{$splitline[1]})) {
                    
                    # Check if the choice is defined.
                    if (! defined ($splitline[1]) || $splitline[1] =~ /^\s*$/) {
                        print "No config setting for option " .
                        $splitline[0] . " in $config_file: - ignoring option.\n";
                    } else {
                        push @{$config_array}, [$splitline[0], $splitline[1]];
                    }
                } else {
                    print "Unknown config setting for option " .
                        $splitline[0] . " in $config_file: " .
                        $splitline[1] . " - ignoring option.\n";
                }
            } else {
                print "Unknown config option in $config_file: " .
                    $splitline[0] . " - ignoring option.\n";
            }
        # If we aren't in the config section yet.            
        } else {
            my @splitline = split /\t/, $line;
            my @details = @splitline[0..3];
            if (scalar @splitline > 4) {
                push @details, $splitline[$#splitline];
            }
            push @{$sample_array}, \@details;
        }
    }
    close($fh);
    return ($sample_array, $config_array);
}

################################################################################
# Subroutine: create_analysis_config_file($config_filename, $config_array)
# 
################################################################################

sub create_analysis_config_file {
    my ($config_filename, $config_array, $sample_list) = @_;
    open (my $fh, ">", $config_filename);
    print {$fh} "#DO NOT EDIT THIS FILE, EVEN AT GUNPOINT! SERIOUSLY!\n";
    
    if (scalar @{$sample_list}) {
        print {$fh} "SAMPLES=" . join(",", @{$sample_list}), "\n";
    }
    
    foreach my $config_pair (@{$config_array}) {
        print {$fh} $config_pair->[0] . "=" . $config_pair->[1] . "\n";
    }

    close($fh);
}

################################################################################
# Subroutine: create_qiime_mapping_file($qiime_mapping_filename, $sample_array)
# 
################################################################################

sub create_qiime_mapping_file {
    my ($qiime_mapping_filename, $sample_array) = @_;
    open (my $fh, ">", $qiime_mapping_filename);
    print {$fh} "$FNB_HEADER\n";
    foreach my $sample_array_ptr (@{$sample_array}) {
        print {$fh} join("\t", @{$sample_array_ptr}), "\n";
    }
    close($fh);
}

sub get_read_counts_from_cd_hit_otu {
    my ($cd_hit_otu_file) = @_;
    open(my $fh, $cd_hit_otu_file);
    <$fh>; #burn headers
    my @sample_counts;
    while (my $line = <$fh>) {
        chomp $line;
        my @splitline = split /\t/, $line;
        push @sample_counts, [$splitline[0], $splitline[$#splitline - 1]];
    }
    close($fh);
    return @sample_counts;
};

sub get_read_counts_from_OTU_table {
    my ($otu_file) = @_;
    open(my $fh, $otu_file);
    my @headers;
    my @counts;
    while (my $line = <$fh>) {
        if ($line =~ /^#/) {
            if ($line =~ /^#OTU\ ID/) {
                my @splitline = split /\t/, $line;
                @headers = @splitline[1..($#splitline-1)];
                @counts = 0 x (scalar @splitline - 2);
            } else {
                next;
            }
        } else {
            my @splitline = split /\t/, $line;
            for (my $i = 0; $i < scalar @splitline - 2; $i++) {
                $counts[$i] += $splitline[$i + 1];
            }
        }
    }
    my @sample_counts;
    for (my $i = 0; $i < scalar @headers; $i++) {
        push @sample_counts, [$headers[$i], $counts[$i]];
    }
    return \@sample_counts;
}

sub convert_hash_to_array {
    my ($config_hash) = @_;
    my @config_array = map {[$_, $config_hash->{$_}]} keys %{$config_hash};
    return \@config_array;
}

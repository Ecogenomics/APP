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
use File::Temp qw(tempfile);
use List::Util qw(min max);
use Statistics::R;

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

########################################################################
# Globals
########################################################################

#
# Here we keep the default acacia config hash
#
our %acacia_config_hash = (
     ANY_DIFF_SIGNIFICANT_FOR_TWO_SEQS => "TRUE",
     AVG_QUALITY_CUTOFF => 30,
     ERROR_MODEL => "Balzer",
     FASTA => "TRUE",
     FASTA_LOCATION => "good.fasta",
     FASTQ => "FALSE",
     FASTQ_LOCATION => "null",
     FILTER_N_BEFORE_POS => 250,
     FLOW_KEY => "TCAG",
     MAXIMUM_MANHATTAN_DISTANCE => 13,
     MAX_RECURSE_DEPTH => 2,
     MAX_STD_DEV_LENGTH => 2,
     MID_FILE => ".dummmy",
     MID_OPTION => "NO_MID",
     MIN_FLOW_TRUNCATION => 150,
     MIN_READ_REP_BEFORE_TRUNCATION => 0.0,
     OUTPUT_DIR => "denoised_acacia",
     OUTPUT_PREFIX => "acacia_out",
     QUAL_LOCATION => "null",
     REPRESENTATIVE_SEQUENCE => "Mode",
     SIGNIFICANCE_LEVEL => -9,
     SPLIT_ON_MID => "FALSE",
     TRIM_TO_LENGTH => 250,
     TRUNCATE_READ_TO_FLOW => ""
);

######################################################################
# Main
######################################################################

my $options = checkParams();

chdir($options->{'d'}) or die "No such directory: " . $options->{'d'} . "\n";
printAtStart();
#qiime_standard_pipeline();
my $config_hash;

open(my $fh, "config.txt");
while (my $line = <$fh>) {
    if ($line =~ /^(\s*)#/) {next};
    if ($line =~ /^(.*)=(.*)$/) {
        $config_hash->{$1} = $2;
    }
}
close($fh);

if (! defined ($config_hash->{PIPELINE})) {
    $config_hash->{PIPELINE} = 'CD_HIT_OTU';
}

my @pipeline = split /,/, $config_hash->{PIPELINE};

if ($pipeline[0] =~ /^(\s*)CD_HIT_OTU(\s*)$/) {
    cd_hit_otu_pipeline($config_hash);
} elsif ($pipeline[0] =~ /^(\s*)QIIME(\s*)$/) {
    qiime_pipeline($config_hash);
} else {
    print "Unknown pipeline: " . $pipeline[0] . "\n";
}

######################################################################
# Globals
######################################################################


######################################################################
# Pipelines
######################################################################


######################################################################
# QIIME Standard Pipeline
#
# This is the pipeline for standard analysis using QIIME.
#
######################################################################

sub qiime_pipeline {
    
    my ($config_hash) = @_;
    
    my $pipeline_modifiers = get_pipeline_modifiers($config_hash->{PIPELINE});
    
    my $GG_OTUS_FASTA = "/srv/whitlam/bio/db/gg/from_www.secongenome.com/2012_10/gg_12_10_otus/rep_set/97_otus.fasta";
    
    my $processing_dir = 'processing';
    my $results_dir = 'results';
    
    mkdir($processing_dir) if (! -e $processing_dir);
    mkdir($results_dir) if (! -e $results_dir);
    
    if (! exists ($config_hash->{TRIM_LENGTH})) {
        $config_hash->{TRIM_LENGTH} = 250;
    }
    
    my @starting_samples = split /,/, $config_hash->{SAMPLES};
    
    my %sample_counts;
    if (-e "sample_counts.txt") {
        %sample_counts = %{get_read_counts_from_sample_counts("sample_counts.txt")}
    }
    
    my %sample_for_analysis_hash;
    
    my $need_to_reanalyse = 0;
    if (! -e "sample_exclusion.txt") {
        $config_hash->{SAMPLES_FOR_ANALYSIS} = $config_hash->{SAMPLES};
        create_analysis_config_file("config.txt", convert_hash_to_array($config_hash));
        %sample_for_analysis_hash = map {$_ => 1} @starting_samples;
        
    } else {
        my @previous_analysis = split /,/, $config_hash->{SAMPLES_FOR_ANALYSIS};
        my %samples_to_use = %{read_sample_exclusion("sample_exclusion.txt")};
        my @sample_for_analysis_array;
        
        foreach my $sample_name (@starting_samples) {
            if ((defined $samples_to_use{$sample_name}) && $samples_to_use{$sample_name}) {
                push @sample_for_analysis_array, $sample_name;
            }
        }
        
        $config_hash->{SAMPLES_FOR_ANALYSIS} = join ",", @sample_for_analysis_array;
        create_analysis_config_file("config.txt", convert_hash_to_array($config_hash));
        
        %sample_for_analysis_hash = map {$_ => 1} @sample_for_analysis_array;
        foreach my $sample_name (@previous_analysis) {
            if ((! exists ($sample_for_analysis_hash{$sample_name})) && ($sample_counts{$sample_name} != 0)) {
                $need_to_reanalyse = 1;
                last;
            }
        }
    }
    
    # If we are running for the first time, or if sample_exclusions.txt has been modified,
    # reset the pipeline. 
    
    if ((! -e $global_qiime_mapping_filename) || ($need_to_reanalyse) || (! exists($config_hash->{PIPELINE_STAGE}))) {
        $config_hash->{PIPELINE_STAGE} = "PREAMBLE";
    }     
        
    if ($config_hash->{PIPELINE_STAGE} eq "PREAMBLE") {
        
    # (Re)create the QIIME mapping file, taking sample_exclusions.txt into account.
        
        my ($sample_array, $config_array) =
            parse_app_config_file("app.config");
        
        my %sample_to_analyse = map {$_ => 1} split /,/, $config_hash->{SAMPLES_FOR_ANALYSIS};
        
        my $new_sample_array;
        foreach my $sample_array_ptr (@{$sample_array}) {
            if (exists $sample_to_analyse{$sample_array_ptr->[0]}) {
                push @{$new_sample_array}, $sample_array_ptr;
            }
        }
        
        create_qiime_mapping_file($global_qiime_mapping_filename, $new_sample_array);
        
        # Begin the pipeline...
        
        my $params = [{-M => 1,
                    -q => 'raw_sequences.qual',
                    -f => 'raw_sequences.fasta',
                    -d => '',
                    -m => $global_qiime_mapping_filename,
                    -b => 'variable_length',
                    -a => 2,
                    -H => 10,
                    -l => $config_hash->{TRIM_LENGTH}}];
        
        if (exists($pipeline_modifiers->{STRICT_FILTERING})) {
            
            $params = [{-M => 1,
                        -q => 'raw_sequences.qual',
                        -f => 'raw_sequences.fasta',
                        -d => '',
                        -m => $global_qiime_mapping_filename,
                        -b => 'variable_length',
                        -a => 2,
                        -H => 10,
                        -l => $config_hash->{TRIM_LENGTH},
                        -s => 25,
                        -w => 10}];
            
        } elsif (exists($pipeline_modifiers->{MEDIUM_FILTERING})) {
            
            $params = [{-M => 1,
                        -q => 'raw_sequences.qual',
                        -f => 'raw_sequences.fasta',
                        -d => '',
                        -m => $global_qiime_mapping_filename,
                        -b => 'variable_length',
                        -a => 2,
                        -H => 10,
                        -l => $config_hash->{TRIM_LENGTH},
                        -s => 20,
                        -w => 10}];
        }
        
        splitLibraries($params);
        
        create_analysis_config_file("config.txt", convert_hash_to_array($config_hash));
        
        if (-z 'seqs.fna') {
          croak "ERROR: No sequences in seqs.fna (no sequences successfully demultiplexed after split_libraries.py).\n" .
          "Check the config file that the barcode/primer sequences are correct.\n"
        }
    
        $config_hash->{PIPELINE_STAGE} = "UCLUST";
        
        create_analysis_config_file("config.txt", convert_hash_to_array($config_hash));
    }
    
    if ($config_hash->{PIPELINE_STAGE} eq "UCLUST") {
    
        my $params = [{'-uchime_ref' => 'seqs.fna'},
                 {'-db' => $QIIME_TAX_blast_file,
                  '-strand' => 'both',
                  '-threads' => 10,
                  '-nonchimeras' => 'good.fasta',
                  '-chimeras' => 'bad.fasta'}];
                  
        uclustRemoveChimeras($params);
        
        $config_hash->{PIPELINE_STAGE} = "ACACIA";
        
        create_analysis_config_file("config.txt", convert_hash_to_array($config_hash));
    }
    
    my $fasta_output_file = 'good.fasta';
    
    if ($config_hash->{PIPELINE_STAGE} eq "ACACIA") {
        
        if (! exists($pipeline_modifiers->{NO_ACACIA})) {
        
            $acacia_config_hash{TRIM_TO_LENGTH} = $config_hash->{TRIM_LENGTH};
        
            my $config_file = create_acacia_config_file(\%acacia_config_hash);
            
            my $params_array = [["-XX:+UseConcMarkSweepGC",
                                 "-Xmx100G"],
                                 {-jar => '$ACACIA'},
                                 {-c => $config_file}];
            
            mkdir('denoised_acacia');
            
            if (run_acacia($params_array)) {
                die "Acacia failed.\n";
            };
            
            unlink($config_file);
            
            $fasta_output_file = "denoised_acacia/acacia_out_all_tags.seqOut";
        
        }
    
        $config_hash->{PIPELINE_STAGE} = "PICK_OTUS";
    
        create_analysis_config_file("config.txt", convert_hash_to_array($config_hash));
        
    }
    
    if ($config_hash->{PIPELINE_STAGE} eq "PICK_OTUS") {
    
        print "Picking OTUs for non normalised data set...\n";
        checkAndRunCommand("pick_otus.py", [{-i => $fasta_output_file,
                                             -s => 0.97,
                                             -o => "$processing_dir/uclust_picked_otus"}], DIE_ON_FAILURE);
        
        $config_hash->{PIPELINE_STAGE} = "REP_SET";
            
        create_analysis_config_file("config.txt", convert_hash_to_array($config_hash));
    }
    
    my @possible_otu_files = glob("$processing_dir/uclust_picked_otus/*_otus.txt");
    if (scalar @possible_otu_files != 1) {
        if (scalar @possible_otu_files) {
            die "Too many possible uclust otu files.\n";
        } else {
            die "No OTU files found.\n";
        }
    }
   
    if ($config_hash->{PIPELINE_STAGE} eq "REP_SET") {
    
        print "Getting a representative set...\n";
        
        
        checkAndRunCommand("pick_rep_set.py", [{-i => $possible_otu_files[0],
                                                -f => $fasta_output_file,
                                                -o => "$processing_dir/rep_set.fa"}], DIE_ON_FAILURE);
        
        $config_hash->{PIPELINE_STAGE} = "ASSIGN_TAXONOMY";
           
        create_analysis_config_file("config.txt", convert_hash_to_array($config_hash));
    }
 
    if ($config_hash->{PIPELINE_STAGE} eq "ASSIGN_TAXONOMY") {
    
        my $TAX_tax_file = $QIIME_TAX_tax_file;
        my $TAX_blast_file = $QIIME_TAX_blast_file;
        my $TAX_aligned_blast_file = $QIIME_TAX_aligned_blast_file;
        my $imputed_file = $QIIME_imputed_file;
        
        checkAndRunCommand("assign_taxonomy.py", [{-i => "$processing_dir/rep_set.fa",
                                                   -t => $TAX_tax_file,
                                                   -b => $TAX_blast_file,
                                                   -m => "blast",
                                                   -a => 10,
                                                   -e => 0.001,
                                                   -o => $processing_dir}], DIE_ON_FAILURE);
        
        $config_hash->{PIPELINE_STAGE} = "GENERATE_OTU_TABLE";
        create_analysis_config_file("config.txt", convert_hash_to_array($config_hash));
    
    }
    
    if ($config_hash->{PIPELINE_STAGE} eq "GENERATE_OTU_TABLE") {
        
        checkAndRunCommand("make_otu_table.py",  [{-i => $possible_otu_files[0],
                                                   -t => "$processing_dir/rep_set_tax_assignments.txt",
                                                   -o => "$results_dir/non_normalised_otu_table.txt"}], DIE_ON_FAILURE);
        
        
        checkAndRunCommand("reformat_otu_table.py",  [{-i => "$results_dir/non_normalised_otu_table.txt",
                                                       -t => "$processing_dir/rep_set_tax_assignments.txt",
                                                       -o => "$results_dir/non_normalised_otu_table_expanded.tsv"}], IGNORE_FAILURE);
        
        my $otu_sample_counts = get_read_counts_from_OTU_table("$results_dir/non_normalised_otu_table.txt");
        
        use Data::Dumper;
        print Dumper($otu_sample_counts);
        
        # If this is the first run through, create a file to record the sample_counts.
        if (! -e "sample_counts.txt") {
            open($fh, ">sample_counts.txt");
            print {$fh} "#NAME\tREAD COUNT\n";
            foreach my $sample_name (@starting_samples) {
                if (! defined $otu_sample_counts->{$sample_name}) {
                    $otu_sample_counts->{$sample_name} = 0; 
                }
                print {$fh} $sample_name. "\t" . $otu_sample_counts->{$sample_name} ."\n";
                $sample_counts{$sample_name} = $otu_sample_counts->{$sample_name};
            }
            close($fh);
        }
                
        my %previous_analysis_hash = map {$_ => 1} split /,/, $config_hash->{SAMPLES_FOR_ANALYSIS};
       
        open($fh, ">sample_exclusion.txt");
        print {$fh} "#NAME\tREAD COUNT\tUSE\n";
        foreach my $sample_name (@starting_samples) {
            my $use = 0;
            if ($previous_analysis_hash{$sample_name} && $sample_counts{$sample_name}) {
                $use = 1;
            }
            print {$fh} $sample_name. "\t" . $sample_counts{$sample_name} . "\t" . $use . "\n";
        }
        close($fh); 
                
        print "Non-normalized OTU table has now been generated.\n";
        print "\nCheck the " . $options->{'d'} . " sample_exclusion.txt and choose which samples you wish to exclude.\n";
        print "When you have chosen which samples to exclude, run the following command:\n";
        print "      app_run_analysis.pl -d " . $options->{'d'} . "\n\n";
    
        $config_hash->{PIPELINE_STAGE} = "NORMALISATION";
            
        create_analysis_config_file("config.txt", convert_hash_to_array($config_hash));
        
        exit();
    
    }
    
    if ($config_hash->{PIPELINE_STAGE} eq "NORMALISATION") {

        use Data::Dumper;
        
        print Dumper(\%sample_counts);
        
        my $min_samples = min(map {$sample_counts{$_}} keys %sample_for_analysis_hash);
        
        my $global_norm_sample_size = (int($min_samples / 50) - 1) * 50;
        
        if ($global_norm_sample_size <= 0) {
            print "Unable to normalise, some samples have too few reads. Review " . $options->{'d'} . "/sample_exclusion.txt " .
                  "and rerun app_run_analysis.pl -d " . $options->{'d'} . "\n\n";
        }
        
        print "Correcting normalisation depth to $global_norm_sample_size sequences...\n";
        
        my $global_norm_num_reps = 1000;
        
        print "Normalizing non normalised table at $global_norm_sample_size sequences... [$global_norm_sample_size, $global_norm_num_reps]\n";
                
        my $params = [{-i => "$results_dir/non_normalised_otu_table.txt",
                       -o => "$processing_dir/rare_tables/",
                       -d => $global_norm_sample_size,
                       -n => $global_norm_num_reps},
                       ["--lineages_included"],
                       ["--k"]];
        
        my $normalised_OTU_file = 
            normalise_OTU_table($params, $processing_dir, $global_norm_sample_size,
                                scalar keys %sample_counts);
            
        `cp $normalised_OTU_file $results_dir/normalised_otu_table.txt`;
        
        checkAndRunCommand("reformat_otu_table.py",  [{-i => "$results_dir/normalised_otu_table.txt",
                                                       -t => "$processing_dir/rep_set_tax_assignments.txt",
                                                       -o => "$results_dir/normalised_otu_table_expanded.tsv"}], IGNORE_FAILURE);

    }
    
}

######################################################################
# CD-HIT-OTU Pipeline
#
# This is the pipeline for analysis using CD-HIT-OTU.
#
######################################################################

sub cd_hit_otu_pipeline {
    
    my ($config_hash) = @_;
    my $pipeline_modifiers = get_pipeline_modifiers($config_hash->{PIPELINE});
    
    my $GG_OTUS_FASTA = "/srv/whitlam/bio/db/gg/from_www.secongenome.com/2012_10/gg_12_10_otus/rep_set/97_otus.fasta";
    
    my $processing_dir = 'processing';
    my $results_dir = 'results';
    
    mkdir($processing_dir) if (! -e $processing_dir);
    mkdir($results_dir) if (! -e $results_dir);
    
    if (! exists ($config_hash->{TRIM_LENGTH})) {
        $config_hash->{TRIM_LENGTH} = 250;
    }
    
    my @starting_samples = split /,/, $config_hash->{SAMPLES};
    
    my %sample_counts;
    if (-e "sample_counts.txt") {
        %sample_counts = %{get_read_counts_from_sample_counts("sample_counts.txt")}
    }
    
    my %sample_for_analysis_hash;
    
    my $need_to_reanalyse = 0;
    if (! -e "sample_exclusion.txt") {
        $config_hash->{SAMPLES_FOR_ANALYSIS} = $config_hash->{SAMPLES};
        create_analysis_config_file("config.txt", convert_hash_to_array($config_hash));
        %sample_for_analysis_hash = map {$_ => 1} @starting_samples;
        
    } else {
        my @previous_analysis = split /,/, $config_hash->{SAMPLES_FOR_ANALYSIS};
        my %samples_to_use = %{read_sample_exclusion("sample_exclusion.txt")};
        my @sample_for_analysis_array;
        
        foreach my $sample_name (@starting_samples) {
            if ((defined $samples_to_use{$sample_name}) && $samples_to_use{$sample_name}) {
                push @sample_for_analysis_array, $sample_name;
            }
        }
        
        $config_hash->{SAMPLES_FOR_ANALYSIS} = join ",", @sample_for_analysis_array;
        create_analysis_config_file("config.txt", convert_hash_to_array($config_hash));
        
        %sample_for_analysis_hash = map {$_ => 1} @sample_for_analysis_array;
        foreach my $sample_name (@previous_analysis) {
            if ((! exists ($sample_for_analysis_hash{$sample_name})) && ($sample_counts{$sample_name} != 0)) {
                $need_to_reanalyse = 1;
                last;
            }
        }
    }
    
    # If we are running for the first time, or if sample_exclusions.txt has been modified,
    # reset the pipeline. 
    
    if ((! -e $global_qiime_mapping_filename) || ($need_to_reanalyse) || (! exists($config_hash->{PIPELINE_STAGE}))) {
        $config_hash->{PIPELINE_STAGE} = "PREAMBLE";
    } 
    
    if ($config_hash->{PIPELINE_STAGE} eq "PREAMBLE") {
        
    # (Re)create the QIIME mapping file, taking sample_exclusions.txt into account.
        
        my ($sample_array, $config_array) =
            parse_app_config_file("app.config");
        
        my %sample_to_analyse = map {$_ => 1} split /,/, $config_hash->{SAMPLES_FOR_ANALYSIS};
        
        my $new_sample_array;
        foreach my $sample_array_ptr (@{$sample_array}) {
            if (exists $sample_to_analyse{$sample_array_ptr->[0]}) {
                push @{$new_sample_array}, $sample_array_ptr;
            }
        }
        
        create_qiime_mapping_file($global_qiime_mapping_filename, $new_sample_array);
        
        # Begin the pipeline...
        
        my $params = [{-M => 1,
                    -q => 'raw_sequences.qual',
                    -f => 'raw_sequences.fasta',
                    -d => '',
                    -m => $global_qiime_mapping_filename,
                    -b => 'variable_length',
                    -a => 20,
                    -H => 40,
                    -k => '',
                    -s => 20,
                    -l => $config_hash->{TRIM_LENGTH},
                    -t => ''}];
        
        
        splitLibraries($params);
        
        if (-z 'seqs.fna') {
            croak "ERROR: No sequences in seqs.fna (no sequences successfully demultiplexed after split_libraries.py).\n" .
            "Check the config file that the barcode/primer sequences are correct.\n"
        }
        
        $config_hash->{PIPELINE_STAGE} = "TRIM_SEQS";
        
        create_analysis_config_file("config.txt", convert_hash_to_array($config_hash));
    }
    
    if ($config_hash->{PIPELINE_STAGE} eq "TRIM_SEQS") {
    
        my $params = [{-f => 'seqs.fna',
                    -q => 'seqs_filtered.qual',
                    -b => $config_hash->{TRIM_LENGTH}}];
        
        truncateFastaAndQual($params);
        
        if (-z 'seqs_filtered.fasta') {
          croak "ERROR: No sequences in seqs_filtered.fna (no sequences successfully trimmed after truncate_fasta_qual_files.py).\n" .
          "Check the config file that the barcode/primer sequences are correct.\n"
        }
        
        $config_hash->{PIPELINE_STAGE} = "ACACIA";
        
        create_analysis_config_file("config.txt", convert_hash_to_array($config_hash));
    
    }
    
    my $fasta_output_file = 'seqs_filtered.fasta';
    
    if ($config_hash->{PIPELINE_STAGE} eq "ACACIA") {
    
        my $fasta_output_file = 'seqs_filtered.fasta';
        
        if (! exists($pipeline_modifiers->{NO_ACACIA})) {
        
            $acacia_config_hash{'FASTA_LOCATION'} = 'seqs_filtered.fasta';
            
            $acacia_config_hash{TRIM_TO_LENGTH} = $config_hash->{TRIM_LENGTH};
            
            my $acacia_config_file = create_acacia_config_file(\%acacia_config_hash);
            
            my $params_array = [["-XX:+UseConcMarkSweepGC",
                                 "-Xmx100G"],
                                 {-jar => '$ACACIA'},
                                 {-c => $acacia_config_file}];
            
            mkdir('denoised_acacia');
            
            if (run_acacia($params_array)) {
                die "Acacia failed.\n";
            };
            
            unlink($acacia_config_file);
            
            $fasta_output_file = "denoised_acacia/acacia_out_all_tags.seqOut";
        }
        
        $config_hash->{PIPELINE_STAGE} = "CD_HIT_OTU";
        
        create_analysis_config_file("config.txt", convert_hash_to_array($config_hash));
    
    }
    
    my $cd_hit_otu_dir = dirname(`which cd-hit-otu-all.pl`);
    
    if ($config_hash->{PIPELINE_STAGE} eq "CD_HIT_OTU") {
    
        print "----------------------------------------------------------------\n";
        print "Start TABLE BASED NORMALISATION data set processing...\n";
        print "----------------------------------------------------------------\n";
        print "Copying reads for analysis...\n";
        
        print "Running CD-HIT-OTU\n";
        checkAndRunCommand("$cd_hit_otu_dir/cd-hit-otu-all.pl",
                           [{-i => $fasta_output_file,
                             -m => "false",
                             -o => "$processing_dir/cd_hit_otu"}], DIE_ON_FAILURE);
        
        my $rep_set_otu_array =
            reformat_CDHIT_repset("$processing_dir/cd_hit_otu/OTU",
                                  "$processing_dir/cd_hit_otu/OTU_numbered");
       
        $config_hash->{PIPELINE_STAGE} = "ASSIGN_TAXONOMY";
            
        create_analysis_config_file("config.txt", convert_hash_to_array($config_hash));
    
    }
    
    if ($config_hash->{PIPELINE_STAGE} eq "ASSIGN_TAXONOMY") {
    
        print "Assigning taxonomy for non normalised data set...\n";
        
        # update our databases (GG by default)
        my $TAX_tax_file = $QIIME_TAX_tax_file;
        my $TAX_blast_file = $QIIME_TAX_blast_file;
        my $TAX_aligned_blast_file = $QIIME_TAX_aligned_blast_file;
        my $imputed_file = $QIIME_imputed_file;
        
        #print "Assign taxonomy method: $assign_taxonomy_method\n";
        #if ($assign_taxonomy_method eq 'blast') {
            checkAndRunCommand("assign_taxonomy.py", [{-i => "$processing_dir/cd_hit_otu/OTU_numbered",
                                                        -t => $TAX_tax_file,
                                                        -b => $TAX_blast_file,
                                                        -m => "blast",
                                                        -a => 10,
                                                        -e => 0.001,
                                                        -o => $processing_dir}], DIE_ON_FAILURE);
        #} elsif ($assign_taxonomy_method eq 'bwasw') {
        #    checkAndRunCommand("assign_taxonomy.py", [{-i => "$processing_dir/cd_hit_otu/OTU_numbered",
        #                                                -t => $TAX_tax_file,
        #                                                -d => $TAX_blast_file,
        #                                                -m => "bwasw",
        #                                                -a => $num_threads,
        #                                                -o => $global_TB_processing_dir}], DIE_ON_FAILURE);
        #} else {
        #    die "Unrecognised assign_taxonomy method: '$assign_taxonomy_method'";
        #}
        
        $config_hash->{PIPELINE_STAGE} = "GENERATE_OTU_TABLE";
            
        create_analysis_config_file("config.txt", convert_hash_to_array($config_hash));
    
    }
    
    if ($config_hash->{PIPELINE_STAGE} eq "GENERATE_OTU_TABLE") {
    
        print "Assigning sample counts to clusters: \n";
        checkAndRunCommand("$cd_hit_otu_dir/clstr_sample_count_matrix.pl",
                           [["_", "$processing_dir/cd_hit_otu/OTU.nr2nd.clstr"]]);
        
        print "Making NON NORMALISED otu table...\n";
        
        create_QIIME_OTU_from_CDHIT("$processing_dir/OTU_numbered_tax_assignments.txt",
                                    "$processing_dir/cd_hit_otu/OTU.nr2nd.clstr.otu.txt",
                                    "$results_dir/non_normalised_otu_table.txt");
        
        print "Reformating OTU table...\n";
        checkAndRunCommand("reformat_otu_table.py",  [{-i => "$results_dir/non_normalised_otu_table.txt",
                                                       -t => "$processing_dir/OTU_numbered_tax_assignments.txt",
                                                       -o => "$results_dir/non_normalised_otu_table_expanded.tsv"}], IGNORE_FAILURE);
    
        my $otu_sample_counts = get_read_counts_from_cd_hit_otu("$processing_dir/cd_hit_otu/OTU.nr2nd.clstr.sample.txt");
            
        # If this is the first run through, create a file to record the sample_counts.
        if (! -e "sample_counts.txt") {
            open($fh, ">sample_counts.txt");
            print {$fh} "#NAME\tREAD COUNT\n";
            foreach my $sample_name (@starting_samples) {
                if (! defined $otu_sample_counts->{$sample_name}) {
                    $otu_sample_counts->{$sample_name} = 0; 
                }
                print {$fh} $sample_name. "\t" . $otu_sample_counts->{$sample_name} ."\n";
                $sample_counts{$sample_name} = $otu_sample_counts->{$sample_name};
            }
            close($fh);
        }
                
        my %previous_analysis_hash = map {$_ => 1} split /,/, $config_hash->{SAMPLES_FOR_ANALYSIS};
        
        
        open($fh, ">sample_exclusion.txt");
        print {$fh} "#NAME\tREAD COUNT\tUSE\n";
        foreach my $sample_name (@starting_samples) {
            my $use = 0;
            if ($previous_analysis_hash{$sample_name} && $sample_counts{$sample_name}) {
                $use = 1;
            }
            print {$fh} $sample_name. "\t" . $sample_counts{$sample_name} . "\t" . $use . "\n";
        }
        close($fh);
        
        print Dumper(\%sample_counts);
        
        print "##############################################\n";
        print "################ READ THIS ###################\n";
        print "##############################################\n";
        
        print "Non-normalized OTU table has now been generated.\n";
        print "Check the " . $options->{'d'} . "/sample_exclusion.txt and choose which samples you wish to exclude.\n";
        print "      vim " . $options->{'d'} . "/sample_exclusion.txt" . "\n\n";
        print "When you have chosen which samples to exclude, run the following command:\n";
        print "      app_run_analysis.pl -d " . $options->{'d'} . "\n\n";
        
        $config_hash->{PIPELINE_STAGE} = "NORMALISATION";
            
        create_analysis_config_file("config.txt", convert_hash_to_array($config_hash));
        
        exit();
    
    }
    
    if ($config_hash->{PIPELINE_STAGE} eq "NORMALISATION") {
                  
        my $min_samples = min(map {$sample_counts{$_}} keys %sample_for_analysis_hash);
        
        my $global_norm_sample_size = (int($min_samples / 50) - 1) * 50;
        
        if ($global_norm_sample_size <= 0) {
            print "Unable to normalise, some samples have too few reads. Review " . $options->{'d'} . "/sample_exclusion.txt " .
                  "and rerun app_run_analysis.pl -d " . $options->{'d'} . "\n\n";
        }
        
        print "Correcting normalisation depth to $global_norm_sample_size sequences...\n";
        
        my $global_norm_num_reps = 1000;
        
        print "Normalizing non normalised table at $global_norm_sample_size sequences... [$global_norm_sample_size, $global_norm_num_reps]\n";
        
        my $params = [{-i => "$results_dir/non_normalised_otu_table.txt",
                    -o => "$processing_dir/rare_tables/",
                    -d => $global_norm_sample_size,
                    -n => $global_norm_num_reps},
                    ["--lineages_included"],
                    ["--k"]];
        
        my $normalised_OTU_file = 
            normalise_OTU_table($params, $processing_dir, $global_norm_sample_size,
                                scalar keys %sample_for_analysis_hash);
            
        `cp $normalised_OTU_file $results_dir/normalised_otu_table.txt`;
        
        checkAndRunCommand("reformat_otu_table.py",  [{-i => "$results_dir/normalised_otu_table.txt",
                                                       -t => "$processing_dir/OTU_numbered_tax_assignments.txt",
                                                       -o => "$results_dir/normalised_otu_table_expanded.tsv"}], IGNORE_FAILURE);
        
        print "Summarizing by taxa.....\n";

        checkAndRunCommand("summarize_taxa.py", [{-i => "$results_dir/normalised_otu_table.txt",
                                                  -o => "$results_dir/breakdown_by_taxonomy/"}], DIE_ON_FAILURE); 
         
        $config_hash->{PIPELINE_STAGE} = "RAREFACTION";
            
        create_analysis_config_file("config.txt", convert_hash_to_array($config_hash));
    }
    
    if (($config_hash->{OTU_TABLES_ONLY}) && (! exists($pipeline_modifiers->{OTU_TABLES_ONLY}))) {
    
        print "APP configured to only create OTU tables. Stopping here.";
        
        return;
    
    }
    
    if ($config_hash->{PIPELINE_STAGE} eq "RAREFACTION") {
        
        my $max_samples = max(map {$sample_counts{$_}} keys %sample_for_analysis_hash);
        
        my $global_rare_X = int($max_samples / 50) * 50;
        my $global_rare_N = 50;
        
        # Do rarefaction in stages (10-100 in 10s), (100-1000 in 50s), 1000-5000 in 100s, 5000-10000 in 500s.
        
        print "Doing rarefaction stages....\n";
        
        my @rarefaction_stages = ({min => 10, max => 100, step => 10},
                                  {min => 100, max => 1000, step => 50},
                                  {min => 1000, max => 10000, step => 100},
                                  {min => 10000, max => 50000, step => 500},
                                  {min => 50000, max => 100000, step => 1000});
    
        foreach my $stage (@rarefaction_stages) {
            if (($global_rare_X < $stage->{max})) {
                checkAndRunCommand("multiple_rarefactions.py", [{-i => "$results_dir/non_normalised_otu_table.txt",
                                                                 -o => "$processing_dir/rarefied_otu_tables/",
                                                                 -m => $stage->{min},
                                                                 -x => $global_rare_X,
                                                                 -s => $stage->{step},
                                                                 -n => $global_rare_N}], DIE_ON_FAILURE);
                last;
            } else {
                checkAndRunCommand("multiple_rarefactions.py", [{-i => "$results_dir/non_normalised_otu_table.txt",
                                                                 -o => "$processing_dir/rarefied_otu_tables/",
                                                                 -m => $stage->{min},
                                                                 -x => $stage->{max} - $stage->{step},
                                                                 -s => $stage->{step},
                                                                 -n => $global_rare_N}], DIE_ON_FAILURE);
            }
        }
        
        $config_hash->{PIPELINE_STAGE} = "ALPHA_DIVERSITY";
            
        create_analysis_config_file("config.txt", convert_hash_to_array($config_hash));
        
    }
    
    if ($config_hash->{PIPELINE_STAGE} eq "ALPHA_DIVERSITY") {
        
        print "Calculating (non-phylogeny dependent) alpha diversity metrics....\n";
        
        my $methods_str = join(",", qw(chao1
                                       chao1_confidence
                                       observed_species
                                       simpson
                                       shannon
                                       fisher_alpha));
        
        checkAndRunCommand("alpha_diversity.py", [{-i => "$processing_dir/rarefied_otu_tables/",
                                                   -o => "$processing_dir/alpha_div/",
                                                   -m => $methods_str}], DIE_ON_FAILURE);
        
        $config_hash->{PIPELINE_STAGE} = "ALPHA_DIVERSITY_COLLATION";
            
        create_analysis_config_file("config.txt", convert_hash_to_array($config_hash));
        
    }
    
    if ($config_hash->{PIPELINE_STAGE} eq "ALPHA_DIVERSITY_COLLATION") {
        
        checkAndRunCommand("collate_alpha.py", [{-i => "$processing_dir/alpha_div/",
                                                 -o => "$processing_dir/alpha_div_collated/"}], DIE_ON_FAILURE);
     
        foreach my $format (("png", "svg")) {
            checkAndRunCommand("make_rarefaction_plots.py", [{-i => "$processing_dir/alpha_div_collated/",
                                                              -m => "qiime_mapping.txt",
                                                              -o => "$results_dir/alpha_diversity/",
                                                              "--resolution" => 300,
                                                              "--imagetype" => $format}], DIE_ON_FAILURE);
        }
        
        $config_hash->{PIPELINE_STAGE} = "TREE_CREATION";
            
        create_analysis_config_file("config.txt", convert_hash_to_array($config_hash));
        
    }
    
    if ($config_hash->{PIPELINE_STAGE} eq "TREE_CREATION") {
    
        my $imputed_file = $QIIME_imputed_file;
        
        print "Treeing non normalised data set...\n";
        checkAndRunCommand("align_seqs.py", [{-i => "$processing_dir/cd_hit_otu/OTU_numbered",
                                              -t => $imputed_file,
                                              -p => 0.6,
                                              -o => "$processing_dir/pynast_aligned"}], DIE_ON_FAILURE);
                
        checkAndRunCommand("filter_alignment.py", [{-i => "$processing_dir/pynast_aligned/OTU_numbered_aligned",
                                                    -o => "$processing_dir"}], DIE_ON_FAILURE);
    }    
}

######################################################################
# Pipeline modules
######################################################################

sub splitLibraries
{
    #-----
    # Wrapper for Qiime split libraries
    #
    my ($params) = @_;
    print "Splitting libraries...\n";
    checkAndRunCommand("split_libraries.py", $params, DIE_ON_FAILURE);
}

sub truncateFastaAndQual
{
    my ($params) = @_;
    print "Trimming reads...\n";
    checkAndRunCommand("truncate_fasta_qual_files.py", $params, DIE_ON_FAILURE);
}

sub uclustRemoveChimeras
{
    #-----
    # Remove chimeras using uclust
    #
    my ($params) = @_;

    
    print "Removing chimeras...\n";
        
    checkAndRunCommand("usearch", $params, DIE_ON_FAILURE);
}

sub run_acacia
{
    #-----
    # run acacia on the data
    #
    my ($params) = @_;
       
    print "Denoising using acacia...\n";
    
    # Because acacia's return value is 1 on success (instead of the traditional
    # zero), we need to ignore on failure and test if the stats file was written
    # to.

    
    checkAndRunCommand("java", $params, IGNORE_FAILURE);
    
    
    # TODO: This is lazy. Shouldn't need to call out to global hash.
    my $stats_file = $acacia_config_hash{OUTPUT_DIR} . "/" .
        $acacia_config_hash{OUTPUT_PREFIX} . "_all_tags.stats";
    
    # If the stats file doesn't exist or is zero size, die and raise an error.
    if ((! -e $stats_file) || (-z $stats_file)) {
        print "\n################# WARNING!!!!!!!!!! #################\n" .
              "The ACACIA stats file was not written to!!!\n" .
              "You should check the acacia_standard_error.txt and\n" .
              "acacia_standard_debug.txt in the QA/denoised_acacia/\n" .
              "directory before proceeding to app_make_results.pl to\n" .
              "ensure ACACIA completed successfully.\n" .
              "#####################################################\n\n";   
        return 1;
    }
    #`java -XX:+UseConcMarkSweepGC -Xmx10G -jar \$ACACIA -c $config_filename`;
    #`sed -i -e "s/all_tags_[^ ]* //" $global_acacia_output_dir/$ACACIA_out_file`;
    return 0;
    
}


sub normalise_OTU_table {
    my ($params, $processing_dir, $global_norm_sample_size, $num_samples) = @_;
    
    checkAndRunCommand("multiple_rarefactions_even_depth.py", $params, DIE_ON_FAILURE);
        
    print "Calculating centroid subsampled table...\n";
    my $centroid_index = find_centroid_table("$processing_dir/rare_tables/",
                                             $global_norm_sample_size,
                                             $num_samples,
                                             "$processing_dir/$tn_dist_file",
                                             "$processing_dir/$tn_log_file",
                                             "$processing_dir/$tn_R_log_file");
    
    return "$processing_dir/rare_tables/rarefaction_$global_norm_sample_size"."_$centroid_index".".txt";

}

###############################################################################
# Subs
###############################################################################

sub get_pipeline_modifiers {
    my %pipeline_modifiers;
    my ($pipeline_string) = @_;
    my @splitline = split /,/, $pipeline_string;
    for (my $i = 1; $i < scalar @splitline; $i++) {
        if ($splitline[$i] =~ /^(\s*)(\S*)(\s*)$/) {
            $pipeline_modifiers{$2}++;
        }
    }
    return \%pipeline_modifiers;
}

sub updateAcaciaConfigHash {
    my ($config_file) = @_;
    open(my $fh, "<", $config_file);
    while (my $line = <$fh>) {
        if ($line =~ /^(.*)=(.*)$/) {
            $acacia_config_hash{$1} = $2;
        }
    }
    close($fh);
}

sub create_acacia_config_file {
    my ($acacia_config_hash) = @_;
    
    my $acacia_config_string =
        join("\n", map({$_ ."=" .$acacia_config_hash->{$_}} sort {$a cmp $b} keys %{$acacia_config_hash}));
    
    # Create a temporary file to dump the config
    my ($tmp_fh, $config_filename) = tempfile("acacia_XXXXXXXX", UNLINK => 0);
    
    # Unbuffer the output of $tmp_fh
    select($tmp_fh);
    $| = 1;
    select(STDOUT);
    
    # Write to the temporary acacia config file
    print {$tmp_fh} $acacia_config_string;
    
    close($tmp_fh);
    
    return $config_filename;
}

sub reformat_CDHIT_repset {
    my ($fasta_file, $output_file) = @_;
    open(my $in_fh, $fasta_file);
    open(my $out_fh, ">$output_file");
    my $return_array;
    my $count = 0;
    while(my $line = <$in_fh>) {
        if ($line =~ /^>(([^\s]*).*)/) {
            print {$out_fh} ">$count $1\n";
            push @{$return_array}, $2;
            $count++;
            
        } else {
            print {$out_fh} $line;
        }
    }
    close($out_fh);
    close($in_fh);
    return $return_array;
}

sub create_QIIME_OTU_from_CDHIT {
    my ($nn_rep_set_tax_assign, $cd_hit_otu_file, $outfile) = @_;
    
    # Read the CD_HIT_OTU OTU table.
    open(my $cd_hit_otu_fh, $cd_hit_otu_file);
    my $header_line = <$cd_hit_otu_fh>;
    my @splitline = split /\t/, $header_line;
    my @headers = @splitline[1..($#splitline-1)];
    
    my @cd_hit_otu_array;
    my $max_sample_count = 0;
    while (my $line = <$cd_hit_otu_fh>) {
        my @splitline = split /\t/, $line;
        #Push everything but the first and last columns;
        my @slice = @splitline[1..($#splitline-1)];
        my $sample_count = scalar @slice;
        $max_sample_count = max($max_sample_count, $sample_count);
        push (@cd_hit_otu_array, \@slice);
    }
    close($cd_hit_otu_fh);
    
    # Create Taxonomy correlations
    my %tax_hash;
    open(my $tax_fh, $nn_rep_set_tax_assign);
    while (my $line = <$tax_fh>) {
        my @splitline = split /\t/, $line;
        $tax_hash{$splitline[0]} = $splitline[1];
    }
    close($tax_fh);
    
    my @sample_counts = (0) x $max_sample_count;
    open(my $out_fh, ">$outfile");
    # Print OTU headers   
    print {$out_fh} "# QIIME v1.3.0 OTU table\n";
    print {$out_fh} join("\t", ("#OTU ID",
                      join("\t", @headers),
                      "Consensus Lineage")), "\n";
    # Print OTU table guts
    for(my $i = 0; $i < scalar @cd_hit_otu_array; $i++) {
        print {$out_fh} join("\t", ($i,
                                    join("\t", @{$cd_hit_otu_array[$i]}),
                                    $tax_hash{$i})), "\n";
        for(my $j = 0; $j < scalar @{$cd_hit_otu_array[$i]}; $j++) {
            $sample_counts[$j] += $cd_hit_otu_array[$i]->[$j];
        }
    }
    close($out_fh);
    return \@sample_counts;
    
}

sub run_R_cmd
{
    #-----
    # Wrapper for running R commands
    #
    my ($R_instance, $cmd, $log_fh) = @_;
    if ($log_fh) {
        print {$log_fh} $cmd, "\n";
    }
    $R_instance->run($cmd);
}

sub find_centroid_table
{
    #-----
    # Find a representative set of $global_norm_sample_size sequences (do this in RRRRRR...)
    #
    my($path, $norm_size, $num_samples, $dist_file, $log_file, $R_cmd_log) = @_;

    my $R_instance = Statistics::R->new();
    
    open (my $R_cmd_log_fh, ">", $R_cmd_log)
        or warn("Unable to open $R_cmd_log for logging R commands. Proceeding without logging.");

    print "Calculating centroid OTU table from tables in $path...\n";

    $R_instance->start();

    my $sampl_p1 = $num_samples + 1;

    # read in the list of distance matricies
    run_R_cmd($R_instance, qq`library(foreign);`, $R_cmd_log_fh);
    run_R_cmd($R_instance, qq`a<-list.files("$path", "rarefaction_$norm_size` . '_[0-9]*.txt");', $R_cmd_log_fh);

    # work out how many there are and allocate an array
    run_R_cmd($R_instance, qq`len_a <- length(a);`, $R_cmd_log_fh);
    run_R_cmd($R_instance, qq`big_frame <- array(0,dim=c($num_samples,$num_samples,len_a));`, $R_cmd_log_fh);

    print "  --start loading data...\n";

    # load each file individually into a big frame
    my $r_str = "for (i in c(1:len_a)) { j <- i - 1; name <- paste(\"$path\",\"rarefaction_$norm_size\",\"_\",j,\".txt\",sep=\"\"); u<-read.table(name,sep=\"\\t\",row.names=1); u[,$sampl_p1]<-NULL; big_frame[,,i]<-as.matrix(dist(t(u), upper=TRUE, diag=TRUE)); i<-i+1; }";
    run_R_cmd($R_instance, $r_str, $R_cmd_log_fh);

    print "  --data loaded, calculating centroid...\n";

    # find the average matrix
    run_R_cmd($R_instance, qq`ave <- big_frame[,,1];`, $R_cmd_log_fh);
    run_R_cmd($R_instance, qq`for (i in c(2:len_a)) { ave <- ave + big_frame[,,i]; }`, $R_cmd_log_fh);
    run_R_cmd($R_instance, qq`ave <- ave/len_a;`, $R_cmd_log_fh);

    print "  --calculating distances of tables to centroid...\n";

    # find the euclidean distance of each matrix from the average
    run_R_cmd($R_instance, qq`dist<-array(0,dim=c(len_a));`, $R_cmd_log_fh);
    run_R_cmd($R_instance, qq`for (i in c(1:len_a)) { dist[i] <- sqrt(sum(big_frame[,,i]-ave)^2); }`, $R_cmd_log_fh);

    # find the min value
    run_R_cmd($R_instance, qq`min_index <- which.min(dist);`, $R_cmd_log_fh);
    my $centroid_otu_index = $R_instance->get('min_index');
    # R indexes from 0
    $centroid_otu_index--;
    print "  --table: $centroid_otu_index is the centroid table\n";

    # make stats on the distances
    # and log what we did
    open my $log_fh, ">", $log_file or die "Could not open log file: $log_file : $!\n";
    run_R_cmd($R_instance, qq`max_dist <- max(dist);`, $R_cmd_log_fh);
    run_R_cmd($R_instance, qq`min_dist <- min(dist);`, $R_cmd_log_fh);
    run_R_cmd($R_instance, qq`range_dist <- max_dist - min_dist;`, $R_cmd_log_fh);
    run_R_cmd($R_instance, qq`mean_dist <- mean(dist);`, $R_cmd_log_fh);
    run_R_cmd($R_instance, qq`median_dist <- median(dist);`, $R_cmd_log_fh);

    print $log_fh "---------------------------------------------------\n";
    print $log_fh "  Centroid OTU table based normalised statistics\n";
    print $log_fh "---------------------------------------------------\n";
    print $log_fh "Max dist:\t".$R_instance->get('max_dist')."\n";
    print $log_fh "Min dist:\t".$R_instance->get('min_dist')."\n";
    print $log_fh "Range:\t".$R_instance->get('range_dist')."\n";
    print $log_fh "Mean:\t".$R_instance->get('mean_dist')."\n";
    print $log_fh "Median:\t".$R_instance->get('median_dist')."\n";

    if(2 < $num_samples)
    {
        run_R_cmd($R_instance, qq`library(permute);`, $R_cmd_log_fh);
        run_R_cmd($R_instance, qq`library(vegan);`, $R_cmd_log_fh);
        run_R_cmd($R_instance, qq`mantel.otu <- mantel(ave,big_frame[,,min_index]);`, $R_cmd_log_fh);
        run_R_cmd($R_instance, qq`m_stat <- mantel.otu\$statistic;`, $R_cmd_log_fh);
        run_R_cmd($R_instance, qq`m_sig <- mantel.otu\$signif;`, $R_cmd_log_fh);
        print $log_fh "Mantel P stat:\t".$R_instance->get('m_sig')."\n";
        print $log_fh "Mantel R stat:\t".$R_instance->get('m_stat')."\n";
    }
    else
    {
        print "Too few samples to perform a mantel test.\n";
    }
    close $log_fh;

    # print all the distances to a file so we can make purdy pictures from them later
    open my $dist_fh, ">", $dist_file or die "Could not open distance file: $dist_file : $!\n";
    my $num_tables = $R_instance->get('len_a');
    foreach my $counter (1..$num_tables)
    {
        print $dist_fh $R_instance->get("dist[$counter]")."\n"
    }
    close $dist_fh;

    close($R_cmd_log_fh);
    # let the user know the result    
    return $centroid_otu_index;
}

sub checkParams {
    my @standard_options = ( "help|h+", "d:s");
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
    if(!exists $options{'d'} ) { print "**ERROR: you MUST give a directory\n"; exec("pod2usage $0"); }

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

    app_run_analysis.pl

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

   Runs an APP analysis in a folder created by app_create_analysis.
   
=head1 SYNOPSIS

    app_run_analysis.pl -d directory [-help|h]

      -d directory                 Directory created by app_create_analysis.pl to run.
      [-help -h]                   Displays basic usage information
         
=cut

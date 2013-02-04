package AppCommon;
require Exporter;
use File::Basename;
use File::Spec;
use File::Temp qw(tempfile);
use List::Util qw(min max);
use Carp;

use AppConfig;

our @ISA = qw(Exporter);
our @EXPORT=qw(
    %known_config_options
    create_unique_analysis_dir
    create_data_file_links
    parse_app_config_files
    read_sample_exclusion
    create_analysis_config_file
    create_qiime_mapping_file
    get_parameter_modifiers
    convert_hash_to_array
    get_read_counts_from_sample_counts
    get_read_counts_from_cd_hit_otu
    get_read_counts_from_OTU_table
    normalise_otu_table
    setup_db_paths
    jack_knifing
);

our %known_config_options = (
    NORMALISE => ['SEQ'],
    MUL_RARE_M => [],
    MUL_RARE_X => [],
    MUL_RARE_S => [],
    MUL_RARE_N => [],
    PIPELINE => [],
    OTU_TABLES_ONLY => [],
    DB => ['SILVA','GG','MERGED']
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
    my ($output_dir, $job_name_array) = @_;
    
    mkdir "$output_dir/raw_files";
    
    foreach my $job (@{$job_name_array}) {
        my $job_dir = $job->[0];
        my $job_name = $job->[1];
        
        if (! $job_dir) {
            $job_dir = '.';
        }
    
        my @fasta_files;
        foreach my $suffix ((".fa", ".fna", ".fasta")) {
            if (-e "$job_dir/$job_name$suffix") {
                push @fasta_files, "$job_dir/$job_name$suffix";
            }
        }
        if (scalar @fasta_files) {
            if (scalar @fasta_files > 1) {
                croak "ERROR: Too many possible FASTA files to choose, remove or rename " .
                "to stop ambiguity. Offending files: " . join(", ", @fasta_files) . "\n";
            }
        } else {
            croak "ERROR: Unable to find fasta file: $job_dir/$job_name.fna\n";
        }
        
        if (! -e "$job_dir/$job_name.qual") {
             croak "ERROR: Unable to find qual file: $job_dir/$job_name.qual\n";
        }
        symlink File::Spec->rel2abs("$job_dir/$job_name.qual"), "$output_dir/raw_files/$job_name.qual"
            or die "Unable to symlink $job_dir/$job_name.qual\n";
        symlink File::Spec->rel2abs("$fasta_files[0]"), "$output_dir/raw_files/$job_name.fasta"
            or die "Unable to symlink $job_dir/$fasta_files[0]\n";
    }
}

################################################################################
# Subroutine: parse_app_config_files($config_files_str)
# 
################################################################################

sub parse_app_config_files {
    my $config_files_str = shift;
    my @config_files = split /,/, $config_files_str;
    
    my $config_array = [];
    my $sample_array = [];
    
    my $config_file_counter = 0;
    
    foreach my $config_file (@config_files) {
        open(my $fh, $config_file);
               
        $config_file_counter++;
        
        my $in_config_section = 0;
        
        while (my $line = <$fh>) {
            chomp $line;
            # Skip the line if starts with a hash(#)
            if ($line =~ /^(\s)*#/) {next};
            
            # Are we in the config section yet?
            if ($line =~ /@@/) {
                $in_config_section = 1;
                # If we are in the config section, and this is the not the first config
                # file, then go to the next config file.
                if ($config_file_counter > 1) {
                    last;
                }
                next;
            }
            
            # If we are in the config section...
            if ($in_config_section) {
                my @splitline = split /=/, $line;
                my @split_without_modifiers = split /,/, $splitline[1];
                my $without_modifers = $split_without_modifiers[0];
                
                # Check if the option is known/allowed.
                if (exists($known_config_options{$splitline[0]})) {
                    
                    my @possible_choices = @{$known_config_options{$splitline[0]}};
                    my %options = map { $_ => 1 } @possible_choices;
                    
                    # Check if the choices for this option are limited.
                    if ((! scalar @possible_choices) || exists ($options{$without_modifers})) {
                        
                        # Check if the choice is defined.
                        if (! defined ($without_modifers) || $splitline[1] =~ /^\s*$/) {
                            print "No config setting for option " .
                            $without_modifers . " in $config_file: - ignoring option.\n";
                        } else {
                            push @{$config_array}, [$splitline[0], $splitline[1]];
                        }
                    } else {
                        print "Unknown config setting for option " .
                            $without_modifers . " in $config_file: " .
                            $without_modifers . " - ignoring option.\n";
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
    }
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


sub get_parameter_modifiers {
    my %parameter_modifiers;
    my ($parameter_string) = @_;
    my @splitline = split /,/, $parameter_string;
    for (my $i = 1; $i < scalar @splitline; $i++) {
        if ($splitline[$i] =~ /^(\s*)(\S*)(\s*)$/) {
            $parameter_modifiers{$2}++;
        }
    }
    return \%parameter_modifiers;
}


sub setup_db_paths {
    my ($extra_param_hash, $db_string) = @_;
    if (! defined($db_string)) {
        $db_string = "GG";
    }
    my $db_modifiers = get_parameter_modifiers($db_string);
    if ($db_string =~ /^GG/) {
        if (exists($db_modifiers->{'11_02'})) {
            $extra_param_hash->{DB}->{OTUS} = "/srv/whitlam/bio/db/gg/qiime_default/gg_otus_4feb2011/rep_set/gg_99_otus_4feb2011.fasta";
            $extra_param_hash->{DB}->{TAXONOMIES} = "/srv/whitlam/bio/db/gg/qiime_default/gg_otus_4feb2011/taxonomies/greengenes_tax.txt";
            $extra_param_hash->{DB}->{IMPUTED} ="/srv/whitlam/bio/db/gg/qiime_default/core_set_aligned.fasta.imputed";
            $extra_param_hash->{DB}->{DESCRIPTION} = "Using Feb 2011 Greengenes OTUS/Taxonomy";
        } else {
            $extra_param_hash->{DB}->{OTUS} = "/srv/whitlam/bio/db/gg/from_www.secongenome.com/2012_10/gg_12_10_otus/rep_set/99_otus.fasta";
            $extra_param_hash->{DB}->{TAXONOMIES} = "/srv/whitlam/bio/db/gg/from_www.secongenome.com/2012_10/gg_12_10_otus/taxonomy/99_otu_taxonomy.txt";
            $extra_param_hash->{DB}->{IMPUTED} ="/srv/whitlam/bio/db/gg/from_www.secongenome.com/2012_10/gg_12_10_otu_99_aligned.fasta";
            $extra_param_hash->{DB}->{DESCRIPTION} = "Using Oct 2012 Greengenes OTUS/Taxonomy";
        }
    } elsif ($db_string =~ /^SILVA/) {
        $extra_param_hash->{DB}->{OTUS} = "/srv/whitlam/bio/db/Silva/QIIME_files_r108/rep_set/Silva_108_rep_set.fna";
        $extra_param_hash->{DB}->{TAXONOMIES} = "/srv/whitlam/bio/db/Silva/QIIME_files_r108/taxa_mapping/Silva_108_taxa_mapping.txt";
        $extra_param_hash->{DB}->{IMPUTED} = "/srv/whitlam/bio/db/Silva/QIIME_files_r108/core_aligned/Silva_108_core_aligned_seqs.fasta";
        $extra_param_hash->{DB}->{DESCRIPTION} = "Using Silva r108 database OTUS/Taxonomy";
    } elsif ($db_string =~ /^MERGED/) {
        $extra_param_hash->{DB}->{OTUS} = "/srv/whitlam/bio/db/merged_gg_silva/201210/merged_gg_silva_seqs.fasta";
        $extra_param_hash->{DB}->{TAXONOMIES} = "/srv/whitlam/bio/db/merged_gg_silva/201210/merged_gg_silva_taxo.txt";
        $extra_param_hash->{DB}->{DESCRIPTION} = "Using Merged Silva r108 and Oct 2012 greengenes database OTUS/Taxonomy";
    } else {
        die "Taxonomic assignment database unknown - $db_string";
    }
    use Data::Dumper;
}

sub read_sample_exclusion {
    my ($sample_exclusion_file) = @_;
    open(my $fh, $sample_exclusion_file);
    <$fh>; #burn headers
    my %sample_use;
    while (my $line = <$fh>) {
        chomp $line;
        my @splitline = split /\t/, $line;
        $sample_use{$splitline[0]} = $splitline[$#splitline ];
    }
    close($fh);
    return \%sample_use;
}

sub get_read_counts_from_sample_counts {
    my ($sample_counts_file) = @_;
    open(my $fh, $sample_counts_file);
    <$fh>; #burn headers
    my %sample_counts;
    while (my $line = <$fh>) {
        chomp $line;
        my @splitline = split /\t/, $line;
        $sample_counts{$splitline[0]} = $splitline[1];
    }
    close($fh);
    return \%sample_counts;
}


sub get_read_counts_from_cd_hit_otu {
    my ($cd_hit_otu_file) = @_;
    open(my $fh, $cd_hit_otu_file);
    <$fh>; #burn headers
    my %sample_counts;
    while (my $line = <$fh>) {
        chomp $line;
        my @splitline = split /\t/, $line;
        $sample_counts{$splitline[0]} = $splitline[$#splitline - 1];
    }
    close($fh);
    return \%sample_counts;
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
    my %sample_counts;
    for (my $i = 0; $i < scalar @headers; $i++) {
        $sample_counts{$headers[$i]} = $counts[$i];
    }
    return \%sample_counts;
}

sub convert_hash_to_array {
    my ($config_hash) = @_;
    my @config_array = map {[$_, $config_hash->{$_}]} keys %{$config_hash};
    return \@config_array;
}

sub normalise_otu_table {
        my $options = shift;
        my $sample_for_analysis_hash_ref = shift;
        my $sample_counts_ref = shift;
        my $results_dir = shift;
        my $processing_dir = shift;
        
        
        my %sample_for_analysis_hash = %{$sample_for_analysis_hash_ref};
        my %sample_counts = %{$sample_counts_ref};
        
        
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
        
        checkAndRunCommand("multiple_rarefactions_even_depth.py", $params, DIE_ON_FAILURE);
        
        print "Calculating centroid subsampled table...\n";
        my $centroid_index = find_centroid_table("$processing_dir/rare_tables/",
                                             $global_norm_sample_size,
                                             scalar keys %sample_for_analysis_hash,
                                             "$processing_dir/$tn_dist_file",
                                             "$processing_dir/$tn_log_file",
                                             "$processing_dir/$tn_R_log_file");
        
        my $normalised_OTU_file = "$processing_dir/rare_tables/rarefaction_$global_norm_sample_size"."_$centroid_index".".txt";
            
        `cp $normalised_OTU_file $results_dir/normalised_otu_table.txt`;
        
        return "$results_dir/normalised_otu_table.txt";
        
}

sub jack_knifing
{
    #-----
    # Do jackknifing for a specific type of matrix
    #
    my ($matrix_type, $raw_otu_table, $raw_tree, $jack_knife_folder, $output_folder) = @_;

    # Produce distance matrix reflecting beta diversity in non-normalised OTU table
    checkAndRunCommand("beta_diversity.py", [{-i => $raw_otu_table,
                                              -o => "$output_folder/$matrix_type",
                                              -m => $matrix_type,
                                              -t => $raw_tree}], DIE_ON_FAILURE);

    my $basename = basename($raw_otu_table);
    $basename =~ s/\.txt//;
    my $beta_otu_table = $matrix_type."_$basename.txt";
    my $upgma_cluster_tree = $matrix_type."_" . $basename . "_upgma.tre";

    checkAndRunCommand("upgma_cluster.py", [{-i => "$output_folder/$matrix_type/$beta_otu_table",
                                             -o => "$output_folder/$matrix_type/$upgma_cluster_tree"}], DIE_ON_FAILURE);

    # Produce distance matrices reflecting beta diversity in the rarefied OTU tables
    checkAndRunCommand("beta_diversity.py", [{-i => $jack_knife_folder,
                                              -o => "$output_folder/$matrix_type/rare_dm/",
                                              -m => $matrix_type,
                                              -t => $raw_tree}], DIE_ON_FAILURE);

    # Perform UPGMA clustering on rarefied distance matrices
    checkAndRunCommand("upgma_cluster.py", [{-i => "$output_folder/$matrix_type/rare_dm/",
                                             -o => "$output_folder/$matrix_type/rare_upgma/"}], DIE_ON_FAILURE);

    # Compare the consensus tree to the beta-derived trees
    checkAndRunCommand("tree_compare.py", [{-s => "$output_folder/$matrix_type/rare_upgma/",
                                            -m => "$output_folder/$matrix_type/$upgma_cluster_tree",
                                            -o => "$output_folder/$matrix_type/upgma_cmp/"}], DIE_ON_FAILURE);

    # Compute principal coordinates
    checkAndRunCommand("principal_coordinates.py", [{-i => "$output_folder/$matrix_type/rare_dm/",
                                                     -o => "$output_folder/$matrix_type/pcoa/"}], DIE_ON_FAILURE);

    # Make PDF of Jackknife tree with labeled support: weighted unifrac command
    my $output_pdf = "$output_folder/$matrix_type"."_betadiv_jackknife_tree.pdf";
    checkAndRunCommand("make_bootstrapped_tree.py", [{-m => "$output_folder/$matrix_type/upgma_cmp/master_tree.tre",
                                                      -s => "$output_folder/$matrix_type/upgma_cmp/jackknife_support.txt",
                                                      -o => $output_pdf}], DIE_ON_FAILURE);

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

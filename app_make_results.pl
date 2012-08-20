#!/usr/bin/perl
###############################################################################
#
#    app_make_results.pl
#    
#    Normalise and complete the QIIME pieline
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
use File::Copy;
use Statistics::R;
use Data::Dumper;

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

#### Get the conf base
my $global_conf_base = basename($options->{'config'});
my @cb_1 = split /_/, $global_conf_base;
my @cb_2 = split /\./, $cb_1[1];
$global_conf_base = $cb_2[0];


### Create a unique output dir for this analysis.
my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst)
    = localtime(time);
    
my $output_dir = sprintf("app_analysis_%04i%02i%02i", $year+1900, $mon+1, $mday);

if (-e $output_dir) {
    my $folder_suffix = 1;
    while (-e $output_dir . "_" . $folder_suffix) {
        $folder_suffix++;
    }
    $output_dir .= "_" . $folder_suffix;
    mkdir($output_dir) or die "Unable to create $output_dir: $!";
}

#### GET THE WORKING DIR
getWorkingDirs(dirname($options->{'config'}), $output_dir);

### Delete any previous analyses
if (-d $global_TB_processing_dir) {
    `rm -rf $global_TB_processing_dir`;
}
if (-d $global_SB_processing_dir) {
    `rm -rf $global_SB_processing_dir`;
}

#### make these dirs
makeResultsDirs(0);

# sample IDs
my %global_samp_ID_list = ();
my $global_num_samples = 0;

# we can compare sequences to the greengenes or the SILVA dbs
my $global_comp_DB_type = "GG";

my $assign_taxonomy_method = 'blast';

# there are a number of different ways to normalise
# by default don't normalise
my $global_norm_style = "TABLE";
my $global_norm_sample_size = 0;

# how many times are we going to resample sequences (by default)
# in second round normalisation
my $global_norm_seq_norm_rounds = 7;

# defaults for rarefication
my $global_rare_M = 50;
# This is made dynamic in parse below
my $global_rare_X = -1;
my $global_rare_S = 50;
my $global_rare_N = 50;

my $global_min_sample_size = 100000000000;
my $global_max_sample_size = -1;

# defaults for otu similarity
my $global_similarity_setting = 0.97;

#defaults for assign taxonomy
my $global_e_value = 0.001;

# number of threads to us (currently this only affects assign taxonomy)
my $num_threads = 5;

# how many reps to do when normalising the OTU table
my $global_norm_num_reps = 1000;

# we only need a subset of these guys to do jacknifing
my $global_JN_file_count = 100;
if($global_JN_file_count > $global_norm_num_reps) { $global_JN_file_count = $global_norm_num_reps - 1; }

#### Override defaults from config or user
if(exists $options->{'identity'}) { $global_similarity_setting = $options->{'identity'}; }
if(exists $options->{'e'}) { $global_e_value = $options->{'e'}; }
if(exists $options->{'a'}) { $assign_taxonomy_method = $options->{'a'}; }
if(exists $options->{'threads'}) { $num_threads = $options->{'threads'}; }

print "Checking if all the config checks out...\t\t";
parse_config_results();

# update our databases (GG by default)
my $TAX_tax_file = $QIIME_TAX_tax_file;
my $TAX_blast_file = $QIIME_TAX_blast_file;
my $TAX_aligned_blast_file = $QIIME_TAX_aligned_blast_file;
my $imputed_file = $QIIME_imputed_file;

if($global_comp_DB_type eq "SILVA")
{
    $TAX_tax_file = $SILVA_TAX_tax_file;
    $TAX_blast_file = $SILVA_TAX_blast_file;
    $TAX_aligned_blast_file = $SILVA_TAX_aligned_blast_file;
    $imputed_file = $SILVA_imputed_file;
} elsif ($global_comp_DB_type eq "MERGED") {
    print "Using the merged database, APP can only generate and normalise OTU tables.\n";
    print "Alpha and Beta diversities will need to be performed manually.\n";
    $TAX_tax_file = $MERGED_TAX_tax_file;
    $TAX_blast_file = $MERGED_TAX_blast_file;
} elsif ($global_comp_DB_type eq 'GG') {
    ;
} else {
    print STDERR "Invalid database defined '$global_comp_DB_type' !!!\n\n";
    exit 1;
}

#### Check the user options and override if required

$TAX_tax_file = overrideDefault($TAX_tax_file, "taxonomy");
$TAX_blast_file = overrideDefault($TAX_blast_file, "blast");
$imputed_file = overrideDefault($imputed_file, "imputed");

print "Using taxonomy file: $TAX_tax_file\n";
print "Using blast file: $TAX_blast_file\n";
print "Using aligned blast file: $TAX_aligned_blast_file\n";
print "Using imputed file: $QIIME_imputed_file\n";

#### Sanity checks for the input files
checkFileExists($TAX_tax_file);
if ($assign_taxonomy_method eq 'blast'){
    checkFileExists("$TAX_blast_file.nsq");
    checkFileExists("$TAX_blast_file.nin");
    checkFileExists("$TAX_blast_file.nhr");
} elsif ($assign_taxonomy_method eq 'bwasw'){
    my @extensions = ('amb','ann','bwt','pac','sa');
    foreach $_ (@extensions)
    {
        checkFileExists($TAX_blast_file.'.'.$_);
    }
}

#### Start the results pipeline!
print "All good!\n";

#### Create a communication bridge with R and start R
my $global_R_instance = Statistics::R->new();
$global_R_instance->start();

####
#### NN DATA SET PROCESSING!
####
print "----------------------------------------------------------------\n";
print "Start TABLE BASED NORMALISATION data set processing...\n";
print "----------------------------------------------------------------\n";
print "Copying reads for analysis...\n";
copy_read_subset("$global_QA_dir/$global_acacia_output_dir/$ACACIA_out_file","$global_TB_processing_dir/$nn_fasta_file");

print "Beginning OTU and Taxonomy module...\n";
print "Picking OTUs for non normalised data set...\n";
checkFileExists("$global_TB_processing_dir/$nn_fasta_file");
checkAndRunCommand("pick_otus.py", [{-i => "$global_TB_processing_dir/$nn_fasta_file",
                                        -s => "$global_similarity_setting",
                                        -o => "$global_TB_processing_dir/uclust_picked_otus"}], DIE_ON_FAILURE);

print "Getting a representative set...\n";
checkFileExists("$global_TB_processing_dir/uclust_picked_otus/$nn_otus_file");
checkFileExists("$global_TB_processing_dir/$nn_fasta_file");

checkAndRunCommand("pick_rep_set.py", [{-i => "$global_TB_processing_dir/uclust_picked_otus/$nn_otus_file",
                                        -f => "$global_TB_processing_dir/$nn_fasta_file"}], DIE_ON_FAILURE);
                      

# if we are doing OTU_AVERAGE (or if we've been asked to) then we need to assign taxonomy here
print "Assigning taxonomy for non normalised data set...\n";
my $nn_rep_set_fasta = "$global_TB_processing_dir/".$nn_fasta_file."_rep_set.fasta";
checkFileExists($nn_rep_set_fasta);

print "Assign taxonomy method: $assign_taxonomy_method\n";
if ($assign_taxonomy_method eq 'blast') {
    checkAndRunCommand("assign_taxonomy.py", [{-i => $nn_rep_set_fasta,
                                                -t => $TAX_tax_file,
                                                -b => $TAX_blast_file,
                                                -m => "blast",
                                                -a => $num_threads,
                                                -e => $global_e_value,
                                                -o => $global_TB_processing_dir}], DIE_ON_FAILURE);
} elsif ($assign_taxonomy_method eq 'bwasw') {
    checkAndRunCommand("assign_taxonomy.py", [{-i => $nn_rep_set_fasta,
                                                -t => $TAX_tax_file,
                                                -d => $TAX_blast_file,
                                                -m => "bwasw",
                                                -a => $num_threads,
                                                -o => $global_TB_processing_dir}], DIE_ON_FAILURE);
} else {
    die "Unrecognised assign_taxonomy method: '$assign_taxonomy_method'";
}

print "Making NON NORMALISED otu table...\n";
my $nn_rep_set_tax_assign = "$global_TB_processing_dir/".$nn_fasta_file."_rep_set_tax_assignments.txt";
checkFileExists($nn_rep_set_tax_assign);
checkFileExists("$global_TB_processing_dir/uclust_picked_otus/$nn_otus_file");
checkAndRunCommand("make_otu_table.py",  [{-i => "$global_TB_processing_dir/uclust_picked_otus/$nn_otus_file",
                                           -t => $nn_rep_set_tax_assign,
                                           -o => $nn_otu_table_file}], DIE_ON_FAILURE);

print "Making NON NORMALISED otu table (Extended format)...\n";


checkAndRunCommand("reformat_otu_table.py",  [{-i => "$nn_otu_table_file",
                                               -t => $nn_rep_set_tax_assign,
                                               -o => "$nn_expanded_otu_table_file"}], IGNORE_FAILURE);


# do rarefaction for unnormalised data
print "Rarefaction...\n";
checkAndRunCommand("multiple_rarefactions.py", [{-i => $nn_otu_table_file,
                                                 -o => "$global_TB_processing_dir/rarefied_otu_tables/",
                                                 -m => $global_rare_M,
                                                 -x => $global_rare_X,
                                                 -s => $global_rare_S,
                                                 -n => $global_rare_N}], DIE_ON_FAILURE);

# normalise the non normalised OTU table
print "Normalizing non normalised table at $global_norm_sample_size sequences... [$global_norm_sample_size, $global_norm_num_reps]\n";
checkAndRunCommand("multiple_rarefactions_even_depth.py", [{-i => $nn_otu_table_file,
                                                            -o => "$global_TB_processing_dir/rare_tables/",
                                                            -d => $global_norm_sample_size,
                                                            -n => $global_norm_num_reps},
                                                            ["--lineages_included"],
                                                            ["--k"]], DIE_ON_FAILURE);

print "Calculating centroid subsampled table...\n";
my $centroid_index = find_centroid_table("$global_TB_processing_dir/rare_tables/",
                                         $global_norm_sample_size,
                                         $tn_dist_file,
                                         $tn_log_file,
                                         $tn_R_log_file);
    
copy("$global_TB_processing_dir/rare_tables/rarefaction_$global_norm_sample_size"."_$centroid_index".".txt",
      "$tn_otu_table_file") or die "Copy Failed $!";

print "Reformatting normalised OTU table...\n";
checkAndRunCommand("reformat_otu_table.py",  [{-i => "$tn_otu_table_file",
                                               -t => $nn_rep_set_tax_assign,
                                               -o => "$tn_expanded_otu_table_file"}], IGNORE_FAILURE);

# Merged database 
if ($global_comp_DB_type eq "MERGED") {
    print "APP only supports OTU table generation and normalisation for the SILVA/GG merged database.\n";
    print "Stopping here.\n";
    exit(0);   
}


print "Summarizing by taxa.....\n";

checkAndRunCommand("summarize_taxa.py", [{-i => "$tn_otu_table_file",
                                          -o => "$global_TB_results_dir/breakdown_by_taxonomy/"}], DIE_ON_FAILURE);



# move 100 of the 100 tables just produced into a new folder for jacknifing
print "Jackknifing in preparation for beta diversity\n";
my $jack_knife_folder = "$global_TB_processing_dir/rare_tables/JN";
`mkdir -p $jack_knife_folder`;
foreach my $jn_file_counter (0..$global_JN_file_count)
{
    my $jn_from_file = "rarefaction_".$global_norm_sample_size."_".$jn_file_counter.".txt";
    `cp $global_TB_processing_dir/rare_tables/$jn_from_file $jack_knife_folder/`;
}

mkdir ("$global_TB_processing_dir/blast_substituted_phylogeny/") or die "Unable to make blast_substituted_phylogeny directory.";

checkAndRunCommand("generate_fasta_from_taxonomy.py",[{-t => $nn_rep_set_tax_assign,
                                                       -r => $nn_rep_set_fasta,
                                                       -b => $TAX_aligned_blast_file,
                                                       -o => "$global_TB_processing_dir/blast_substituted_phylogeny/non_normalised.fa_rep_set_aligned_pfiltered.fasta"}], DIE_ON_FAILURE);

# We only have to align and filter the de novo phylogeny sequences because the 
# BLAST substitution will have taken the sequences from the aligned taxonomy file
# to produce an aligned FASTA

print "Treeing non normalised data set...\n";
checkAndRunCommand("align_seqs.py", [{-i => $nn_rep_set_fasta,
                                      -t => $imputed_file,
                                      -p => 0.6,
                                      -o => "$global_TB_processing_dir/de_novo_phylogeny/pynast_aligned"}], DIE_ON_FAILURE);

my $nn_rep_set_aligned = "$global_TB_processing_dir/de_novo_phylogeny/pynast_aligned/".$nn_fasta_file."_rep_set_aligned.fasta";
&checkFileExists($nn_rep_set_aligned);
checkAndRunCommand("filter_alignment.py", [{-i => $nn_rep_set_aligned,
                                            -o => "$global_TB_processing_dir/de_novo_phylogeny/"}], DIE_ON_FAILURE);

# From here on in, processing will occur twice, once for the de-novo alignment
# and another for the blast substituted phylogeny
#
my @phylogeny_dirs = ("de_novo_phylogeny", "blast_substituted_phylogeny");

foreach my $phylogeny_dir (@phylogeny_dirs) {  
    
    # Define processing and results directories
    my $phylogeny_proc_dir = "$global_TB_processing_dir/$phylogeny_dir";
    my $phylogeny_results_dir = "$global_TB_results_dir/$phylogeny_dir";
    
    my $nn_rep_set_aligned_filtered = "$phylogeny_proc_dir/".$nn_fasta_file."_rep_set_aligned_pfiltered.fasta";
    &checkFileExists($nn_rep_set_aligned_filtered);
    checkAndRunCommand("make_phylogeny.py", [{-i => $nn_rep_set_aligned_filtered,
                                              -r => "midpoint"}], DIE_ON_FAILURE);

    my $nn_rep_set_aligned_filtered_tree = "$phylogeny_proc_dir/$nn_fasta_file"."_rep_set_aligned_pfiltered.tre";
    &checkFileExists($nn_rep_set_aligned_filtered_tree);
    move( $nn_rep_set_aligned_filtered_tree, "$phylogeny_proc_dir/$nn_tree_file") or die "Move failed $!";

    print "Calculating alpha-diversity metrics...\n";
    checkAndRunCommand("alpha_diversity.py", [{-i => "$global_TB_processing_dir/rarefied_otu_tables/",
                                               -t => "$phylogeny_proc_dir/$nn_tree_file",
                                               -o => "$phylogeny_proc_dir/alpha_div/",
                                               -m => join(",", qw(chao1
                                                                  chao1_confidence
                                                                  PD_whole_tree
                                                                  observed_species
                                                                  simpson
                                                                  shannon
                                                                  fisher_alpha))}], DIE_ON_FAILURE);                                                       
                                                           
    checkAndRunCommand("collate_alpha.py", [{-i => "$phylogeny_proc_dir/alpha_div/",
                                             -o => "$phylogeny_proc_dir/alpha_div_collated/"}], DIE_ON_FAILURE);
    
    # make the same image twice (two different formats)
    foreach my $format (("png", "svg")) {
        checkAndRunCommand("make_rarefaction_plots.py", [{-i => "$phylogeny_proc_dir/alpha_div_collated/",
                                                          -m => "$global_QA_dir/qiime_mapping.txt",
                                                          -o => "$phylogeny_results_dir/alpha_diversity/",
                                                          "--resolution" => 300,
                                                          "--imagetype" => $format}], DIE_ON_FAILURE);
    }
    
    print "Jacknifed beta diversity....\n";
    jack_knifing("weighted_unifrac", $nn_otu_table_file, "$phylogeny_proc_dir/$nn_tree_file", $jack_knife_folder, $phylogeny_results_dir);
    jack_knifing("unweighted_unifrac", $nn_otu_table_file, "$phylogeny_proc_dir/$nn_tree_file", $jack_knife_folder, $phylogeny_results_dir);
    jack_knifing("euclidean", $nn_otu_table_file, "$phylogeny_proc_dir/$nn_tree_file", $jack_knife_folder, $phylogeny_results_dir);
    jack_knifing("hellinger", $nn_otu_table_file, "$phylogeny_proc_dir/$nn_tree_file", $jack_knife_folder, $phylogeny_results_dir);


    # beta for the normalised table
    my @beta_methods = ('weighted_unifrac','unweighted_unifrac','euclidean','hellinger','binary_euclidean','chord');
    foreach my $matrix_type (@beta_methods)
    {    
        my $params = [{-i => "$tn_otu_table_file",
                       -o => "$phylogeny_results_dir/beta_diversity/$tn_prefix/$matrix_type",
                       -m => $matrix_type}];
        
        # qiime spews when you give it a tree for some methods 
        if(!($matrix_type =~ /euclid/))
        {
            push @{$params}, {-t => "$phylogeny_proc_dir/$nn_tree_file"};
        }
        
        checkAndRunCommand("beta_diversity.py", $params, DIE_ON_FAILURE);
    
        
        my $beta_file = $matrix_type."_".$tn_prefix."_otu_table.txt";
        my $upgma_file = $matrix_type."_".$tn_prefix."_otu_table_upgma.tre";
        my $pcoa_file = $matrix_type."_".$tn_prefix."_pcoa.txt";
    
        # Perform UPGMA clustering on rarefied distance matrices 
        checkAndRunCommand("upgma_cluster.py", [{-i => "$phylogeny_results_dir/beta_diversity/$tn_prefix/$matrix_type/$beta_file",
                                                    -o => "$phylogeny_results_dir/beta_diversity/$tn_prefix/$matrix_type/$upgma_file"}], DIE_ON_FAILURE);    
        
        # Compute principal coordinates
        checkAndRunCommand("principal_coordinates.py", [{-i => "$phylogeny_results_dir/beta_diversity/$tn_prefix/$matrix_type/$beta_file",
                                                            -o => "$phylogeny_results_dir/beta_diversity/$tn_prefix/$matrix_type/$pcoa_file"}], DIE_ON_FAILURE);
    }
}

# shuffle the results around

####
#### END NN DATA SET PROCESSING
####

####
#### SEQ NORMALISED DATA SET PROCESSING
####
# these file will exist if we have used the seq centroid method
if($global_norm_style eq "SEQ")
{
    print "----------------------------------------------------------------\n";
    print "Start SEQ normalised data processing...\n";
    print "----------------------------------------------------------------\n";

    # make the seq processing dir if we need to
    makeResultsDirs(1);

    # find centroid sequences
    find_centroid_sequences($global_SB_processing_dir."/norm_tables/", 7,
                            "$global_TB_processing_dir/$nn_fasta_file",
                            "$global_SB_processing_dir/$sn_fasta_file",
                            $sn_dist_file, $sn_log_file);

    # once we've found our guy, lets make an otu table with assigned taxonomies   
    print "Picking OTUs for SEQ normalised data set...\n";
    checkAndRunCommand("pick_otus.py", [{-i => "$global_SB_processing_dir/$sn_fasta_file",
                                         -s => $global_similarity_setting,
                                         -o => "$global_SB_processing_dir/uclust_picked_otus"}], DIE_ON_FAILURE);

    print "Gettting a representitive set...\n";
    checkAndRunCommand("pick_rep_set.py", [{-i => "$global_SB_processing_dir/uclust_picked_otus/$sn_otus_file",
                                            -f => "$global_SB_processing_dir/$sn_fasta_file"}], DIE_ON_FAILURE);

    # if we are doing OTU_AVERAGE (or if we've been asked to) then we need to assign taxonomy here
    print "Assigning taxonomy for SEQ normalised data set...\n";
    my $sn_rep_set_fasta = "$global_SB_processing_dir/".$sn_fasta_file."_rep_set.fasta";
    print "Assign taxonomy method: $assign_taxonomy_method\n";
        if ($assign_taxonomy_method eq 'blast') {
            checkAndRunCommand("assign_taxonomy.py", [{-i => $sn_rep_set_fasta,
                                                        -t => $TAX_tax_file,
                                                        -b => $TAX_blast_file,
                                                        -m => "blast",
                                                        -a => $num_threads,
                                                        -e => $global_e_value,
                                                        -o => $global_SB_processing_dir}], DIE_ON_FAILURE);
        } elsif ($assign_taxonomy_method eq 'bwasw') {
            checkAndRunCommand("assign_taxonomy.py", [{-i => $sn_rep_set_fasta,
                                                        -t => $TAX_tax_file,
                                                        -d => $TAX_blast_file,
                                                        -m => "bwasw",
                                                        -a => $num_threads,
                                                        -o => $global_SB_processing_dir}], DIE_ON_FAILURE);
        } else {
            die "Unrecognised assign_taxonomy method: '$assign_taxonomy_method'";
        }
    print "Making SEQ normalised otu table...\n";
    my $sn_rep_set_tax_assign = "$global_SB_processing_dir/".$sn_fasta_file."_rep_set_tax_assignments.txt";
    checkAndRunCommand("make_otu_table.py", [{-i => "$global_SB_processing_dir/uclust_picked_otus/$sn_otus_file",
                                              -t => $sn_rep_set_tax_assign,
                                              -o => $sn_otu_table_file}], DIE_ON_FAILURE);
    
    checkAndRunCommand("reformat_otu_table.py",  [{-i => $sn_otu_table_file,
                                                   -t => $sn_rep_set_tax_assign,
                                                   -o => "$sn_expanded_otu_table_file"}], IGNORE_FAILURE);

    print "Summarizing by taxa.....\n";
    checkAndRunCommand("summarize_taxa.py", [{-i => $sn_otu_table_file,
                                              -o => "$global_SB_results_dir/breakdown_by_taxonomy/"}], DIE_ON_FAILURE);
    
    mkdir ("$global_SB_processing_dir/blast_substituted_phylogeny/") or die "Unable to make blast_substituted_phylogeny directory.";

    checkAndRunCommand("generate_fasta_from_taxonomy.py",[{-t => $sn_rep_set_tax_assign,
                                                           -r => $sn_rep_set_fasta,
                                                           -b => $TAX_aligned_blast_file,
                                                           -o => "$global_SB_processing_dir/blast_substituted_phylogeny/$sn_fasta_file"."_rep_set_aligned_pfiltered.fasta"}], DIE_ON_FAILURE);
    
    print "Treeing SEQ normalised data set...\n";
    checkAndRunCommand("align_seqs.py", [{-i => $sn_rep_set_fasta,
                                          -t => $imputed_file,
                                          -p => 0.6,
                                          -o => "$global_SB_processing_dir/de_novo_phylogeny/pynast_aligned"}], DIE_ON_FAILURE);
    
    my $sn_rep_set_aligned = "$global_SB_processing_dir/de_novo_phylogeny/pynast_aligned/".$sn_fasta_file."_rep_set_aligned.fasta";
    &checkFileExists($sn_rep_set_aligned);
    checkAndRunCommand("filter_alignment.py", [{-i => $sn_rep_set_aligned,
                                                -o => "$global_SB_processing_dir/de_novo_phylogeny/"}], DIE_ON_FAILURE);
    
    # From here on in, processing will occur twice, once for the de-novo alignment
    # and another for the blast substituted phylogeny
    #
    my @phylogeny_dirs = ("de_novo_phylogeny", "blast_substituted_phylogeny");
    
    foreach my $phylogeny_dir (@phylogeny_dirs) {  
    
        # Define processing and results directories
        my $phylogeny_proc_dir = "$global_SB_processing_dir/$phylogeny_dir";
        my $phylogeny_results_dir = "$global_SB_results_dir/$phylogeny_dir";
        if (! -e $phylogeny_results_dir) {
             mkdir ($phylogeny_results_dir)
        }
    
        my $sn_rep_set_aligned_filtered = "$phylogeny_proc_dir/".$sn_fasta_file."_rep_set_aligned_pfiltered.fasta";
        &checkFileExists($sn_rep_set_aligned_filtered);
        checkAndRunCommand("make_phylogeny.py", [{-i => $sn_rep_set_aligned_filtered,
                                                  -r => "midpoint"}], DIE_ON_FAILURE);
        
        my $sn_rep_set_aligned_filtered_tree = "$phylogeny_proc_dir/".$sn_fasta_file."_rep_set_aligned_pfiltered.tre";
        move( $sn_rep_set_aligned_filtered_tree, "$phylogeny_proc_dir/$sn_tree_file");
    
        print "Alpha and beta diversity for SEQ normalised otu table...\n";
        # alpha
        
        checkAndRunCommand("alpha_diversity.py", [{-i => $sn_otu_table_file,
                                           -t => "$phylogeny_proc_dir/$sn_tree_file",
                                           -o => "$phylogeny_results_dir/alpha_diversity.txt",
                                           -m => join(",", qw(chao1
                                                              chao1_confidence
                                                              PD_whole_tree
                                                              observed_species
                                                              simpson
                                                              shannon
                                                              fisher_alpha))}], DIE_ON_FAILURE);
        #beta
        my @beta_methods = ('weighted_unifrac','unweighted_unifrac','euclidean','hellinger','binary_euclidean','chord');
        foreach my $matrix_type (@beta_methods)
        {
            my $params = {-i => $sn_otu_table_file,
                          -o => "$phylogeny_results_dir/beta_diversity/$matrix_type",
                          -m => $matrix_type};
        
            # qiime spews when you give it a tree for some methods 
            if(!($matrix_type =~ /euclid/))
            {
                $params->{"-t"} = "$phylogeny_proc_dir/$sn_tree_file";
            }
            
            checkAndRunCommand("beta_diversity.py", [$params], DIE_ON_FAILURE);
            
            my $beta_file = $matrix_type."_".$sn_prefix."_otu_table.txt";
            my $upgma_file = $matrix_type."_".$sn_prefix."_otu_table_upgma.tre";
            my $pcoa_file = $matrix_type."_".$sn_prefix."_pcoa.txt";
            
            # Perform UPGMA clustering on rarefied distance matrices 
            checkAndRunCommand("upgma_cluster.py", [{-i => "$phylogeny_results_dir/beta_diversity/$matrix_type/$beta_file",
                                                     -o => "$phylogeny_results_dir/beta_diversity/$matrix_type/$upgma_file"}], DIE_ON_FAILURE);
    
            # Compute principal coordinates
            checkAndRunCommand("principal_coordinates.py", [{-i => "$phylogeny_results_dir/beta_diversity/$matrix_type/$beta_file",
                                                             -o => "$phylogeny_results_dir/beta_diversity/$matrix_type/$pcoa_file"}], DIE_ON_FAILURE);
        }

    }
}
    
print "Results are located in: $output_dir/results/\n";

####
#### END SEQ NORMALISED DATA SET PROCESSING
####

# stop the interpreter
$global_R_instance->stop();

######################################################################
# CUSTOM SUBS
######################################################################
sub jack_knifing
{
    #-----
    # Do jackknifing for a specific type of matrix
    #
    my ($matrix_type, $raw_otu_table, $raw_tree, $jack_knife_folder, $phylogeny_folder) = @_;

    # Produce distance matrix reflecting beta diversity in non-normalised OTU table
    checkAndRunCommand("beta_diversity.py", [{-i => $raw_otu_table,
                                              -o => "$phylogeny_folder/beta_diversity/$nn_prefix/$matrix_type",
                                              -m => $matrix_type,
                                              -t => $raw_tree}], DIE_ON_FAILURE);

    # Perform UPGMA clustering on non-normalised distance matrix
    my $beta_otu_table = $matrix_type."_".$nn_prefix."_otu_table.txt";
    my $upgma_cluster_tree = $matrix_type."_".$nn_prefix."_otu_table_upgma.tre";

    checkAndRunCommand("upgma_cluster.py", [{-i => "$phylogeny_folder/beta_diversity/$nn_prefix/$matrix_type/$beta_otu_table",
                                             -o => "$phylogeny_folder/beta_diversity/$nn_prefix/$matrix_type/$upgma_cluster_tree"}], DIE_ON_FAILURE);

    # Produce distance matrices reflecting beta diversity in the rarefied OTU tables
    checkAndRunCommand("beta_diversity.py", [{-i => $jack_knife_folder,
                                              -o => "$phylogeny_folder/beta_diversity/$nn_prefix/$matrix_type/rare_dm/",
                                              -m => $matrix_type,
                                              -t => $raw_tree}], DIE_ON_FAILURE);

    # Perform UPGMA clustering on rarefied distance matrices 
    checkAndRunCommand("upgma_cluster.py", [{-i => "$phylogeny_folder/beta_diversity/$nn_prefix/$matrix_type/rare_dm/",
                                             -o => "$phylogeny_folder/beta_diversity/$nn_prefix/$matrix_type/rare_upgma/"}], DIE_ON_FAILURE);

    # Compare the consensus tree to the beta-derived trees
    checkAndRunCommand("tree_compare.py", [{-s => "$phylogeny_folder/beta_diversity/$nn_prefix/$matrix_type/rare_upgma/",
                                            -m => "$phylogeny_folder/beta_diversity/$nn_prefix/$matrix_type/$upgma_cluster_tree",
                                            -o => "$phylogeny_folder/beta_diversity/$nn_prefix/$matrix_type/upgma_cmp/"}], DIE_ON_FAILURE);

    # Compute principal coordinates
    checkAndRunCommand("principal_coordinates.py", [{-i => "$phylogeny_folder/beta_diversity/$nn_prefix/$matrix_type/rare_dm/",
                                                     -o => "$phylogeny_folder/beta_diversity/$nn_prefix/$matrix_type/pcoa/"}], DIE_ON_FAILURE);

    # Make PDF of Jackknife tree with labeled support: weighted unifrac command
    my $output_pdf = "$phylogeny_folder/$matrix_type"."_betadiv_jackknife_tree.pdf";
    checkAndRunCommand("make_bootstrapped_tree.py", [{-m => "$phylogeny_folder/beta_diversity/$nn_prefix/$matrix_type/upgma_cmp/master_tree.tre",
                                                      -s => "$phylogeny_folder/beta_diversity/$nn_prefix/$matrix_type/upgma_cmp/jackknife_support.txt",
                                                      -o => $output_pdf}], DIE_ON_FAILURE);

}

sub run_R_cmd
{
    #-----
    # Wrapper for running R commands
    #
    my ($cmd, $log_fh) = @_;
    if ($log_fh) {
        print {$log_fh} $cmd, "\n";
    }
    $global_R_instance->run($cmd);
}

sub find_centroid_table
{
    #-----
    # Find a representative set of $global_norm_sample_size sequences (do this in RRRRRR...)
    #
    my($path, $norm_size, $dist_file, $log_file, $R_cmd_log) = @_;

    open (my $R_cmd_log_fh, ">", $R_cmd_log)
        or warn("Unable to open $R_cmd_log for logging R commands. Proceeding without logging.");

    print "Calculating centroid OTU table from tables in $path...\n";

    $global_R_instance->start();

    my $sampl_p1 = $global_num_samples + 1;

    # read in the list of distance matricies
    run_R_cmd(qq`library(foreign);`, $R_cmd_log_fh);
    run_R_cmd(qq`a<-list.files("$path", "*.txt");`, $R_cmd_log_fh);

    # work out how many there are and allocate an array
    run_R_cmd(qq`len_a <- length(a);`, $R_cmd_log_fh);
    run_R_cmd(qq`big_frame <- array(0,dim=c($global_num_samples,$global_num_samples,len_a));`, $R_cmd_log_fh);

    print "  --start loading data...\n";

    # load each file individually into a big frame
    my $r_str = "for (i in c(1:len_a)) { j <- i - 1; name <- paste(\"$path\",\"rarefaction_$norm_size\",\"_\",j,\".txt\",sep=\"\"); u<-read.table(name,sep=\"\\t\",row.names=1); u[,$sampl_p1]<-NULL; big_frame[,,i]<-as.matrix(dist(t(u), upper=TRUE, diag=TRUE)); i<-i+1; }";
    run_R_cmd($r_str, $R_cmd_log_fh);

    print "  --data loaded, calculating centroid...\n";

    # find the average matrix
    run_R_cmd(qq`ave <- big_frame[,,1];`, $R_cmd_log_fh);
    run_R_cmd(qq`for (i in c(2:len_a)) { ave <- ave + big_frame[,,i]; }`, $R_cmd_log_fh);
    run_R_cmd(qq`ave <- ave/len_a;`, $R_cmd_log_fh);

    print "  --calculating distances of tables to centroid...\n";

    # find the euclidean distance of each matrix from the average
    run_R_cmd(qq`dist<-array(0,dim=c(len_a));`, $R_cmd_log_fh);
    run_R_cmd(qq`for (i in c(1:len_a)) { dist[i] <- sqrt(sum(big_frame[,,i]-ave)^2); }`, $R_cmd_log_fh);

    # find the min value
    run_R_cmd(qq`min_index <- which.min(dist);`, $R_cmd_log_fh);
    my $centroid_otu_index = $global_R_instance->get('min_index');
    # R indexes from 0
    $centroid_otu_index--;
    print "  --table: $centroid_otu_index is the centroid table\n";

    # make stats on the distances
    # and log what we did
    open my $log_fh, ">", $log_file or die "Could not open log file: $log_file : $!\n";
    run_R_cmd(qq`max_dist <- max(dist);`, $R_cmd_log_fh);
    run_R_cmd(qq`min_dist <- min(dist);`, $R_cmd_log_fh);
    run_R_cmd(qq`range_dist <- max_dist - min_dist;`, $R_cmd_log_fh);
    run_R_cmd(qq`mean_dist <- mean(dist);`, $R_cmd_log_fh);
    run_R_cmd(qq`median_dist <- median(dist);`, $R_cmd_log_fh);

    print $log_fh "---------------------------------------------------\n";
    print $log_fh "  Centroid OTU table based normalised statistics\n";
    print $log_fh "---------------------------------------------------\n";
    print $log_fh "Max dist:\t".$global_R_instance->get('max_dist')."\n";
    print $log_fh "Min dist:\t".$global_R_instance->get('min_dist')."\n";
    print $log_fh "Range:\t".$global_R_instance->get('range_dist')."\n";
    print $log_fh "Mean:\t".$global_R_instance->get('mean_dist')."\n";
    print $log_fh "Median:\t".$global_R_instance->get('median_dist')."\n";

    if(2 < $global_num_samples)
    {
        run_R_cmd(qq`library(permute);`, $R_cmd_log_fh);
        run_R_cmd(qq`library(vegan);`, $R_cmd_log_fh);
        run_R_cmd(qq`mantel.otu <- mantel(ave,big_frame[,,min_index]);`, $R_cmd_log_fh);
        run_R_cmd(qq`m_stat <- mantel.otu\$statistic;`, $R_cmd_log_fh);
        run_R_cmd(qq`m_sig <- mantel.otu\$signif;`, $R_cmd_log_fh);
        print $log_fh "Mantel P stat:\t".$global_R_instance->get('m_sig')."\n";
        print $log_fh "Mantel R stat:\t".$global_R_instance->get('m_stat')."\n";
    }
    else
    {
        print "Too few samples to perform a mantel test.\n";
    }
    close $log_fh;

    # print all the distances to a file so we can make purdy pictures from them later
    open my $dist_fh, ">", $dist_file or die "Could not open distance file: $dist_file : $!\n";
    my $num_tables = $global_R_instance->get('len_a');
    foreach my $counter (1..$num_tables)
    {
        print $dist_fh $global_R_instance->get("dist[$counter]")."\n"
    }
    close $dist_fh;

    close($R_cmd_log_fh);
    # let the user know the result    
    return $centroid_otu_index;
}

sub find_centroid_sequences
{
    #-----
    # Find the centroid of a set of otu tables
    #
    my($path, $num_reps, $seqs, $final_seq_file, $dist_file, $log_file) = @_;
    
    print "Normalising read counts...\n";
    print "Input: $seqs\n";
    print "Num reps: $num_reps\n";
    
    # make a bunch of normalised read files
    `mkdir -p $path`;
    foreach my $file_counter (1..$num_reps)
    {
        print ".";
        # make a normalised file
        my $norm_seq_file_root = $path."norm_".$file_counter;
        my $norm_seq_file = $norm_seq_file_root.".fa";
        my $norm_otus_file = $norm_seq_file_root."_otus.txt";
        my $norm_otu_table_file = $norm_seq_file_root."_otu_table.txt";
        checkAndRunCommand("app_normalise_reads.pl", [{-in => $seqs,
                                                       -n => $global_norm_sample_size,
                                                       -o => $norm_seq_file}], DIE_ON_FAILURE); 
    
        # make an otu table for each guy
        checkAndRunCommand("pick_otus.py", [{-i => $norm_seq_file,
                                             -s => $global_similarity_setting,
                                             -o => $path}], DIE_ON_FAILURE);
        
        checkAndRunCommand("make_otu_table.py", [{-i => $norm_otus_file,
                                                  -o => $norm_otu_table_file}], DIE_ON_FAILURE);    
    }

    print "\n";
    
    # find out which table is closest to the ave found for centroid OTU table found previously. Use R
    # variables should still be in memory
    my $sampl_p1 = $global_num_samples + 1;
    
    # make a 3d data frame again
    $global_R_instance->run(qq`new_frame <- array(0,dim=c($global_num_samples,$global_num_samples,$num_reps));`);
    
    # load each file individually into a big frame
    my $r_str = "for (i in c(1:$num_reps)) { name <- paste(\"$path"."norm_\",i,\"_otu_table.txt\",sep=\"\"); u<-read.table(name,sep=\"\\t\",row.names=1); u[,$sampl_p1]<-NULL; new_frame[,,i]<-as.matrix(dist(t(u), upper=TRUE, diag=TRUE)); i<-i+1; }";
    $global_R_instance->run($r_str);    
    
    # find the euclidean distance of each matrix from the average
    $global_R_instance->run(qq`dist_seq<-array(0,dim=c($num_reps));`);
    $global_R_instance->run(qq`for (i in c(1:$num_reps)) { dist_seq[i] <- sqrt(sum(new_frame[,,i]-ave)^2); }`);
    
    # find the min value
    $global_R_instance->run(qq`closest_norm <- which.min(dist_seq);`);
    my $closest_norm_index = $global_R_instance->get('closest_norm');
    
    # copy the sequence file over that we wish to use
    my $cp_cmd = "cp  $path"."norm_$closest_norm_index".".fa ".$final_seq_file;
    `$cp_cmd`; 
    
    # make stats on the distances
    # and log what we did
    open my $log_fh, ">", $log_file or die "Could not open log file: $log_file : $!\n";
    $global_R_instance->run(qq`max_dist <- max(dist_seq);`);
    $global_R_instance->run(qq`min_dist <- min(dist_seq);`);
    $global_R_instance->run(qq`range_dist <- max_dist - min_dist;`);
    $global_R_instance->run(qq`mean_dist <- mean(dist_seq);`);
    $global_R_instance->run(qq`median_dist <- median(dist_seq);`);
    
    print $log_fh "-----------------------------------------------\n";
    print $log_fh "  Read based normalised statistics\n";
    print $log_fh "-----------------------------------------------\n";
    print $log_fh "Max dist:\t".$global_R_instance->get('max_dist')."\n";
    print $log_fh "Min dist:\t".$global_R_instance->get('min_dist')."\n";
    print $log_fh "Range:\t".$global_R_instance->get('range_dist')."\n";
    print $log_fh "Mean:\t".$global_R_instance->get('mean_dist')."\n";
    print $log_fh "Median:\t".$global_R_instance->get('median_dist')."\n";
  
    if(2 < $global_num_samples)
    {
        $global_R_instance->run(qq`mantel.otu <- mantel(ave,new_frame[,,closest_norm]);`);
        $global_R_instance->run(qq`m_stat <- mantel.otu\$statistic;`);
        $global_R_instance->run(qq`m_sig <- mantel.otu\$signif;`);
        print $log_fh "Mantel P stat:\t".$global_R_instance->get('m_sig')."\n";
        print $log_fh "Mantel R stat:\t".$global_R_instance->get('m_stat')."\n";
    }
    else
    {
        print "Too few samples to perform a mantel test.\n";
    }
    
    # print all the distances to a file so we can make purdy pictures from them later
    open my $dist_fh, ">", $dist_file or die "Could not open distance file: $dist_file : $!\n";
    foreach my $counter (1..$num_reps)
    {
        print $dist_fh $global_R_instance->get("dist_seq[$counter]")."\n"         
    }
    close $dist_fh;

    print "Done\n";
}

sub find_global_norm_sample_size
{
    #-----
    # pick a sane mimium number of reads
    # Take nearest multiple of 50 under the minimum size... ...more or less
    #
    my $twenny_under = $global_min_sample_size - 20;
    my $nearest_fifty = int($global_min_sample_size / 50) * 50;
    $global_norm_sample_size = $nearest_fifty;
    if($nearest_fifty > $twenny_under)
    {
        $global_norm_sample_size -= 50;
    }
    if($global_norm_sample_size <= 0)
    {
        die "Your least abundant sample with $global_min_sample_size sequences is too small to continue!\nPlease remove this sample (and any with similar numbers) and try again.\n";
    }
    print "Normalised sample size calculated at: $global_norm_sample_size reads\n";
    return $global_norm_sample_size;
}

sub copy_read_subset
{
    #-----
    # copy over the denoised reads into the processing dir
    # only take reads whose IDs are in the $global_samp_ID_list
    #
    my($source_fasta, $target_fasta) = @_;
    open my $s_fh, "<", $source_fasta or die "**ERROR: could not open file: $source_fasta $!\n";
    open my $t_fh, ">", $target_fasta or die "**ERROR: could not open file: $target_fasta $!\n";
    my $good_seq = 0;
    my %seen_seqs = ();
    while(<$s_fh>)
    {
        if($_ =~ /^>/)
        {
            # header
            my @components = split / /, $_;
            # we need to ignore this guy if he's a duplicate
            if(!exists $seen_seqs{$components[0]})
            {
                my $header = $components[0];
                $header =~ s/>([^_]*)_.*/$1/;
                if(1 == $global_samp_ID_list{$header})
                {
                    $good_seq = 1;
                    print $t_fh $_;
                }
                $seen_seqs{$components[0]} = 1;
            }
        }
        elsif(1 == $good_seq)
        {
            print $t_fh $_;
            $good_seq = 0;
        }
    }
    close $s_fh;
    close $t_fh;
}

sub parse_config_results
{
    #-----
    # parse the app config file and produce a qiime mappings file
    #
    open my $conf_fh, "<", $options->{'config'} or die $!;
    print $global_mapping_file;
    open my $mapping, ">", $global_mapping_file or die $!;
    print $mapping "$FNB_HEADER\n";
    
    my $full_rarefaction = 0;
    
    while(<$conf_fh>)
    {
        next if($_ =~ /^#/);
        last if($_ =~ /^@/);
        chomp $_;
        my @fields = split /\t/, $_;
        
        # save the MID for later
        $global_samp_ID_list{$fields[$FNA{'SampleID'}]} = $fields[$FNA{'USE'}];
        
        # we need to find the globally maximal number of sequences for any USED sample
        if("1" eq int $fields[$FNA{'USE'}])
        {
            $global_num_samples++;
            my $sample_size = int $fields[$FNA{'ACC'}];
            if($sample_size > $global_rare_X)
            {
                $global_rare_X = $sample_size;
            }
            if($sample_size < $global_min_sample_size)
            {
                $global_min_sample_size = $sample_size;
            }
            print $mapping "$fields[$FNA{'SampleID'}]\t$fields[$FNA{'BarcodeSequence'}]\t$fields[$FNA{'LinkerPrimerSequence'}]\t$fields[$FNA{'Description'}]\n";
        }
    }
    close $mapping;
    
    print "\t...Processing $global_num_samples samples\n";
    # user options section
    while(<$conf_fh>)
    {
        chomp $_;
        my @fields = split /=/, $_;
        if($#fields > 0)
        {
            if($fields[0] eq "DB")
            {
               if($fields[1] eq "SILVA")
                {
                    $global_comp_DB_type = "SILVA";
                } elsif($fields[1] eq "MERGED")
                {
                    $global_comp_DB_type = "MERGED";
                } elsif($fields[1] eq '0') {
                    # this is just the default
                } else {
                    print STDERR "Invalid database defined '$fields[1]' !!!\n\n";
                    exit 1;
                }
            }
            elsif($fields[0] eq "NORMALISE")
            {
                # is this guy set?
                if($fields[1] ne "")
                {
                    chomp $fields[1];
                    my @norm_fields = split /,/, $fields[1];
                    $global_norm_style = $norm_fields[0];

                    # check the normailisation style is sane
                    if($global_norm_style ne "SEQ" and $global_norm_style ne "TABLE")
                    {
                        die "You must specify 'SEQ' or 'TABLE' as normalisation methods (if you specify anything)\n";
                    }

                    # see if the user decided on how many sequences to normalise to!
                    if($#norm_fields == 0)
                    {
                        # user did not specify an amount to normalise by
                        # select automatically
                        print "Finding normalisation size automatically\n";
                        find_global_norm_sample_size();
                    }
                    else
                    {
                        # user selected an amount, use this
                        $global_norm_sample_size = int $norm_fields[1];
                        print "Using user defined sample size of: $global_norm_sample_size\n";
                        if($global_norm_sample_size <= 0)
                        {
                            die "You need to specify a number greater than or equal to zero (or none at all).\nHint: try 'NORMALISE=$global_norm_style' or 'NORMALISE=$global_norm_style,XXX' where XXX is the number of sequences to normalise to\n";
                        }
                    }
                }
            }
            elsif($fields[0] eq "MUL_RARE_M")
            {
                # is this guy set?
                if($fields[1] ne "")
                {
                    my $user_min = int($fields[1]);
                    if($user_min > 0)
                    {
                        $global_rare_M = $user_min;
                    }
                }
            }
            elsif($fields[0] eq "MUL_RARE_X")
            {
                # is this guy set?
                if($fields[1] ne "")
                {
                    my $user_max =  int($fields[1]);
                    if($user_max < $global_rare_X)
                    {
                        $global_rare_X = $user_max;
                    }
                }
            }
            elsif($fields[0] eq "MUL_RARE_S")
            {
                # is this guy set?
                if($fields[1] ne "")
                { 
                    my $user_step = int($fields[1]);
                    if($user_step > 0)
                    {
                        $global_rare_S = $user_step;
                    }
                }
            }
            elsif($fields[0] eq "MUL_RARE_N")
            {
                # is this guy set?
                if($fields[1] ne "")
                { 
                    my $user_num = int($fields[1]);
                    if($user_num > 0)
                    {
                        $global_rare_N = $user_num;
                    }
                }
            }
            elsif($fields[0] eq "NUM_THREADS")
            {
                if($fields[1] ne "")
                { 
                    my $user_num = int($fields[1]);
                    if($user_num > 0)
                    {
                        $num_threads = $user_num;
                    }
                }
            }
            elsif($fields[0] eq "TRUNCATE_RAREFACTION_DEPTH")
            {
                if((uc($fields[1]) eq "FALSE") || ($fields[1] == 0))
                {
                    $full_rarefaction = 1;
                }
            }
        }
    }    
    close $conf_fh;
    
    # finally, check to see if we've picked a normailisation value
    if(0 == $global_norm_sample_size)
    {
        # user did not specify an amount to normalise by
        # select automatically
        print "Finding normalisation size automatically\n";
        find_global_norm_sample_size();
    }
    if (! $full_rarefaction) {
        if ($global_rare_X > 20000) {
            print "#####################################################\n" .
                  "WARNING: A sample has more than 20000 reads, this\n".
                  "can cause delays in calculating alpha diversity and\n" .
                  "normalized OTU tables due to the number of OTU table\n" .
                  "rarefactions that need to be calcuated.\n\n" .
                  "APP will only calculate rarefactions to 20000 reads\n" .
                  "and the normalisation table size will also be reduced\n" .
                  "to this size (if necessary).\n\n" .
                  "To overwrite this behaviour, add the following line\n" .
                  "to the bottom of the app config file:\n\n" .
                  "TRUNCATE_RAREFACTION_DEPTH=FALSE\n\n" .
                  "#####################################################\n";
            $global_rare_X = 20000;
            if ($global_norm_sample_size > 20000) {
                $global_norm_sample_size = 20000;
            }
        }
    } 
}

######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    my @standard_options = ( "help|h+", "config|c:s", "identity|i:i", "e:f",
                             "b|blast:s", "t|taxonomy:s", "i|imputed:s",
                             "a|assign-taxonomy-method:s",'threads:i');
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
    if(! -e $options{'config'} ) {
        print "**ERROR: " . $options{'config'} ." does not exist\n";
        exec("pod2usage $0");
    }
    if(exists $options{'identity'})
    {
        if(($options{'identity'} <= 0) || ($options{'identity'} > 1))
        {
            die "Identity must be an integer greater than 0 and no greater than 1\n";
        }
    }
    my @taxonomy_asignment_methods = ('blast','bwasw');
    if(exists $options{'assign-taxonomy-method'})
    {
        if(!(grep {$_ eq $options{'assign-taxonomy-method'}} @taxonomy_asignment_methods))
        {
            die "Taxonomy assignment method '$options{'assign-taxonomy-method'}' is not acceptable.";
        }
    }
    
    
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

sub overrideDefault {
    #-----
    # Set and override default values for parameters
    #
    my ($default_value, $option_name) = @_;
    if(exists $options->{$option_name}) 
    {
        return $options->{$option_name};
    }
    return $default_value;
}


__DATA__

=head1 NAME

    app_make_results.pl

=head1 COPYRIGHT

   copyright (C) 2011 Michael Imelfort and Paul Dennis, 2012 Connor Skennerton,
        Adam Skarshewski and Ben Woodcroft.

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

   Insert detailed description here

=head1 SYNOPSIS

    app_make_results.pl -c|config CONFIG_FILE [-help|h]

      -c CONFIG_FILE               App config file to be processed
      [-i identity VALUE]          Set blast identity for OTU clustering (pick_otus.py) [default: 97%]
      [-e EVALUE]                  Set e-value for blast (assign_taxonomy.py) [default 0.001]      
      [-b FILE]                    Path to a custom blast database / bwa database
      [-t FILE]                    Path to a custom taxonomy for otus
      [-i FILE]                    Path to a custom imputed file
      [-a FILE]                    assign_taxonomy method [default blast, alternative bwasw (for BWA)]
      [--threads NUM_THREADS]      Use this many threads where possible [default 5]
      [-help -h]                   Displays basic usage information
         
=cut


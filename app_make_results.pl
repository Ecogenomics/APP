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

#### GET THE WORKING DIR
getWorkingDirs($options->{'config'});

#### make these dirs
makeResultsDirs(0);

# sample IDs
my %global_samp_ID_list = ();
my $global_num_samples = 0;

# we can compare sequences to the greengenes or the SILVA dbs
my $global_comp_DB_type = "GG";

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

# how many reps to do when normalising the OTU table
my $global_norm_num_reps = 1000;

# we only need a subset of these guys to do jacknifing
my $global_JN_file_count = 100;
if($global_JN_file_count > $global_norm_num_reps) { $global_JN_file_count = $global_norm_num_reps - 1; } 

#### Override defaults from config or user
if(exists $options->{'identity'}) { $global_similarity_setting = $options->{'identity'}; }
if(exists $options->{'e'}) { $global_e_value = $options->{'e'}; }

print "Checking is all the config checks out...\t\t";
parse_config_results();

# update our databases (GG by default)
my $TAX_tax_file = $QIIME_TAX_tax_file;
my $TAX_blast_file = $QIIME_TAX_blast_file;
my $imputed_file = $QIIME_imputed_file;

if($global_comp_DB_type eq "SILVA")
{
    $TAX_tax_file = $SILVA_TAX_tax_file;
    $TAX_blast_file = $SILVA_TAX_blast_file;
    $imputed_file = $SILVA_imputed_file;
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

print "Picking OTUs for non normalised data set...\n";
`pick_otus.py -i $global_TB_processing_dir/$nn_fasta_file -s $global_similarity_setting -o $global_TB_processing_dir/uclust_picked_otus`;

print "Gettting a representitive set...\n";
`pick_rep_set.py -i $global_TB_processing_dir/uclust_picked_otus/$nn_otus_file -f $global_TB_processing_dir/$nn_fasta_file`;

# if we are doing OTU_AVERAGE (or if we've ben asked to) then we need to assign taxonomy here
print "Assigning taxonomy for non normalised data set...\n";
my $nn_rep_set_fasta = "$global_TB_processing_dir/".$nn_fasta_file."_rep_set.fasta";
`assign_taxonomy.py -i $nn_rep_set_fasta -t $TAX_tax_file -b $TAX_blast_file -m blast -e $global_e_value -o $global_TB_processing_dir`;

print "Treeing non normalised data set...\n";
`align_seqs.py -i $nn_rep_set_fasta -t $imputed_file -p 0.6 -o $global_TB_processing_dir/pynast_aligned`;
my $nn_rep_set_aligned = "$global_TB_processing_dir/pynast_aligned/".$nn_fasta_file."_rep_set_aligned.fasta";
`filter_alignment.py -i $nn_rep_set_aligned -o $global_TB_processing_dir`;
my $nn_rep_set_aligned_filtered = "$global_TB_processing_dir/".$nn_fasta_file."_rep_set_aligned_pfiltered.fasta";
`make_phylogeny.py -i $nn_rep_set_aligned_filtered -r midpoint`; 
my $nn_rep_set_aligned_filtered_tree = "$global_TB_processing_dir/".$nn_fasta_file."_rep_set_aligned_pfiltered.tre";
`mv $nn_rep_set_aligned_filtered_tree $nn_tree_file`;

print "Making NON NORMALISED otu table...\n";
my $nn_rep_set_tax_assign = "$global_TB_processing_dir/".$nn_fasta_file."_rep_set_tax_assignments.txt";
`make_otu_table.py -i $global_TB_processing_dir/uclust_picked_otus/$nn_otus_file -t $nn_rep_set_tax_assign -o $nn_otu_table_file`;

# do rarefaction for unnormalised data
print "Rarefaction and diversity...\n";
`multiple_rarefactions.py -i $nn_otu_table_file -o $global_TB_processing_dir/rarefied_otu_tables/ -m $global_rare_M -x $global_rare_X -s $global_rare_S -n $global_rare_N`;
`alpha_diversity.py -i $global_TB_processing_dir/rarefied_otu_tables/ -t $nn_tree_file -o $global_TB_processing_dir/alpha_div/ -m chao1,chao1_confidence,PD_whole_tree,observed_species,simpson,shannon,fisher_alpha`;
`collate_alpha.py -i $global_TB_processing_dir/alpha_div/ -o $global_TB_processing_dir/alpha_div_collated/`;

# make the same image twice (two different formats)
`make_rarefaction_plots.py -i $global_TB_processing_dir/alpha_div_collated/ -m $global_QA_dir/qiime_mapping.txt -o $global_TB_results_dir/alpha_diversity/ --resolution 300 --imagetype svg`;
`make_rarefaction_plots.py -i $global_TB_processing_dir/alpha_div_collated/ -m $global_QA_dir/qiime_mapping.txt -o $global_TB_results_dir/alpha_diversity/ --resolution 300 --imagetype png`;

# normalise the non normalised OTU table
print "Normalizing non normalised table at $global_norm_sample_size sequences... [$global_norm_sample_size, $global_norm_num_reps]\n";
`multiple_rarefactions_even_depth.py -i $nn_otu_table_file -o $global_TB_processing_dir/rare_tables/ -d $global_norm_sample_size -n $global_norm_num_reps --lineages_included --k`;

my $centroid_index = find_centroid_table("$global_TB_processing_dir/rare_tables/", $global_norm_sample_size, $tn_dist_file, $tn_log_file);

my $cp_str = "cp $global_TB_processing_dir/rare_tables/rarefaction_$global_norm_sample_size"."_$centroid_index".".txt $tn_otu_table_file";
`$cp_str`;

print "Summarizing by taxa.....\n";
`summarize_taxa.py -i $tn_otu_table_file -o $global_TB_results_dir/breakdown_by_taxonomy/`;

# move 100 of the 100 tables just produced into a new folder for jacknifing
my $jack_knife_folder = "$global_TB_processing_dir/rare_tables/JN";
`mkdir -p $jack_knife_folder`;
foreach my $jn_file_counter (0..$global_JN_file_count)
{
    my $jn_from_file = "rarefaction_".$global_norm_sample_size."_".$jn_file_counter.".txt";
    `cp $global_TB_processing_dir/rare_tables/$jn_from_file $jack_knife_folder/`;
}

print "Jacknifed beta diversity....\n";
jack_knifing("weighted_unifrac", $nn_otu_table_file, $nn_tree_file, $jack_knife_folder);
jack_knifing("unweighted_unifrac", $nn_otu_table_file, $nn_tree_file, $jack_knife_folder);
jack_knifing("euclidean", $nn_otu_table_file, $nn_tree_file, $jack_knife_folder);
jack_knifing("hellinger", $nn_otu_table_file, $nn_tree_file, $jack_knife_folder);

# beta for the normalised table
my @beta_methods = ('weighted_unifrac','unweighted_unifrac','euclidean','hellinger','binary_euclidean','chord');
foreach my $matrix_type (@beta_methods)
{
    my $beta_str = "beta_diversity.py -i $tn_otu_table_file -o $global_TB_results_dir/beta_diversity/$tn_prefix/$matrix_type -m $matrix_type";
    # qiime spews when you give it a tree for some methods
    if(!($matrix_type =~ /euclid/))
    {
        $beta_str .= " -t $nn_tree_file";
    }
    `$beta_str`;
    
    my $beta_file = $matrix_type."_".$tn_prefix."_otu_table.txt";
    my $upgma_file = $matrix_type."_".$tn_prefix."_otu_table_upgma.tre";
    my $pcoa_file = $matrix_type."_".$tn_prefix."_pcoa.txt";

    # Perform UPGMA clustering on rarefied distance matrices 
    `upgma_cluster.py -i $global_TB_results_dir/beta_diversity/$tn_prefix/$matrix_type/$beta_file -o $global_TB_results_dir/beta_diversity/$tn_prefix/$matrix_type/$upgma_file`;
    
    # Compute principal coordinates
    `principal_coordinates.py -i $global_TB_results_dir/beta_diversity/$tn_prefix/$matrix_type/$beta_file -o $global_TB_results_dir/beta_diversity/$tn_prefix/$matrix_type/$pcoa_file`;
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
    find_centroid_sequences($global_SB_processing_dir."/norm_tables/", 7, "$global_TB_processing_dir/$nn_fasta_file", "$global_SB_processing_dir/$sn_fasta_file", $sn_dist_file, $sn_log_file);
    
    # once we've found our guy, lets make an otu table with assigned taxonomies   
    print "Picking OTUs for SEQ normalised data set...\n";
    `pick_otus.py -i "$global_SB_processing_dir/$sn_fasta_file" -s $global_similarity_setting -o $global_SB_processing_dir/uclust_picked_otus`;

    print "Gettting a representitive set...\n";
    `pick_rep_set.py -i $global_SB_processing_dir/uclust_picked_otus/$sn_otus_file -f "$global_SB_processing_dir/$sn_fasta_file"`;

    # if we are doing OTU_AVERAGE (or if we've ben asked to) then we need to assign taxonomy here
    print "Assigning taxonomy for SEQ normalised data set...\n";
    my $sn_rep_set_fasta = "$global_SB_processing_dir/".$sn_fasta_file."_rep_set.fasta";
    `assign_taxonomy.py -i $sn_rep_set_fasta -t $TAX_tax_file -b $TAX_blast_file -m blast -e $global_e_value -o $global_SB_processing_dir`;
    
    print "Treeing SEQ normalised data set...\n";
    `align_seqs.py -i $sn_rep_set_fasta -t $imputed_file -p 0.6 -o $global_SB_processing_dir/pynast_aligned`;
    my $sn_rep_set_aligned = "$global_SB_processing_dir/pynast_aligned/".$sn_fasta_file."_rep_set_aligned.fasta";
    `filter_alignment.py -i $sn_rep_set_aligned -o $global_SB_processing_dir`;
    my $sn_rep_set_aligned_filtered = "$global_SB_processing_dir/".$sn_fasta_file."_rep_set_aligned_pfiltered.fasta";
    `make_phylogeny.py -i $sn_rep_set_aligned_filtered -r midpoint`; 
    my $sn_rep_set_aligned_filtered_tree = "$global_SB_processing_dir/".$sn_fasta_file."_rep_set_aligned_pfiltered.tre";    
    `mv $sn_rep_set_aligned_filtered_tree $sn_tree_file`;
    
    print "Making SEQ SEQ normalised otu table...\n";
    my $sn_rep_set_tax_assign = "$global_SB_processing_dir/".$sn_fasta_file."_rep_set_tax_assignments.txt";
    `make_otu_table.py -i $global_SB_processing_dir/uclust_picked_otus/$sn_otus_file -t $sn_rep_set_tax_assign -o $sn_otu_table_file`;

    print "Alpha and beta diversity for SEQ normalised otu table...\n";
    # beta
    my @beta_methods = ('weighted_unifrac','unweighted_unifrac','euclidean','hellinger','binary_euclidean','chord');
    foreach my $matrix_type (@beta_methods)
    {
        # basic beta diversity call
        my $beta_str = "beta_diversity.py -i $sn_otu_table_file -o $global_SB_results_dir/beta_diversity/$matrix_type -m $matrix_type";
        # qiime spews when you give it a tree fro some methods
        if(!($matrix_type =~ /euclid/))
        {
            $beta_str .= " -t $sn_tree_file";
        }
        `$beta_str`;
        
        my $beta_file = $matrix_type."_".$sn_prefix."_otu_table.txt";
        my $upgma_file = $matrix_type."_".$sn_prefix."_otu_table_upgma.tre";
        my $pcoa_file = $matrix_type."_".$sn_prefix."_pcoa.txt";
        # Perform UPGMA clustering on rarefied distance matrices 
        `upgma_cluster.py -i $global_SB_results_dir/beta_diversity/$matrix_type/$beta_file -o $global_SB_results_dir/beta_diversity/$matrix_type/$upgma_file`;
        
        # Compute principal coordinates
        `principal_coordinates.py -i $global_SB_results_dir/beta_diversity/$matrix_type/$beta_file -o $global_SB_results_dir/beta_diversity/$matrix_type/$pcoa_file`;
        
    }

    # alpha
    `alpha_diversity.py -i $sn_otu_table_file -t $sn_tree_file -o $global_SB_results_dir/alpha_diversity.txt -m chao1,chao1_confidence,PD_whole_tree,observed_species,simpson,shannon,fisher_alpha`;

    print "Summarizing by taxa.....\n";
    `summarize_taxa.py -i $sn_otu_table_file -o $global_SB_results_dir/breakdown_by_taxonomy/`;
}


print "Results are located in: $global_working_dir/results/\n";

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
    my ($matrix_type, $raw_otu_table, $raw_tree, $jack_knife_folder ) = @_;
    
    # Produce distance matrix reflecting beta diversity in non-normalised OTU table
    `beta_diversity.py -i $raw_otu_table -o $global_TB_results_dir/beta_diversity/$nn_prefix/$matrix_type -m $matrix_type -t $raw_tree`;
        
    # Perform UPGMA clustering on non-normalised distance matrix
    my $beta_otu_table = $matrix_type."_".$nn_prefix."_otu_table.txt";
    my $upgma_cluster_tree = $matrix_type."_".$nn_prefix."_otu_table_upgma.tre";

    `upgma_cluster.py -i $global_TB_results_dir/beta_diversity/$nn_prefix/$matrix_type/$beta_otu_table -o $global_TB_results_dir/beta_diversity/$nn_prefix/$matrix_type/$upgma_cluster_tree`;
    
    # Produce distance matrices reflecting beta diversity in the rarefied OTU tables
    `beta_diversity.py -i $jack_knife_folder -o $global_TB_results_dir/beta_diversity/$nn_prefix/$matrix_type/rare_dm/ -m $matrix_type -t $raw_tree`;
    
    # Perform UPGMA clustering on rarefied distance matrices 
    `upgma_cluster.py -i $global_TB_results_dir/beta_diversity/$nn_prefix/$matrix_type/rare_dm/ -o $global_TB_results_dir/beta_diversity/$nn_prefix/$matrix_type/rare_upgma/`;
    
    # Compare the consensus tree to the beta-derived trees
    `tree_compare.py -s $global_TB_results_dir/beta_diversity/$nn_prefix/$matrix_type/rare_upgma/ -m $global_TB_results_dir/beta_diversity/$nn_prefix/$matrix_type/$upgma_cluster_tree -o $global_TB_results_dir/beta_diversity/$nn_prefix/$matrix_type/upgma_cmp/`;
    
    # Compute principal coordinates
    `principal_coordinates.py -i $global_TB_results_dir/beta_diversity/$nn_prefix/$matrix_type/rare_dm/ -o $global_TB_results_dir/beta_diversity/$nn_prefix/$matrix_type/pcoa/`;
        
    # Make PDF of Jackknife tree with labeled support: weighted unifrac command
    my $output_pdf = "$global_TB_results_dir/$matrix_type"."_betadiv_jackknife_tree.pdf"; 
    `make_bootstrapped_tree.py -m $global_TB_results_dir/beta_diversity/$nn_prefix/$matrix_type/upgma_cmp/master_tree.tre -s $global_TB_results_dir/beta_diversity/$nn_prefix/$matrix_type/upgma_cmp/jackknife_support.txt -o $output_pdf`;

}

sub run_R_cmd
{
    #-----
    # Wrapper for running R commands
    #
    my ($cmd) = @_;
    #print "$cmd\n";
    $global_R_instance->run($cmd);
}

sub find_centroid_table
{
    #-----
    # Find a representative set of $global_norm_sample_size sequences (do this in RRRRRR...)
    #
    my($path, $norm_size, $dist_file, $log_file) = @_;
    
    print "Calculating centroid OTU table from tables in $path...\n";
    
    $global_R_instance->start();

    my $sampl_p1 = $global_num_samples + 1;

    # read in the list of distance matricies
    run_R_cmd(qq`library(foreign);`);
    run_R_cmd(qq`a<-list.files("$path", "*.txt");`);
    
    # work out how many there are and allocate an array
    run_R_cmd(qq`len_a <- length(a);`);
    run_R_cmd(qq`big_frame <- array(0,dim=c($global_num_samples,$global_num_samples,len_a));`);
    
    print "  --start loading data...\n";
    
    # load each file individually into a big frame
    my $r_str = "for (i in c(1:len_a)) { j <- i - 1; name <- paste(\"$path\",\"rarefaction_$norm_size\",\"_\",j,\".txt\",sep=\"\"); u<-read.table(name,sep=\"\\t\",row.names=1); u[,$sampl_p1]<-NULL; big_frame[,,i]<-as.matrix(dist(t(u), upper=TRUE, diag=TRUE)); i<-i+1; }";
    run_R_cmd($r_str);    

    print "  --data loaded, calculating centroid...\n";
    
    # find the average matrix
    run_R_cmd(qq`ave <- big_frame[,,1];`);
    run_R_cmd(qq`for (i in c(2:len_a)) { ave <- ave + big_frame[,,i]; }`);
    run_R_cmd(qq`ave <- ave/len_a;`);
    
    print "  --calculating didtances of tables to centroid...\n";
    
    # find the euclidean distance of each matrix from the average
    run_R_cmd(qq`dist<-array(0,dim=c(len_a));`);
    run_R_cmd(qq`for (i in c(1:len_a)) { dist[i] <- sqrt(sum(big_frame[,,i]-ave)^2); }`);
    
    # find the min value
    run_R_cmd(qq`min_index <- which.min(dist);`);
    my $centroid_otu_index = $global_R_instance->get('min_index');
    # R indexes from 0
    $centroid_otu_index--;
    print "  --table: $centroid_otu_index is the centroid table\n";
    
    # make stats on the distances
    # and log what we did
    open my $log_fh, ">", $log_file or die "Could not open log file: $log_file : $!\n";
    run_R_cmd(qq`max_dist <- max(dist);`);
    run_R_cmd(qq`min_dist <- min(dist);`);
    run_R_cmd(qq`range_dist <- max_dist - min_dist;`);
    run_R_cmd(qq`mean_dist <- mean(dist);`);
    run_R_cmd(qq`median_dist <- median(dist);`);
    
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
        run_R_cmd(qq`library(permute);`);
        run_R_cmd(qq`library(vegan);`);
        run_R_cmd(qq`mantel.otu <- mantel(ave,big_frame[,,min_index]);`);
        run_R_cmd(qq`m_stat <- mantel.otu\$statistic;`);
        run_R_cmd(qq`m_sig <- mantel.otu\$signif;`);
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
        `app_normalise_reads.pl -in $seqs -n $global_norm_sample_size -o $norm_seq_file`; 
    
        # make an otu table for each guy
        `pick_otus.py -i $norm_seq_file -s $global_similarity_setting -o $path`;
        `make_otu_table.py -i $norm_otus_file -o $norm_otu_table_file`;    
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
    open my $mapping, ">", $global_mapping_file or die $!;
    print $mapping "$FNB_HEADER\n";

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
}

######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    my @standard_options = ( "help|h+", "config|c:s", "identity|i:i", "e:i");
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
    if(exists $options{'identity'})
    {
        if(($options{'identity'} <= 0) || ($options{'identity'} > 1))
        {
            die "Identity must be an integer greater than 0 and  no greater than 1\n";
        }
    }
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

    app_make_results.pl

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

   Insert detailed description here

=head1 SYNOPSIS

    app_make_results.pl -c|config CONFIG_FILE [-help|h]

      -c CONFIG_FILE               App config file to be processed
      [-i identity VALUE]          Set blast identity [default: 97%]
      [-e EVALUE]                  Set e-value for blast (assign taxonomy) [default 0.001]
      [-help -h]                   Displays basic usage information
         
=cut


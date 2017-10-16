#!/bin/perl

#########################################################################################
#                                                                                       #
#    Input: List of Entrez UIDs (integer identifiers, e.g. PMID, GI, Gene ID)           #
#    EFetch Output: Formatted data records (e.g. abstracts, FASTA                       #
#                                                                                       #
#########################################################################################


use LWP::Simple;
use Term::ProgressBar;
use Sys::CPU;
use Parallel::ForkManager;
use strict;
use warnings;
use diagnostics;
use Data::Dumper;


#Check required arguments
if (scalar(@ARGV) != 2)
{
    print "Usage: perl ncbiAccessionsFromTaxid.pl <query> <accession.list>\n\n";
    #perl ncbiAccessionsFromTaxid.pl "txid1639[Organism:exp]" ~/Desktop/listeria_mono.list
    print "Example: perl ncbiAccessionsFromTaxid.pl \"txid1639[Organism:exp]\" ~/Desktop/listeria_mono.list";
    exit;
}

#print(Dumper \@ARGV);

#Command line arguments
my $query = $ARGV[0];
my $outfile = $ARGV[1];

#select database
my $db = "protein";

#assemble the esearch URL
my $base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
my $url = $base . "esearch.fcgi?db=$db&term=$query&usehistory=y";

# https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term=txid1639[Organism:exp]&usehistory=y
#post the esearch URL
my $output = get($url);
die "Could not complete query search on NCBI server." unless $output;

#parse WebEnv, QueryKey and Count (# records retrieved)
my $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
my $key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
my $count = $1 if ($output =~ /<Count>(\d+)<\/Count>/);

print "There is $count $db accessions for $query";

#open output file for writing
open(OUT, ">>", $outfile) || die "Can't open file!\n";

#progress bar
print "\n\nDownloading list of accessions\n\n";
#my $progress = Term::ProgressBar->new($count);

# https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&WebEnv=NCID_1_2591324_130.14.22.215_9001_1469634077_1008305355_0MetA0_S_MegaStore_F_1&query_key=1&&retmax=500&&rettype=acc
#create all the URLs to retrieve and store them in an array
#retrieve data in batches of 200
my $retmax = 500;
my @URLs;
for (my $retstart = 0; $retstart < $count; $retstart += $retmax) {
        my $efetch_url = $base ."efetch.fcgi?db=$db&WebEnv=$web";
        $efetch_url .= "&query_key=$key&retstart=$retstart";
        $efetch_url .= "&retmax=$retmax&rettype=acc&retmode=text";
        #my $efetch_out = get($efetch_url);
        #print OUT "$efetch_out";
        
        push(@URLs, $efetch_url);
        
        #update progress bar
        #$progress->update($_);
}

#Use parallel loop to download multple URLs at the same time

#Prepare for parallel loop
my $number_of_cpus = Sys::CPU::cpu_count();
my $MAX_PROCESSES = $number_of_cpus;
my $pm = Parallel::ForkManager->new($MAX_PROCESSES);

foreach my $link (@URLs)
{
    $pm-> start and next; # Forks and returns the pid for the child
        
    my $efetch_out = get($link);
    if ($efetch_out) { print OUT $efetch_out; }
    else { print "Could not download accessions from $link\n"; }
    
    $pm->finish; # terminate the child process
}
$pm->wait_all_children;

close OUT;

print ("\n\ndone\n");

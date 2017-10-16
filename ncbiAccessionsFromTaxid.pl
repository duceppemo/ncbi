#!/bin/perl

#########################################################################################
#                                                                                       #
#    Input: List of Entrez UIDs (integer identifiers, e.g. PMID, GI, Gene ID)           #
#    EFetch Output: Formatted data records (e.g. abstracts, FASTA                       #
#                                                                                       #
#########################################################################################


use LWP::Simple;
use Term::ProgressBar;

use strict;
use warnings;
use diagnostics;
use Data::Dumper;

use Getopt::Long;
use Pod::Usage;


####################
#                  #
#   Option flags   #
#                  #
####################


# Default values
my $db;
my $query;
my $outfile;
my $help;
my $man;


GetOptions(
    'db|d=s' => \$db,  # string
    'query|q=s' => \$query,  # string
    'output|o=s' => \$outfile,  # string
    'help|h' => \$help,
    'man' => \$man) or pod2usage(2);

pod2usage(1) if $help;
pod2usage(-exitval => 0, -verbose => 2) if $man;


##############
#            #
#   Checks   #
#            #
##############


#check if arguments are not empty
if (!defined($db) || !defined($query) || !defined($outfile))
{
    pod2usage(1)
}


############
#          #
#   Main   #
#          #
############


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
open(OUT, ">", $outfile) || die "Can't open file!\n";

#progress bar
print "\n\nDownloading list of accessions\n\n";

# https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&WebEnv=NCID_1_2591324_130.14.22.215_9001_1469634077_1008305355_0MetA0_S_MegaStore_F_1&query_key=1&&retmax=500&&rettype=acc
#create all the URLs to retrieve and store them in an array
#retrieve data in batches of 200
my $retmax = 500;
my $progress = Term::ProgressBar->new($count/$retmax);
for (my $retstart = 0; $retstart < $count; $retstart += $retmax) {
        my $efetch_url = $base ."efetch.fcgi?db=$db&WebEnv=$web";
        $efetch_url .= "&query_key=$key&retstart=$retstart";
        $efetch_url .= "&retmax=$retmax&rettype=acc&retmode=text";
        my $efetch_out = get($efetch_url);
        print OUT "$efetch_out";
        
        #update progress bar
        $progress->update($_);
}

close OUT;

print ("\ndone\n");


__END__
=head1 NAME

ncbiAccessionsFromTaxid.pl - Get the accession numbers related to the result of the text query

=head1 SYNOPSIS

perl ncbiAccessionsFromTaxid.pl --db <protein|nucleotide> --query 'my query' --output <accession.list>

=head1 OPTIONS

=over 4

=item B<--help [-h]>

Print this help

=item B<--man>

Print script manual

=item B<--db [-d]>

NCBI database. Allowed values are 'protein' or 'nucleotide'

=item B<--query [-q]>

String to look for. E.g. 'txid33090[Organism] AND "internal transcribed spacer"'

=item B<--output [-o]>

Output file. This is the list of accession numbers.

=back

=head1 DESCRIPTION

B<ncbiAccessionsFromTaxid.pl> will search a query string on NCBI and return a list of corresponding nucleotide or protein accession numbers,
depending which database is selected.

=cut

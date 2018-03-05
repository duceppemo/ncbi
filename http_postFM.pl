#!/bin/perl


#########################################################################################
#                                           #
#   Input: List of Entrez UIDs (integer identifiers, e.g. PMID, GI, Gene ID)    #
#   EFetch Output: Formatted data records (e.g. abstracts, FASTA)           #
#                                           #
#########################################################################################


# TODO ->  add option for format


use LWP::Simple;
use strict;
use warnings;
use diagnostics;
use Parallel::ForkManager;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use Sys::CPU;
use File::Path qw(make_path);

my $cpu = Sys::CPU::cpu_count();


####################
#                  #
#   Option flags   #
#                  #
####################

# https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly

# Default values
my $input;
my $output;
my $type;
my $split;
my $nproc = $cpu;  # all cpu
my $retmax = 300;  # size of batch to download at once
my $help;
my $man;

GetOptions(
    'input|i=s' => \$input,  # string
    'output|o=s' => \$output,  # string
    'type|t=s' => \$type,  # string
    'split|s' => \$split, # flag
    'nproc|n=i' => \$nproc,  # integer
    'block|b=i' => \$retmax,
    'help|h' => \$help,
    'man' => \$man) or pod2usage(2);

pod2usage(1) if $help;
pod2usage(-exitval => 0, -verbose => 2) if $man;


##############
#            #
#   Checks   #
#            #
##############


#check if mandatory arguments are not empty
if ( !(defined($input)) or !(defined($output)) or !(defined($type)) )
{
    print("Missing madatory arguments!\n");
    pod2usage(1);
}

#check if database type is good
if ($type ne ("nucleotide" || "protein"))
{
    print "Sequence type must be either 'protein' or 'nucleotide'\n";
    pod2usage(1);
}

# check is correct ouptut
if ($split)
{
    eval { make_path($output) };
    if ($@) {
      print "Couldn't create $output: $@";
    }
    # opendir(my $dir_fh, $output) or die "Output directory does not exist: $!\n";
    # closedir($dir_fh);
}

#Set $db
my $db;
if ($type eq "protein") {$db = $type;} else {$db = "nuccore";}


##################
#                #
#   Parse list   #
#                #
##################


#Open the input file, in read-only mode
open(my $infilehandle, "<", $input) or die "Error opening $input: $!\n";
my @query;

if ($split)
{
    while(my $line = <$infilehandle>)
    {
        chomp($line);
        $line =~ s/^\s+|\s+$//g;  # trim remove white space from both ends of a string
        next unless length($line);  # skip empty lines or lines with only whitespace characters
        
        push(@query, $line);
    }
}
else
{
    my $counter = 0;
    my $string;

   #Retrieve data in batches of $retmax
    #create array with file
    while(my $line = <$infilehandle>)
    {
        chomp($line);
        $line =~ s/^\s+|\s+$//g;  # trim remove white space from both ends of a string
        next unless length($line);  # skip empty lines or lines with only whitespace characters

        $counter += 1;

        if ($counter < $retmax)
        {
            $string .= $line . ",";
        }
        else
        {
            $string .= $line;
            push(@query, $string);
            $counter = 0;
            $string = "";
        }

    }
    #remainig IDs
    #delete the last comma that has no ID following
    chop($string);
    push(@query, $string);
}

# print Dumper(\@query);  # Debug

#Close input file
close ($infilehandle);


####################
#                  #
#   NCBI e-utils   #
#                  #
####################


print "Downloading using $nproc slots:\n";

my $pm = new Parallel::ForkManager($nproc);

#Retrieve data in batches of $retmax
my $size = @query;
my $i = 0;

my $base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi";

# Check for split option
if($split)
{
    foreach my $entry (@query)
    {
        select(undef, undef, undef, 0.3); # sleep 0.3s
        # sleep(1);  # max 3 requests per seconds for ncbi
        $i += 1;
        print "\r$i/$size";  # Display which block of 500 sequences is being downloaded

        my $pid = $pm->start and next;

        #open output file for writing
        open(my $output_fh, ">", "$output/$entry.fasta") || die "Can't write to $output/$entry.fasta: $!\n";

        # my $efetch_url = $base_url . "?db=$db&id=$entry&rettype=gp&retmode=text";  # GenPept
        my $efetch_url = $base_url . "?db=$db&id=$entry&rettype=fasta&retmode=text";  # fasta
        # my $efetch_url = $base_url . "?db=$db&id=$entry";#&rettype=fasta&retmode=text";  # ASN.1

        my $efetch_out = get($efetch_url);
        die "Could not get $efetch_url" unless defined $efetch_out;
        print {$output_fh} ("$efetch_out");
        close ($output_fh);

        $pm->finish; # Terminates the child process
    }
}
else
{
    #open output file for writing
    open(my $output_fh, ">>", $output) || die "Can't write to $output: $!\n";

    foreach my $entry (@query)
    {
        select(undef, undef, undef, 0.3); # sleep 0.3s
        # sleep(1);  # max 3 requests per seconds for ncbi
        $i += 1;
        print "\r$i/$size";  # Display which block of $retmax sequences is being downloaded
        my $pid = $pm->start and next;

        # my $efetch_url = $base_url . "?db=$db&id=$entry&rettype=gp&retmode=text";  # GenPept
        my $efetch_url = $base_url . "?db=$db&id=$entry&rettype=fasta&retmode=text";  # fasta
        # my $efetch_url = $base_url . "?db=$db&id=$entry";#&rettype=fasta&retmode=text";  # ASN.1

        my $efetch_out = get($efetch_url);
        die "Could not get $efetch_url" unless defined $efetch_out;
        print {$output_fh} ("$efetch_out");

        $pm->finish; # Terminates the child process
    }
    close ($output_fh);
}

print ("\nDone!\n");


__END__
=head1 NAME

http_postFM.pl - Download protein or nucleotide sequences from genbank from a list of IDs.

=head1 SYNOPSIS

perl http_postFM.pl [options] -i proteinID.list -o out.fasta -t protein

=head1 OPTIONS

=over 4

=item B<--help [-h]>

Print this help

=item B<--man>

Print script manual

=item B<--input [-i]>

List of Entrez IDs (e.g. Accession numbers, Gene IDs, etc.). Mandatory

=item B<--output [-o]>

Output fasta file in normal mode.
Output folder in split mode. Mandatory

=item B<--type [-t]>

Type of sequences in the list. "protein" or "nucleotide".
Only one type at the time is supported. Mandatory

=item B<--split [-s]>

Write an individual file for each accession number.
Default writes all sequences in a single file. Optional

=item B<--block [-b]>

Size of block (number of sequences).
Default is 300. Optional

=item B<--nproc [-n]>

Maximum number of block of sequences to downloaded in parallel.
A new block is downloaded every 300 millisecond to respect NCBI's terms.
Default is equal to number of CPUs available. Optional

=back

=head1 DESCRIPTION

B<http_postFM.pl> will parallel download blocks of 200 sequences in fasta format from a list of Entrez UIDs (e.g. accession numbers).
There is a 0.3 second delay before each block starts to follow NCBI's rule of maximum three requests per second.

The ID list must contain one ID per line. You have to make sure that the sequence type matches the ID type in the list.
Only one sequence type is allowed, i.e. nucleotide or protein.

Perl dependencies: Parallel::ForkManager, Sys::CPU    

=cut

#!/usr/bin/perl -w

=head1 NAME

simulate_reads

=head1 SYNOPSIS

	 perl simulate_reads.pl [--help] [--verbose] [--include] 
                            [--readlen=INT]
                            --fwd=STRING --rev=STRING 
                            --out=STRING --db=STRING

	 Outputs only the part of a sequence found between the given primers 
	
        --include: keep primer sequences in the output; default, remove
        --readlen: length of the reads to simmulate. Default 0, extract full amplicon
		--fwd, rev: primers to generate amplicons
        --out: prefix for the output file, will be prefixed as .fasta if length = 0
                        and as _r1.fasta,_r2.fasta if length>0
		--db: fasta file (gapped or not); output will be unaligned
		--help: This info.
	    

=head1 AUTHOR

luisa.hugerth@scilifelab.se

=cut

use warnings;
use strict;

use Bio::SeqIO;
use Getopt::Long;
use Pod::Usage;

my $help = 0;
my $verbose = 0;
my $include = 0;
my $readlen = 0;
my $fwd;
my $rev;
my $outputstring;
my $dbfile;

GetOptions(
  "help!" => \$help,
  "verbose!" => \$verbose,
  "fwd=s" => \$fwd,
  "rev=s" => \$rev,
  "include!" => \$include,
  "readlen=i" => \$readlen,
  "out=s" => \$outputstring,
  "db=s" => \$dbfile
);

pod2usage(0) if $help;

pod2usage(-msg => "Need a dbfile", -exitval => 1) unless $dbfile;
pod2usage(-msg => "Please provide an output path", -exitval => 1) unless $outputstring;
pod2usage(-msg => "Read length must be positive", -exitval => 1) if $readlen < 0;
my $infile = Bio::SeqIO->new(-file=>"$dbfile", -format => "fasta") or die "Can't open $dbfile";

sub regenerate{
	my $seq=$_[0];
	$seq=~s/U/T/g;
	$seq=~s/R/[AG]/g;
	$seq=~s/Y/[CT]/g;
	$seq=~s/S/[GC]/g;
	$seq=~s/W/[AT]/g;
	$seq=~s/K/[GT]/g;
	$seq=~s/M/[AC]/g;
	$seq=~s/B/[CGT]/g;
	$seq=~s/D/[AGT]/g;
	$seq=~s/H/[ACT]/g;
	$seq=~s/V/[ACG]/g;
	$seq=~s/N/[ACTG]/g;
	return $seq;
}

$fwd = regenerate($fwd);
$rev = regenerate($rev);

my $fwdlen = $readlen;
my $revlen = $readlen;
my $ff;
my $rf;
if ($readlen > 0){
    $fwdlen = $readlen - length($fwd);
    $revlen = $readlen - length($rev);
    my $ffile = $outputstring."_r1.fasta";
    my $rfile = $outputstring."_r2.fasta";
    open ($ff, ">", $ffile) or die "Could not open file '$ffile' $!";
    open ($rf, ">", $rfile) or die "Could not open file '$rfile' $!";
}
else{
    my $filename = $outputstring.".fasta";
    open ($ff, ">", $filename) or die "Could not open file '$filename' $!";
}    

while (my $record = $infile->next_seq){
	my $seq = $record->seq;
	$seq=~s/_//g;
	$seq=~s/U/T/g;
	if ($seq=~/$fwd/ and $seq=~/$rev/){
        if ($readlen == 0){
            my $amplicon;
    		if ($include){
	    		$seq=~/($fwd.+$rev)/;
		    	$amplicon = $1;
    		}
	    	 else{
		    	$seq=~/$fwd(.+)$rev/;
			    $amplicon = $1;
    		 }
	    	print $ff ">".$record->display_id."\n";
		    print $ff "$amplicon\n";
        }
        else{
            my $fwd_read;
            my $rev_read;
            if ($include){
                $seq=~/($fwd.+$rev)/;
                $fwd_read = substr($1, 0, $readlen);
                $rev_read = substr($1, -$readlen);
            }
             else{
                $seq=~/$fwd(.+)$rev/;
                $fwd_read = substr($1, 0, $fwdlen);
                $rev_read = substr($1, -$revlen);
             }
            print $ff ">".$record->display_id."\n";
            print $rf ">".$record->display_id."\n";
            print $ff "$fwd_read\n";
            print $rf "$rev_read\n";
        }
	}
}



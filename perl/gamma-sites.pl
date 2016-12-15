#!/usr/bin/env perl

use warnings;
use strict;

use File::Slurp;
use JSON;
use Math::CDF qw(qgamma);
use Clone qw(clone);

use Data::Dumper;

BEGIN {
    use FindBin;
    use lib $FindBin::Bin;
    push @INC, $FindBin::Bin;
}

use Getopt::Long;
use File::Basename;

use SequenceIterator qw(iterseq printseq revcomp);
use Stockholm;

my ($progname) = fileparse($0);

my ($shape, $scale, $bins) = (1, 1, 2);
my ($pretty, $verbose) = (0, 0);

my $usage = "";
$usage .= "$progname -- add discretized-gamma prior to historian rate matrix\n";
$usage .= "\n";
$usage .= "Usage: $progname [opts] params.json\n";
$usage .= "\n";
$usage .= "Options:\n";
$usage .= " -shape <shape>  Set shape parameter (default $shape)\n";
$usage .= " -scale <scale>  Set scale parameter (default $scale)\n";
$usage .= " -bins <bins>    Set number of bins (default $bins)\n";
$usage .= " -pretty         Pretty-print JSON output\n";
$usage .= " -verbose        Print debugging info\n";
$usage .= "\n";

GetOptions ("shape=f" => \$shape,
	    "scale=f" => \$scale,
	    "bins=i"  => \$bins,
	    "pretty" => \$pretty,
	    "verbose" => \$verbose)
    or die $usage;

@ARGV == 1 or die $usage;
my ($paramsFilename) = @ARGV;

my @quartiles = map (qgamma ($_ / ($bins+1), $shape, $scale), 1..$bins);
warn "Quartiles: (@quartiles)\n" if $verbose;

my $paramsFile = read_file ($paramsFilename);
my $params = decode_json ($paramsFile);

my $mixture = clone ($params);
delete $mixture->{'subrate'};
delete $mixture->{'rootprob'};

my @cpt = map ({ 'subrate' => clone ($params->{'subrate'}),
		 'rootprob' => clone ($params->{'rootprob'}),
		 'weight' => 1/$bins }, 1..$bins);

for my $bin (0..$bins-1) {
    while (my ($src, $destHash) = each %{$cpt[$bin]->{'subrate'}}) {
	for my $dest (keys %$destHash) {
	    $destHash->{$dest} *= $quartiles[$bin];
	}
    }
}

$mixture->{'mixture'} = \@cpt;

print to_json ($mixture, { 'pretty' => $pretty });

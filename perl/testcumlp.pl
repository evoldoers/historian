#!/usr/bin/env perl

use warnings;

die "Usage: $0 <file> <nCells>" unless @ARGV == 2;

my $file = $ARGV[0];
my $nCells = $ARGV[1];

local *FILE;
open FILE, "<$file" or die $!;

my $clp;
my $line = 0;
my $n = 0;
while (<FILE>) {
    ++$line;
    if (/\"cumLogProb\": \"([^\"]+)\"/) {
	$clp = $1;
    } else {
	if (/\"fwdLogProb\": \"([^\"]+)\"/) {
	    my $flp = $1;
	    if (defined($clp)) {
		if ($clp eq $flp) {
		    ++$n;
		} else {
		    print "not ok: on lines ", $line-1, "-$line of $file: cumLogProb = $clp, fwdLogProb = $flp\n";
		    exit;
		}
	    }
	}
	$clp = undef;
    }
}

close FILE;

if ($n != $nCells) {
    print "not ok: in $file, found $n adjacent cumLogProb/fwdLogProb lines; expected $nCells\n";
    exit;
}

print "ok: in $file, $n cumLogProb's match succeeding fwdLogProb's\n";

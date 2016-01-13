#!/usr/bin/env perl

use warnings;
use File::Temp;

die "Usage: $0 <prog> <args...> <expected>" unless @ARGV >= 3;
my $expected = pop @ARGV;
my ($prog, @args) = @ARGV;

my $fh = File::Temp->new();
my $fname = $fh->filename;

system "$prog @args >$fname";
my $diff = `diff $fname $expected`;

if (length $diff) {
    print "`$prog @args` does not match $expected:\n";
    print `diff -y $fname $expected`;
    print "not ok\n";
} else {
    print "ok: `$prog @args` matches $expected\n";
}

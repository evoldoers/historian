#!/usr/bin/env perl -w

use JSON;

my (@alph, $initial, $mutate, @cpt);
while (<>) {
    if (/\(chain/) {
	if (@alph) {
	    push @cpt, { 'subrate' => $mutate, 'rootprob' => $initial };
	}
	@alph = ();
	$initial = {};
	$mutate = {};
    } elsif (/initial..state..(.)....prob.([0-9\.e\+\-]+)/) {
	$initial->{$1} = $2 + 0;
	push @alph, $1;
    } elsif (/mutate..from..(.)....to..(.)....rate.([0-9\.e\+\-]+)/) {
	$mutate->{$1}->{$2} = $3 + 0;
    }
}

my $model = { 'alphabet' => join("",@alph),
	      'mixture' => \@cpt,
	      'insrate' => .01,
	      'delrate' => .01,
	      'insextprob' => .66,
	      'delextprob' => .66 };

print to_json ($model, { 'pretty' => 1 });

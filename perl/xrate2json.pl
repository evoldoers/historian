#!/usr/bin/env perl -w

my (@alph, %initial, %mutate);
while (<>) {
    if (/initial..state..(.)....prob.([0-9\.e\+\-]+)/) {
	$initial{$1} = $2;
	push @alph, $1;
    } elsif (/mutate..from..(.)....to..(.)....rate.([0-9\.e\+\-]+)/) {
	$mutate{$1}->{$2} = $3;
    }
}

print
    "{\n",
    "  \"alphabet\": \"", join("",@alph), "\",\n",
    "  \"rootprob\":\n",
    "  {\n",
    join (",\n", map ("    \"$_\": $initial{$_}", @alph)), "\n",
    "  },\n",
    "  \"subrate\":\n",
    "  {\n";
for my $i (@alph) {
    print "    \"$i\": { ", join (", ", map ("\"$_\": ".$mutate{$i}->{$_}, grep(exists $mutate{$i}->{$_},@alph))), " }", ($i eq $alph[$#alph] ? "" : ","), "\n";
}
print
    "  },\n",
    "  \"insrate\": 0.01,\n",
    "  \"delrate\": 0.01,\n",
    "  \"insextprob\": 0.66,\n",
    "  \"delextprob\": 0.66\n",
    "}\n";

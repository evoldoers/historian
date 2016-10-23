#!/usr/bin/env perl -w

die "Usage: $0 <model>" unless @ARGV == 1;
my ($model) = @ARGV;
my $ucmodel = uc($model);

open HDR, ">src/${model}.h";
print HDR map ("$_\n",
	       "#ifndef ${ucmodel}_MODEL_INCLUDED",
	       "#define ${ucmodel}_MODEL_INCLUDED",
	       "",
	       "#include \"model.h\"",
	       "",
	       "RateModel ${model}Model();",
	       "",
	       "extern const char* ${model}ModelText;",
	       "",
	       "#endif /* ${ucmodel}_MODEL_INCLUDED */");
close HDR;

open CPP, ">src/${model}.cpp";
print CPP map ("$_\n",
	       "#include \"${model}.h\"",
	       "#include \"jsonutil.h\"",
	       "",
	       "RateModel ${model}Model() {",
	       "  RateModel m;",
	       "  ParsedJson pj (${model}ModelText);",
	       "  m.read (pj.value);",
	       "  return m;",
	       "}",
	       "",
	       "const char* ${model}ModelText =");

open JSON, "model/${model}.json" or die "model/${model}.json: $!";
my $q = chr(34);
while (<JSON>) {
    chomp;
    s/$q/\\$q/g;
    print CPP chr(34), $_, "\\n", chr(34), "\n";
}
print CPP ";\n";

close CPP;

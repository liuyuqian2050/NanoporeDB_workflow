use strict;
use warnings;

my $group=$ARGV[0];
my $path=$ARGV[1];
if (-e "$path/Group") {
	`rm -rf $path/Group/`;
}
`mkdir $path/Group`;

my @g=();
open IN,"3_nanopore_allminingF_id.txt";
while (<IN>) {
	chomp($_);
	push(@g,$_);
}
close IN;
my $count=@g;my $len=1;
if ($count<=$group) {
	$group=$count;
}

my $base_len = int($count / $group);
my $remainder = $count % $group;

my $t = 1; 
my $current_count = 0; 
my $target = $base_len + ($remainder > 0 ? 1 : 0);

for (my $j = 0; $j < @g; $j++) {
    my $t = ($j % $group) + 1; 
    open OUT, ">>", "$path/Group/Group_$t.txt";
    print OUT $g[$j], "\n";
    close OUT;
}

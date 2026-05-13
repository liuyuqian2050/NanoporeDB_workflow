use strict;
use warnings;

my $g=$ARGV[0];
my $path=$ARGV[1];
my $file="$path/Group/Group_$g\.txt";
my @pp=split(/\//,$path);
my $basename=pop(@pp);
my $path1=join("/",@pp);

my $pdbpath="$path1/1nanopore_query/PDB_nanopore";
my $afpath="$path/nanopore_AF3";
my $out="$path/USalign_AF3log";
if (-e "$out") {
	`rm -rf $out`;
}
`mkdir $out`;

opendir DIR,"$pdbpath";
my @dots=readdir(DIR);
close DIR;

my @p=();
for (my $j=0;$j<@dots ;$j++) {
	if ($dots[$j]=~ /(.*?)\.cif/) {
		my $pdb=$1;
		push(@p,$pdb);
	}
}
my $t=0;
open IN,"$file";
while (<IN>) {
	chomp($_);
	my $id=$_;
	my $cif="$afpath/$id\.cif";
	for (my $k=0;$k<@p ;$k++) {
		my $pdb="$pdbpath/$p[$k]\.cif";
		my $log="$out/".$id."_".$p[$k].".log";
		`USalign -mol prot -mm 1 -ter 1 $cif $pdb -outfmt 2 >$log`;
		$t++;
		print $t,"\t",$id,"\t",$p[$k],"\n";
	}
}
close IN;

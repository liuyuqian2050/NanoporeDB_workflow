use strict;
use warnings;
use Bio::SeqIO;

my $group=$ARGV[0];
die "Usage: perl $0 <groups>\n" unless defined $group;
my $path=`pwd`;
chomp($path);
my $paths4="$path/4Multimer_prediction";

my @g=();
open IN,"3_nanopore_allminingF_id.txt";
while (<IN>) {
	chomp($_);
	push(@g,$_);
}
close IN;
my $countall=@g;
if ($countall<=$group) {
	$group=$countall;
}

my $struc_filter="$path/5nanopore_struc_merge";
if (-e "$struc_filter") {
	`rm -rf $struc_filter`;
}
`mkdir $struc_filter`;
my $pathAFM="$paths4/nanopore_AFM";
my $pathAF3="$paths4/nanopore_AF3";

`perl bin/score_merge.pl $pathAFM $pathAF3`;
`perl bin/group_split.pl $group $paths4`;

######### US-align AFM,AF3-> PDB
my @childs=();
my $count=0;
for (my $j=0;$j<$group ;$j++) {
	$count++;                          #   print "fork process $count\n";
	my $pid = fork();			
	if ($pid) {                        #   parent
		sleep 1;
		push(@childs, $pid);
	}
	elsif ($pid == 0) {
		print "       ","fork process $count\n";	
		`perl bin/runTMscoreAFMmulti.pl $count $paths4`;
		`perl bin/runTMscoreAF3multi.pl $count $paths4`;
		exit(0);
	}
	else {
		die "couldn't fork: $!\n";
	}
}
foreach (@childs) {
	waitpid($_,0);
}

######### USalign merge
my $path3="$paths4/USalign_AF3log";
opendir DIR,"$path3";
my @dots=readdir(DIR);
close DIR;
open OUT,">$paths4/AF3_USalign_TMscore.txt";
for (my $j=0;$j<@dots ;$j++) {
	if ($dots[$j]=~ /(.*?)\.log/) {
		my @a=split(/\_/,$1);
		my $query=shift(@a);
		my $target=join("_",@a);
		open IN,"$path3/$dots[$j]";
		while (<IN>) {
			chomp($_);
			my @a=split(/\s+/,$_);
			if ($a[0]!~ /^#/) {
				my $tm1=$a[2];
				my $tm2=$a[3];
				print OUT $query,"\t",$target,"\t",$tm1,"\t",$tm2,"\n";
			}
		}
		close IN;
	}
}
close OUT;

my $path2="$paths4/USalign_AFMlog";
opendir DIR,"$path2";
my @dot=readdir(DIR);
close DIR;
open OUT,">$paths4/AFM_USalign_TMscore.txt";
for (my $j=0;$j<@dot ;$j++) {
	if ($dot[$j]=~ /(.*?)\.log/) {
		my @a=split(/\_/,$1);
		my $query=shift(@a);
		my $target=join("_",@a);
		open IN,"$path2/$dot[$j]";
		while (<IN>) {
			chomp($_);
			my @a=split(/\s+/,$_);
			if ($a[0]!~ /^#/) {
				my $tm1=$a[2];
				my $tm2=$a[3];
				print OUT $query,"\t",$target,"\t",$tm1,"\t",$tm2,"\n";
			}
		}
		close IN;
	}
}
close OUT;


my %AList=();my %BList=();my $t1=0;my $t2=0;
open IN,"$paths4/AFM_USalign_TMscore.txt";
while (<IN>) {
	chomp($_);
	my @a=split(/\t/,$_);
	if ($a[2]>=0.8 or $a[3]>=0.8) {
		if (! $AList{$a[0]}) {
			$AList{$a[0]}=1;
			$t1++;
		}
	}
}
close IN;
print "AFM:",$t1,"\n";

open IN,"$paths4/AF3_USalign_TMscore.txt";
while (<IN>) {
	chomp($_);
	my @a=split(/\t/,$_);
	if ($a[2]>=0.8 or $a[3]>=0.8) {
		if (! $BList{$a[0]}) {
			$BList{$a[0]}=1;
			$t2++;
		}
	}
}
close IN;
print "AF3:",$t2,"\n";

my $m1=0;my $m2=0;my $m3=0;
open OUT,">5_nanopore_ALLmerge.txt";
open IN,"$paths4/nanopore_ALLpred.record";
while (<IN>) {
	chomp($_);
	my @a=split(/\t/,$_);
	if ($AList{$a[0]} or $BList{$a[0]}) {
		if ($AList{$a[0]} and $BList{$a[0]}) {
			if ($a[3]>=$a[$#a] and $a[3]>=0.7) {
				print OUT $a[0],"\t","AFM\t",$a[3],"\t","AFM:",$a[3],";AF3:",$a[$#a],"\n";
			}
			if ($a[3]<$a[$#a] and $a[$#a]>=0.7) {
				print OUT $a[0],"\t","AF3\t",$a[$#a],"\t","AFM:",$a[3],";AF3:",$a[$#a],"\n";
			}
		}
		elsif ($AList{$a[0]} and $a[3]>=0.7) {
			print OUT $a[0],"\t","AFM\t",$a[3],"\n";
		}
		elsif ($BList{$a[0]} and $a[$#a]>=0.7) {
			print OUT $a[0],"\t","AF3\t",$a[$#a],"\n";
		}
	}
}
close IN;
close OUT;


######### Structure merge
my $h=0;my $fcp=0;
open IN,"5_nanopore_ALLmerge.txt";
while (<IN>) {
	chomp($_);
	my @a=split(/\t/,$_);
	if ($a[1]=~ /AFM/) {
		my $pdb=$pathAFM."/".$a[0].".pdb";
		my $json=$pathAFM."/".$a[0].".json";
		if (-e "$pdb") {
			`cp $pdb $struc_filter`;
			`cp $json $struc_filter`;
			$fcp++;
		}
		else {
			print $pdb," not found!\n"
		}
		$h++;
	}
	else {
		my $pdb=$pathAF3."/".$a[0].".cif";
		my $json=$pathAF3."/".$a[0].".json";
		if (-e "$pdb") {
			`cp $pdb $struc_filter`;
			`cp $json $struc_filter`;
			$fcp++;
		}
		else {
			print $pdb," not found!\n"
		}
		$h++;
	}
}
close IN;
if ($h != $fcp) {
	print "Please check if the structure/json name is consistent with 3_nanopore_allminingF_id.txt!";
}
print "Structure merge:",$h,"\n";

use strict;
use warnings;

my $path_AFM=$ARGV[0];
my $path_AF3=$ARGV[1];
my @p=split(/\//,$path_AFM);
my $basename=pop(@p);
my $path=join("/",@p);

`python bin/countplddt_chainCA.py $path_AFM $path/AFM_pLDDT.txt`;
`python bin/countplddt_chainCA.py $path_AF3 $path/AF3_pLDDT.txt`;

my %CList=();my $tm=0;
open IN,"$path/AFM_pLDDT.txt";
while (<IN>) {
	chomp($_);
	my @a=split(/\t/,$_);
	my $plddt=$a[1];
	if ($a[1]>1) {
		$plddt=$a[1]/100;
	} 
	if ($a[0]=~ /(.*?)\.pdb/) {
		my $id=$1;
		$CList{$id}=$plddt;
		$tm++;
	}
}
close IN;
print "AFM_plddt_record: ",$tm,"\n";

opendir DIR,"$path_AFM";
my @dots=readdir(DIR);
close DIR;
my %AList=();my %BList=();
for (my $j=0;$j<@dots ;$j++) {
	if ($dots[$j]=~ /(.*?)\.json/) {
		my $id=$1;
		my $iptm=0;my $ptm=0;my $plddt=0;
		open IN,"$path_AFM/$dots[$j]";
		while (<IN>) {
			chomp($_);
			if ($_=~ /\"iptm\"\: (.*?)\,/) {
				$iptm=$1;
			}
			elsif ($_=~ /\"ptm\"\: (.*?)\,/) {
				$ptm=$1;
			}
		}
		close IN;
		if ($CList{$id}) {
			$plddt=$CList{$id};
			$AList{$id}=$iptm."\t".$ptm."\t".$plddt;
		}
		else {
			print "$id\_plddt not record in file!\n";
		}
	}
}

my %DList=();my $t=0;
open IN,"$path/AF3_pLDDT.txt";
while (<IN>) {
	chomp($_);
	my @a=split(/\t/,$_);
	my $plddt=$a[1]/100;
	if ($a[0]=~ /(.*?)\.cif/) {
		my $id=$1;
		$DList{$id}=$plddt;
		$t++;
	}
}
close IN;
print "AF3_plddt_record: ",$t,"\n";

opendir DIR,"$path_AF3";
my @dot=readdir(DIR);
close DIR;
for (my $j=0;$j<@dot ;$j++) {
	if ($dot[$j]=~ /(.*?)\.json/) {
		my $id=$1;
		my $iptm=0;my $ptm=0;my $plddt=0;
		open IN,"$path_AF3/$dot[$j]";
		while (<IN>) {
			chomp($_);
			if ($_=~ /\"iptm\"\: (.*?)\,/) {
				$iptm=$1;
			}
			elsif ($_=~ /\"ptm\"\: (.*?)\,/) {
				$ptm=$1;
			}
		}
		close IN;
		if ($DList{$id}) {
			$plddt=$DList{$id};
			$BList{$id}=$iptm."\t".$ptm."\t".$plddt;
		}
		else {
			print "$id\_plddt not record in file!\n";
		}
	}
}

my $w=0;
open OUT,">$path/nanopore_ALLpred.record";
print OUT "ID\tiptm_AFM\tptm_AFM\tplddt_AFM\tiptm_AF3\tptm_AF3\tplddt_AF3\n";
foreach  (sort keys %AList) {
	if (! $BList{$_}) {
		print $_,"\n";
	}
	else {
		print OUT $_,"\t",$AList{$_},"\t",$BList{$_},"\n";
		$w++;
	}
}
close OUT;
print "nanopore_all:",$w,"\n";

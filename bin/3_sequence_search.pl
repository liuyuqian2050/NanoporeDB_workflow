use strict;
use warnings;
use Bio::SeqIO;

my $Path_mmseqs2_uniref90=$ARGV[0];#cretedb with mmseqs2 (mmseqs createdb MGnify90FL.fa MGnify90FL)
my $Path_mmseqs2_MGYP90FL=$ARGV[1];#cretedb with mmseqs2 (mmseqs createdb uniref90.fasta uniref90)
my $Path_uniref90_fasta=$ARGV[2];#fasta DB of uniref90
my $Path_MGYP90FL_fasta=$ARGV[3];#fasta DB of MGYP90FL(full length sequences)
my $threads=$ARGV[4];

`mmseqs easy-search 2_nanopore_AFDBt70_ALLF.pep $Path_mmseqs2_uniref90 3_nanopore_uniref90_0.7.aln tmp -c 0.9 -e 1e-4 --threads $threads --format-output query,target,fident,alnlen,qcov,tcov,qstart,qend,tstart,tend,evalue,bits`;
`mmseqs easy-search 2_nanopore_AFDBt70_ALLF.pep $Path_mmseqs2_MGYP90FL 3_nanopore_MGYP90FL_0.7.aln tmp -c 0.9 -e 1e-4 --threads $threads --format-output query,target,fident,alnlen,qcov,tcov,qstart,qend,tstart,tend,evalue,bits`;

my %DList=();my %EList=();
open OUT,">3_nanopore_allmining.pep";
my $m1=0;my $m2=0;
my $seqio2=new Bio::SeqIO(-file => "2_nanopore_AFDBt70_ALLF.pep", -format => 'fasta');
while (my $seq=$seqio2->next_seq) {
	if (! $DList{$seq->seq()}) {
		$DList{$seq->seq()}=$seq->display_id;
		print OUT ">",$seq->display_id," ",$seq->description,"\n",$seq->seq(),"\n";
		$m1++;
	}
	else {
		my @b=split(/\_/,$seq->display_id);
		my $rt=pop(@b);
		my $id=join("_",@b);
		if (! $EList{$id}) {
			$EList{$id}=1;
			$m2++;
		}
	}
}
print "afdb_fetch: ",$m1,"\t","dup: ",$m2,"\n";

my $uniref90fetch="3_nanopore_uniref90_0.7.aln";
my $MGYP90FLfetch="3_nanopore_MGYP90FL_0.7.aln";

my %AList=();my %BList=();my $t1=0;my $t2=0;
open IN,"$uniref90fetch";
while (<IN>) {
	chomp($_);
	my @a=split(/\s+/,$_);
	if (! $AList{$a[1]}) {
		$AList{$a[1]}=1;
		$t1++;
	}
}
close IN;
print $uniref90fetch,": ",$t1,"\n";

open IN,"$MGYP90FLfetch";
while (<IN>) {
	chomp($_);
	my @a=split(/\s+/,$_);
	if (! $BList{$a[1]}) {
		$BList{$a[1]}=1;
		$t2++;
	}
}
close IN;
print $MGYP90FLfetch,": ",$t2,"\n";

my %SList=();
my $r1=0;my $r2=0;
my $seqio=new Bio::SeqIO(-file => "$Path_uniref90_fasta", -format => 'fasta');
while (my $seq=$seqio->next_seq) {
	if ($AList{$seq->display_id}) {
		if (! $SList{$seq->seq()}) {
			$SList{$seq->seq()}=1;
			print OUT ">",$seq->display_id," ",$seq->description,"\n",$seq->seq(),"\n";
			$r1++;
		}
		else {
			$r2++;
		}
	}
}
print "uniref90_fetch: ",$r1,"\tdup: ",$r2,"\n";

my %GList=();
my $r11=0;my $r12=0;
my $seqio1=new Bio::SeqIO(-file => "$Path_MGYP90FL_fasta", -format => 'fasta');
while (my $seq=$seqio1->next_seq) {
	if ($BList{$seq->display_id}) {
		if (! $GList{$seq->seq()}) {
			$GList{$seq->seq()}=1;
			print OUT ">",$seq->display_id," ",$seq->description,"\n",$seq->seq(),"\n";
			$r11++;
		}
		
		else {
			$r12++;
		}
	}
}
close OUT;
print "MGYP90FL_fetch: ",$r11,"\tdup: ",$r12,"\n";

my $af=0;my $uni=0;my $mgy=0;my $pd=0;my %HList=();
open OUT,">3_nanopore_allminingF.fasta";
open OUT1,">3_nanopore_allminingF_id.txt";
my $seqio3=new Bio::SeqIO(-file => "3_nanopore_allmining.pep", -format => 'fasta');
while (my $seq=$seqio3->next_seq) {
	if (! $HList{$seq->seq()}) {
		$HList{$seq->seq()}=1;
		if ($seq->display_id=~ /^AF-/) {
			$af++;
			my @d=split(/\-F1/,$seq->display_id);
			my $idout=$d[0];
			print OUT ">",$idout," ",$seq->description,"\n",$seq->seq(),"\n";
			print OUT1 $idout,"\n";
		}
		elsif ($seq->display_id=~ /^UniRef/) {
			$uni++;
			print OUT ">",$seq->display_id," ",$seq->description,"\n",$seq->seq(),"\n";
			print OUT1 $seq->display_id,"\n";
		}
		elsif ($seq->display_id=~ /^MGYP/) {
			$mgy++;
			print OUT ">",$seq->display_id," ",$seq->description,"\n",$seq->seq(),"\n";
			print OUT1 $seq->display_id,"\n";
		}
		else {
			$pd++;
		}
		
	}
}
close OUT;
close OUT1;
print "AFDB_fetch:",$af,"\tPDB_fetch:",$pd,"\tuniref_fetch:",$uni,"\tMGYP_fetch:",$mgy,"\n";

use strict;
use warnings;
use Bio::SeqIO;

my $path_pdb=$ARGV[0];#1nanopore_query/PDB_nanopore
my $path_afdb=$ARGV[1];#1nanopore_query/AFDB_nanopore
my $afdb_foldseek_path=$ARGV[2];#download the Foldseek format afdb
my $threads=40;
if ($ARGV[3]) {
	$threads=$ARGV[3];
}
my @p=split(/\//,$path_afdb);
my $basename=pop(@p);
my $path=join("/",@p);

#####1. retain the high-confidence monomeric structures (pLDDT>=70)
print "1. retain the high-confidence monomeric structures (pLDDT>=70).\n";
`python bin/1_clean_plddt70.py $path_afdb`;
my $path_afdb70=$path_afdb;
$path_afdb70 =~ s/\/?$/_t70/;

#####2.monomer search
print "2.monomer search...\n";
print "	AFDB search step1...\n";
`foldseek easy-search $path_pdb $path_afdb70 2_nanoporeF.aln tmp -c 0.7 -e 0.01 --alignment-type 1 --tmscore-threshold 0.8 --threads $threads --format-output query,target,alntmscore,qtmscore,ttmscore,qcov,tcov,lddt,prob`;
my $outpathf=$path_afdb70."F";
if (! -e "$outpathf") {
	`mkdir $outpathf`;
}

my $tmaln="2_nanoporeF.aln";
my %FList=();
open IN,"$tmaln";
while (<IN>) {
	chomp($_);
	my @a=split(/\t/,$_);
	if (! $FList{$a[1]}) {
		$FList{$a[1]}=1;
		my $pdb=$a[1].".pdb";
		if (-e "$path_afdb70/$pdb") {
			`cp $path_afdb70/$pdb $outpathf`;
		}
	}
}
close IN;

###################################
print "	AFDB search step2...\n";
`foldseek easy-search $outpathf $afdb_foldseek_path/afdb 2_nanopore_AFDBall.aln tmp -c 0.7 -e 0.01 --alignment-type 1 --tmscore-threshold 0.8 --threads $threads --format-output query,target,alntmscore,qtmscore,ttmscore,qcov,tcov,lddt,prob`;

my $path1="$path/AFDB_nanopore_t70F_AFDBadd";
if (! -e "$path1") {
	`mkdir $path1`;
}

opendir DIR,"$outpathf";
my @dots=readdir(DIR);
close DIR;

my $t0=0;my %EList=();
for (my $j=0;$j<@dots ;$j++) {
	if ($dots[$j]=~ /(.*?)\.pdb/) {
		my $id=$1;
		$t0++;
		$EList{$id}=1;
	}
}

my $t3=0;
open IN,"2_nanopore_AFDBall.aln";
while (<IN>) {
	chomp($_);
	my @a=split(/\s+/,$_);
	my $id=$a[1];
	if (! $EList{$id}) {
		$EList{$id}=1;
		my $pdb=$id.".pdb";
		$t3++;
		my $url="https://alphafold.ebi.ac.uk/files/$pdb";
		if ((! -e "$outpathf/$pdb") and (! -e "$path1/$pdb")) {
			`wget -P $path1/ --no-check-certificate $url`;
		}
	}
}
close IN;
`python bin/1_clean_plddt70.py $path1`;

#####3. fetch sequences
print "3. fetch sequences...\n";
my $input_folder1=$path1."_t70";
my $pep_t70F="3AFDB_nanopore_t70F.fasta0";
my $pep_t70Fafdb="3AFDB_nanopore_t70F_AFDBall_t70.fasta0";
my $out="2_nanopore_AFDBt70_ALL.pep";
my $out1="2_nanopore_PDB.pep";

`python bin/fetch_seqs0.py $outpathf $pep_t70F`;
`python bin/fetch_seqs0.py $input_folder1 $pep_t70Fafdb`;
`python bin/fetch_seqs0.py $path_pdb PDB.fas`;

my %AList=();
open OUT,">$out";
my $seqio   = new Bio::SeqIO(-file => "$pep_t70F", -format => 'fasta');
while ( my $seq = $seqio->next_seq ) {
	if ($seq->display_id !~ /^AF-/) {
		my @a=split(/\_/,$seq->display_id);
		my $id=$a[0];
		if (! $AList{$id}) {
			$AList{$id}=$seq->seq();
			print OUT ">",$seq->display_id,"\n",$seq->seq(),"\n";
		}
		else {
			if ($AList{$id} ne $seq->seq()) {
				print OUT ">",$seq->display_id,"\n",$seq->seq(),"\n";
			}
		}
	}
	else {
		print OUT ">",$seq->display_id,"\n",$seq->seq(),"\n";
	}
}
my $seqio2   = new Bio::SeqIO(-file => "$pep_t70Fafdb", -format => 'fasta');
while ( my $seq = $seqio2->next_seq ) {
	print OUT ">",$seq->display_id,"\n",$seq->seq(),"\n";
}
close OUT;

###################################
my %GList=();my %HList=();my $p1=0;my$p2=0;
opendir DIR,"$path_pdb";
my @dot=readdir(DIR);
close DIR;
for (my $j=0;$j<@dot ;$j++) {
	if ($dot[$j]=~ /(.*?)\.cif/ or $dot[$j]=~ /(.*?)\.pdb/) {
		$p1++;
	}
}

open OUT,">2_nanopore_AFDBt70_ALLF.pep";
open OUT1,">$out1";
my $seqio3=new Bio::SeqIO(-file => "PDB.fas", -format => 'fasta');
while (my $seq=$seqio3->next_seq) {
	my $id=$seq->display_id;
	if (! $HList{$seq->seq()}) {
		$HList{$seq->seq()}=1;
		print OUT1 ">",$id," ",$seq->description,"\n",$seq->seq(),"\n";
		print OUT ">",$id," ",$seq->description,"\n",$seq->seq(),"\n";
	}
}

my $seqio4=new Bio::SeqIO(-file => "$out", -format => 'fasta');
while (my $seq=$seqio4->next_seq) {
	my $id=$seq->display_id;
	if (! $GList{$seq->seq()}) {
		$GList{$seq->seq()}=1;
		$p2++;
		print OUT ">",$id," ",$seq->description,"\n",$seq->seq(),"\n";
	}
}
close OUT;
close OUT1;
print "nanopore_count: PDB:",$p1,"\tAFDB_filter:",$p2,"\n";

unlink("$pep_t70F");
unlink("$pep_t70Fafdb");
unlink("$out");
unlink("PDB.fas");

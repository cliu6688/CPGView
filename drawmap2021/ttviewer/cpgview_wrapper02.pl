#!/usr/bin/perl
#
use FindBin qw($Bin);

my $inputfile = shift;
my $projectid = shift;

print STDERR "$Bin $inputfile $projectid\n";

if (! -e $inputfile) {
	print STDERR "\ninputfile does not exist\n\nUsage: $0 inputfile projectid (optional\n\n";
	exit;
}

if ($projectid eq "") {
	$projectid = time();
	$projectid =~ s/\s+//g;
}

my $outdir = "/var/www/html/upload/d".$projectid;
#my $outdir = "/tmp/d".$projectid;
my $inputfile1 = $projectid.".gbf";

my $cmd = "mkdir $outdir";

`$cmd` unless (-e $outdir);
`cp $inputfile $outdir/$inputfile1`;

my $python = "/biodata4/home/cliu/cpgview/bin/python";
my $data_prep = $Bin."/data_process.py";

my $cmd = "(cd $outdir; $python $data_prep $inputfile1 dc dt)";
print STDERR "$cmd\n"; `$cmd`;

my $outfile_c = $projectid."_cis.pdf";
my $cmd = "(cd $Bin; Rscript viewCSGene.R $outdir/dc/cis_splicing_gene.csv $outdir/dc/cis_splicing_subgene.csv $outdir/$outfile_c)";
print STDERR "$cmd\n"; `$cmd`;

my $outfile_t = $projectid."_trans.pdf";
my $cmd = "(cd $Bin; Rscript viewTSGene.R $outdir/dt/trans_splicing_gene.csv $outdir/dt/trans_splicing_subgene.csv $outdir/$outfile_t)";
print STDERR "$cmd\n"; `$cmd`;

my $outfile_cc = $projectid."_circle";
my $cmd = "(cd $Bin/Chloroplot/; Rscript chloroplot_Genes.R $outdir/$inputfile1 $outdir/$outfile_cc)";
print STDERR "$cmd\n"; `$cmd`;

my $msg = "\n\nThe run is completed successfully! You file can be found in $outdir\n\n";
print STDERR $msg;

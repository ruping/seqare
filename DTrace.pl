#!/usr/bin/perl -w

use Getopt::Long;
use Data::Dumper;
use strict;
use File::Glob ':glob';
use File::Basename;
use FindBin qw($RealBin);
use Parallel::ForkManager;
use lib "$RealBin/lib";
use bwaMapping;
use snvCalling;
use seqStats;

my %options;
my %runlevel;
my %runTask;

$options{'noexecute'}   = 0;
$options{'quiet'}       = 0;
$options{'runlevels'}   = undef;
$options{'runTask'}     = undef;
$options{'readlen'}     = 0;
$options{'mapper'}      = "bwa";
$options{'sampleName'}  = 'SRP';
$options{'FASTQ1'}      = 'SRP';
$options{'FASTQ2'}      = 'SRP';
$options{'fastqFiles1'} = 'SRP';
$options{'fastqFiles2'} = 'SRP';
$options{'bams'}        = 'SRP';
$options{'bamID'}       = 1;
$options{'lanepath'}    = 'SRP';
$options{'threads'}     = 1;
$options{'splitChr'}    = undef;
$options{'help'}        = undef;
$options{'qcOFF'}       = undef;
$options{'root'}        = "$RealBin/../PIPELINE";
$options{'readpool'}    = 'SRP';
$options{'gf'}          = "png";                 #the format used in html report
$options{'bzip'}        = undef;                 #to allow bzip compressed fastq files
$options{'Rbinary'}     = 'R';
$options{'seqType'}     = 'WXS,paired-end';      #experimental types
$options{'tmpDir'}      = '';
$options{'bin'}         = "$RealBin/";
$options{'configure'}   = "SRP";

$options{'somaticInfo'} = "SRP";


if (@ARGV == 0) {
  helpm();
} else {
  printf STDERR "\n# $0 %s\n",join(" ",@ARGV);
}


GetOptions(
           "task|t=s"     => \$options{'task'},
           "sampleName=s" => \$options{'sampleName'},
           "FASTQ1=s"     => \$options{'FASTQ1'},
           "FASTQ2=s"     => \$options{'FASTQ2'},
           "fastqFiles1=s"=> \$options{'fastqFiles1'},
           "fastqFiles2=s"=> \$options{'fastqFiles2'},
           "bams=s"       => \$options{'bams'},
           "bamID=s"      => \$options{'bamID'},
           "qcOFF"        => \$options{'qcOFF'},
           "runID=s"      => \$options{'runID'},
           "runlevel=s"   => \$options{'runlevels'},
           "runTask=s"    => \$options{'runMethod'},
           "seqType=s"    => \$options{'seqType'},
           "noexecute"    => \$options{'noexecute'},
           "quiet"        => \$options{'quiet'},
           "splitChr"     => \$options{'splitChr'},
           "readlen=i"    => \$options{'readlen'},
           "mapper=s"     => \$options{'mapper'},
           "threads=i"    => \$options{'threads'},
           "gf=s"         => \$options{'gf'},
           "root=s"       => \$options{'root'},
           "readpool=s"   => \$options{'readpool'},
           "bzip"         => \$options{'bzip'},
           "Rbinary=s"    => \$options{'Rbinary'},
           "help|h"       => \$options{'help'},
           "configure=s"  => \$options{'configure'},
           "somaticInfo=s"=> \$options{'somaticInfo'},
           "tmpDir=s"     => \$options{'tmpDir'},
          );

#print help
helpm() if ($options{'help'});


### Read configuration and set all paths----------------------------------
my %confs;
open IN, "$options{'configure'}";
while ( <IN> ) {
  chomp;
  next if /^#/;
  my @cols = split /\t/;
  $confs{$cols[0]} = $cols[1];
}
close IN;

#translate environment variable
foreach my $confele (keys %confs){
  while ($confs{$confele} =~ /\$([A-Za-z0-9]+)/g) {
    my $eleName = $1;
    my $eleTranslate;
    if (exists ($confs{$eleName})) {
      $eleTranslate = $confs{$eleName};
      $confs{$confele} =~ s/\$$eleName/$eleTranslate/;
    } else {
      die("can't translate eleName: $eleName\n");
    }
  }
}
print STDERR Dumper (\%confs);
#-------------------------------------------------------------------------

### Frequently used names-------------------------------------------------
my @chrs = split(/\,/, $confs{'chrs'});
#-------------------------------------------------------------------------

#decompression option-----------------------------------------------------
$options{'decompress'} = "gzip -d -c";
$options{'compress'}   = "gzip";
$options{'zipSuffix'}  = "gz";
if ( $options{'bzip'} ) {
  $options{'decompress'} = "bzip2 -d -c";
  $options{'compress'}   = "bzip2";
  $options{'zipSuffix'}  = "bz2";
}
#-------------------------------------------------------------------------

### Already specified full path fastq files-------------------------------
if ($options{'fastqFiles1'} ne 'SRP'){
  $options{'fastqFiles1'} =~ s/\,/ /g;
}
if ($options{'fastqFiles1'} ne 'SRP'){
  $options{'fastqFiles2'} =~ s/\,/ /g;
}
#-------------------------------------------------------------------------

### Runlevel/Task check up------------------------------------------------
if ($options{'runlevels'}) { #true runlevels
  foreach my $r (split /\,/,$options{'runlevels'}) {
    my $from=1;
    my $to=20;
    if ($r=~/^(\d+)/) {
      $from=$1;
    }
    if ($r=~/\-(\d+)$/) {
      $to=$1;
    } elsif ($r!~/\-/) {
      $to=$from;
    }
    for (my $i=$from;$i<=$to;$i++) {
      $runlevel{$i}=1;
    }
  }
} elsif ($options{'runTask'}) {
  foreach my $task (split(/\,/, $options{'runTask'})) {
    $runTask{$task} = '';
  }
} else {
  print STDERR "no runlevel or runTask has been set, exit.\n";
  helpm();
}
#-------------------------------------------------------------------------

if ($options{'root'} eq "$RealBin/../PIPELINE") {
  if (-e "$RealBin/../PIPELINE") {
    print STDERR "no root dir given, analysis will be run under $options{'root'}.\n";
  }
  else {
    print STDERR "no root dir given, $options{'root'} does not exist, please do -h or --help to check how to set root dir.\n";
    helpm();
  }
} else {
  $options{'readpool'} = $options{'root'} if $options{'readpool'} eq 'SRP';
}


#store somatic information------------------------------------------------
my %somatic;
my %germline;                   #may have multiple tumors
if (-s "$options{'somaticInfo'}") {

  open IN, "$options{'somaticInfo'}";
  while ( <IN> ) {
    chomp;
    s/[\s\n]$//;
    my @columns = split /\t/;
    my $tumor = $columns[0];
    my $normal = $columns[1];

    $somatic{$tumor} = $normal;
    push(@{$germline{$normal}}, $tumor) if $normal ne 'undef';
  }
  close IN;
  #print STDERR Dumper (\%somatic);
  #print STDERR Dumper (\%germline);
}
#-------------------------------------------------------------------------


###
###preparation the lane and read path enviroment
###

if ($options{'lanepath'} eq 'SRP') {
  printtime();
  $options{'lanepath'} = "$options{'root'}/$options{'sampleName'}";   #define lane path
  print STDERR "####### lane name is set to $options{'sampleName'} #######\n\n";
  unless (-e "$options{'lanepath'}") {
    my $cmd = "mkdir -p $options{'lanepath'}";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }
}

if ($options{'bams'} ne 'SRP') {
  unless (-e "$options{'lanepath'}/02_MAPPING") {
    my $cmd = "mkdir -p $options{'lanepath'}/02_MAPPING";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }
  my $finalBam = "$options{'lanepath'}/02_MAPPING/$options{'sampleName'}\.sorted\.ir\.rmDup\.bam";
  unless (-s "$finalBam") {
    my $cmd = "ln -s $options{'bams'} $finalBam";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }
  unless (-s "$finalBam\.bai"){
    my $cmd = "ln -s $options{'bams'}\.bai $finalBam\.bai";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }
  goto REALSTEPS;
}

if ($options{'readpool'} ne 'SRP' and $options{'FASTQ1'} ne 'SRP' and $options{'fastqFiles1'} eq 'SRP') {

  printtime();
  print STDERR "####### preparing directories #######\n\n";

  unless (-e "$options{'lanepath'}/01_READS/") {
    my $cmd = "mkdir -p $options{'lanepath'}/01_READS/";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  my @fastqFile1 = split(/\,/, $options{'FASTQ1'});
  $options{'fastqFiles1'} =~ s/^SRP//;
  foreach my $fastqFile1 (@fastqFile1) {
    if ($fastqFile1 !~ /\.[bg]z2?$/){
      die "\[error\]: $fastqFile1 must be gzip or bzipped!\n";
    }
    $fastqFile1 = $options{'readpool'}.'/'.$fastqFile1;
    $options{'fastqFiles1'} .= $fastqFile1." ";
  }
  $options{'fastqFiles1'} =~ s/\s$//;

  if ($options{'FASTQ2'} ne 'SRP') {
    my @fastqFile2 = split(/\,/, $options{'FASTQ2'});
    $options{'fastqFiles2'} =~ s/^SRP//;
    foreach my $fastqFile2 (@fastqFile2) {
      $fastqFile2 = $options{'readpool'}.'/'.$fastqFile2;
      $options{'fastqFiles2'} .= $fastqFile2." ";
    }
    $options{'fastqFiles2'} =~ s/\s$//;
  }

  print STDERR "lanefile1:\t$options{'fastqFiles1'}\n";
  print STDERR "lanefile2:\t$options{'fastqFiles2'}\n"; #if paired end

}


foreach my $fastqFile1 (split(" ", $options{'fastqFiles1'})){
  my $cmd = "ln -s $fastqFile1 $options{'lanepath'}/01_READS/";
  my $fastqFile1Basename = basename($fastqFile1);
  RunCommand($cmd,$options{'noexecute'},$options{'quiet'}) unless (-s "$options{'lanepath'}/01_READS/$fastqFile1Basename");
}


if ($options{'fastqFiles2'} ne 'SRP'){
  foreach my $fastqFile2 (split(" ", $options{'fastqFiles2'})){
    my $cmd = "ln -s $fastqFile2 $options{'lanepath'}/01_READS/";
    my $fastqFile2Basename = basename($fastqFile2);
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'}) unless (-s "$options{'lanepath'}/01_READS/$fastqFile2Basename");
  }
}


if ( $options{'readlen'} == 0 ) { #read length not set
  my @fastqFiles1Temp = split(/\s/, $options{'fastqFiles1'});
  my $first_second_line = `$options{'decompress'} "$fastqFiles1Temp[0]" | head -2 | grep -v "^@"`;
  $options{'readlen'} = length($first_second_line) - 1;
  print STDERR "read length is not set, will take the original read length ($options{'readlen'} bp)\n";
}

#######################################################################################################
REALSTEPS:

###
###runlevel1: QC
###

my $runlevels = 1;
if (exists($runlevel{$runlevels}) or exists($runTask{'QC'})) {

  printtime();
  print STDERR "####### runlevel $runlevels now #######\n\n";

  my $qc_out1 = "$options{'lanepath'}/01_READS/mate1.qc";    ######
  my $qc_out2;

  if ($options{'fastqFiles2'} ne 'SRP') {
    $qc_out2 = "$options{'lanepath'}/01_READS/mate2.qc";
  }

  unless (-e "$qc_out1") {
    my $cmd = "$options{'decompress'} $options{'fastqFiles1'} | $options{'bin'}/fastx_quality_stats -Q33 -o $qc_out1";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'}) unless ($options{'qcOFF'});
  }
  unless (($options{'fastqFiles2'} eq 'SRP') or (-e "$qc_out2")) {
    my $cmd = "$options{'decompress'} $options{'fastqFiles2'} | $options{'bin'}/fastx_quality_stats -Q33 -o $qc_out2";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'}) unless ($options{'qcOFF'});
  }

  printtime();
  print STDERR "####### runlevel $runlevels done #######\n\n";

}


###
###runlevel2: do the mapping and generate the statistics
###


$runlevels = 2;
if (exists($runlevel{$runlevels}) or exists($runTask{'mapping'}) or exists($runTask{'indelRealignment'}) or exists($runTask{'MarkDuplicates'})) {

  my $rawBam = "$options{'lanepath'}/02_MAPPING/$options{'sampleName'}\.bam";
  my $sortedBam = "$options{'lanepath'}/02_MAPPING/$options{'sampleName'}\.sorted\.bam";
  my $irBam = ($options{'splitChr'})?"$options{'lanepath'}/02_MAPPING/$options{'sampleName'}\.sorted\.ir\.$chrs[0]\.bam":"$options{'lanepath'}/02_MAPPING/$options{'sampleName'}\.sorted\.ir\.bam";
  my $finalBam = ($options{'splitChr'})?"$options{'lanepath'}/02_MAPPING/$options{'sampleName'}\.sorted\.ir\.$chrs[0]\.rmDup\.bam":"$options{'lanepath'}/02_MAPPING/$options{'sampleName'}\.sorted\.ir\.rmDup\.bam";

  printtime();
  print STDERR "####### runlevel $runlevels now #######\n\n";

  unless (-e "$options{'lanepath'}/02_MAPPING") {
    my $cmd = "mkdir -p $options{'lanepath'}/02_MAPPING";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  my @allBams = bsd_glob("$options{'lanepath'}/02_MAPPING/$options{'sampleName'}*\.bam");
  if ($#allBams == -1 or exists($runTask{'mapping'})) {
    my $ReadGroup = '@RG'."\tID:".$options{'bamID'}."\tSM\:".$options{'sampleName'};
    my $cmd;
    if ($options{'fastqFiles2'} eq 'SRP') { #single end
      $cmd = bwaMapping->bwaPairMapping($ReadGroup, $options{'threads'}, $confs{'BWAINDEX'}, $rawBam, $options{'fastqFiles1'});
    } else {                    #paired-end
      $cmd = bwaMapping->bwaPairMapping($ReadGroup, $options{'threads'}, $confs{'BWAINDEX'}, $rawBam, $options{'fastqFiles1'}, $options{'fastqFiles2'});
    }
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  unless (-s "$finalBam" and !(exists($runTask{'indelRealignment'}) or exists($runTask{'MarkDuplicates'})) ) {    #processing bam
    if (-s "$rawBam" and !(-s "$sortedBam")) {     #must sort
      my $cmd = bwaMapping->bamSort($options{'threads'}, $rawBam."\.tmp", $sortedBam, $rawBam);
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      $cmd = bwaMapping->bamIndex($sortedBam);     #index it
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    if (-s "$rawBam" and -s "$sortedBam") {
      my $cmd = "rm $rawBam -f";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }

    if ((-s "$sortedBam" and !(-s "$irBam")) or exists($runTask{'indelRealignment'})) { #indel realignment
      my $indelTargetList = $sortedBam."\.target_intervals.list";
      my $CHR = 'ALL';
      if ($options{'splitChr'}) {     #if split chr, folk it up
        my $chrBatches = partitionArray(\@chrs, $options{'threads'});
        foreach my $chrBatch (@{$chrBatches}) {
          my $manager = Parallel::ForkManager->new($options{'threads'});
          my $processedChroms = "chromosome ";
          foreach my $chrom (@{$chrBatch}) {
            $manager->start and next;
            $CHR = $chrom;
            $indelTargetList =~ s/\.target_intervals.list/\.$CHR\.target_intervals.list/;
            $irBam =~ s/\.bam/\.$CHR\.bam/;
            my $cmd = bwaMapping->indelRealignment1($confs{'gatkBin'}, $sortedBam, $confs{'GFASTA'}, $confs{'KNOWNINDEL1'}, $confs{'KNOWNINDEL2'}, $CHR, $indelTargetList);
            RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
            $cmd = bwaMapping->indelRealignment2($confs{'gatkBin'}, $sortedBam, $confs{'GFASTA'}, $indelTargetList, $confs{'KNOWNINDEL1'}, $confs{'KNOWNINDEL2'}, $CHR, $irBam);
            RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
            $cmd = bwaMapping->bamIndex($irBam);     #index it
            RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
            $manager->finish;
            $processedChroms .= $chrom.',';
          }
          $manager->wait_all_children;
          print STDERR "$processedChroms have been processed!\n";
        }
      } else {
        my $cmd = bwaMapping->indelRealignment1($confs{'gatkBin'}, $sortedBam, $confs{'GFASTA'}, $confs{'KNOWNINDEL1'}, $confs{'KNOWNINDEL2'}, $CHR, $indelTargetList);
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
        $cmd = bwaMapping->indelRealignment2($confs{'gatkBin'}, $sortedBam, $confs{'GFASTA'}, $indelTargetList, $confs{'KNOWNINDEL1'}, $confs{'KNOWNINDEL2'}, $CHR, $irBam);
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
        $cmd = bwaMapping->bamIndex($irBam);     #index it
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }
    }
    if (-s "$sortedBam" and -s "$irBam") {
      my $cmd = "rm $sortedBam $sortedBam\.bai -f";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }

    if ((-s "$irBam" and !(-s "$finalBam")) or exists($runTask{'MarkDuplicates'})) {  #rmDup
      my $rmDupMetric = $irBam.".rmDupMetric";
      my $cmd = bwaMapping->MarkDuplicates($confs{'MarkDuplicatesBin'}, $irBam, $finalBam, $rmDupMetric);
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      $cmd = bwaMapping->bamIndex($finalBam);     #index it
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    if (-s "$irBam" and -s "$finalBam") {
      my $cmd = "rm $irBam $irBam\.bai -f";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
  }

  printtime();
  print STDERR "####### runlevel $runlevels done #######\n\n";
}


###
###runlevel3: STATS
###

$runlevels = 3;
if (exists $runlevel{$runlevels}) {

  unless (-e "$options{'lanepath'}/03_STATS") {
    my $cmd = "mkdir -p $options{'lanepath'}/03_STATS";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  #basic read counting stats:
  my $mappingStats = "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.mapping.stats";
  my $finalBam = "$options{'lanepath'}/02_MAPPING/$options{'sampleName'}\.sorted\.ir\.rmDup\.bam";
  unless (-s "$mappingStats"){
    my $cmd = seqStats->mappingStats("$options{'bin'}/Rseq_bam_stats", $finalBam, $options{'readlen'}, $mappingStats);
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  #loren curve
  my $bedCover = "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.bedcoverNoDup";
  my $lorenzCover = "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.lorenzNoDup";
  unless (-s "$lorenzCover") {
    unless (-s "$bedCover") {
      my $cmd = seqStats->grepStarts("$options{'bin'}/grep_starts", $confs{'targetRegion'}, $finalBam, $bedCover);
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    my $cmd = seqStats->getLorenz("$options{'bin'}/lorenzCurveNGS.pl", $bedCover, $lorenzCover);
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }
  if (-s "$bedCover" and -s "$lorenzCover") {
    my $cmd = "rm $bedCover -f";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  printtime();
  print STDERR "####### runlevel $runlevels now #######\n\n";

  printtime();
  print STDERR "####### runlevel $runlevels done #######\n\n";

}


###
###runlevel4: SNV calling
###

$runlevels = 4;
if (exists $runlevel{$runlevels}) {

  unless (-e "$options{'lanepath'}/04_SNV") {
    my $cmd = "mkdir -p $options{'lanepath'}/04_SNV";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  printtime();
  print STDERR "####### runlevel $runlevels now #######\n\n";

  #my $finalBam = ($options{'splitChr'})?"$options{'lanepath'}/02_MAPPING/$options{'sampleName'}\.sorted\.ir\.$chrs[0]\.rmDup\.bam":"$options{'lanepath'}/02_MAPPING/$options{'sampleName'}\.sorted\.ir\.rmDup\.bam";
  my $finalBam = "$options{'lanepath'}/02_MAPPING/$options{'sampleName'}\.sorted\.ir\.rmDup\.bam";
  my $normalBam;
  if ($options{'somaticInfo'} eq "SRP"){
    print STDERR "ERROR: somaticInfo is not provided! Must set for somatic calling!\n";
    exit 22;
  } elsif ( !exists( $somatic{$options{'sampleName'}} ) ){
    print STDERR "ERROR: $options{'sampleName'} is not in the somatic hash table!\n";
  } else { #get normal bam
    my $normalSampleName = $somatic{$options{'sampleName'}};
    $normalBam = "$options{'root'}/$normalSampleName/02_MAPPING/$normalSampleName\.sorted\.ir\.rmDup\.bam";
    unless (-s "$normalBam"){
      print STDERR "ERROR: $normalBam is not found, please run mapping and processing for $normalSampleName!!\n";
      exit 22;
    }
  }

  my $muTectOut = "$options{'lanepath'}/04_SNV/$options{'sampleName'}\.mutect";
  my $vcfOutTmp = "$options{'lanepath'}/04_SNV/$options{'sampleName'}\.mutect.vcf";
  my $vcfOut = "$options{'lanepath'}/04_SNV/$options{'sampleName'}\.mutect.genome.vcf";
  my $vcfOutSorted = "$options{'lanepath'}/04_SNV/$options{'sampleName'}\.mutect.genome.sorted.vcf";
  my $vcfMultiAnno = "$options{'lanepath'}/04_SNV/$options{'sampleName'}\.mutect.genome.sorted.vcf.$confs{'species'}_multianno.txt";
  my $vcfMultiAnnoVCF = "$options{'lanepath'}/04_SNV/$options{'sampleName'}\.mutect.genome.sorted.vcf.$confs{'species'}_multianno.vcf";
  my $vcfMultiAnnoMod = "$options{'lanepath'}/04_SNV/$options{'sampleName'}\.mutect.genome.sorted.vcf.$confs{'species'}_multianno.mod.vcf";
  unless (-s "$muTectOut") {
    my $cmd = snvCalling->muTectCalling($confs{'muTectBin'}, $finalBam, $normalBam, $confs{'GFASTA'}, $confs{'muTectCOSMIC'}, $confs{'muTectDBSNP'}, $muTectOut, $vcfOutTmp);
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  #annoVar annotate---------------------------------------------------------------------
  if (-s "$muTectOut" and !-s "$vcfMultiAnnoMod") {

    my $cmd = snvCalling->muTect2vcf("$options{'bin'}/mutect2vcf.pl", $muTectOut, $vcfOut);                                                #convert mutect 2 vcf
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});

    $cmd = snvCalling->vcfSort($vcfOut, $vcfOutSorted);                                                                                    #sort vcf
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});

    $cmd = snvCalling->runAnnovar("$confs{'ANNOVARDIR'}/table_annovar.pl", $vcfOutSorted, $confs{'ANNOVARDB'}, $confs{'species'}, );       #table annovar
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});

    $cmd = snvCalling->convertVCFannovar("$options{'bin'}/convert_annovar_vcf.pl", $vcfMultiAnno, $vcfMultiAnnoVCF);
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});

  }

  #rm temporary files
  if (-s "$vcfMultiAnnoMod" and -s "$vcfOutTmp") {
    my $cmd = "rm -rf $vcfOutTmp $vcfOutTmp\.idx";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }
  if (-s "$vcfMultiAnnoMod" and -s "$vcfOut") {
    my $cmd = "rm -rf $vcfOut";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }
  if (-s "$vcfMultiAnnoMod" and -s "$vcfOutSorted\.avinput") {
    my $cmd = "rm -rf $vcfOutSorted\.avinput";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }
  if (-s "$vcfMultiAnnoMod" and -s "$vcfOutSorted\.invalid_input") {
    my $cmd = "rm -rf $vcfOutSorted\.invalid_input";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }
  if (-s "$vcfMultiAnnoMod" and -s "$vcfOutSorted\.refGene.invalid_input") {
    my $cmd = "rm -rf $vcfOutSorted\.refGene.invalid_input";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }
  if (-s "$vcfMultiAnnoMod" and -s "$vcfMultiAnno") {
    my $cmd = "rm -rf $vcfMultiAnno $vcfMultiAnnoVCF";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }
  #------------------------------------------------------------------------------------

  printtime();
  print STDERR "####### runlevel $runlevels done #######\n\n";

}



###------------###################################################################################################
##  sub-region  #################################################################################################
###------------###################################################################################################

sub RunCommand {
  my ($command,$noexecute,$quiet) = @_ ;
  unless ($quiet){
    printtime();
    print STDERR "$command\n\n";
  }
  unless ($noexecute) {
    system($command);
  }
}

sub helpm {
  print STDERR "\nGENERAL OPTIONS:\n\t--runlevel\tthe steps of runlevel, from 1-3, either rl1-rl2 or rl. See below for options for each runlevel.\n\tor\n";
  print STDERR "\t--runTask\tthe specific task. e.g., \'QC\', \'indelRealignment\', \'MarkDuplicates\', etc\n\n";
  print STDERR "\t--configure\tthe tab delimited file containing conf info for annotations\n";
  print STDERR "\t--sampleName\tthe name of the lane needed to be processed (must set for runlevel 1-5)\n";
  print STDERR "\t--seqType\tcomma separated, possible arguments \'paired-end\', \'single-end\', \'WXS\' and \'WGS\' (default).\n";
  print STDERR "\t--root\t\tthe root directory of the pipeline (default is \$bin/../PIPELINE/, MUST set using other dir)\n";
  print STDERR "\t--species\tspecify the reference version of the species, such as hg19 (default), mm10.\n";
  print STDERR "\t--patient\tthe patient id, which will be written into the target file for edgeR\n";
  print STDERR "\t--tissue\tthe tissue type name (like \'normal\', \'cancer\'), for the target file for running edgeR and cuffdiff\n";
  print STDERR "\t--Rbinary\tthe name of R executable, default is \'R\'. Set if your R binary name is different.\n\n";

  #print STDERR "CONTROL OPTIONS FOR EACH RUNLEVEL:\n";
  #print STDERR "runlevel 1: quality checking and insert size estimatiion using part of reads\n";
  print STDERR "\t--readpool\tthe directory where all the read files with names ending with \.f(ast)?q\.[gb]z2\? located.\n";
  print STDERR "\t--qcOFF\t\tturn off quality check.\n";
  print STDERR "\t--FASTQ1\t\tcomma separated zipped fastq file names for mate 1 (without dir name).\n";
  print STDERR "\t--FASTQ2\t\tcomma separated zipped fastq file names for mate 2 (without dir name).\n";
  print STDERR "\t--fastqFiles1\tcomma separated zipped fastq file names for mate 1 (with dir name).\n";
  print STDERR "\t--fastqFiles2\tcomma separated zipped fastq file names for mate 2 (with dir name).\n";
  print STDERR "\t--readlen\tthe sequenced read length (default the length of the first read in the fastq file)\n";
  print STDERR "\t--bamID\t\tthe ID for read group record in bam file, default is 1. set to the id needed to differentiate seq-experiments\n";
  print STDERR "\t--splitChr\tsplit Chromosomes. (default not set)\n";

  print STDERR "\nrunlevel 2: mapping and report of mapping statistics\n";
  print STDERR "\t--mapper\tthe mapper for read-alignment, now support \'bwa\' (default).\n";
  print STDERR "\t--gf\t\tthe graphical format in mapping report, \'png\' (default) or \'pdf\' (when a x11 window is not available)\n";

  print STDERR "\nrunlevel 3: STATS\n";

  print STDERR "\nrunlevel 4: SNV calling (muTect)\n";
  print STDERR "\t--somaticInfo\tsample information for tumor and normal pair (tab delimited)\n";

  print STDERR "\nOTHER OPTIONS\n";
  print STDERR "\t--noexecute\tdo not execute the command, for testing purpose\n";
  print STDERR "\t--quiet\t\tdo not print the command line calls and time information\n";
  print STDERR "\t--threads\tthe number of threads used for the mapping (default 1)\n";
  print STDERR "\t--help\t\tprint this help message\n";
  print STDERR "\t--tmpDir\ttmp dir for generating large tmp files\n";
  print STDERR "\t--bzip\t\tthe read fastq files are bziped, rather than gziped (default).\n";
  print STDERR "\t--bin\t\tbinary directories (default is where this script is executed.').\n\n";

  #print STDERR "\nSynopsis: RTrace.pl --runlevel 1 --sampleName <sample1> --runID <ID> --root <dir_root> --anno <dir_anno> 2>>run.log\n";
  #print STDERR "Synopsis: RTrace.pl --runlevel 2 --sampleName <sample1> --runID <ID> --root <dir_root> --anno <dir_anno> --patient <ID> --tissue <type> --threads <N> 2>>run.log\n";
  #print STDERR "Synopsis: RTrace.pl --runlevel 3-4 --sampleName <sample1> --runID <ID> --root <dir_root> --anno <dir_anno> --RA 1 --threads <N> 2>>run.log\n";
  #print STDERR "Synopsis: RTrace.pl --runlevel 5 --sampleName <sample1> --runID <ID> --root <dir_root> --anno <dir_anno> --threads <N> 2>>run.log\n";
  #print STDERR "Synopsis: RTrace.pl --runlevel 7 --root <dir_root> --anno <dir_anno> --priordf 1 2>>run.log\n";
  #print STDERR "remember to set --species option to a string, such as mm10, if it is not a human sample!!!\n\n";

  #print STDERR "Runlevel dependencies (->): 4->3->2->1, 6->5->2->1, 7->2->1\n\n";
  exit 0;
}

sub printtime {
  my @time = localtime(time);
  printf STDERR "\n[".($time[5]+1900)."\/".($time[4]+1)."\/".$time[3]." ".$time[2].":".$time[1].":".$time[0]."]\t";
}

sub round {
    my $number = shift;
    my $tmp = int($number);
    if ($number >= ($tmp+0.5)){
      $tmp++;
    }
    return $tmp;
}

sub uniqueArray {
   my $array = shift;
   my %arraytmp;
   foreach my $item (@{$array}){
     $arraytmp{$item} = '';
   }
   my @arraytmp = keys %arraytmp;
   return @arraytmp;
}

sub ceiling {
  my ($num) = @_;
  return int($num) + ($num > int($num));
}


sub partitionArray {
    my ($arr, $N) = @_;

    my @res;
    my $i = 0;

    while ($i + $N-1 <= $#$arr) {
        push @res, [@$arr[$i .. $i+$N-1]];
        $i += $N;
    }

    if ($i <= $#$arr) {
        push @res, [@$arr[$i .. $#$arr]];
    }
    return \@res;
}

package bwaMapping;

use strict;

#
# BWA Mapping
#

sub bwaPairMapping {

  my ($class, $ReadGroup, $threads, $bwaindex, $outBam, $readfiles1, $readfiles2) = @_;

  my $cmd = "bwa mem -r 1.2 -t $threads -R \'$ReadGroup\' $bwaindex \'\<zcat $readfiles1\' \'\<zcat $readfiles2\' | samtools view -bS - >$outBam";

  return $cmd;

}

sub bwaSingleMapping {

  my ($class, $ReadGroup, $threads, $bwaindex, $outBam, $readfiles1) = @_;

  my $cmd = "bwa mem -r 1.2 -t $threads -R \'$ReadGroup\' $bwaindex \'\<zcat $readfiles1\' | samtools view -bS - >$outBam";

  return $cmd;

}

sub bamSort {

  my ($class, $threads, $bamTmp, $outBam, $inBam) = @_;

  my $cmd = "samtools sort -@ $threads -O bam -T $bamTmp -o $outBam $inBam";

  return $cmd;

}

sub bamIndex {

  my ($class, $inBam) = @_;

  my $cmd = "samtools index $inBam";

  return $cmd;
}

sub indelRealignment1 {

  my ($class, $gatkBin, $inBam, $gfasta, $knownindel1, $knownindel2, $CHR, $outList) = @_;

  my $cmd = "java -Xmx2g -jar $gatkBin -T RealignerTargetCreator -R $gfasta -I $inBam -known $knownindel1 -known $knownindel2 -L $CHR -o $outList";
  if ($CHR eq 'ALL') {
    $cmd = "java -Xmx2g -jar $gatkBin -T RealignerTargetCreator -R $gfasta -I $inBam -known $knownindel1 -known $knownindel2 -o $outList";
  }

  return $cmd;

}

sub indelRealignment2 {

  my ($class, $gatkBin, $inBam, $gfasta, $targetList, $knownindel1, $knownindel2, $CHR, $outBam) = @_;

  my $cmd = "java -Xmx2g -jar $gatkBin -T IndelRealigner -R $gfasta -I $inBam -targetIntervals $targetList -known $knownindel1 -known $knownindel2 -compress 5 -L $CHR -o $outBam";
  if ($CHR eq 'ALL') {
    $cmd = "java -Xmx2g -jar $gatkBin -T IndelRealigner -R $gfasta -I $inBam -targetIntervals $targetList -known $knownindel1 -known $knownindel2 -compress 5 -o $outBam";
  }

  return $cmd;

}

sub MarkDuplicates {

  my ($class, $MarkDuplicatesBin, $inBam, $outBam, $metric) = @_;

  my $cmd = "java -Xmx2g -jar $MarkDuplicatesBin I=$inBam O=$outBam M=$metric REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=LENIENT";

  return $cmd;

}



1;



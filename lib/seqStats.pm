package seqStats;

use strict;

#
# sequencing Stats
#

sub mappingStats {

  my ($class, $bamStatsBin, $BAM, $readlen, $mappingStatsOut) = @_;

  my $cmd = "$bamStatsBin --mapping $BAM --readlength $readlen --maxIntron 23000 --type multiMis >$mappingStatsOut";

  return $cmd;

}

sub grepStarts {

  my  ($class, $grepStartsBin, $targetRegion, $BAM, $bedCover) = @_;

  my $cmd = "$grepStartsBin --region $targetRegion --mapping $BAM >$bedCover";

  return $cmd;

}


sub getLorenz {

  my  ($class, $lorenzCurveBin, $bedCover, $lorenzCover) = @_;

  my $cmd = "perl $lorenzCurveBin $bedCover >$lorenzCover";

  return $cmd;

}

sub bed2wig {

  my ($class, $bed2wigBin, $bedCount, $wigOut) = @_;

  my $cmd = "perl $bed2wigBin $bedCount >$wigOut";

  return $cmd;

}


sub insertSize {

  my ($class, $BAM, $outINS) = @_;

  my $cmd = "samtools view -f 0x2 -F 0x400 $BAM | cut -f 9 | awk \'\$1\>0 \&\& \$1\<1000\' >$outINS";

  return $cmd;

}


sub plotInsertSize {

  my ($class, $insertSizeRbin, $path, $sampleName, $insGZ, $outPDF) = @_;

  my $cmd = "$insertSizeRbin $path $sampleName $insGZ $outPDF";

  return $cmd;

}


1;



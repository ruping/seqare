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


1;



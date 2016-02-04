package cnaCalling;

use strict;

#
# cna Calling
#

sub runTitan {

  my ($class, $titanRBin, $PATH, $sampleName, $alleleCount, $tumorWig, $normalWig, $gcWig, $mapWig, $exons) = @_;

  my $cmd = "$titanRBin $PATH $sampleName $alleleCount $tumorWig $normalWig $gcWig $mapWig $exons";

  return $cmd;

}


1;



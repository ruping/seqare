package cnaCalling;

use strict;

#
# cna Calling
#

sub runTitan {

  my ($class, $titanRBin, $PATH, $sampleName, $alleleCount, $tumorWig, $normalWig, $gcWig, $mapWig, $plp, $plpe, $nc, $ncm, $exons) = @_;

  my $cmd = "$titanRBin $PATH $sampleName $alleleCount $tumorWig $normalWig $gcWig $mapWig $plp, $plpe, $nc, $ncm, $exons";

  return $cmd;

}


1;



package cnaCalling;

use strict;

#
# cna Calling
#

sub runTitan {

  my ($class, $RscriptBin, $titanRBin, $PATH, $sampleName, $alleleCount, $tumorWig, $normalWig, $gcWig, $mapWig, $plp, $plpe, $nc, $ncm, $exons) = @_;

  my $cmd = "$RscriptBin, $titanRBin $PATH $sampleName $alleleCount $tumorWig $normalWig $gcWig $mapWig $plp $plpe $nc $ncm $exons";

  return $cmd;

}


1;



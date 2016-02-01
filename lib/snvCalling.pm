package snvCalling;

use strict;

#
# snv Calling
#

sub muTectCalling {

  my ($class, $muTectBin, $BAM, $NORMALBAM, $gfasta, $COSMIC, $DBSNP, $muTectOut, $vcfOut) = @_;

  $cmd = "java -Xmx2g -jar $muTectBin -rf BadCigar --analysis_type MuTect --reference_sequence $gfasta --cosmic $COSMIC --dbsnp $DBSNP --input_file:normal $NORMALBAM --input_file:tumor $BAM --enable_extended_output --out $muTectOut -vcf $vcfOut";

  return $cmd;

}


sub muTect2vcf {

  my ($class, $mutect2vcfBin, $inMutect, $outVCF) = @_;

  my $cmd = "perl $mutect2vcfBin $inMutect >$outVCF";

  return $cmd;

}


sub vcfSort {

  my ($class, $inVCF, $outVCF) = @_;

  my $cmd = "cat $inVCF | vcf-sort >$outVCF";

  return $cmd;

}


sub runAnnovar {

  my ($class, $annovarBin, $inVCF, $ANNOVARDB, $species, ) = @_;

  my $cmd = "perl $annovarBin $inVCF $ANNOVARDB -buildver $species -remove -protocol refGene,cytoBand,genomicSuperDups,phastConsElements46way,tfbsConsSites,gwasCatalog,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,exac03,popfreq_max_20150413,ljb26_all,clinvar_20150629 -operation g,r,r,r,r,r,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput --outfile $inVCF";

  return $cmd;
}


sub convertVCFannovar {

  my ($class, $convertVCFannovarBin, $vcfMultiAnno, $vcfMultiAnnoVCF) = @_;

  my $cmd = "perl $convertVCFannovarBin $vcfMultiAnno $vcfMultiAnnoVCF";

  return $cmd;

}



1;



use strict;

my $files = shift;
my $headern = shift;

if ($headern eq ''){
  print STDERR "no header number provided, now assume \^\# as default header lines.\n";
  $headern = -1;
}

my @files = split(',', $files);

my $fc = 1;
foreach my $file (sort @files) {
  if ($fc == 1) {
    open IN, "$file";
    while ( <IN> ) {
      print "$_";
    }
    close IN;
  } else {
    my $hc = 0;
    open IN, "$file";
    while ( <IN> ) {
      $hc ++;
      if ($headern >= 0) {
        if ($hc <= $headern){
          next;
        } else {
          print "$_" if ($_ ne '');
        }
      } elsif ($headern == -1) {
        if ($_ =~ /^#/){
          next;
        } else {
          print "$_" if ($_ ne '');
        }
      }
    }
    close IN;
  }
  $fc ++;
}

exit;

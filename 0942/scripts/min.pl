my $min = 999999.9999;
my $line = 0;
my $minline = 0;
my $linetxt;
while (<>) {
  $line++;
  chomp($_);
  my @line = split(/\s+/, $_);
  foreach my $item (@line) {
    if (($item =~ /^-?\d+\.\d*$/) && ($item < $min)) {
      $min = $item;
      $minline = $line;
      $linetxt = $_;
    }
  }
}
print "At line $minline: $min\n";
print "$linetxt\n";


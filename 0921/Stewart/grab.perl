#  This parses the output of phc.  Initially is make a long list with all
#  components of all solutions shoved in, and then it accesses them.
#
open(PHC,"Dietmaier.phc");
my $copy = 0;
my $nsols = -1;
my @sols = ();

open (ALPHA, ">", "RealSols.maple");

print ALPHA "# Maple file containing the real solutions to Dietmaier's system\n";
print ALPHA "Points:=[\n";

while (<PHC>) {  #parse each line
  my $line = $_;
  chomp($line);
  if ($line =~ m/solution/) {
    $nsols++;
    if ($copy==0){ print ALPHA " [\n";}
    $copy = 1;
  }
  if ($line =~ m/failure/) {
    $copy = 0;
  }
  if ($copy == 1) {
    if ($line =~ m/^\s(n|h).*:\s+([0-9\.E\-]+)\s+([0-9\.E\-]+)/ ) {
       print ALPHA " $2 + ($3*I),\n";
    }
    if ($line =~ m/^==/ ) {
      print ALPHA " NULL],\n";
      $copy = 0;
    }
  }
}

print ALPHA "NULL]:\n";

close (PHC);

close (ALPHA); 

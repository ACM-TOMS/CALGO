#  This parses the output of realDistinctSolns, creating a maple file.
#
open(DATA,"realDistinctSolns");
my $copy = 0;
my $nsols = -1;
my @sols = ();

open (ALPHA, ">", "RealSols.maple");

print ALPHA "# Maple file containing the real distinct solutions to the nine-point problem\n";
print ALPHA "Points:=[\n";
$iter = 0;

while (<DATA>) {  #parse each line
  my $line = $_;
  chomp($line);
  if ($line =~ m/([0-9\.e\-]+)/) {
    $iter++;
    if ($iter==1) {
      print ALPHA "[$1,\n";
    }
    if ($iter<12 and $iter>1) {
      print ALPHA " $1,\n";
    }
    if ($iter==12) {
      print ALPHA " $1],\n";
      $iter=0;
    }
  }
}

print ALPHA "NULL]:\n";

close (DATA);

close (ALPHA);

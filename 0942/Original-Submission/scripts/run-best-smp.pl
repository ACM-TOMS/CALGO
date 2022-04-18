#!/usr/bin/perl

my $params = "hwc/results.ps701n1.params.cut";
my $outfile = "runs/results.ps701n1.smp";

my @threads = (1, 2, 4, 8, 16, 32);


# Execute best configurations with SMP enabled
sub run_best_smp
{
  open(my $fhi, '<', $params) or die $!;
  open(my $fho, '>', $outfile) or die $!;
  while (<$fhi>) {
    $line++;
    chomp($_);

    # Run test
    my @line = split(/\s+/, $_);
    @line[0] = "@line[0].smp";

    print $fho "#@line\n";

    # FOR EACH THREAD SIZE OF THREADS GROUP
    foreach $thread (@threads) {
      print "Running: OMP_NUM_THREADS=$thread @line[0..8]\n";
      my @output = split (/\s+/, `OMP_NUM_THREADS=$thread @line[0..8]`);
      print $fho "$thread \t @output[9..11]\n";
    }

    print $fho "\n";
  }

  close($fhi);
  close($fho);
}

# Main function
&run_best_smp();


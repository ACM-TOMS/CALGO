#!/usr/bin/perl

#my $papiex = "papiex -q -n -e";

my $params = "hwc/results.casegpu.params";
my $outfile = "hwc/results.casegpu.hwc";

my @sets = ("PAPI_L1_DCA,PAPI_L1_DCM",
            "PAPI_L2_DCA,PAPI_L2_DCM",
            "PAPI_TOT_INS,PAPI_TOT_CYC",
            "PAPI_FP_OPS");


# Execute best configurations with hwc
sub run_best_hwc
{
  open(my $fhi, '<', $params) or die $!;
  open(my $fho, '>', $outfile) or die $!;
  while (<$fhi>) {
    $line++;
    chomp($_);

    # Run test
    my @line = split(/\s+/, $_);
    @line[0] = "@line[0].hwc";

    print $fho "#@line\n";

    # FOR EACH SET OF HWC
    foreach $set (@sets) {
      print "Running: PAPI_COUNTERS=$set @line[0..8]\n";
      my $output = `PAPI_COUNTERS=$set @line[0..8] 2>&1| grep PAPI`;
      print $fho "$output";
    }

    print $fho "\n";
  }

  close($fhi);
  close($fho);
}

# Main function
&run_best_hwc();


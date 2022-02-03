#!/usr/bin/perl

my $tmp = "tmp.hwc";
my $command = "likwid-perfctr -o $tmp -C 0 -g";

my $params = "hwc/results.juropa.params48";
my $outfile = "hwc/results.juropa.hwc48";

my @sets = ("FLOPS_DP",
            "CACHE",
            "L2CACHE",
            "L3CACHE",
            "MEM",
            "DATA",
            "VIEW");


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
      print "Running: Group=$set @line[0..8]\n";
      `$command $set ./bin/@line[0..8]`;
      my $output = `cat tmp.hwc`;
      print $fho "$output";
    }

    print $fho "\n";
  }

  close($fhi);
  close($fho);
}

# Main function
&run_best_hwc();


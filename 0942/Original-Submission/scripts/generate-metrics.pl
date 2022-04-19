#!/usr/bin/perl
# Generate HWC metrics for further plotting

my @plat   = ("IBM Power6", "Intel Nehalem", "AMD Opteron", "IBM BlueGene/P", "Intel Nehalem", "IBM Power7");
my @pname  = ("huygens", "juropa", "louhi", "jugene", "casegpu", "ps701n1");
my @size   = ("512", "512", "512", "256", "512", "512");

my @algo   = ("Naive", "Naive+Semi", "Rivera", "Rivera+Semi",
              "Time-skewing", "Time-skewing+Semi", "Cache oblivious", "Cache oblivious+Semi");

my @mname  = ("cyc", "fpmem", "oi", "bw", "ipc", "l1", "l2", "l3", "gflops", "traffic", "mps", "time");
my @metric = ("Cycles * 10^9", "Ratio FP/Mem", "Operational Intensity", "Useful bandwidth (in GB/sec)", "IPC",
              "L1 Misses * 10^6", "L2 Misses * 10^6", "L3 Misses * 10^6", "GFlops", "Memory traffic R+W (GBytes)", "Million points per second", "Time in seconds");
my @operat = ("PAPI_TOT_CYC * 10^-9", "PAPI_FP_OPS/PAPI_L1_DCA", "FP_OPS/((UNC_QMC_NORMAL_READS_ANY+UNC_QMC_WRITES)*64)", "(((N+2*stencil)^3+N^3)*8bytes*10^-9)/sec", "PAPI_TOT_INS/PAPI_TOT_CYC",
              "PAPI_L1_DCM * 10^-6", "PAPI_L2_DCM * 10^-6", "UNC_L3_MISS_ANY", "PAPI_FP_OPS*10^-9/sec", "Memory GBytes*10^-9", "N^3/sec", "sec");

# Set how many timesteps and lengths we have in .hwc files
my @timestep = (1, 2, 4, 8);
my @lengths  = (1, 2, 4, 7, 14);

# By default timestep = 4
my @gentstep = (1, 4, 8); # Add those timesteps to be generated

# Depending on the performance tool used to extract hwc
# values are in the first column or in the second one
# PAPI/PAPIEX values first (0), LIKWID values second column (1)
my $colvals = 0;


# Compute metrics for files
sub cmetric
{
  my $val;
  my ($al, $ts, $platt, $set, $vals) = @_;
  print "Set for $al with timestep $ts in $platt platform:\n";
  print map { "$_ => $set->{$_}\n" } keys %{$set};
  print "\n";

  my @line = split(/\s+/, $set->{'run'});
  my $size = $line[1];
  my $tstep = $line[7];
  my $length = $line[8];
  my $time = pop(@line);

  # For Naive or Rivera multiply by number of timesteps
  # These runs have not been run with several steps
  my $delta = $ts/$tstep;

  # Multiply by delta those metrics which may be affected

  # Cycles
  if ($platt eq "juropa" && defined($set->{"CPU_CLK_UNHALTED_CORE"})) {$val = $set->{"CPU_CLK_UNHALTED_CORE"} * 10**-9 * $delta;}
  elsif ($platt eq "louhi" && defined($set->{"CPU_CLOCKS_UNHALTED"})) {
    $val = $set->{"CPU_CLOCKS_UNHALTED"} * 10**-9 * $delta;
  }
  else                    {$val = $set->{"PAPI_TOT_CYC"} * 10**-9 * $delta;}
  push(@{$vals->{'cyc'}}, $val);

  # FPMEM (fpmem)
  # Synthetic flops are: N^3 * 2 (Multiply-Add) * 2 (directions) * 3 (axis) * length (stencil order) + 1 (central point)
  # Real flops are: $set->{"PAPI_FP_OPS"}/$set->{"PAPI_L1_DCA"};
  #                 $set->{"PAPI_FP_OPS"}/$set->{"PAPI_LST_INS"};
  #if ($set->{"PAPI_LST_INS"} != 0) { $val = (($size**3)*(2*2*3*$length+1)*$tstep)/$set->{"PAPI_LST_INS"}; }
  #else {                             $val = (($size**3)*(2*2*3*$length+1)*$tstep)/$set->{"PAPI_L1_DCA"};  }
  if ($platt eq "juropa" && defined($set->{"MEM_INST_RETIRED_LOADS"}) && defined($set->{"MEM_INST_RETIRED_STORES"})) {
    $val = (($size**3)*(2*2*3*$length+1)*$tstep)/($set->{"MEM_INST_RETIRED_LOADS"}+$set->{"MEM_INST_RETIRED_STORES"});
  }
  elsif ($platt eq "louhi" && defined($set->{"DATA_CACHE_ACCESSES"})) {
    $val = (($size**3)*(2*2*3*$length+1)*$tstep)/($set->{"DATA_CACHE_ACCESSES"});
  }
  else {$val = (($size**3)*(2*2*3*$length+1)*$tstep)/$set->{"PAPI_L1_DCA"};}
  push(@{$vals->{'fpmem'}}, $val);

  # Operational Intensity (oi)
  #if    ($platt eq "juropa") {$val = $set->{"FP_COMP_OPS_EXE_SSE_DOUBLE_PRECISION"}/(($set->{"UNC_QMC_NORMAL_READS_ANY"}+$set->{"UNC_QMC_WRITES_FULL_ANY"})*64);}
  if ($platt eq "juropa" && defined($set->{"FP_COMP_OPS_EXE_SSE_DOUBLE_PRECISION"}) && defined($set->{"UNC_L3_MISS_ANY"})) {
    $val = $set->{"FP_COMP_OPS_EXE_SSE_DOUBLE_PRECISION"}/(($set->{"UNC_L3_MISS_ANY"})*64);
  }
  elsif ($platt eq "louhi" && defined($set->{"SSE_RETIRED_ADD_DOUBLE_FLOPS"}) && defined($set->{"SSE_RETIRED_MULT_DOUBLE_FLOPS"})) {
    $val = ($set->{"SSE_RETIRED_ADD_DOUBLE_FLOPS"}+$set->{"SSE_RETIRED_MULT_DOUBLE_FLOPS"})/
                     (($set->{"CPU_TO_DRAM_LOCAL_TO_0"}+$set->{"CPU_TO_DRAM_LOCAL_TO_1"}+$set->{"CPU_TO_DRAM_LOCAL_TO_2"}+$set->{"CPU_TO_DRAM_LOCAL_TO_3"})*64);
  }
  #elsif ($platt eq "huygens" && defined($set->{"PM_L3SA_MISS"}) && defined($set->{"PM_L3SB_MISS"})) {
  elsif ($platt eq "huygens" && defined($set->{"PM_L2_MISS"})) {
    #$val = $set->{"PAPI_FP_OPS"}/(($set->{"PM_L3SA_MISS"}+$set->{"PM_L3SB_MISS"})*128);
    $val = $set->{"PAPI_FP_OPS"}/($set->{"PM_L2_MISS"}*128);
  }
  else {$val = 0.0;}
  push(@{$vals->{'oi'}}, $val);

  # Bandwidth (bw)
  $val = ((($size+2*$length)**3+$size**3)*8*10**-9*$tstep)/($time);
  push(@{$vals->{'bw'}}, $val);

  # IPC (ipc)
  if ($platt eq "juropa" && defined($set->{"INSTR_RETIRED_ANY"}) && defined($set->{"CPU_CLK_UNHALTED_CORE"})) {
    $val = $set->{"INSTR_RETIRED_ANY"}/($set->{"CPU_CLK_UNHALTED_CORE"});
  }
  elsif ($platt eq "louhi" && defined($set->{"INSTRUCTIONS_RETIRED"}) && defined($set->{"CPU_CLOCKS_UNHALTED"})) {
    $val = $set->{"INSTRUCTIONS_RETIRED"}/($set->{"CPU_CLOCKS_UNHALTED"});
  }
  elsif ($platt eq "jugene" && defined($set->{"PAPI_FP_INS"}) && defined($set->{"PAPI_LD_INS"}) && defined($set->{"PAPI_SR_INS"})) {
    $val = ($set->{"PAPI_FP_INS"}+$set->{"PAPI_LD_INS"}+$set->{"PAPI_SR_INS"})/$set->{"PAPI_TOT_CYC"};
  }
  else {$val = $set->{"PAPI_TOT_INS"}/$set->{"PAPI_TOT_CYC"};}
  push(@{$vals->{'ipc'}}, $val);

  # L1 misses (l1)
  if ($platt eq "juropa" && defined($set->{"L1D_REPL"})) {
    $val = $set->{"L1D_REPL"} * 10**-6 * $delta;
  }
  elsif ($platt eq "louhi" && defined($set->{"DATA_CACHE_REFILLS_L2_ALL"}) && defined($set->{"DATA_CACHE_REFILLS_NORTHBRIDGE_ALL"})) {
    $val = ($set->{"DATA_CACHE_REFILLS_L2_ALL"}+$set->{"DATA_CACHE_REFILLS_NORTHBRIDGE_ALL"}) * 10**-6 * $delta;
  }
  elsif ($platt eq "huygens" && defined($set->{"PM_LD_MISS_L1"}) && defined($set->{"PM_ST_MISS_L1"})) {
    $val = ($set->{"PM_LD_MISS_L1"}+$set->{"PM_ST_MISS_L1"}) * 10**-6 * $delta;
  }
  else {$val = $set->{"PAPI_L1_DCM"} * 10**-6 * $delta;}
  push(@{$vals->{'l1'}}, $val);

  # L2 misses (l2)
  if ($platt eq "juropa" && defined($set->{"L2_RQSTS_MISS"})) {
    $val = $set->{"L2_RQSTS_MISS"} * 10**-6 * $delta;
  }
  elsif ($platt eq "louhi" && defined($set->{"L2_MISSES_ALL"})) {
    $val = $set->{"L2_MISSES_ALL"} * 10**-6 * $delta;
  }
  elsif ($platt eq "huygens" && defined($set->{"PM_L2_MISS"})) {
    $val = $set->{"PM_L2_MISS"} * 10**-6 * $delta;
  }
  else {$val = $set->{"PAPI_L2_DCM"} * 10**-6 * $delta;}
  push(@{$vals->{'l2'}}, $val);

  # L3 misses (l3)
  if ($platt eq "juropa" && defined($set->{"UNC_L3_MISS_ANY"})) {
    $val = $set->{"UNC_L3_MISS_ANY"} * 10**-6 * $delta;
  }
  elsif ($platt eq "louhi" && defined($set->{"L3_MISSES_ALL_ALL_CORES"})) {
    $val = $set->{"L3_MISSES_ALL_ALL_CORES"} * 10**-6 * $delta;
  }
  elsif ($platt eq "huygens" && defined($set->{"PM_L3SA_MISS"})) {
    #$val = ($set->{"PM_L3SA_MISS"}+$set->{"PM_L3SB_MISS"}) * 10**-6 * $delta;
    $val = ($set->{"PM_L3SA_MISS"}) * 10**-6 * $delta;
  }
  elsif (defined($set->{"PAPI_L3_TCM"})) {
    $val = $set->{"PAPI_L3_TCM"} * 10**-6 * $delta;
  }
  else {$val = $set->{"PAPI_L3_DCM"} * 10**-6 * $delta;}
  push(@{$vals->{'l3'}}, $val);

  # GFlops (gflops)
  $val = (($size**3)*(2*2*3*$length+1)*$tstep*10**-9)/($time);
  push(@{$vals->{'gflops'}}, $val);

  # Memory traffic (traffic)
  if ($platt eq "juropa" && defined($set->{"UNC_QMC_NORMAL_READS_ANY"}) && defined($set->{"UNC_QMC_WRITES_FULL_ANY"})) {
    $val = (($set->{"UNC_QMC_NORMAL_READS_ANY"}+$set->{"UNC_QMC_WRITES_FULL_ANY"})*64) * 10**-9 * $delta;
  }
  elsif ($platt eq "louhi" && defined($set->{"CPU_TO_DRAM_LOCAL_TO_0"})) {
    $val = (($set->{"CPU_TO_DRAM_LOCAL_TO_0"}+$set->{"CPU_TO_DRAM_LOCAL_TO_1"}+$set->{"CPU_TO_DRAM_LOCAL_TO_2"}+$set->{"CPU_TO_DRAM_LOCAL_TO_3"})*64) * 10**-9 * $delta;
  }
  else {$val = 0.0;} #($set->{"PAPI_L1_DCA"}*8)/(($size**3)*$tstep);}
  push(@{$vals->{'traffic'}}, $val);

  # Million points per second (mps)
  $val = (($size**3)*$tstep*10**-6)/$time;
  push(@{$vals->{'mps'}}, $val);

  # Time in seconds (time)
  $val = $time*$delta;
  push(@{$vals->{'time'}}, $val);

  return $vals;
}


# Build platforms table for every metric with best execution on .hwc
sub plotmetrics {

  my $indplat = 0;

  # FOR EACH PLATFORM
  foreach $p (@pname) {

    # Open input file
    print "Scanning results.$p.hwc metric: $m\n";
    open(my $fhi, '<', "results.$p.hwc") or die $!;

    # FOR EACH TIMESTEP
    foreach $t (@timestep) {

      my $indmet = 0;
      my %handles;

      # FOR EACH METRIC
      foreach $m (@mname) {

        # If $t timestep must be generated then create output file
        if (grep {$_ eq $t} @gentstep) {
          # Open output file
          open(my $fho, '>', "results.$p.t$t.$m") or die $!;
          $handles{$m} = $fho;

          print {$handles{$m}} "# Platform: $p\n";
          print {$handles{$m}} "# Timesteps: $t\n";
          print {$handles{$m}} "# Size: $size[$indplat]\n";
          print {$handles{$m}} "# Metric: $m ($metric[$indmet] - $operat[$indmet])\n";
          my @title = map {sprintf("%-20s", $_)} map {qq|"$_"|} @algo;
          print {$handles{$m}} sprintf("%-20s", "Algorithms") . "@title\n";

          $indmet++;
        }
      }

      # FOR EACH LENGTH
      foreach $l (@lengths) {

        my %set;
        my %vals;
        %vals = ();
        my @keys = keys %vals;
        while (scalar @{$vals{$keys[0]}} != scalar @algo) {

          $_ = <$fhi>;

          chomp($_);

          # New set
          if (my ($found) = m/(#)/o) {
#           print "Found: $found\n";
            %set = ();
            $set{'run'} = $_;
          }
          # End of the set
          elsif (length $_ == 0) {
            # Generate metrics for current platform, timestep and algorithm
            &cmetric($algo[scalar @{$vals{$keys[0]}}], $t, $p, \%set, \%vals);

            # If parsed all algorithmic metric values and the current timestep must be generated dump into file
            @keys = keys %vals;
            if ((scalar @{$vals{$keys[0]}} == scalar @algo) && (grep {$_ eq $t} @gentstep)) {

              print "Metrics for timestep $t in $p platform:\n";
              print map { "$_ => @{$vals{$_}}\n" } keys %vals;
              print "\n";

              foreach $m (@mname) {
#               print {$handles{$m}} sprintf("%-20s", qq|"$l"|), map {sprintf("%-21.6f", $_)} @{$vals{$m}};
                print {$handles{$m}} sprintf("%-20s", $l), map {sprintf("%-21.6f", $_)} @{$vals{$m}};
                print {$handles{$m}} "\n";
              }
            }
          }
          # Continue reading set
          else {
            # Get PAPI code and value
            my @line = split(/\s+/, $_);
            if ($colvals == 0) {$set{$line[1]} = $line[0];}
            else               {$set{$line[0]} = $line[1];}
          }
        }
      }

      # Close all opened files
      foreach $m (@mname) {
        close($handles{$m});
      }
    }

    $indplat++;
    close($fhi);
  }
}


# Main function
&plotmetrics( );


#!/usr/bin/perl
# Generate scalability metrics for further plotting

my @plat   = ("IBM Power7");
my @pname  = ("ps701n1");
my @size   = ("512x1024x1024");

my @algo   = ("Naive", "Naive+Semi",
              "Rivera", "Rivera+Semi",
              "Time-skewing", "Time-skewing+Semi",
              "Cache oblivious", "Cache oblivious+Semi");

#my @algo   = ("Naive Time", "Naive Speedup", "Naive+Semi Time", "Naive+Semi Speedup",
#              "Rivera Time", "Rivera Speedup", "Rivera+Semi Time", "Rivera+Semi Speedup",
#              "Time-skewing Time", "Time-skewing Speedup", "Time-skewing+Semi Time", "Time-skewing+Semi Speedup",
#              "Cache oblivious Time", "Cache oblivious Speedup", "Cache oblivious+Semi Time", "Cache oblivious+Semi Speedup");

# Set how many timesteps and lengths we have in .hwc files
my @timestep = (2, 4); #(1, 2, 4, 8);
#my @lengths  = ("Length 1", "Length 2", "Length 4", "Length 7", "Length 14");
my @lengths  = (1, 2, 4, 7, 14);

# By default timestep = 4
my @gentstep = (2, 4); # Add those timesteps to be generated


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
  # These runs have been not run with several steps
  my $delta = $ts/$tstep;

  # Generate lines for each number of cores
  foreach $c (keys %{$set}) {
    if ($c != 'run') {
      if (length $set->{$c} == 0) {
        print "WARNING!!!: Hash $c void ($set->{$c}) set to 1\n";
        $set->{$c} = 1;
      }
      # Push time and speedup
      $val = $set->{$c}*$delta;
      push(@{$vals->{$c}}, $val);
      $val = $set->{'1'}/$set->{$c};
      push(@{$vals->{$c}}, $val);
    }
  }

  return $vals;
}


# Build platforms table for every metric with best execution on .hwc
sub plotscala {

  my $indplat = 0;

  # FOR EACH PLATFORM
  foreach $p (@pname) {

    # Open input file
    print "Scanning results.$p.smp metric: $m\n";
    open(my $fhi, '<', "results.$p.smp") or die $!;

    # FOR EACH TIMESTEP
    foreach $t (@timestep) {

      my $indlen = 0;
      my %handles;

      # FOR EACH LENGTH
      foreach $l (@lengths) {

        # If $t timestep must be generated then create output file
        if (grep {$_ eq $t} @gentstep) {
          # Open output file
          open(my $fho, '>', "results.$p.t$t.l$l.smp.dat") or die $!;
          $handles{$l} = $fho;

          print {$handles{$l}} "# Platform: $p\n";
          print {$handles{$l}} "# Timesteps: $t\n";
          print {$handles{$l}} "# Length: $l\n";
          print {$handles{$l}} "# Size: $size[$indplat]\n";
          print {$handles{$l}} "# Metric: Scalability - Time\n";
          my @title = map {sprintf("%-19s", $_)} map {qq|"$_"|} @algo;
          print {$handles{$l}} sprintf("%-10s", "Algos") . "@title\n";

          $indlen++;
        }
      }

      # FOR EACH LENGTH
      foreach $l (@lengths) {

        my %set;
        my %vals;
        %vals = ();
        my @keys = keys %vals;
        while (scalar @{$vals{'1'}} != (scalar @algo)*2) {

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
            # Generate metrics for current platform, timestep, length and algorithm
            &cmetric($algo[(scalar @{$vals{'1'}})/2], $t, $p, \%set, \%vals);

            # If parsed all algorithmic metric values and the current timestep must be generated then dump it into file
            @keys = keys %vals;
            if ((scalar @{$vals{'1'}} == (scalar @algo)*2) && (grep {$_ eq $t} @gentstep)) {

              print "Metrics for timestep $t length $l in $p platform:\n";
              print map { "$_ => @{$vals{$_}}\n" } keys %vals;
              print "\n";

              foreach $c (sort { $a <=> $b } keys %vals) {
#               print {$handles{$l}} sprintf("%-10s", qq|"$l"|), map {sprintf("%-21.6f", $_)} @{$vals{$m}};
                print {$handles{$l}} sprintf("%-10s", $c), map {sprintf("%-10.5f", $_)} @{$vals{$c}};
                print {sprintf("%-10.5f", $_)} @{$vals{$c}};
                print {$handles{$l}} "\n";
              }
            }
          }
          # Continue reading set
          else {
            # Get number of cores and avg time (min, max, avg)
            my @line = split(/\s+/, $_);
            $set{$line[0]} = $line[3];
          }
        }
      }

      # Close all opened files
      foreach $l (@lengths) {
        close($handles{$l});
      }
    }

    $indplat++;
    close($fhi);
  }
}


# Main function
&plotscala( );


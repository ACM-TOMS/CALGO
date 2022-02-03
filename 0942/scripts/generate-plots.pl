#!/usr/bin/perl

my @dirs = ("huygens-runs", "juropa-runs", "louhi-runs", "jugene-runs", "casegpu-runs", "curie-runs", "ps701n1-runs");
my @plat = ("IBM Power6", "Intel Nehalem", "AMD Opteron", "IBM BlueGene/P", "Intel Nehalem", "Intel Nehalem-EX", "IBM Power7");
my @pname = ("huygens", "juropa", "louhi", "jugene", "casegpu", "curie", "ps701n1");
my @size = ("512", "512", "512", "256", "512", "512", "512");
my @algo = ("Naive", "Naive+Semi", "Rivera", "Rivera+Semi",
            "Time-skewing", "Time-skewing+Semi", "Cache oblivious", "Cache oblivious+Semi");
my @timesteps = (1, 2, 4, 8);
my @lengths   = (1, 2, 4, 7, 14);

%algo_map = ( 
          "Naive" => "probe",
          "Naive+Semi" => "semi_probe",
          "Rivera" => "blocked_probe",
          "Rivera+Semi" => "blocked_semi_probe",
          "Time-skewing" => "timeskew_probe",
          "Time-skewing+Semi" => "timeskew_semi_probe",
          "Cache oblivious" => "oblivious_probe",
          "Cache oblivious+Semi" => "oblivious_semi_probe"
);

# reading params
my @params = ( dirs, plat, pname, size, algo, timesteps, lengths );
open(my $fh, '<', $ARGV[0]) or die $!;
while (<$fh>) {
  chomp($_);
  my @line = split(/\s+/, $_);
  foreach $param (@params) {
    if( @line[0] =~ /^$param/) {
      for( $i = 2; $i <= $#line; $i++ ) {
        @$param[$i - 2] = @line[$i];
      }
    }
  }
}
close($fh);

my $indplat = 0;

# Compute minimum between several files
# Returns a vector composed of (value, linenum, filename, linetext)
# of the best performance file
sub min
{
  my $minval = 999999.9999;
  my $minline = 0;
  my $minfile;
  my $linetxt;
  foreach $file (@_) {
    my $line = 0;
    print "Scanning $file\n";
    open(my $fh, '<', $file) or die $!;
    while (<$fh>) {
      $line++;
      chomp($_);
      # Get last element
      my @line = split(/\s+/, $_);
      #@line = pop(@line);
      foreach my $item (@line) {
        # If best time for SIZExSIZExSIZE, it stores as best run
        if (($item =~ /^-?\d+\.\d*$/) && ($item < $minval) &&
            ((@line[1] == @size[$indplat]) && (@line[2] == @size[$indplat]) && (@line[3] == @size[$indplat]))) {
          $minval = $item;
          $minline = $line;
          $minfile = $file;
          $linetxt = $_;
        }
      }
    }
  }
  close($fh);
  print "Minimum at line $minline of $minfile: $minval\n";
  print "$linetxt\n\n";
  my @min = ($minval, $minline, $minfile, $linetxt);
}


# Build platforms table for a specific t and l with every optimization
sub plotall {

  # FOR EACH LENGTH
  foreach $l (@lengths) {

    # FOR EACH TIMESTEP
    foreach $t (@timesteps) {

      my @filters;
      foreach $a (@algo) {
        push ( @filters, "\\.$algo_map{ $a }.*\\.t$t.l$l\$" );
      }   

      print "dir: ",@dirs[0],"\n";
      open(my $fh, '>', @dirs[0]."/results.t$t.l$l.dat") or die $!;
      print $fh "# Platforms: @plat\n";
      print $fh "# Timesteps: $t, Length: $l\n";
      print $fh "# Sizes: @size\n";
      my @title = map {sprintf("%-20s", $_)} map {qq|"$_"|} @algo;
      print $fh sprintf("%-20s", "Algorithms") . "@title\n";

      $indplat = 0;

      # FOR EACH PLATFORM DIRECTORY
      foreach $dir (@dirs) {
 
        my $indfilter = 0;
 
        # FOR EACH FILTER ALGORITHM
        foreach $filter (@filters) 
        {
          print "Filter: $filter\n";
          opendir(my $dh, $dir) || die ("Cannot open directory $dir");
          my @files = grep { /$filter/ } readdir($dh);
          my $indfile = 0;
  
          # FOR EACH FILE IN DIRECTORY AND FILTER ADD PATH
          foreach my $file (@files) {
            @files[$indfile] = $dir."/".$file;
            $indfile++;
          }

          #print "@files" . "\n";
          my @minn = &min( @files );

          ## If non time-step execution multiply by timestep
          #if (($indfilter == 0) || ($indfilter == 1) || ($indfilter == 2) || ($indfilter == 3)) {
          #  @minn[0] = @minn[0] * $t;
          #}
          push(@m, sprintf("%-20s", @minn[0]));

          $indfilter++;
          closedir($dh);

        }
        print $fh sprintf("%-20s", map {qq|"$_"|} @plat[$indplat]) . "@m \n";
        @m = ();

        $indplat++;
      }

      close($fh);
    }
  }
}


# Build a platform table for a specific t with every optimization
sub plotplat {

  $indplat = 0;

  # FOR EACH PLATFORM DIRECTORY
  foreach $dir (@dirs) {

    # File that contains the best executions parameters
    open(my $fhline, '>', @dirs[0]."/results.@pname[$indplat].params") or die $!;

    # FOR EACH TIMESTEP
    foreach $t (@timesteps) {

      open(my $fh, '>', @dirs[0]."/results.@pname[$indplat].t$t.dat") or die $!;
      print $fh "# Platform: @plat[$indplat]\n";
      print $fh "# Timestep: $t\n";
      print $fh "# Size: @size[$indplat]\n";
      my @title = map {sprintf("%-20s", $_)} map {qq|"$_"|} @algo;
      print $fh sprintf("%-20s", "Algorithms") . "@title\n";

      # FOR EACH LENGTH
      foreach $l (@lengths) {

        my @filters;
        foreach $a (@algo) {
          push ( @filters, "\\.$algo_map{ $a }.*\\.t$t.l$l\$" );
        }   

        my $indfilter = 0;
    
        # FOR EACH FILTER ALGORITHM
        foreach $filter (@filters) 
        {
          print "Filter: $filter\n";
          opendir(my $dh, $dir) || die ("Cannot open directory $dir");
          my @files = grep { /$filter/ } readdir($dh);
          my $indfile = 0;
    
          # FOR EACH FILE IN DIRECTORY AND FILTER ADD PATH
          foreach my $file (@files) {
            @files[$indfile] = $dir."/".$file;
            $indfile++;
          }

          #print "@files" . "\n";
          my @minn = &min( @files );

          ## If non time-step execution multiply by timestep
          #if (($indfilter == 0) || ($indfilter == 1) || ($indfilter == 2) || ($indfilter == 3)) {
          #  @minn[0] = @minn[0] * $t;
          #}
          push(@m, sprintf("%-20s", @minn[0]));

          # Push the execution line into the best executions
          print $fhline @minn[3] . "\n";

          $indfilter++;
          closedir($dh);
        }

        print $fh sprintf("%-20s", map {qq|"$_"|} "length $l") . "@m \n";
        @m = ();

      }

      close( $fh );
    }

    close ( $fhline );

    $indplat++;
  }
}


# Build a blocking platform table for a specific t and l with every optimization
sub plottile {

  $indplat = 0;

  # FOR EACH PLATFORM DIRECTORY
  foreach $dir (@dirs) {

    # FOR EACH TIMESTEP
    foreach $t (@timesteps) {

      # FOR EACH LENGTH
      foreach $l (@lengths) {
      
        my @filters;
        foreach $a (@algo) {
          push ( @filters, "\\.$algo_map{ $a }.*\\.t$t.l$l\$" );
        }   

        my $indfilter = 0;
        my @minfiles;

        # FOR EACH FILTER ALGORITHM (GET THE BETTER FILES)
        foreach $filter (@filters) 
        {
          print "Filter: $filter\n";
          opendir(my $dh, $dir) || die ("Cannot open directory $dir");
          my @files = grep { /$filter/ } readdir($dh);
          my $indfile = 0;
    
          # FOR EACH FILE IN DIRECTORY AND FILTER ADD PATH
          foreach my $file (@files) {
            @files[$indfile] = $dir."/".$file;
            $indfile++;
          }
 
          #print "@files" . "\n";
          my @minn = &min( @files );
          push( @minfiles, [@minn] );

          $indfilter++;
          closedir($dh);
        }

        my $indfilter = 0;
        my @data;
        my @blocks;

        # FOR EACH FILTER ALGORITHM (GET DATA FROM BETTER FILES)
        foreach $filter (@filters)
        {
          print "Getting data from: $minfiles[$indfilter][2]\n";
          open(my $fh, '<', $minfiles[$indfilter][2]) or die $!;

          my @vals;
          while (<$fh>) {
            $line++;
            chomp($_);
#           $_ =~ s/^\s+//; #remove leading spaces
#           $_ =~ s/\s+$//; #remove trailing spaces
            # Get last element
            my @line = split(/\s+/, $_);
            my $item = pop(@line);
#            print "Processing Item: $item\n";
            if ($item =~ /^-?\d+\.\d*$/) {
              push(@vals, $item);
              if ($indfilter == 2) {
                push( @blocks, "@line[4]x@line[5]x@line[6]" );
              }
            }
          }
          push( @data, [@vals] );

          $indfilter++;
          close( $fh );
        }

        # FOR EACH FILTER ALGORITHM (GENERATE OUTPUT MERGED)
        open(my $fh, '>', @dirs[0]."/results.tile.@pname[$indplat].t$t.l$l.dat") or die $!;
        print $fh "# Platform: @plat[$indplat]\n";
        print $fh "# Timestep: $t, Length: $l\n";
        print $fh "# Size: @size[$indplat]\n";
        my @title = map {sprintf("%-19s", $_)} map {qq|"$_"|} @algo;
        print $fh sprintf("%-20s", "Algorithms") . "@title\n";

        my $val;
        my $indelem = 0;

        # FOR EACH ELEMENT IN RIVERA AND TIMESKEW RUNS GENERATE OUTPUT
        foreach $elem (@{$data[2]}) {

#          print "Elem: $elem\n";
          print $fh sprintf("%-20s", map {qq|"$_"|} @blocks[$indelem]);

          my $indfilter = 0;
          foreach $filter (@filters)
          {
            # If non time-step execution multiply by timestep
            # Naive
            if (($indfilter == 0) || ($indfilter == 1)) {
              $val = ${$minfiles[$indfilter]}[0]; # * $t;
            }
            # Rivera
            elsif (($indfilter == 2) || ($indfilter == 3)) {
              $val = ${$data[$indfilter]}[$indelem]; # * $t;
            }
            # Timeskew
            elsif (($indfilter == 4) || ($indfilter == 5)) {
              $val = ${$data[$indfilter]}[$indelem];
            }
            # Cache-oblivious
            elsif (($indfilter == 6) || ($indfilter == 7)) {
              $val = ${$minfiles[$indfilter]}[0];
            }
            print $fh sprintf("%-20s", $val);

            $indfilter++;
          }
          print $fh "\n";

          $indelem++;
        }
        close( $fh );
      }
    }
    $indplat++;
  }
}



# Main function
&plotall( );
&plotplat( );
&plottile( );


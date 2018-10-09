#!/usr/bin/perl

=pod

=head1 NAME

 mintia_assembly.pl

=head1 SYNOPSIS

 mintia_assembly.pl --input FILE[S] --vectorSeq FILE --dirOutputs STR

=head1 DESCRIPTION

 This tool is the first step of the fosmid assembly and annotation pipeline.
 It assembles raw reads, looks for and removes the cloning vector, and extracts
 the  longest and  the most  covered contigs.  It has  been build to handle two
 types of raw  reads as inputs:  single (454, ion torrent reads, ...) or paired
 (Illumina,...) reads.
 This tools is not able to process PacBio or Oxford Nanopore reads.

=head1 OPTIONS

=over 8

=item B<-i, --input> FILE[S]

 Fastq(.gz) file(s)
 For each sample one OR two fastq file must be provided:
 - Paired data must contain R1/R2: R[12].f[ast]q[.gz]
   Ex: sampleName1_R1.fastq sampleName1_R2.fastq
       sampleName2_R1.fq.gz sampleName2_R2.fq.gz
 - Single data
   Ex: sampleName.f[ast]q[.gz]

=item B<-v, --vectorSeq> FILE

 Path to the vector fasta file

=item B<-l, --length> INT

 Fosmid's expected length [40000]

=item B<--minimalContigLength> INT

 Contig's minimum length [1000]

=item B<--minimalContigDepth> INT

 Contig's minimum depth [8]

=item B<-c, --maxDepth> INT

 Coverage, maximum depth use to filter input reads [300]

=item B<-d, --dirOutputs> STR

 Path to the outputs directory

=item B<-Z, --zipOutput> STR

 Zip output name [mintia_assembly.zip]

=item B<-H, --htmlOutput> STR

 HTML output name [mintia_assembly.html]

=item B<-L, --logOutput> STR

 Log output file name [mintia.log]

=item B<-t, --threads>

 number of threads for SPADES [8]

=item B<-h, --help>

 Print help

=back

=head1 VERSION

 Mintia_v0.2

=head1 AUTHORS

 Philippe   Bardou - INRA Toulouse - support.sigenae@inra.fr
 Christophe  Klopp - INRA Toulouse
 Sandrine Laguerre - INRA Toulouse
 Sarah       Maman - INRA Toulouse
 Sabrina Rodriguez - INRA Toulouse

=head1 COPYRIGHT

 2018 INRA

=head1 LICENSE

 GNU GPLv3

=cut

use strict;
use warnings;
use List::Util ;
use File::Spec ;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use Pod::Usage;
use Getopt::Long;

use POSIX;
use 5.010;  #for filesize

my $SPADES = "/usr/local/bioinfo/src/SPAdes/SPAdes-3.11.0-Linux/bin/spades.py";
my $CROSS_MATCH = "/usr/local/bioinfo/bin/cross_match";

=head2 function locate_regex

 Usage    : locate_regex( $regex, $string )
 Function : Find regex location. Check if there is a match and
            return the start and end location of match(s)
 Return   : [array] of start/stop
 Args     : [str] regex
            [str] string
=cut
sub locate_regex {
    my ($regex, $string) = @_;
    my @a_res;
    while($string =~ /$regex/g) {
        push(@a_res, [ $-[0], $+[0]-1 ]);
    }
    return @a_res;
}

=head2 function read_crossmatch

 Usage    : read_crossmatch( $file )
 Function : Parse crossmatch output file
 Return   : [array] of start/stop crossmatch hit(s)
 Args     : [str] crossmatch output file name

=cut
sub read_crossmatch {
  my ($file) = @_;
  my @a_res;
  open(FILE, "$file") || die "Error: Enabled to open $file\n";
  while(my $line=<FILE>) {
    if($line=~/^\s*\d+\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s+(\S+)\s+(\d+)\s+(\d+)/) {
      push(@a_res, [ $2, $3 ]);
    }
    elsif($line=~//) { }
  }
  close(FILE);
  return @a_res;
}

=head2 function run_crossmatch

 Usage    : run_crossmatch( $seq1, $seq2 )
 Function : Run crossmatch on seq1 VS seq2 and return new sequence
 Return   : [str] new sequence
 Args     : [str] seq1
            [str] seq2

=cut
sub run_crossmatch {
  my ($seq1, $seq2, $outputDir) = @_;

  my $tmp1 = File::Temp->new(TEMPLATE=>'tmp1XXXXXXX', DIR=>$outputDir, SUFFIX=>'.fa', UNLINK=>0);
  my $tmp2 = File::Temp->new(TEMPLATE=>'tmp2XXXXXXX', DIR=>$outputDir, SUFFIX=>'.fa', UNLINK=>0);
  my $resC = File::Temp->new(TEMPLATE=>'resCXXXXXXX', DIR=>$outputDir, SUFFIX=>'.fa', UNLINK=>0);

  open(TMP1, ">$tmp1") || die "Error: Enabled to create $tmp1\n";
  open(TMP2, ">$tmp2") || die "Error: Enabled to create $tmp2\n";
  print TMP1 ">seq1\n$seq1\n";
  print TMP2 ">seq2\n$seq2\n";
  close(TMP1);
  close(TMP2);

  # Processing the sequences
  `($CROSS_MATCH $tmp1 $tmp2 -minmatch 10 -minscore 20 > $resC) >& /dev/null`;
  my @a_res = read_crossmatch($resC);

  # Processing crossmatch result
  my $newcoord = 0;
  for(my $i=0; $i<=$#a_res; $i++) {
    # Check if the overlap is at the end of seq1 and at the start of seq2
	if ((length($seq1)-$a_res[$i][1]) < 10) {
      $newcoord = $a_res[$i][0];
    }
  }
  if($newcoord > 0) { return substr($seq1,0,length($seq1)-$newcoord+1).$seq2; }
  else              { return $seq1.$seq2; }
}

=head2 function remove_vector

 Usage    : remove_vector( $regex, $string )
 Function : Find regex location. Check if there is a match and
            return the start and end location of match(s)
 Return   : [str] new sequence without XXed block (vector)
 Args     : [str] regex
            [str] string

=cut
sub remove_vector {
  my ($ref_a_coord, $seq, $outputDir) = @_;

  # Processing the different vector sequence locations
  # If the vector sequence is close to the beginning or the end of the contig
  # then it is removed => (case 1 and/or case2)
  # If the vector sequence is in the middle of the contig the the contig is
  # modified (split and join) => (case4)
  my ($case, $start, $stop, $start1, $stop1) = (0, 0, 0, 0, 0);
  my $len = length($seq);
  for(my $i=0; $i<=$#{$ref_a_coord}; $i++) {
    # Case 1 XXed block at the beginning of the scaffold
    # Seq= ^[ATGC]{,10}XXXXX...XXXX[ATGC]+$
    #                              ^=startPos
    if($$ref_a_coord[$i][0]<10) {
      $start = $$ref_a_coord[$i][1]+1;
      $case += 1;
    }
    # Case 2 XXed block at the end of the scaffold
    # Seq= ^[ATGC]+XXXXX...XXXX[ATGC]{,10}$
    #             ^=stopPos
    elsif(($len-$$ref_a_coord[$i][1])<10) {
      $stop  = $$ref_a_coord[$i][0];
      $case += 2;
    }
    # Case 4 XXed block in the middle of the scaffold
    # Seq= ^[ATGC]{10,}XXXXX...XXXX[ATGC]{10,}$
    #                 ^=stopPos    ^=startPos
    else {
      $start1 = $$ref_a_coord[$i][1]+1;
      $stop1  = $$ref_a_coord[$i][0];
	  $case  += 4;
    }
    print LOG "AAA:$case\n";
  }
  # XXXXX----------------
  if   ($case == 1) { return substr($seq,$start,$len-$start+1); }
  # ----------------XXXXX
  elsif($case == 2) { return substr($seq,0,$stop);              }
  # XXXXX-----------XXXXX
  elsif($case == 3) { return substr($seq,$start,$stop-$start);  }
  # --------XXXXX--------
  elsif($case == 4) { return run_crossmatch(substr($seq,$start1,$len-$start1+1),
                                            substr($seq,      0,$stop1),
                                            $outputDir); }
  # XXXXX-----XXXXX------
  elsif($case == 5) { return run_crossmatch(substr($seq,$start1,$len-$start1+1),
                                            substr($seq, $start,$stop1-$start),
                                            $outputDir); }
  # ------XXXXX-----XXXXX
  elsif($case == 6) { return run_crossmatch(substr($seq,$start1,$stop-$start),
                                            substr($seq,      0,      $stop1),
                                            $outputDir); }
  # XXXXX---XXXXX---XXXXX
  elsif($case == 7) { return run_crossmatch(substr($seq,$start1,$stop-$start),
                                            substr($seq, $start,$stop1-$start),
                                            $outputDir); }
  else              { return "Error" }
}



sub xx_bargraph {
	my (%h_param) = @_;

	my $res = "";
	foreach my $k (keys (%h_param)) {
		$res .= "\nvar chart_".$k." = new Highcharts.Chart({\n";
		$res .= "chart: {renderTo:\'$k\',type:\'bar\',spacing:[15,5,5,5]},\n";
		$res .= '
	    title: {text:null},
	    exporting: {enabled:false},
	    credits: \'enabled\',
	    xAxis: {visible:false},
	    yAxis: {title: {text:null},reversedStacks: false,endOnTick:false,max:';
		$res .= ($h_param{$k}{"len"}+1).",min:0,labels:{y:10}},";
		$res .= '
	    legend: {enabled:false},
	    tooltip: {padding: 5, positioner: function () {return { x:2, y:2 };},
								shared:true,
								shadow: false,
        				borderWidth: 0,
								formatter: function() {
									var s = [];
									var pos = 0;
									var nb  = 1;
									$.each(this.points, function(i, point) {
										if(point.series.color.stops[0][1] == "#ba2020") {
											s.push( "<span><b>BlockX_"+nb+":</b> "+
															pos+"-"+(pos+point.y-1)+
															" ("+((pos+point.y-1)-pos+1)+"bp)</span>");
											nb++;
										}
										pos += point.y;
									});
									if(s.length == 0) {s.push("No XXed block found")}
									return s.join(\'<br>\');
								},
			},
	    plotOptions: {
	      series: {
	      		stacking: \'normal\',
	          dataLabels: {
	              enabled: true,
	              y:+20,
	              color:\'#666666\'
	          }
	      }
	    },
	    series: [';
		my $c1 = "{radialGradient:{cx:0.5,cy:0.5,r:0.5},stops:[[0,'#159fc0'],[1,'#067495']]}";
		my $c2 = "{radialGradient:{cx:0.5,cy:0.5,r:0.5},stops:[[0,'#ba2020'],[1,'#940e0e']]}";
		# No XX found
		if($h_param{$k}{"XX"} !~ /^\[/) {
			$res .= "{color: $c1, data: [".$h_param{$k}{"len"}."]}";
	  }
		else {
			$h_param{$k}{"XX"} =~s/\[//g;
			my @a_coord = split(']',$h_param{$k}{"XX"});
			my $nb   = 0;
			my $prev = 0;
			foreach my $coord (@a_coord) {
				$res .= "," if($nb);
				my ($start, $stop) = split('-',$coord);
				if($start-$prev != 0) {
					$res .= "{color: $c1,	data: [".($start-$prev)."]},";
				}
				$res .= "{color: $c2, data: [".($stop-$start+1)."]}";
				$nb++;
				$prev = $stop+1;
			}
			if(($h_param{$k}{"len"}-$prev-1) != -1) {
				$res .= ",{color: $c1, data: [".($h_param{$k}{"len"}-$prev-1)."]}";
			}
		}
		$res .= ']
		});';
	}
	return $res;
}


MAIN:
{
	# Parameters
	my @a_inputSeq  = ();
	my $fosmidLen   = 40000;
	my $minCtgLen   = 1000;
	my $minCtgDepth = 8;
	my $vectorSeq   = undef;
	my $coverage    = 300;
	my $outputDir   = undef;
	my $outputZip   = "mintia_assembly.zip";
	my $outputHtml  = "mintia_assembly.html";
	my $outputLog   = "mintia.log";
	my $threads     = 8;
	my $version     = 0;
	my $help        = 0;

	# Data store in hash like:
	# %h_sample = (
	#    samplename1 => {
	#          R1 => 'PATH_TO_R1_INPUT_FQ_FILE',
	#          R2 => 'PATH_TO_R2_INPUT_FQ_FILE',      # if paired
	#          R1_F => 'PATH_TO_R1_FILTERED_FQ_FILE',
	#          R2_F => 'PATH_TO_R2_FILTERED_FQ_FILE', # if paired
	#          nbScaffold => int,
	#          scaffold => {
	#            ID1 => {
	#              len => int,
	#              cov => int,
	#              seq => str,
	#              newseq => str
	#            }
	#            [...]
	#          }
	#    }
	#    [...]
	# )
  	my %h_sample  = ();
	my %h_sampleR = ();  ## Same hash use for store removed scaffold

	GetOptions(
		'i|input=s{1,}'         => \@a_inputSeq,
		'l|length=i'            => \$fosmidLen,
		'minimalContigLength=i' => \$minCtgLen,
		'minimalContigDepth=i'  => \$minCtgDepth,
		'v|vectorSeq=s'         => \$vectorSeq,
		'c|maxDepth=i'          => \$coverage,
		'd|dirOutputs=s'        => \$outputDir,
		'Z|zipOutput=s'         => \$outputZip,
		'H|htmlOutput=s'        => \$outputHtml,
		'L|logOutput=s'         => \$outputLog,
		't|threads=i'           => \$threads,
		'version'               => \$version,
		'h|help'                => \$help
	);
	pod2usage(
		-verbose => 99,
		-sections => "SYNOPSIS|DESCRIPTION|OPTIONS|VERSION"
	) if($help);
	pod2usage(
		-verbose => 99,
		-sections => "VERSION"
	) if($version);
	$version = `$0 --version`;
	$version =~ s/^Ver.*:\s+(\S+)\s+/$1/;
	pod2usage("$0: '-i, --input' is required.")      if($#a_inputSeq == -1);
	foreach my $i (@a_inputSeq) {
		pod2usage("$0: $i doesn't exist.")   if !(-e $i);
		pod2usage("$0: $i is not readable.") if !(-r $i);
	}
	pod2usage("$0: '-v, --vectorSeq' is required.")  if !defined($vectorSeq);
	pod2usage("$0: $vectorSeq doesn't exist.")       if !(-e $vectorSeq);
	pod2usage("$0: $vectorSeq is not readable.")     if !(-r $vectorSeq);
	pod2usage("$0: '-d, --dirOutputs' is required.") if !defined($outputDir);

	if(! -e $outputDir) { mkdir $outputDir || die "Error: Enabled to create output dir $outputDir."; }
	open(LOG, ">$outputDir/$outputLog")    || die "Error: Enabled to create $outputDir/$outputLog";

	#Check input fastq file(s)
	print LOG "## Inputs...";
	foreach my $i (@a_inputSeq) {
    	#Find sample name
		my $sn = basename($i);
		if($sn =~ /^(.*)R[12].f[ast]*q[.gz]*$/) { $sn = $1; $sn =~s/[\.\_]$//; }
    	elsif($sn =~/^(.*).f[ast]*q[.gz]*$/)    { $sn = $1; }
		else {
    	  pod2usage("Error: enable to extract sample name from $i (must contain .f[ast]q[.gz]).")
	    }
		#Build $h_sample
		if(exists($h_sample{$sn}{"R1"})) {
			if(exists($h_sample{$sn}{"R2"})) {
			     pod2usage("Error: more than 2 fastq files related to $sn.");
			}
		else { $h_sample{$sn}{"R2"} = $i; }
		}
		else { $h_sample{$sn}{"R1"} = $i; }
	}
	print LOG (keys %h_sample) . " sample(s) found:\n";
	foreach my $k (sort keys(%h_sample)) {
    	print LOG " - $k: ";
    	if(exists($h_sample{$k}{"R2"})) { print LOG "paired-end\n" }
    	else                            { print LOG "single-read\n" }
	}

	#Filter using Max depth by sample
	print LOG "\n## Read filter using max depth ($coverage)\n";
	my $maxLen = $coverage*$fosmidLen;
	foreach my $k (sort keys(%h_sample)) {
		my $zip = "cat ";
		my $filetype = `file -bsiL $h_sample{$k}{"R1"}`;
	    if($filetype =~ /application\/x\-gzip/) { $zip = "gunzip -c "; }
    	open(READ1, "$zip $h_sample{$k}{'R1'} |") || die "Error: Enabled to open $h_sample{$k}{'R1'}.";

		if(exists($h_sample{$k}{"R2"})) {
		  $zip = "cat ";
		  $filetype = `file -bsiL $h_sample{$k}{'R2'}`;
		  if($filetype =~ /application\/x\-gzip/) { $zip = "gunzip -c "; }
		  open(READ2, "$zip $h_sample{$k}{'R2'} |") || die "Error: Enabled to open $h_sample{$k}{'R2'}.";
		}

		#outputDir
		open(READ1FILTER, ">$outputDir/$k"."_R1.fq") || die "Error: Enabled to create $outputDir/$k"."_R1.fq";
		$h_sample{$k}{"R1_F"} = "$outputDir/$k"."_R1.fq";
		if(exists($h_sample{$k}{"R2"})) {
		  open(READ2FILTER, ">$outputDir/$k"."_R2.fq") || die "Error: Enabled to create $outputDir/$k"."_R2.fq";
		  $h_sample{$k}{"R2_F"} = "$outputDir/$k"."_R2.fq";
		}

    	my $curLen          = 0;
		my $totalLen        = 0;
	    my $nbFragment      = 0;
		my $nbFragmentTotal = 0;
		my ($r1head, $r1seq, $r1sep, $r1qual);
		my ($r2head, $r2seq, $r2sep, $r2qual);
		while(! eof READ1) {
			$r1head = <READ1>; $r1seq = <READ1>; $r1sep = <READ1>; $r1qual = <READ1>;
			#Validate Fastq R1
			if($r1head!~/^@/ || $r1seq!~/^[ATGCNatgcn]*$/ || $r1sep!~/^\+/) {
				print "$r1head\n$r1seq\n$r1sep\n";
				pod2usage("Error: $h_sample{$k}{'R1'} is not a valid fastq file.")
			}
			if(exists($h_sample{$k}{'R2'})) {
				$r2head = <READ2>; $r2seq = <READ2>; $r2sep = <READ2>; $r2qual = <READ2>;
				#Validate Fastq R2
				if($r2head!~/^@/ || $r2seq!~/^[ATGCNatgcn]+$/ || $r2sep!~/^\+/) {
					pod2usage("Error: $h_sample{$k}{'R2'} is not a valid fastq file.")
				}
			}
			if($curLen<$maxLen) {
				$nbFragment++;
				$curLen+=length($r1seq)-1;
	    		print READ1FILTER "$r1head$r1seq$r1sep$r1qual";
				if(exists($h_sample{$k}{'R2'})) {
					$curLen+=length($r2seq)-1;
					print READ2FILTER "$r2head$r2seq$r2sep$r2qual";
				}
			}
			$nbFragmentTotal++;
			$totalLen+=length($r1seq)-1;
			$totalLen+=length($r2seq)-1 if(exists($h_sample{$k}{'R2'}));
		}
    	close READ1;
    	close READ1FILTER ;
    	if(exists($h_sample{$k}{'R2'})) {
      		close READ2;
     		close READ2FILTER;
    	}
    	print LOG " - $k: $nbFragment/$nbFragmentTotal reads take into account\n";
		$h_sample{$k}{'nbFragment'} = $nbFragmentTotal;
		$h_sample{$k}{'totalLen'}   = $totalLen;
		if($curLen != $totalLen) {
			$h_sample{$k}{'nbFragment_F'} = $nbFragment;
			$h_sample{$k}{'totalLen_F'}   = $curLen;
		}
  	}

	#SPADES
	print LOG "\n## Run SPADES\n";
	foreach my $k (sort keys(%h_sample)) {
		print LOG " - $k...";
		if(exists($h_sample{$k}{'R2_F'})) {
		  `$SPADES -1 $h_sample{$k}{'R1_F'} -2 $h_sample{$k}{'R2_F'} -t $threads --careful -o $outputDir/$k`;
		}
		else {
		  `$SPADES -s $h_sample{$k}{'R1_F'} -t $threads --careful -o $outputDir/$k`;
		}
		print LOG "done\n";
	}

  #CROSS MATCH (on the scaffolds.fasta generated by SPADES)
  print LOG "\n## Run Crossmatch\n";
  foreach my $k (sort keys(%h_sample)) {
    print LOG " - $k...";
    `($CROSS_MATCH $outputDir/$k/scaffolds.fasta $vectorSeq -minmatch 9 -minscore 30 -screen > $outputDir/$k/crossmatch-scaffolds.stdout) >& $outputDir/$k/crossmatch-scaffolds.stderr`;
    print LOG "done\n";
  }

  #Find longuest scaffold and
  #Filter scaffold using --minimalContigLength and Contig's minimum length [1000]
  print LOG "\n## Filter scaffolds\n";
  foreach my $k (sort keys(%h_sample)) {
    print LOG " - $k...";
    open(SCAF, "$outputDir/$k/scaffolds.fasta.screen") || die "Error: Enabled to open $outputDir/$k/scaffolds.fasta.screen.";
    my ($id, $len, $cov, $seq) = (0,0,0,"");
    my $nbscaffold = 0;
    while(my $line = <SCAF>) {
      next if($line=~/^\s*$/);
      if($line=~/^>/ && $len>$minCtgLen && $cov>$minCtgDepth) {
        $h_sample{$k}{"scaffolds"}{$id}{"len"} = $len;
        $h_sample{$k}{"scaffolds"}{$id}{"cov"} = sprintf("%.2f", $cov);
        $h_sample{$k}{"scaffolds"}{$id}{"seq"} = $seq;
      }
			elsif($line=~/^>/ && $nbscaffold!=0) {
				$h_sampleR{$k}{"scaffolds"}{$id}{"len"} = $len;
      	$h_sampleR{$k}{"scaffolds"}{$id}{"cov"} = sprintf("%.2f", $cov);
			}
      if($line=~/^>NODE\_(\d+)\_length\_(\d+)\_cov\_([\d\.]+)\s+$/) {
        $id  = $1;
        $len = $2;
        $cov = $3;
				$seq = "";
        $nbscaffold++;
      }
      elsif($line=~/^[ATGCNXatgcnx]+$/) {
        chomp $line;
        $seq.=$line;
      }
      else {
        print "Error: $outputDir/$k/scaffolds.fasta.screen unexpected line\n$line";
        exit();
      }
    }
		# Last scaffold
    if($len>$minCtgLen && $cov>$minCtgDepth) {
      $h_sample{$k}{"scaffolds"}{$id}{"len"} = $len;
      $h_sample{$k}{"scaffolds"}{$id}{"cov"} = sprintf("%.2f", $cov);
      $h_sample{$k}{"scaffolds"}{$id}{"seq"} = $seq;
    }
		else {
			$h_sampleR{$k}{"scaffolds"}{$id}{"len"} = $len;
      $h_sampleR{$k}{"scaffolds"}{$id}{"cov"} = sprintf("%.2f", $cov);
		}
    close SCAF;
    print LOG (keys %{$h_sample{$k}{"scaffolds"}})."/$nbscaffold selected\n";
    $h_sample{$k}{"nbScaffold"} = $nbscaffold;
  }

  #Find vector (XX regions) from crossmatch output (scaffolds.fasta.screen)
  print LOG "\n## Find and remove XXed block(s):\n";
  foreach my $k1 (sort keys(%h_sample)) {
    print LOG " - $k1:\n";
    foreach my $k2 (sort {$a <=> $b} keys(%{$h_sample{$k1}{"scaffolds"}})) {
      print LOG "   - Scaffold $k2:\n";
      print LOG "     - Length: ".$h_sample{$k1}{"scaffolds"}{$k2}{"len"}." bp\n";
      print LOG "     - XXed block coordinate: ";
      my @a_xxpos = locate_regex("X+", $h_sample{$k1}{"scaffolds"}{$k2}{"seq"});
      if(scalar @a_xxpos) {
        my $xxlen = 0;
        for(my $i=0; $i<=$#a_xxpos; $i++) {
          print LOG "[$a_xxpos[$i][0]-$a_xxpos[$i][1]]\t";
          $h_sample{$k1}{"scaffolds"}{$k2}{"XX"} .= "[$a_xxpos[$i][0]-$a_xxpos[$i][1]]";
          $xxlen += $a_xxpos[$i][1]-$a_xxpos[$i][0]+1;
        }
        print LOG "\n     - XXed block length: $xxlen bp\n";
        my $newseq = remove_vector(\@a_xxpos, $h_sample{$k1}{"scaffolds"}{$k2}{"seq"}, $outputDir);
        $h_sample{$k1}{"scaffolds"}{$k2}{"newseq"} = $newseq;
				print LOG "     - New scaffold length: ";
				if($newseq =~ /Error/) { print LOG "Error - Unexpected number of XXed blocks\n"; }
				else                   { print LOG length($newseq)." bp\n";	}
      }
      else {
        print LOG "No vector found!\n";
        $h_sample{$k1}{"scaffolds"}{$k2}{"XX"} = "-";
				$h_sample{$k1}{"scaffolds"}{$k2}{"newseq"} = $h_sample{$k1}{"scaffolds"}{$k2}{"seq"};
      }
    }
  }

  # Create output files:
  #  - scaffold filtered vector cleaned fasta file by sample
  #  - info table
  #  - report html:
  #    - summary info/stat
  #    - histogramme of length (scaffold with/without vector, vector)
  #    - plot length vs coverage

  # Table
  print "#Sample\tPaired\tNb_scaffold\tNb_Filtered\t(ID\tLength\tCoverage\tXXed_Block\tClean_length)xN\n";
  foreach my $k1 (sort keys(%h_sample)) {
    print "$k1\t";                                       #Sample
    if(exists($h_sample{$k1}{"R2"})) { print "Yes\t"; }  #Paired
    else                             { print "No\t";  }
    print $h_sample{$k1}{"nbScaffold"}."\t";             #Nb_scaffold
    print "".(keys %{$h_sample{$k1}{"scaffolds"}})."\t"; #Nb_Filtered

    #By filtered scaffold
    foreach my $k2 (sort {$a <=> $b} keys(%{$h_sample{$k1}{"scaffolds"}})) {
      print "$k2\t";                                              #ID
      print $h_sample{$k1}{"scaffolds"}{$k2}{"len"}."\t";         #Length
      print $h_sample{$k1}{"scaffolds"}{$k2}{"cov"}."\t";         #Coverage
      print $h_sample{$k1}{"scaffolds"}{$k2}{"XX"}."\t";          #XXed_Block
			if($h_sample{$k1}{"scaffolds"}{$k2}{"newseq"} =~ /Error/) { #Clean_length
				print "Error - Unexpected number of XXed blocks\t";
			}
			else {
				print length($h_sample{$k1}{"scaffolds"}{$k2}{"newseq"})."\t";
			}
    }
    print "\n";
  }



	##############################
	## HTML Report
	##############################
	open (HTML, ">$outputDir/mintia.html") || die "Error: Enabled to create $outputDir/mintia.html\n";
  print HTML '<!doctype html>
<html lang="en">
  <head>
    <!-- Required meta tags -->
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
    <meta name="description" content="Mintia assembly module">
    <meta name="author" content="Philippe Bardou">

    <!-- Bootstrap CSS -->
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css" integrity="sha384-Gn5384xqQ1aoWXA+058RXPxPg6fy4IWvTNh0E263XmFcJlSAwiGgFAW/dAiS6JXm" crossorigin="anonymous">
    <!-- Icons -->
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/open-iconic/1.1.1/font/css/open-iconic-bootstrap.min.css">
    <!-- Datatable -->
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.16/css/dataTables.bootstrap4.min.css">

    <!-- Personalised css -->
    <style type="text/css">
      body { font-size: .875rem; }
      .sidebar {
        position: fixed;
        top: 0;
        bottom: 0;
        left: 0;
        z-index: 100; /* Behind the navbar */
        padding: 0;
        box-shadow: inset -1px 0 0 rgba(0, 0, 0, .1);
      }
      .sidebar-sticky {
        position: -webkit-sticky;
        position: sticky;
        top: 48px; /* Height of navbar */
        height: calc(100vh - 48px);
        padding-top: 4rem;
        overflow-x: hidden;
        overflow-y: auto; /* Scrollable contents if viewport is shorter than content. */
      }
      .sidebar .nav-link {
        margin-right: 4px;
        color: #212529;
      }
      .sidebar .nav-link:hover,
      .sidebar .nav-link.active { color: #999; }
      .sidebar-heading {
        font-size: .75rem;
        text-transform: uppercase;
      }
      .navbar-brand {
        padding-top: .75rem;
        padding-bottom: .75rem;
        font-size: 1rem;
        background-color: rgba(0, 0, 0, .7);
        box-shadow: inset -1px 0 0 rgba(0, 0, 0, .7);
				z-index: 200;
      }
      .navbar .form-control {
        padding: .75rem 1rem;
        border-width: 0;
        border-radius: 0;
				z-index: 200;
      }
      .form-control-dark {
        color: #fff;
        background-color: rgba(255, 255, 255, .1);
        border-color: rgba(255, 255, 255, .1);
      }
      .form-control-dark:focus {
        border-color: transparent;
        box-shadow: 0 0 0 3px rgba(255, 255, 255, .25);
      }
      .border-top { border-top: 1px solid #e5e5e5; }
      .border-bottom { border-bottom: 1px solid #e5e5e5; }
			.valn { vertical-align: middle !important; }
			.anchor{
  			display: block;
  			height: 83px; /*same height as header*/
  			margin-top: -83px; /*same height as header*/
  			visibility: hidden;
			}
		}
    </style>

    <title>Mintia assembly report</title>
  </head>
  <body>
    <nav class="navbar navbar-dark sticky-top bg-dark flex-md-nowrap p-0">
      <a class="navbar-brand col-sm-3 col-md-2 mr-0" href="#">Mintia assembly report</a>
    </nav>

    <div class="container-fluid">
      <div class="row">
        <nav class="col-md-2 d-none d-md-block bg-light sidebar">
          <div class="sidebar-sticky">
            <ul class="nav flex-column">
              <li class="nav-item">
                <a class="nav-link" href="#inputs-parameters">
                  <span class="oi oi-list-rich" aria-hidden="true"></span>
                  Inputs and parameters
                </a>
              </li>
							<li class="nav-item" style="padding-left:10px">
								<a class="nav-link" href="#parameters">
									<span class="oi oi-list" aria-hidden="true"></span>
									Parameters
								</a>
							</li>
							<li class="nav-item" style="padding-left:10px">
								<a class="nav-link" href="#cloning-vector">
									<span class="oi oi-pie-chart" aria-hidden="true"></span>
									Cloning vector
								</a>
							</li>
							<li class="nav-item" style="padding-left:10px">
              	<a class="nav-link" href="#fastq-files">
              		<span class="oi oi-file" aria-hidden="true"></span>
              		Fastq files
              	</a>
            	</li>
							<li class="nav-item">
								<a class="nav-link" href="#assembly-detection">
									<span class="oi oi-eye aria-hidden="true"></span>
									Assembly and vector detection
								</a>
							</li>
              <li class="nav-item">
                <a class="nav-link" href="#cov-len">
                  <span class="oi oi-graph aria-hidden="true"></span>
                  Scaffold coverage VS length
                </a>
              </li>
							<li class="nav-item">
                <a class="nav-link" href="#download">
									<span class="oi oi-data-transfer-download"></span>
									Download file
								</a>
							</li>
            </ul>
          </div>
';
	print HTML "<div style=\"text-align:center;font-size:smaller;color:darkgrey;margin-top:-25px\">
		Produced by $version<br>
		Copyright © 2018, <img src=\"http://www.inra.fr/extension/itkinra/design/inra/images/favicon.ico\">
		<a style=\"color:#212529;\" href=\"http://inra.fr\" target=\"_blank\">INRA</a><br>
		Designed by the <a style=\"color:#212529;\" href=\"http://sigenae.org\" target=\"_blank\">Sigenae</a> team.</div>";
	print HTML '
        </nav>

				<main role="main" class="col-md-9 ml-sm-auto col-lg-10 pt-3 px-4">
          <div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center pb-2 border-bottom">
						<h1 class="h4">Inputs and parameters</h1>
						<span class="anchor" id="inputs-parameters"></span>
					</div>
';
	print HTML '
					<div class="d-flex">
					  <div class="mt-4 mr-4 pl-0 col-md-4">
							<h5>Parameters</h5>
							<span class="anchor" id="parameters"></span>
				    	<ul class="list-group">
';
  print HTML "<li class=\"list-group-item d-flex justify-content-between align-items-center\">Fosmid\'s expected length <span class=\"badge badge-danger badge-pill ml-4\">$fosmidLen"."bp</span></li>";
	print HTML "<li class=\"list-group-item d-flex justify-content-between align-items-center\">Maximum coverage <span class=\"badge badge-danger badge-pill ml-4\">$coverage</span></li>";
	print HTML "<li class=\"list-group-item d-flex justify-content-between align-items-center\">Contig\'s minimum length <span class=\"badge badge-danger badge-pill ml-4\">$minCtgLen"."bp</span></li>";
	print HTML "<li class=\"list-group-item d-flex justify-content-between align-items-center\">Contig\'s minimum depth <span class=\"badge badge-danger badge-pill ml-4\">$minCtgDepth</span></li>";
	print HTML	'
					    </ul>
						</div>

						<div class="mt-4">
							<h5>Cloning vector</h5>
							<span class="anchor" id="cloning-vector"></span>
';
	open(VECTOR, "$vectorSeq") || die "Error: Enabled to open $vectorSeq\n";
	my ($id,$len,$gc,$A,$T,$C,$G, $o) = ("", 0, 0, 0, 0, 0, 0, 0);
	while(my $line=<VECTOR>) {
		chomp $line;
		if($line =~/^>(\S+)/) { $id = $1; }
		else {
			$len += length($line);
			$A   += ($line =~ tr/aA/aA/);
			$T   += ($line =~ tr/tT/tT/);
			$C   += ($line =~ tr/cC/cC/);
			$G   += ($line =~ tr/gG/gG/);
		}
	}
	$gc = sprintf("%.2f", ($G+$C)/$len*100);
	$o  = $len-($A+$T+$C+$G);
	close(VECTOR);

	print HTML '
							<div class="d-flex">
							  <ul class="list-group d-flex flex-row flex-wrap">
';
	print HTML "<li class=\"list-group-item w-100\"><b> Name:</b> $id</li>\n";
	print HTML "<li class=\"list-group-item w-50\"><b> Length:</b> $len</li>\n";
	print HTML "<li class=\"list-group-item w-50\"><b> GC%:</b> $gc</li>\n";
	print HTML "<li class=\"list-group-item w-50\"><b> A:</b> $A</li>\n";
	print HTML "<li class=\"list-group-item w-50\"><b> T:</b> $T</li>\n";
	print HTML "<li class=\"list-group-item w-50 mb-0\" style=\"border-bottom-left-radius:.25rem;\"><b> G:</b> $C</li>\n";
	print HTML "<li class=\"list-group-item w-50 mb-0\" style=\"border-bottom-left-radius:0;\"><b> C:</b> $G</li>\n";

	print HTML '
								</ul>
								<div id="compseq-graph"></div>
							</div>
						</div>
					</div>

					<h5 class="mt-4">Fastq files</h5>
					<span class="anchor" id="fastq-files"></span>
					<table class="table table-striped table-bordered mb-0" style="width:100%;">
        		<thead>
            	<tr>
                <th>Sample</th>
								<th style="text-align:center">Paired</th>
                <th nowrap style="text-align:center">Nb. of seq.</th>
								<th nowrap style="text-align:center">Full seq. length</th>
								<th style="text-align:center">Coverage</th>
								<th nowrap style="text-align:center">Nb. of filtered seq.</th>
								<th nowrap style="text-align:center">Filtered seq. length</th>
								<th nowrap style="text-align:center">Coverage (filtered)</th>
            	</tr>
        		</thead>
        		<tbody>
';
	foreach my $k1 (sort keys(%h_sample)) {
		my $class = "class=\"valn text-right\"";
		print HTML "<tr><td class=\'valn\'>$k1</td>";
		if(exists($h_sample{$k1}{"R2"})) { print HTML "<td class=\'valn\'>Yes</td>"; }
		else                             { print HTML "<td class=\'valn\'>No</td>";  }
		print HTML "<td $class>".$h_sample{$k1}{"nbFragment"}."</td>";
		print HTML "<td $class>".$h_sample{$k1}{"totalLen"}."</td>";
		printf HTML "<td $class>%d</td>", ($h_sample{$k1}{"totalLen"}/$fosmidLen);
		if(exists($h_sample{$k1}{"nbFragment_F"})) {
			print HTML "<td $class>".$h_sample{$k1}{"nbFragment_F"}."</td>";
			print HTML "<td $class>".$h_sample{$k1}{"totalLen_F"}."</td>";
			printf HTML "<td $class>%d</td>", ($h_sample{$k1}{"totalLen_F"}/$fosmidLen);
		}
		else {
			print HTML "<td $class>-</td><td $class>-</td><td $class>-</td>";
		}
		print HTML "</tr>";
	}
	print HTML '
						</tbody>
					</table>

					<div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center mt-5 pb-2 border-bottom">
						<h1 class="h4">Assembly and vector detection</h1>
						<span class="anchor" id="assembly-detection"></span>
          </div>
					<table class="table table-striped table-bordered mt-4 mb-0" style="width:100%">
        		<thead>
            	<tr>
                <th>Sample</th>
                <th nowrap style="text-align:center">Scaffold</th>
								<th nowrap style="text-align:center">Filtered</th>
                <th nowrap style="text-align:center">ID</th>
';
	print HTML "<th nowrap style='text-align:center'>Length (>$minCtgLen)</th>\n";
	print HTML "<th nowrap style='text-align:center'>Coverage (>$minCtgDepth)</th>\n";
	print HTML '
								<th nowrap style="text-align:center;width:275px;">Vector position</th>
								<th nowrap style="text-align:center">Clean length</th>
            	</tr>
        		</thead>
        		<tbody>
';

	my %h_xxGraph;
	my $sm = 0;
	foreach my $k1 (sort keys(%h_sample)) {
		my $class = "class=\"valn text-right\"";
		my $rowspan = "";
		if((keys(%{$h_sample{$k1}{"scaffolds"}})) > 1) {
			$rowspan = " rowspan=\"".(keys(%{$h_sample{$k1}{"scaffolds"}}))."\"";
		}
		print HTML "<tr><td class='valn' $rowspan>$k1</td>\n";
		print HTML "<td $class $rowspan>".$h_sample{$k1}{"nbScaffold"}."</td>\n";
		print HTML "<td $class $rowspan>".(keys %{$h_sample{$k1}{"scaffolds"}})."</td>\n";
		#By filtered scaffold
		my $scaff = 0;
		foreach my $k2 (sort {$a <=> $b} keys(%{$h_sample{$k1}{"scaffolds"}})) {
			print HTML "<tr>\n" if($scaff != 0);
			print HTML "<td $class>$k2</td>\n";
			print HTML "<td $class>".$h_sample{$k1}{"scaffolds"}{$k2}{"len"}."</td>\n";
			print HTML "<td $class>".$h_sample{$k1}{"scaffolds"}{$k2}{"cov"}."</td>\n";
			$h_xxGraph{"xxGraph_$sm$scaff"}{"len"} = $h_sample{$k1}{"scaffolds"}{$k2}{"len"};
			$h_xxGraph{"xxGraph_$sm$scaff"}{"XX"}  = $h_sample{$k1}{"scaffolds"}{$k2}{"XX"};
			print HTML "<td $class><div id='xxGraph_$sm$scaff' style='min-width: 250px; max-width: 250px; height:55px; margin: 0 auto'></div></td>\n";
			if($h_sample{$k1}{"scaffolds"}{$k2}{"newseq"} =~ /Error/) {
				print HTML "<td $class><b>Error:</b> Unexpected number of XXed blocks</td></tr>\n";
			}
			else {
				print HTML "<td $class>".length($h_sample{$k1}{"scaffolds"}{$k2}{"newseq"})."</td></tr>\n";
			}
			$scaff++;
		}
		$sm++;
	}
	print HTML '
						</tbody>
					</table>

	        <div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center mt-5 pb-2 border-bottom">
            <h1 class="h4">Scaffold coverage VS length</h1>
						<span class="anchor" id="cov-len"></span>
          </div>
					<div class="mt-4" id="covVSlen-graph"></div>

					<span class="anchor" id="download"></span>
	        <div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center mt-5 pb-2 border-bottom">
            <h1 class="h4">Download file</h1>
          </div>
					<table class="table table-striped table-bordered mt-4 mb-5" style="width:100%">
        		<thead>
            	<tr>
                <th>Sample</th>
                <th nowrap style="text-align:center">Scaffold sequence</th>
								<th nowrap style="text-align:center">Vector position</th>
								<th nowrap style="text-align:center">Clean scaffold sequence</th>
        			</tr>
        		</thead>
        		<tbody>
';
	my $download      = "var dwnl = [";
	my $allfasta      = "var allfasta = \"";
	my $allgff        = "var allgff = \"";
	my $allcleanfasta = "var allclean = \"";
	my $dwnlid=0;
	foreach my $k1 (sort keys(%h_sample)) {
		print HTML "<tr><td class='valn' >$k1</td>\n";

		my $seq = "";
		foreach my $k2 (sort {$a <=> $b} keys(%{$h_sample{$k1}{"scaffolds"}})) {
			$seq .= ">$k1#$k2";
			$seq .= "#".$h_sample{$k1}{"scaffolds"}{$k2}{"len"};
			$seq .= "#".$h_sample{$k1}{"scaffolds"}{$k2}{"cov"}."\\n";
			$seq .= join("\\n", $h_sample{$k1}{"scaffolds"}{$k2}{"seq"}=~/(.{1,60})/g)."\\n";
		}
		if($dwnlid != 0) { $download .= ",";}
		$download .= "\"$seq\"";
		$allfasta .= "$seq";
		print HTML "<td style=\"text-align:right\"><button id=\"dwnl-$dwnlid\" file=\"$k1.fasta\" class=\"btn btn-sm btn-outline-secondary\">
								Fasta (".(keys %{$h_sample{$k1}{"scaffolds"}}).")</button></td>\n";

		$dwnlid++;
		$seq = "";
		foreach my $k2 (sort {$a <=> $b} keys(%{$h_sample{$k1}{"scaffolds"}})) {
			my $xx = $h_sample{$k1}{"scaffolds"}{$k2}{"XX"};
			$xx =~s/\[//g;
			$xx =~s/^-$//g;
			my @a_coord = split(']', $xx);
			my $cmp = 1;
			foreach my $coord (@a_coord) {
				$seq .= "$k1#$k2\\t"; #seqname
				$seq .= "$version\\tBlockXX_$cmp\\t"; #source feature
				my ($start, $stop) = split('-',$coord);
				$seq .= ($start+1)."\\t".($stop+1)."\\t"; # start end
				$seq .= ".\\t.\\t.\\n"; #score strand frame
				$cmp++;
			}
		}
		if($dwnlid != 0) { $download .= ",";}
		$download .= "\"$seq\"";
		$allgff .= "$seq";
		print HTML "<td style=\"text-align:right\"><button id=\"dwnl-$dwnlid\" file=\"$k1.gff\" class=\"btn btn-sm btn-outline-secondary\">
								GFF</button></td>\n";

		$dwnlid++;
		$seq = "";
		foreach my $k2 (sort {$a <=> $b} keys(%{$h_sample{$k1}{"scaffolds"}})) {
			$seq .= ">$k1#k2";
			$seq .= "#".length($h_sample{$k1}{"scaffolds"}{$k2}{"newseq"})."\\n";
			$seq .= join("\\n", $h_sample{$k1}{"scaffolds"}{$k2}{"newseq"}=~/(.{1,60})/g)."\\n";
		}
		if($dwnlid != 0) { $download .= ",";}
		$download .= "\"$seq\"";
		$allcleanfasta .= "$seq";
		print HTML "<td style=\"text-align:right\"><button id=\"dwnl-$dwnlid\" file=\"$k1.clean.fasta\" class=\"btn btn-sm btn-outline-secondary\">
								Fasta (".(keys %{$h_sample{$k1}{"scaffolds"}}).")</button></td>\n";
		print HTML "</tr>\n";
		$dwnlid++;
	}
	$download .= "];";
	print HTML "\n<script>$download</script>\n";
	print HTML "\n<script>$allfasta\"</script>\n";
	print HTML "\n<script>$allgff\"</script>\n";
	print HTML "\n<script>$allcleanfasta\"</script>\n";
	print HTML '		<tr>
							<td></td>';
	print HTML "\n<td style=\"text-align:center\"><button id=\"allfasta\" file=\"all.fasta\" class=\"btn btn-sm btn-outline-secondary\">
								All Fasta</button></td>\n";
	print HTML "<td style=\"text-align:center\"><button id=\"allgff\" file=\"all.gff\" class=\"btn btn-sm btn-outline-secondary\">
								All GFF</button></td>\n";
	print HTML "<td style=\"text-align:center\"><button id=\"allcleanfasta\" file=\"allcleaned.fasta\" class=\"btn btn-sm btn-outline-secondary\">
								All cleaned Fasta</button></td>\n";
	print HTML '								
						</tr>
						</tbody>
					</table>
        </main>
      </div>
    </div>

    <!-- jQuery first, then Popper.js, then Bootstrap JS -->
    <script type="text/javascript" language="javascript" src="https://code.jquery.com/jquery-3.2.1.slim.min.js" integrity="sha384-KJ3o2DKtIkvYIK3UENzmM7KCkRr/rE9/Qpg6aAZGJwFDMVNA/GpGFF93hXpG5KkN" crossorigin="anonymous"></script>
    <script type="text/javascript" language="javascript" src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.12.9/umd/popper.min.js" integrity="sha384-ApNbgh9B+Y1QKtv3Rn7W3mgPxhU9K/ScQsAP7hUibX39j7fakFPskvXusvfa0b4Q" crossorigin="anonymous"></script>
    <script type="text/javascript" language="javascript" src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/js/bootstrap.min.js" integrity="sha384-JZR6Spejh4U02d8jOt6vLEHfe/JQGiRRSQQxSfFWpi1MquVdAyjUar5+76PVCmYl" crossorigin="anonymous"></script>

    <!-- Datatable -->
    <script type="text/javascript" language="javascript" src="https://cdn.datatables.net/1.10.16/js/jquery.dataTables.min.js"></script>
    <script type="text/javascript" language="javascript" src="https://cdn.datatables.net/1.10.16/js/dataTables.bootstrap4.min.js"></script>

		<! Highcharts -->
		<script src="https://code.highcharts.com/highcharts.js"></script>
		<script src="https://code.highcharts.com/modules/exporting.js"></script>
		<script src="https://code.highcharts.com/modules/export-data.js"></script>

    <!-- compseq-graph -->
		<script type="text/javascript" language="javascript">
		// Radialize the colors
		Highcharts.setOptions({
		    colors: Highcharts.map([\'#d72b2b\', \'#FF9655\', \'#64E572\', \'#159fc0\', \'#24CBE5\', \'#6AF9C4\', \'#FFF263\', \'#50B432\', \'#DDDF00\'], function (color) {
		        return {
		            radialGradient: {
		                cx: 0.5,
		                cy: 0.3,
		                r: 0.7
		            },
		            stops: [
		                [0, color],
		                [1, Highcharts.Color(color).brighten(-0.4).get(\'rgb\')] // darken
		            ]
		        };
		    })
		});

		var compseqGraph = new Highcharts.Chart({
		chart: {
			renderTo: "compseq-graph",
			plotBackgroundColor: null,
			plotBorderWidth: 0,
			width: 360,
			height: 180,
			margin: [-60, -25, -70, -25]
		},
		title: {
';
	print HTML "text: \"Seq. composition<br/>(GC%: <b>$gc</b>)\",";
	print HTML '
					align: "center",
					verticalAlign: "middle",
					y: 50,
					style: {"font-size": "13px"}
			},
			credits: false,
			tooltip: {
				formatter: function () {
					return \'<b>\' + this.key + \'</b>: \' + this.y;
				}
			},
			plotOptions: {
					pie: {
							dataLabels: {
									enabled: true,
									distance: -40,
									style: {
											fontWeight: "bold",
											color: "white",
											textShadow: "0px 1px 2px black"
									},
						format: "{point.name}: <b>{point.percentage:.1f}%</b>"
							},
							startAngle: -90,
							endAngle: 90,
							center: ["50%", "75%"]
					}
			},
			series: [{
					type: "pie",
					name: "Sequence composition",
					innerSize: "50%",
					data: [
';
	print HTML "['A', $A],['T', $T], ['C', $C], ['G', $G], ['Other', $o]";
	print HTML '
					]
			}]
		});
';
	print HTML xx_bargraph(%h_xxGraph);
	print HTML '
		var covVSlenGraph = new Highcharts.Chart({
			chart: {
				renderTo: "covVSlen-graph",
		    type: \'scatter\',
		    zoomType: \'xy\'
	    },
	    title: {
	        text: \'Coverage versus length\'
	    },
	    subtitle: {
';
	print HTML "text: \'Produced by $version\'";
	print HTML "
	    },
			credits: false,
     	colors:['rgba(186, 31, 31,.7)','rgba(186, 31, 31,.3)',
							'rgba( 15,129,162,.7)','rgba( 15,129,162,.3)',
							'rgba( 89,218,103,.7)','rgba( 89,218,103,.3)',
							'rgba(245,140, 77,.7)','rgba(245,140, 77,.3)',
							'rgba(148, 45,196,.7)','rgba(148, 45,196,.3)',
							'rgba(162, 91,109,.7)','rgba(162, 91,109,.3)',
							'rgba(215,222, 28,.7)','rgba(215,222, 28,.3)',
							'rgba(  0,  0,  0,.7)','rgba(  0,  0,  0,.3)',
							'rgba(227,140,181,.7)','rgba(227,140,181,.3)'],
	    xAxis: {
	        title: {
	            enabled: true,
	            text: \'Coverage\'
	        },
					min: 0,
	        startOnTick: true,
	        endOnTick: true,
	        showLastLabel: true,

					plotBands: [{
							color: \'rgba(215, 43, 43, .1)\',
							from: 0,
							to: $minCtgDepth,
							zIndex: 3
					}],
					plotLines: [{
											color: \'darkred\',
											dashStyle: \'dot\',
											width: 2,
											value: $minCtgDepth,
											label: {
																	rotation:-90,
																	x:12,
																	verticalAlign: \'top\',
                									textAlign: \'right\',
																	style: {fontStyle: \'italic\',color: \'darkred\'},
																	text: \'Scaffold minimum coverage\'
											},
											zIndex: 3
									}]
	    },
	    yAxis: {
	        title: {
	            text: \'Length\'
	        },
					min:0,
					plotBands: [{
							color: \'rgba(215, 43, 43, .1)\',
							from: 0,
							to: $minCtgLen,
							zIndex: 3
					}],
					plotLines: [{
											color: \'darkred\',
											dashStyle: \'dot\',
											width: 2,
											value: $minCtgLen,
											label: {
												align: \'right\',
												x:-12,
																	style: {fontStyle: \'italic\',color: \'darkred\'},
																	text: \'Scaffold minimum length\'
											},
											zIndex: 3
									}]
	    },
	    plotOptions: {
	        scatter: {
	            marker: {
	                radius: 5,
	                states: {
	                    hover: {
	                        enabled: true,
	                        lineColor: \'rgb(100,100,100)\'
	                    }
	                }
	            },
	            states: {
	                hover: {
	                    marker: {
	                        enabled: false
	                    }
	                }
	            },
	            tooltip: {
	                headerFormat: \'<b>{series.name}</b><br>\',
	                pointFormat: \'<b>Coverage:</b> {point.x}<br><b>Length:</b> {point.y}\'
	            }
	        }
	    },
	    series: [\n";
	my $f1 = 0;
	foreach my $k1 (sort keys(%h_sample)) {
		print HTML "," if($f1);
		print HTML "\n{name:\'$k1\', data:[";
		my $f2 = 0;
		foreach my $k2 (sort {$a <=> $b} keys(%{$h_sample{$k1}{"scaffolds"}})) {
			print HTML "," if($f2);
			print HTML "[".$h_sample{$k1}{"scaffolds"}{$k2}{"cov"};
			print HTML ",".$h_sample{$k1}{"scaffolds"}{$k2}{"len"}."]";
			$f2++;
		}
		print HTML "]}";
		$f1++;
		if(exists($h_sampleR{$k1}{"scaffolds"})) {
			print HTML ",\n{name:\'$k1"."_Removed\', data:[";
			$f2 = 0;
			foreach my $k2 (sort {$a <=> $b} keys(%{$h_sampleR{$k1}{"scaffolds"}})) {
				print HTML "," if($f2);
				print HTML "[".$h_sampleR{$k1}{"scaffolds"}{$k2}{"cov"};
				print HTML ",".$h_sampleR{$k1}{"scaffolds"}{$k2}{"len"}."]";
				$f2++;
			}
			print HTML "]}";
		}
	}
print HTML ']
		});

		<!-- Download -->
		$(\'[id^=dwnl]\').click(function(){
			download($(this).attr(\'file\'), dwnl[$(this).attr(\'id\').split(\'-\')[1]]);
		});
		<!-- Download ALL fasta -->
		$(\'[id=allfasta]\').click(function(){
			download($(this).attr(\'file\'), allfasta);
		});
		<!-- Download ALL gff -->
		$(\'[id=allgff]\').click(function(){
			download($(this).attr(\'file\'), allgff);
		});
		<!-- Download ALL cleaned fasta -->
		$(\'[id=allcleanfasta]\').click(function(){
			download($(this).attr(\'file\'), allcleanfasta);
		});

		function download(filename, text) {
    	var element = document.createElement(\'a\');
    	element.setAttribute(\'href\', \'data:text/plain;charset=utf-8,\' + encodeURIComponent(text));
    	element.setAttribute(\'download\', filename);
    	element.style.display = \'none\';
    	document.body.appendChild(element);
			element.click();
			document.body.removeChild(element);
		}
		</script>
  </body>
</html>
';
  close(HTML);
}

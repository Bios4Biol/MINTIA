#!/usr/bin/perl

=pod

=head1 NAME

 mintia.pl - Fosmid assembly and annotation pipeline.

=head1 SYNOPSIS

 mintia.pl check
 mintia.pl assemble -i FASTQ_FILE[S] -v FASTA_FILE -d STR
 mintia.pl annotate

=head1 CHECK SYNOPSIS

 mintia.pl check

=head1 ASSEMBLE SYNOPSIS

 mintia.pl assemble -i FASTQ_FILE[S] -v FASTA_FILE -d STR

=head1 ANNOTATE SYNOPSIS

 mintia.pl annotate -i FASTA_FILE -n NR_DMND_FILE -u UNIPROT_DMND_FILE -F -d STR

=head1 COMMANDS

 check    - step 0 to check the dependencies
 assemble - step 1 to assemble raw reads...
 annotate - step 2 to annotate filtered and cleaned scaffold(s)

=head1 DESCRIPTION

 Step 0: check the dependencies

 Step 1: assembles raw  reads, looks  for and removes the  cloning vector,  and
 extracts the longest and the most covered contigs. It has been build to handle
 two  types of  raw  reads  as inputs:  single (454, ion torrent reads, ...) or
 paired (Illumina,...) reads.
 This tools is not able to process PacBio or Oxford Nanopore reads.

 Step 2: annotate filtered and cleaned scaffold(s) provided by the step 1.

=head1 MAIN OPTIONS

=over 8

=item B<--version>

 Print version
 
=item B<-h, --help>

 Print help

=back

=head1 CHECK OPTIONS

=over 8

=item B<-h, --help>

 Print help

=back

=head1 ASSEMBLE OPTIONS

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

 Vector fasta file

=item B<--length> INT

 Fosmid's expected length [40000]

=item B<--minimalContigLength> INT

 Contig's minimum length [1000]

=item B<--minimalContigDepth> INT

 Contig's minimum depth [8]

=item B<-c, --maxDepth> INT

 Coverage, maximum depth use to filter input reads [300]

=item B<-d, --dirOutputs> STR

 Path to the outputs directory

=item B<-H, --htmlOutput> STR

 HTML output name [mintia_assemble.html]

=item B<-L, --logOutput> STR

 Log output file name [mintia_assemble.log]

=item B<-t, --threads>

 number of threads for SPADES [8]

=item B<-h, --help>

 Print help

=back

=head1 ANNOTATE OPTIONS

=over 8

=item B<-i, --input> FILE

 Fasta(.gz) file
 
=item B<-s, --separator> CHAR [#]

 Which separator allows retreiving the fosmid name
 This separator will be also use to create ORF id
 Example: >fosmidName1#contig1... |
          >fosmidName1#contig2... | => fosmidName1
          >fosmidName1#contig3... |

=item B<-n, --nrDB> FILE

 Non-redundant proteins database indexed for Diamond (Ex: the nr.dmnd)
 Required by -F.

=item B<-u, --uniprotDB> FILE

 Proteic sequence database indexed for Diamond (Ex: the uniprot_sprot.dmnd)
 Required by -F and -M.

=item B<-F, --FunctionalAndTaxonomic>

 Run functional and taxonomic annotations
 -n, --nrDB and -u, --uniprotDB must be provided
 
=item B<-e, --evalue> FLOAT

 Maximum diamond e-value to report alignments [10e-8]

=item B<-q, --queryCover> INT

 Minimum diamond query cover% to report an alignment [50]
 
=item B<-M, --Megan> FILE

 Run MEGAN - A license file must be provided
 -n, --nrDB must be provided

=item B<-C, --Cog> FILE

 Run annotations with COGs, DB COGs path file

=item B<-c, --cMaxEvalue> FLOAT

 Max Evalue for Rps-Blast COGs filtering [10e-8]

=item B<-S, --SubmissionFiles>

 Build submission files
 
=item B<-D, --DiamondAgainstPrivateDB> FILE

 Run diamond against your own protein reference FASTA file

=item B<-t, --threads> INT

 Number of threads for Blast [8]
 
=item B<-d, --dirOutputs> STR

 Path to the outputs directory

=item B<-H, --htmlOutput> STR

 HTML output name [mintia_annotate.html]
 
=item B<-L, --logOutput> STR

 Log output file name [mintia_annotate.log]
 
=item B<-k, --keepTmpFiles>

 Keep temporary files

=item B<-h, --help>

 Print help

=back

=head1 VERSION

 Mintia_v1.0

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
use File::Copy;
use Pod::Usage;
use Getopt::Long qw(:config no_ignore_case);
use Term::ANSIColor;
use POSIX;
use 5.010;  #for filesize

#my $DIAMOND_NR_DB  = "/bank/diamonddb/nr.dmnd";
#my $DIAMOND_UP_DB  = "/bank/diamonddb/uniprot_sprot.dmnd";
my $MINTIA_VERSION = "Mintia_v1.0";


########################################################################
# HTML templates
########################################################################
my $HTML_HEADER = '<!doctype html>
<html lang="en">
  <head>
    <!-- Required meta tags -->
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
    <meta name="description" content="###TITLE###">
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
	  hr{
		margin: 0;
		border: 0;  
	  }
	  .igv-popover-table {
		margin: 5px;
		font-size: smaller;
	  }
	  .igv-popover-item {
		font-weight: bold;
		vertical-align: top;
		text-transform: capitalize;
		white-space: nowrap;
	  }
	  .igv-popover-value {
	    overflow: auto;
	    white-space: normal;
	  }
	  .btn-sm {
		  padding: .03rem .5rem 0 .5rem
	  }
	  .textoverflow {
		text-overflow: ellipsis;
		overflow : hidden;
		white-space: nowrap;
		max-width:0;
		width:30%;
	  }
	  .textoverflow:hover {
		text-overflow: clip;
		white-space: normal;
		word-break: break-all;
	  }
	  .zoom-in-notice-container div {
		font-size: small;
		text-transform: full-width;
	  }
    </style>

	<! Highcharts -->
	<script src="https://code.highcharts.com/highcharts.js"></script>
	<script src="https://code.highcharts.com/modules/exporting.js"></script>
	<script src="https://code.highcharts.com/modules/export-data.js"></script>

    <title>###TITLE###</title>
  </head>
  <body>
    <nav class="navbar navbar-dark sticky-top bg-dark flex-md-nowrap p-0">
      <a class="navbar-brand col-sm-3 col-md-2 mr-0" href="#">###TITLE###</a>
    </nav>
    <div class="container-fluid">
      <div class="row">
		###MENU###
		<div style="text-align:center;font-size:smaller;color:darkgrey;margin-top:-25px">
			Produced by '.$MINTIA_VERSION.'<br>
			Copyright © 2018, <img src="http://www.inra.fr/extension/itkinra/design/inra/images/favicon.ico">
			<a style="color:#212529;" href="http://inra.fr" target="_blank">INRA</a><br>
			Designed by the <a style="color:#212529;" href="http://sigenae.org" target="_blank">Sigenae</a> team.
		</div>
	</nav>
	<main role="main" class="col-md-9 ml-sm-auto col-lg-10 pt-3 px-4">
';

my $HTML_FOOTER = '
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
	
	<!-- IGV JS-->
    <script src="https://igv.org/web/release/2.0.1/dist/igv.min.js"></script>
	
	<script type="text/javascript" language="javascript">
		$(function() {
			$(".pop").on("click", function() {
				$("#megantitle").text($(this).find("img").attr("data-title"));
				$("#meganimg").attr("src", $(this).find("img").attr("src"));
				$("#meganview").modal("show");   
			});	
			$(".btn").on("click", function() {
				$(".igv-popover").css("display", "none");
			});	
		});
	</script>
  </body>
</html>
';


########################################################################
# Assemble functions 
########################################################################
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
  open(FILE, "$file") || die "Error: Unabled to open $file\n";
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

  my $tmp1 = File::Temp->new(TEMPLATE=>'tmp1XXXXXXX', DIR=>$outputDir, SUFFIX=>'.fa', UNLINK=>1);
  my $tmp2 = File::Temp->new(TEMPLATE=>'tmp2XXXXXXX', DIR=>$outputDir, SUFFIX=>'.fa', UNLINK=>1);
  my $resC = File::Temp->new(TEMPLATE=>'resCXXXXXXX', DIR=>$outputDir, SUFFIX=>'.fa', UNLINK=>1);

  open(TMP1, ">$tmp1") || die "Error: Unabled to create $tmp1\n";
  open(TMP2, ">$tmp2") || die "Error: Unabled to create $tmp2\n";
  print TMP1 ">seq1\n$seq1\n";
  print TMP2 ">seq2\n$seq2\n";
  close(TMP1);
  close(TMP2);

  # Processing the sequences
  `(cross_match $tmp1 $tmp2 -minmatch 10 -minscore 20 > $resC) >& /dev/null`;
  my @a_res = read_crossmatch($resC);

  # Processing crossmatch result
  my $newcoord = 0;
  for(my $i=0; $i<=$#a_res; $i++) {
    # Check if the overlap is at the end of seq1 and at the start of seq2
	if ((length($seq1)-$a_res[$i][1]) < 10) {
      $newcoord = $a_res[$i][0];
    }
  }
  if($newcoord > 0) { return substr($seq1,0,$newcoord+1).$seq2; }
  #return substr($seq1,0,length($seq1)-$newcoord+1).$seq2; }
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

########################################################################
# Check main
########################################################################
=head2 procedure check

 Usage    : check()
 Function : Step 0 check

=cut
sub check {
	my $help = 0;
	
	GetOptions(
		'h|help' => \$help
	) || pod2usage(-verbose => 99,-sections => "NAME|CHECK SYNOPSIS|CHECK OPTIONS");
	pod2usage(
		-verbose => 99,
		-sections => "NAME|CHECK SYNOPSIS|CHECK OPTIONS"
	) if($help);
	
	print "\n##############################################\n";
	print "        $MINTIA_VERSION check dependencies\n";
	print "##############################################\n";
	# Step1 Assemble
	print "\n- Step 1 - assemble:\n";
	my $err = `which spades.py >& /dev/null`;
	print "  => spades...........";
	if($?) { print colored(['bold red'], 'unavailable in the PATH!', "\n"); }
	else   {
		print colored(['bold green'], "ok");
		my $version = `spades.py --version`;
		chomp($version);
		$version =~ s/^.*\s/version:/;
		print "...$version\n";
	}
	
	$err = `which cross_match >& /dev/null`;
	print "  => cross_match......";
	if($?) { print colored(['bold red'], 'unavailable in the PATH!', "\n"); }
	else   {
		print colored(['bold green'], "ok");
		my $version = `cross_match 2>&1 | head -n3 | tail -n1`;
		chomp($version);
		$version =~ s/^cross_match version\s/version:/;
		print "...$version\n";
	}
	
    # Step2 Annotate
    print "\n- Step 2 - annotate:\n";
	$err = `which prokka >& /dev/null`;
	print "  => prokka...........";
	if($?) { print colored(['bold red'], 'unavailable in the PATH!', "\n"); }
	else   {
		print colored(['bold green'], "ok");
		my $version = `prokka --version 2>&1`;
		chomp($version);
		$version =~ s/^\w+\s+/version:/;
		print "...$version\n";
	}
	
	$err = `which diamond >& /dev/null`;
	print "  => diamond..........";
	if($?) { print colored(['bold red'], 'unavailable in the PATH!', "\n"); }
	else   {
		print colored(['bold green'], "ok");
		my $version = `diamond --version 2> /dev/null`;
		$version =~ s/^.*version\s+/version:/;
		chomp($version);
		print "...$version\n";
	}
	
	$err = `which xvfb-run >& /dev/null`;
	print "  => xvfb-run.........";
	if($?) { print colored(['bold red'], 'unavailable in the PATH!', "\n"); }
	else   { print colored(['bold green'], "ok\n"); }
	
	$err = `which MEGAN >& /dev/null`;
	print "  => MEGAN............";
	if($?) { print colored(['bold red'], 'unavailable in the PATH!', "\n"); }
	else   {
		print colored(['bold green'], "ok");
		my $version = `xvfb-run MEGAN -g --version | head -n1`;
		$version =~ s/^.*version\s([\d\.]+).*$/version:$1/;
		chomp($version);
		print "...$version\n";
	}
	
	$err = `which rpsblast >& /dev/null`;
	print "  => rpsblast.........";
	if($?) { print colored(['bold red'], 'unavailable in the PATH!', "\n"); }
	else   {
		print colored(['bold green'], "ok");
		my $version = `rpsblast -version | head -n1`;
		$version =~ s/^.*\s([\d\.]+).*$/version:$1/;
		chomp($version);
		print "...$version\n";
	}
	
	$err = `which samtools >& /dev/null`;
	print "  => samtools.........";
	if($?) { print colored(['bold red'], 'unavailable in the PATH!', "\n"); }
	else   {
		print colored(['bold green'], "ok");
		my $version = `samtools --version | head -n1`;
		$version =~ s/^.*\s([\d\.]+).*$/version:$1/;
		chomp($version);
		print "...$version\n";
	}
	
	$err = `which tabix >& /dev/null`;
	print "  => tabix............";
	if($?) { print colored(['bold red'], 'unavailable in the PATH!', "\n"); }
	else   {
		print colored(['bold green'], "ok");
		my $version = `tabix 2>&1 | head -n2 | tail -n1`;
		$version =~ s/^.*:\s+(.*)$/version:$1/;
		chomp($version);
		print "...$version\n";
	}
	print "\n";
}

########################################################################
# Assemble main 
########################################################################
=head2 procedure assemble

 Usage    : assemble()
 Function : Step 1 assemble

=cut
sub assemble {
	# Parameters
	my @a_inputSeq  = ();
	my $fosmidLen   = 40000;
	my $minCtgLen   = 1000;
	my $minCtgDepth = 8;
	my $vectorSeq   = undef;
	my $coverage    = 300;
	my $outputDir   = undef;
	my $outputHtml  = "mintia_assemble.html";
	my $outputLog   = "mintia_assemble.log";
	my $threads     = 8;
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
		'length=i'              => \$fosmidLen,
		'minimalContigLength=i' => \$minCtgLen,
		'minimalContigDepth=i'  => \$minCtgDepth,
		'v|vectorSeq=s'         => \$vectorSeq,
		'c|maxDepth=i'          => \$coverage,
		'd|dirOutputs=s'        => \$outputDir,
		'H|htmlOutput=s'        => \$outputHtml,
		'L|logOutput=s'         => \$outputLog,
		't|threads=i'           => \$threads,
		'h|help'                => \$help
	) || pod2usage(-verbose => 99,-sections => "NAME|ASSEMBLE SYNOPSIS|ASSEMBLE OPTIONS");
	pod2usage(
		-verbose => 99,
		-sections => "NAME|ASSEMBLE SYNOPSIS|ASSEMBLE OPTIONS"
	) if($help);
	pod2usage(
		-message => "$0 assemble: '-i, --input' is required.\n",
		-verbose => 99,
		-sections => "NAME|ASSEMBLE SYNOPSIS|ASSEMBLE OPTIONS") if($#a_inputSeq == -1);
	foreach my $i (@a_inputSeq) {
		pod2usage(
			-message => "$0 assemble: $i doesn't exist.\n",
			-verbose => 99,
			-sections => "NAME|ASSEMBLE SYNOPSIS|ASSEMBLE OPTIONS") if !(-e $i);
		pod2usage(
			-message => "$0 assemble: $i is not readable.\n"
			-verbose => 99,
			-sections => "NAME|ASSEMBLE SYNOPSIS|ASSEMBLE OPTIONS") if !(-r $i);
	}
	pod2usage(
		-message => "$0 assemble: '-v, --vectorSeq' is required.\n",
		-verbose => 99,
		-sections => "NAME|ASSEMBLE SYNOPSIS|ASSEMBLE OPTIONS") if !defined($vectorSeq);
	pod2usage(
		-message => "$0 assemble: $vectorSeq doesn't exist.\n",
		-verbose => 99,
		-sections => "NAME|ASSEMBLE SYNOPSIS|ASSEMBLE OPTIONS") if !(-e $vectorSeq);
	pod2usage(
		-message => "$0 assemble: $vectorSeq is not readable.\n",
		-verbose => 99,
		-sections => "NAME|ASSEMBLE SYNOPSIS|ASSEMBLE OPTIONS") if !(-r $vectorSeq);
	pod2usage(
		-message => "$0 assemble: '-d, --dirOutputs' is required.\n",
		-verbose => 99,
		-sections => "NAME|ASSEMBLE SYNOPSIS|ASSEMBLE OPTIONS") if !defined($outputDir);
			
	if(! -e $outputDir) { mkdir $outputDir || die "Error: Unabled to create output dir $outputDir."; }
	open(LOG, ">$outputDir/$outputLog")    || die "Error: Unabled to create $outputDir/$outputLog";

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
    	open(READ1, "$zip $h_sample{$k}{'R1'} |") || die "Error: Unabled to open $h_sample{$k}{'R1'}.";

		if(exists($h_sample{$k}{"R2"})) {
		  $zip = "cat ";
		  $filetype = `file -bsiL $h_sample{$k}{'R2'}`;
		  if($filetype =~ /application\/x\-gzip/) { $zip = "gunzip -c "; }
		  open(READ2, "$zip $h_sample{$k}{'R2'} |") || die "Error: Unabled to open $h_sample{$k}{'R2'}.";
		}

		#outputDir
		open(READ1FILTER, ">$outputDir/$k"."_R1.fq") || die "Error: Unabled to create $outputDir/$k"."_R1.fq";
		$h_sample{$k}{"R1_F"} = "$outputDir/$k"."_R1.fq";
		if(exists($h_sample{$k}{"R2"})) {
		  open(READ2FILTER, ">$outputDir/$k"."_R2.fq") || die "Error: Unabled to create $outputDir/$k"."_R2.fq";
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
			`spades.py -1 $h_sample{$k}{'R1_F'} -2 $h_sample{$k}{'R2_F'} -t $threads --careful -o $outputDir/$k`;
		}
		else {
			`spades.py -s $h_sample{$k}{'R1_F'} -t $threads --careful -o $outputDir/$k`;
		}
		print LOG "done\n";
	}

  #CROSS MATCH (on the scaffolds.fasta generated by SPADES)
  print LOG "\n## Run Crossmatch\n";
  foreach my $k (sort keys(%h_sample)) {
    print LOG " - $k...";
    `(cross_match $outputDir/$k/scaffolds.fasta $vectorSeq -minmatch 9 -minscore 30 -screen > $outputDir/$k/crossmatch-scaffolds.stdout) >& $outputDir/$k/crossmatch-scaffolds.stderr`;
    print LOG "done\n";
  }

  #Find longuest scaffold and
  #Filter scaffold using --minimalContigLength and Contig's minimum length [1000]
  print LOG "\n## Filter scaffolds\n";
  foreach my $k (sort keys(%h_sample)) {
    print LOG " - $k...";
    open(SCAF, "$outputDir/$k/scaffolds.fasta.screen") || die "Error: Unabled to open $outputDir/$k/scaffolds.fasta.screen.";
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

  # Result on stdout AND create output file (cleaned scaffolds input of annotate module)
  open (CLEAN, ">$outputDir/mintia_assemble.fasta") || die "Error: Unabled to create $outputDir/mintia_assemble.fasta\n";
  my $assembleFasta = "";
  print LOG "##\n#Sample\tPaired\tNb_scaffold\tNb_Filtered\t(ID\tLength\tCoverage\tXXed_Block\tClean_length)xN\n";
  foreach my $k1 (sort keys(%h_sample)) {
    print LOG "$k1\t";                                       #Sample
    if(exists($h_sample{$k1}{"R2"})) { print LOG "Yes\t"; }  #Paired
    else                             { print LOG "No\t";  }
    print LOG $h_sample{$k1}{"nbScaffold"}."\t";             #Nb_scaffold
    print LOG "".(keys %{$h_sample{$k1}{"scaffolds"}})."\t"; #Nb_Filtered

    #By filtered scaffold
    foreach my $k2 (sort {$a <=> $b} keys(%{$h_sample{$k1}{"scaffolds"}})) {
	  #If seq lenght>0
	  if(length($h_sample{$k1}{"scaffolds"}{$k2}{"newseq"})>0) {
		print CLEAN ">$k1#$k2";
		print CLEAN "#".length($h_sample{$k1}{"scaffolds"}{$k2}{"newseq"})."\n";
		print CLEAN join("\n", $h_sample{$k1}{"scaffolds"}{$k2}{"newseq"}=~/(.{1,60})/g)."\n";
	  }
      print LOG "$k2\t";                                              #ID
      print LOG $h_sample{$k1}{"scaffolds"}{$k2}{"len"}."\t";         #Length
      print LOG $h_sample{$k1}{"scaffolds"}{$k2}{"cov"}."\t";         #Coverage
      print LOG $h_sample{$k1}{"scaffolds"}{$k2}{"XX"}."\t";          #XXed_Block
			if($h_sample{$k1}{"scaffolds"}{$k2}{"newseq"} =~ /Error/) { #Clean_length
				print LOG "Error - Unexpected number of XXed blocks\t";
			}
			else {
				print LOG length($h_sample{$k1}{"scaffolds"}{$k2}{"newseq"})."\t";
			}
    }
    print LOG "\n";
  }
  close CLEAN;
  close(LOG);


	##############################
	## HTML Report
	##############################
	open (HTML, ">$outputDir/$outputHtml") || die "Error: Unabled to create $outputDir/$outputHtml\n";
	print HTML '<!doctype html>
<html lang="en">
  <head>
    <!-- Required meta tags -->
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
    <meta name="description" content="Mintia assemble module">
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
    </style>

    <title>Mintia assemble report</title>
  </head>
  <body>
    <nav class="navbar navbar-dark sticky-top bg-dark flex-md-nowrap p-0">
      <a class="navbar-brand col-sm-3 col-md-2 mr-0" href="#">Mintia assemble report</a>
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
									<span class="oi oi-eye" aria-hidden="true"></span>
									Assemble and vector detection
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
          </div>';
	print HTML "<div style=\"text-align:center;font-size:smaller;color:darkgrey;margin-top:-25px\">
		Produced by $MINTIA_VERSION<br>
		Copyright © 2018, <img src=\"http://www.inra.fr/extension/itkinra/design/inra/images/favicon.ico\">
		<a style=\"color:#212529;\" href=\"http://inra.fr\" target=\"_blank\">INRA</a><br>
		Designed by the <a style=\"color:#212529;\" href=\"http://sigenae.org\" target=\"_blank\">Sigenae</a> team.</div>";
	print HTML '
        </nav>

				<main role="main" class="col-md-9 ml-sm-auto col-lg-10 pt-3 px-4">
          <div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center pb-2 border-bottom">
						<h1 class="h4">Inputs and parameters</h1>
						<span class="anchor" id="inputs-parameters"></span>
					</div>';
	print HTML '
					<div class="d-flex">
					  <div class="mt-4 mr-4 pl-0 col-md-4">
							<h5>Parameters</h5>
							<span class="anchor" id="parameters"></span>
				    	<ul class="list-group">';
	print HTML "<li class=\"list-group-item d-flex justify-content-between align-items-center\">Fosmid\'s expected length <span class=\"badge badge-danger badge-pill ml-4\">$fosmidLen"."bp</span></li>";
	print HTML "<li class=\"list-group-item d-flex justify-content-between align-items-center\">Maximum coverage <span class=\"badge badge-danger badge-pill ml-4\">$coverage</span></li>";
	print HTML "<li class=\"list-group-item d-flex justify-content-between align-items-center\">Contig\'s minimum length <span class=\"badge badge-danger badge-pill ml-4\">$minCtgLen"."bp</span></li>";
	print HTML "<li class=\"list-group-item d-flex justify-content-between align-items-center\">Contig\'s minimum depth <span class=\"badge badge-danger badge-pill ml-4\">$minCtgDepth</span></li>";
	print HTML	'
					    </ul>
						</div>

						<div class="mt-4">
							<h5>Cloning vector</h5>
							<span class="anchor" id="cloning-vector"></span>';
	open(VECTOR, "$vectorSeq") || die "Error: Unabled to open $vectorSeq\n";
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
							  <ul class="list-group d-flex flex-row flex-wrap" style="width:50%;">';
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
        		<tbody>';
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
						<h1 class="h4">Assemble and vector detection</h1>
						<span class="anchor" id="assembly-detection"></span>
          </div>
					<table class="table table-striped table-bordered mt-4 mb-0" style="width:100%">
        		<thead>
            	<tr>
                <th>Sample</th>
                <th nowrap style="text-align:center">Scaffold</th>
								<th nowrap style="text-align:center">Filtered</th>
                <th nowrap style="text-align:center">ID</th>';
	print HTML "<th nowrap style='text-align:center'>Length (>$minCtgLen)</th>\n";
	print HTML "<th nowrap style='text-align:center'>Coverage (>$minCtgDepth)</th>\n";
	print HTML '
								<th nowrap style="text-align:center;width:275px;">Vector position</th>
								<th nowrap style="text-align:center">Clean length</th>
            	</tr>
        		</thead>
        		<tbody>';

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
	my $download = "var dwnl = [";
	my $allfasta = "var allfasta = \"";
	my $allgff   = "var allgff = \"";
	my $allclean = "var allclean = \"";
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
				$seq .= "$MINTIA_VERSION\\tBlockXX_$cmp\\t"; #source feature
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
			if(length($h_sample{$k1}{"scaffolds"}{$k2}{"newseq"})>0) {
				$seq .= ">$k1#$k2";
				$seq .= "#".length($h_sample{$k1}{"scaffolds"}{$k2}{"newseq"})."\\n";
				$seq .= join("\\n", $h_sample{$k1}{"scaffolds"}{$k2}{"newseq"}=~/(.{1,60})/g)."\\n";
			}
		}
		if($dwnlid != 0) { $download .= ",";}
		$download .= "\"$seq\"";
		$allclean .= "$seq";
		print HTML "<td style=\"text-align:right\"><button id=\"dwnl-$dwnlid\" file=\"$k1.clean.fasta\" class=\"btn btn-sm btn-outline-secondary\">
								Fasta (".(keys %{$h_sample{$k1}{"scaffolds"}}).")</button></td>\n";
		print HTML "</tr>\n";
		$dwnlid++;
	}
	$download .= "];";
	print HTML "\n<script>$download</script>\n";
	print HTML "\n<script>$allfasta\"</script>\n";
	print HTML "\n<script>$allgff\"</script>\n";
	print HTML "\n<script>$allclean\"</script>\n";
	print HTML '		<tr>
							<td></td>';
	print HTML "\n<td style=\"text-align:center\"><button id=\"allfasta\" file=\"all.fasta\" class=\"btn btn-sm btn-outline-secondary\">
								All Fasta</button></td>\n";
	print HTML "<td style=\"text-align:center\"><button id=\"allgff\" file=\"all.gff\" class=\"btn btn-sm btn-outline-secondary\">
								All GFF</button></td>\n";
	print HTML "<td style=\"text-align:center\"><button id=\"allclean\" file=\"allcleaned.fasta\" class=\"btn btn-sm btn-outline-secondary\">
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
	title: {';
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
	print HTML "text: \'Produced by $MINTIA_VERSION\'";
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
		$(\'[id=allclean]\').click(function(){
			download($(this).attr(\'file\'), allclean);
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









########################################################################
# Annotate functions 
########################################################################

=head2 function split_fasta_by_id

 Usage    : split_fasta_by_id( $fastaFile, $separator, $ouputDir, $fastaName )
 Function : Split a multifasta file by id and create one folder by id and the
            corresponding $fastaName file
			Id will be selected between > and separator or end of line
 Return   : Array of id
 Args     : [str] fastaFile
			[str] separator
            [str] outputDir
            [str] fastaName

=cut
sub split_fasta_by_id {
  my ($fastaFile, $separator, $outputDir, $fastaName) = @_;
  
  my @a_id;
  open(IN, "$fastaFile") || die "Error: Unabled to open $fastaFile\n";
  my $curId = "";
  while(my $line = <IN>) {
	  if($line =~ /^>(.*)$/) {
		my $id = $1;
		if($separator ne "" && index($id, $separator)!=-1) {
			$id = substr($id, 0, index($id, $separator));
		}
		if($curId ne $id) {
			close OUT if($curId ne "");
			$curId = $id;
			if(! -e "$outputDir/$curId") { mkdir "$outputDir/$curId" || die "Error: Unabled to create output dir $outputDir/$curId."; }
			open(OUT, ">$outputDir/$curId/$fastaName") || die "Error: Unabled to create $outputDir/$curId/$fastaName\n";
			push(@a_id, $curId);
		}
	  }
	  print OUT $line;
  }
  close IN;
  close OUT;
  return @a_id;
}


=head2 function split_xml_by_id

 Usage    : split_xml_by_id( $xmlFile, $separator, $outputDir, $xmlName )
 Function : Split a xml file by id and create one folder by id and the
            corresponding $xmlName file
			Id will be selected between > and separator or end of line
 Return   : Array of id
 Args     : [str] xmlFile
			[str] separator
			[str] outputDir
			[str] xmlName

=cut
sub split_xml_by_id {
	my ($xmlFile, $separator, $outputDir, $xmlName) = @_;
  
	my @a_id;
	open(IN, "$xmlFile") || die "Error: Unabled to open $xmlFile\n";
	my $header = "";
	my @a_hits = ();
	my $footer = "";
	my $line = <IN>;
	while($line !~ /<Iteration>/) {
		$header .= $line if($line !~ /<BlastOutput_query-/);
		$line = <IN>;
	}
	while($line !~ /<\/BlastOutput_iterations>/) {
		my $hit = $line;
		$line = <IN>;
		while($line !~ /<Iteration>/ && $line !~ /<\/BlastOutput_iterations>/) {
			$hit .= $line;
			$line = <IN>;
		}
		push(@a_hits, $hit);
	}
	$footer = $line;
	while($line = <IN>) {
		$footer .= $line;
	}
 	close IN;
 	
 	my $curId  = "";
 	foreach my $hit (@a_hits) {
		if($hit =~ /<Iteration_query-def>(.*)<\/Iteration_query-def>/) {
			my $id = $1;
			if($separator ne "" && index($id, $separator)!=-1) {
				$id = substr($id, 0, index($id, $separator));
			}
			if($curId ne $id) {
				if($curId ne "") {
					print OUT $footer;
					close OUT;
				}
				$curId = $id;
				
				if(! -e "$outputDir/$curId") { mkdir "$outputDir/$curId" || die "Error: Unabled to create output dir $outputDir/$curId."; }
				open(OUT, ">$outputDir/$curId/$xmlName") || die "Error: Unabled to create $outputDir/$curId/$xmlName\n";
				print OUT $header;
				push(@a_id, $curId);
			}
			print OUT $hit;
		}
	}
	print OUT $footer;
	close OUT;
	return @a_id;
}

########################################################################
# Annotate main 
########################################################################
=head2 function annotate

 Usage    : annotate()
 Function : Step 2 annotate

=cut
sub annotate {
	# Parameters
	my $inputSeq     = "";
	my $dbNR         = "";
	my $dbUniP       = "";
	my $separator    = "#";
	my $outputDir    = undef;
	my $outputHtml   = "mintia_annotate.html";
	my $outputLog    = "mintia_annotate.log";
	my $help         = 0;
	my $funAndTaxo   = 0;
	my ($diamond_evalue, $diamond_queryCover) = (10e-8, 50);
	my $megan        = 0;
	my ($cog,  $cog_cMaxEvalue) = (0, 10e-8);
	my ($diamond, $ncbi)        = (0, 0);
	my $threads      = 8;
	my $keep         = 0;

	GetOptions(
		'i|input=s{1}'                => \$inputSeq,
		's|separator=s'               => \$separator,
		'n|nrDB=s'                    => \$dbNR,
		'u|uniprotDB=s'               => \$dbUniP,
		'F|FunctionalAndTaxonomic'    => \$funAndTaxo,
		'e|evalue=f'                  => \$diamond_evalue,
		'q|queryCover=i'              => \$diamond_queryCover,
		'M|Megan=s'                   => \$megan,
		'C|Cog=s'                     => \$cog,
		'c|cMaxEvalue=f'              => \$cog_cMaxEvalue,
		'S|SubmissionFiles'           => \$ncbi,
		'D|DiamondAgainstPrivateDB=s' => \$diamond,
		't|threads=i'                 => \$threads,
		'd|dirOutputs=s'              => \$outputDir,
		'H|htmlOutput=s'              => \$outputHtml,
		'L|logOutput=s'               => \$outputLog,
		'k|keepTmpFiles'              => \$keep,
		'h|help'                      => \$help
	) || pod2usage(-verbose => 99,-sections => "NAME|ANNOTATE SYNOPSIS|ANNOTATE OPTIONS");
	pod2usage(
		-verbose => 99,
		-sections => "NAME|ANNOTATE SYNOPSIS|ANNOTATE OPTIONS"
	) if($help);
	pod2usage(
		-message => "$0 annotate: '-i, --input' is required.\n",
		-verbose => 99,
		-sections => "NAME|ANNOTATE SYNOPSIS|ANNOTATE OPTIONS") if ($inputSeq eq "");
	pod2usage(
		-message => "$0 annotate: $inputSeq doesn't exist.\n",
		-verbose => 99,
		-sections => "NAME|ANNOTATE SYNOPSIS|ANNOTATE OPTIONS") if !(-e $inputSeq);
	pod2usage(
			-message => "$0 annotate: $inputSeq is not readable.\n"
			-verbose => 99,
			-sections => "NAME|ANNOTATE SYNOPSIS|ANNOTATE OPTIONS") if !(-r $inputSeq);
	if($dbNR ne "") {
		pod2usage(
			-message => "$0 annotate: $dbNR doesn't exist.\n",
			-verbose => 99,
			-sections => "NAME|ANNOTATE SYNOPSIS|ANNOTATE OPTIONS") if !(-e $dbNR);
		pod2usage(
			-message => "$0 annotate: $dbNR is not readable.\n"
			-verbose => 99,
			-sections => "NAME|ANNOTATE SYNOPSIS|ANNOTATE OPTIONS") if !(-r $dbNR);
	}
	if($dbUniP ne "") {
		pod2usage(
			-message => "$0 annotate: $dbUniP doesn't exist.\n",
			-verbose => 99,
			-sections => "NAME|ANNOTATE SYNOPSIS|ANNOTATE OPTIONS") if !(-e $dbUniP);
		pod2usage(
			-message => "$0 annotate: $dbUniP is not readable.\n"
			-verbose => 99,
			-sections => "NAME|ANNOTATE SYNOPSIS|ANNOTATE OPTIONS") if !(-r $dbUniP);
	}	
	if($funAndTaxo) {
		pod2usage(
			-message => "$0 annotate: -n, --nrDB is required by -F.\n",
			-verbose => 99,
			-sections => "NAME|ANNOTATE SYNOPSIS|ANNOTATE OPTIONS") if ($dbNR eq "");
		pod2usage(
			-message => "$0 annotate: -u, --uniprotDB is required by -F.\n",
			-verbose => 99,
			-sections => "NAME|ANNOTATE SYNOPSIS|ANNOTATE OPTIONS") if ($dbUniP eq "");
	}
	if($megan) {
		pod2usage(
			-message => "$0 annotate: $megan doesn't exist.\n",
			-verbose => 99,
			-sections => "NAME|ANNOTATE SYNOPSIS|ANNOTATE OPTIONS") if !(-e $megan);
		pod2usage(
			-message => "$0 annotate: $megan is not readable.\n"
			-verbose => 99,
			-sections => "NAME|ANNOTATE SYNOPSIS|ANNOTATE OPTIONS") if !(-r $megan);
		pod2usage(
			-message => "$0 annotate: -n, --nrDB is required by -M.\n",
			-verbose => 99,
			-sections => "NAME|ANNOTATE SYNOPSIS|ANNOTATE OPTIONS") if ($dbNR eq "");
	}
	if($cog) {
		pod2usage(
			-message => "$0 annotate: $cog doesn't exist.\n",
			-verbose => 99,
			-sections => "NAME|ANNOTATE SYNOPSIS|ANNOTATE OPTIONS") if !(-e "$cog.aux");
		pod2usage(
				-message => "$0 annotate: $cog is not readable.\n"
				-verbose => 99,
				-sections => "NAME|ANNOTATE SYNOPSIS|ANNOTATE OPTIONS") if !(-r "$cog.aux");
	}
	pod2usage(
		-message => "$0 annotate: '-d, --dirOutputs' is required.\n",
		-verbose => 99,
		-sections => "NAME|ANNOTATE SYNOPSIS|ANNOTATE OPTIONS") if !defined($outputDir);
	pod2usage(
		-message => "$0 annotate: '-s, --separator' must be a CHAR ($separator).\n",
		-verbose => 99,
		-sections => "NAME|ANNOTATE SYNOPSIS|ANNOTATE OPTIONS") if(length($separator)!=1);
	
	if(! -e $outputDir) { mkdir $outputDir || die "Error: Unabled to create output dir $outputDir."; }
	open(LOG, ">$outputDir/$outputLog")    || die "Error: Unabled to create $outputDir/$outputLog";
	
	####
	## Check input fasta file and create output folder tree structure (one dir by fosmid/fastaId)
	####
	print LOG "## Check input fasta file\n";
	my %h_fasta = ();
	my $zip = "cat ";
	my $filetype = `file -bsiL $inputSeq`;
	if($filetype =~ /application\/x\-gzip/) { $zip = "gunzip -c "; }
    open(INPUT, "$zip $inputSeq |") || die "Error: Unabled to open $inputSeq.";
	my $currentName = "";
	while(my $line = <INPUT>) {
		if($line =~ /^>(.*)$/) {
			$currentName = $1;
		}
		elsif($line =~ /^[ATCGNXatgcnx]+$/) {
			chomp $line;
			$h_fasta{$currentName} .= $line;
		}
		else {
			pod2usage("Error: $inputSeq is not a valid fasta file.")
		}
	}
	close(INPUT);
	my @a_fosmid = split_fasta_by_id($inputSeq, $separator, $outputDir, "fosmid.fasta");

	####
	## Run prokka, parse and generate output files (gff, nuc and prot)
	####
	my $prokkaGFF = "$outputDir/prokka.gff";
	my $prokkaTFA = "$outputDir/prokka.fasta";
	my $prokkaPRO = "$outputDir/prokka-prot.fasta";
	
	# %h_orf = (
	#    fastaID => {
	#          ORF      => 'Number of prokka ORF',
	#          uniORF   => 'Number of annotated ORF against Uniprot',
	#          nrORF    => 'Number of annotated ORF against NR',
	#          ORFName1 => {
	#             "Start"  => 'Start coord on fasta',
	#             "Length" => 'length of the orf',
	#             "Strand" => 'Strand coord on fasta',
	
	#             "uniNbHSP"  => 'Number of HSP from diamond VS Uniprot',                    |
	#             "uniMaxCov" => 'Max HSPcov from diamond VS Uniprot: qcovhsp',              |
	#             "uniCov"    => 'HSP Coverage'},    |                                       | => From Uniprot
	#             "uniIden"   => 'HSP % identity'},  | => For the best HSP (max of cov*Iden) |
	#             "uniSpecies"=> 'HSP Species'}      |                                       |
	
	#             "nrNbHSP"  => 'Number of HSP from diamond VS NR',                          |
	#             "nrMaxCov" => 'Max HSPcov from diamond VS Uniprot: qcovhsp',               |
	#             "nrCov"    => 'HSP Coverage'},    |                                        | => From NR
	#             "nrIden"   => 'HSP % identity'},  | => For the best HSP (max of cov*Iden)  |
	#             "nrSpecies"=> 'HSP Species'}      |                                        |
	#          }
	#          ...
	#    }
	# )
	my %h_orf = ();
	
	print LOG "## Run prokka\n";
	`prokka $inputSeq --outDir $outputDir --cpu $threads --prefix prokka_tmp_ --locustag "##TEMPLATE##" --force >& $outputDir/prokka_tmp_.stderr`;
	my $orfcmp      = 1;
	my $name        = "";
	my $prokkagff   = "";
	my %h_prokka2mintia;
	open(GFF, ">$prokkaGFF") || die "Error: Unabled to create $prokkaGFF";
	open(PROKKAGFF, "$outputDir/prokka_tmp_.gff") || die "Error: Unabled to open $outputDir/prokka_tmp_.gff";
	while(my $line=<PROKKAGFF>) {
		if($line !~ /##TEMPLATE##/) {
			print GFF $line if($line =~ /^##/);
		}
		else {
			my @a_line = split(/\t/, $line);
			 # New contig ?
			if($name ne $a_line[0]) {
				if($name ne "") {
					$orfcmp--;
					print LOG " - $name: $orfcmp ORF(s) found\n";
					$h_orf{$name}{"ORF"}    = $orfcmp;
					$h_orf{$name}{"uniORF"} = 0;
					$h_orf{$name}{"nrORF"}  = 0;
				}
				$orfcmp = 1;
				$name = $a_line[0];
			}
			$line =~s/(##TEMPLATE##\_\d+);/ORF$orfcmp;/g;
			$h_prokka2mintia{$1} = "$name$separator"."ORF$orfcmp";
			print GFF $line;
			$h_orf{$name}{"ORF$orfcmp"}{"Start"}      = $a_line[3];
			$h_orf{$name}{"ORF$orfcmp"}{"Length"}     = ($a_line[4]-$a_line[3]+1);
			$h_orf{$name}{"ORF$orfcmp"}{"Strand"}     = 0;
			$h_orf{$name}{"ORF$orfcmp"}{"Strand"}     = 16 if($a_line[6] eq "-");
			$h_orf{$name}{"ORF$orfcmp"}{"uniNbHSP"}   = 0;
			$h_orf{$name}{"ORF$orfcmp"}{"uniMaxCov"}  = 0;
			$h_orf{$name}{"ORF$orfcmp"}{"uniCov"}     = 0;
			$h_orf{$name}{"ORF$orfcmp"}{"uniIden"}    = 0;
			$h_orf{$name}{"ORF$orfcmp"}{"uniSpecies"} = "";
			$h_orf{$name}{"ORF$orfcmp"}{"nrNbHSP"}    = 0;
			$h_orf{$name}{"ORF$orfcmp"}{"nrMaxCov"}   = 0;
			$h_orf{$name}{"ORF$orfcmp"}{"nrCov"}      = 0;
			$h_orf{$name}{"ORF$orfcmp"}{"nrIden"}     = 0;
			$h_orf{$name}{"ORF$orfcmp"}{"nrSpecies"}  = "";
			$orfcmp++;
		}
	}
	close(PROKKAGFF);
	close(GFF);
	$orfcmp--;
	print LOG " - $name: $orfcmp ORF(s) found\n";
	$h_orf{$name}{"ORF"}    = $orfcmp;
	$h_orf{$name}{"uniORF"} = 0;
	$h_orf{$name}{"nrORF"}  = 0;
	print LOG " ==> Create $prokkaGFF..........done\n";
	
	open(NUC, ">$prokkaTFA") || die "Error: Unabled to create $prokkaTFA";
	open(PROKKAFFN, "$outputDir/prokka_tmp_.ffn") || die "Error: Unabled to open $outputDir/prokka_tmp_.ffn";
	while(my $line=<PROKKAFFN>) {
		if($line =~ /^>(\S+)/) {
			$line=~s/>(\S+)/>$h_prokka2mintia{$1}/;
		}
		print NUC $line;
	}
	close(PROKKAFFN);
	close(NUC);
	print LOG " ==> Create $prokkaTFA........done\n";

	open(PRO, ">$prokkaPRO") || die "Error: Unabled to create $prokkaPRO";
	open(PROKKAFAA, "$outputDir/prokka_tmp_.faa") || die "Error: Unabled to open $outputDir/prokka_tmp_.faa";
	while(my $line=<PROKKAFAA>) {
		if($line =~ /^>(\S+)/) {
			$line=~s/>(\S+)/>$h_prokka2mintia{$1}/;
		}
		print PRO $line;
	}
	close(PROKKAFAA);
	close(PRO);
	print LOG " ==> Create $prokkaPRO...done\n";
	
	if($ncbi) {
		print LOG " ==> Create submission file(s) by fosmid.......";
		open(PROKKAGBK, "$outputDir/prokka_tmp_.gbk") || die "Error: Unabled to open $outputDir/prokka_tmp_.gbk";
		while(my $line=<PROKKAGBK>) {
			if($line =~ /^LOCUS/) {
				my @a_line = split(/\s+/, $line);
				
				my $id = "";
				if($separator ne "" && index($a_line[1], $separator)!=-1) {
					$id = substr($a_line[1], 0, index($a_line[1], $separator));
				}
				else {
					$id = $a_line[1];
				}
								
				my $outfile = "$outputDir/$id/prokka.gbk";
				open(NCBI, ">>$outfile")  || die "Error: Unabled to create $outfile";
				print NCBI $line;
			}
			elsif($line =~ /^\/\/$/) {
				print NCBI $line;
				close(NCBI);
			}
			elsif($line =~ /\/locus_tag=/) {
				$line=~s/locus_tag="(\S+)"/locus_tag="$h_prokka2mintia{$1}"/;
				print NCBI $line;
			}
			else { print NCBI $line; }
		}
		close(PROKKAGBK);
		print LOG "done\n";
	}
		
	####
	## Contigs and ORFs diamond blastx against nr and uniprot to identity missed ORFs by prokka
	## Rq: diamond param -k 100000 --max-hsps... 
	####
	if($funAndTaxo) {
		print LOG "## Run functional and taxonomic annotation\n";
		print LOG " - Run diamond-blastx $inputSeq against NR........";
		`diamond blastx --db $dbNR --query $inputSeq --threads $threads --outfmt 6 qseqid sseqid pident nident length mismatch gaps gapopen qstart qend sstart send evalue bitscore stitle qcovhsp -k 1000 -e $diamond_evalue --out $outputDir/contigs-nr.diamond.tsv >& $outputDir/contigs-nr.diamond_tmp_.stderr`;

		print LOG "done\n - Run diamond-blastx $inputSeq against Uniprot...";
		`diamond blastx --db $dbUniP --query $inputSeq --threads $threads --outfmt 6 qseqid sseqid pident nident length mismatch gaps gapopen qstart qend sstart send evalue bitscore stitle qcovhsp -k 1000 -e $diamond_evalue --out $outputDir/contigs-uniprot.diamond.tsv >& $outputDir/contigs-uniprot.diamond_tmp_.stderr`;
		
		print LOG "done\n - Run diamond-blastx $prokkaTFA against NR........";
		`diamond blastx --db $dbNR --query $prokkaTFA --threads $threads --outfmt 6 qseqid sseqid pident nident length mismatch gaps gapopen qstart qend sstart send evalue bitscore stitle qcovhsp -k 1000 --query-cover $diamond_queryCover -e $diamond_evalue --out $outputDir/prokka-nr.diamond.tsv >& $outputDir/prokka-nr.diamond_tmp_.stderr`;
		
		print LOG "done\n - Run diamond-blastx $prokkaTFA against Uniprot...";
		`diamond blastx --db $dbUniP --query $prokkaTFA --threads $threads --outfmt 6 qseqid sseqid pident nident length mismatch gaps gapopen qstart qend sstart send evalue bitscore stitle qcovhsp -k 1000 --query-cover $diamond_queryCover -e $diamond_evalue --out $outputDir/prokka-uniprot.diamond.tsv >& $outputDir/prokka-uniprot.diamond_tmp_.stderr`;
		print LOG "done\n";
	}
	
	# Split fasta by fosmid need for HTML report (and Megan)
	my @a_fosmidFromprokkaTfa   = split_fasta_by_id($prokkaTFA, $separator, $outputDir, "prokka.fasta");
	
	####
	## MEGAN must be run by fosmid and takes as input:
	## - fosmid orf fasta file (extract from prokka.fasta)
	## - fosmid diamond-blastx of ORFs against nr (extract from prokka-nr.diamond.xml)
	## MEGAN produces a .rma file and a jpg tree file 
	####
	if($megan) {
		print LOG "## Run MEGAN\n";
		print LOG " - Run diamond-blastx $prokkaTFA against NR........";
		`diamond blastx --db $dbNR --query $prokkaTFA --threads $threads --outfmt 5 -k 1000 --query-cover $diamond_queryCover -e $diamond_evalue --out $outputDir/prokka-nr.diamond.xml >& $outputDir/prokka-nr.diamond_tmp_.stderr`;
		print LOG "done\n";
				
		# Split fasta and xml by fosmid to run MEGAN by fosmid
		print LOG " - Prepare inputs files: split FASTA and XML files by fosmid...";
		my @a_fosmidFromDiamondXml  = split_xml_by_id("$outputDir/prokka-nr.diamond.xml", $separator, $outputDir, "diamond-prokkaVSnr_tmp_.xml");
			
		### TODO valider que les tableaux sont identiques !!!!
		print LOG "done\n";
		
		for(my $i=0;$i<=$#a_fosmid;$i++) {
			# Build config file
			my $path = "$outputDir/$a_fosmid[$i]";
			open(MEGCFG, ">$path/megan_tmp_.cfg") || die "Error: Unabled to create $path/megan.cfg\n";
			print MEGCFG "import blastFile='$path/diamond-prokkaVSnr_tmp_.xml' fastaFile='$path/prokka.fasta' meganFile='$path/megan_tmp_.rma';\n"
						. "set windowSize=1200 x 100;\n"
						. "collapse rank='Species';\n"
						. "select nodes='all';\n"
						. "nodeLabels assigned='true';\n"
						. "nodeLabels summarized='true';\n"
						. "set fillColor='red';\n"
						. "set font='arial-italic-12';\n"
						. "set magnifier='true';\n"
						. "select nodes='none';\n"
						. "zoom full;\n"
						. "exportimage file='$path/megan.png' format='png' replace='true' textAsShapes='false' title='$a_fosmid[$i]';\n"
						. "export what=tree file='$path/megan.newick' simplify='false' showInternalLabels='true' showUnassigned=true;\n"
						. "quit;\n";
			close MEGCFG;
			print LOG " - Run MEGAN on $a_fosmid[$i]...";
			`(xvfb-run MEGAN -g -E -L $megan -c $path/megan_tmp_.cfg) >& $path/megan_tmp_.log`;
			if(-e "$path/megan.png") {}
			else {`(xvfb-run MEGAN -g -E -L $megan -c $path/megan_tmp_.cfg) >& $path/megan_tmp_.log`;}
			
			print LOG "done\n";
		}
	}
	
	
	# %h_nbOrfByCog = (
	#    fosmidID => {
	#          cogcategory => {
	#             "COGs" => [],
	#             "ORFs" => []
	#          }
	# })
	my %h_nbOrfByCog = ();
	if($cog) {
		print LOG "## Run rpsblast (COGs)...";
		`(rpsblast -query $prokkaPRO -db $cog -num_threads $threads -evalue $cog_cMaxEvalue -out $outputDir/cogs_$cog_cMaxEvalue.rpsblast) >& $outputDir/cogs_$cog_cMaxEvalue.rpsblast_tmp_.stderr`;
		
		my $curFosId = "";
		my $curOrfId = "";
		my $firstHit = 1;
		open(COGOUT, "$outputDir/cogs_$cog_cMaxEvalue.rpsblast") || die "Error: Unabled to open $outputDir/cogs_$cog_cMaxEvalue.rpsblast\n";
		while(my $line=<COGOUT>) {	
			if($line =~ /^Query=\s+(.*)#(ORF\d+)/) {
				$curFosId = $1;
				$curOrfId = $2;
				$firstHit = 1;
			}
			elsif($firstHit && $line =~ /^>.*(COG\d+)*(COG\d+,\s+.*)$/) {
				my $cogdescr = $2;
				# Catch full description
				while($line = <COGOUT>) {
					if($line !~/Length\s*=/) {
						$line =~ s/^\s+/ /;
						$line =~ s/\s+$//;
						$cogdescr .= $line;
					}
					else {
						$cogdescr=~s/\s+/ /g;
						$cogdescr=~s/\.$//;
						last;
					}
				}
				# Find category
				my $cat = "";
				if($cogdescr =~ /\[(.*)\]/) {
					$cat = $1;
					$cat =~s/\s+\/.*$//;
					$cat=~s/\s+/ /g;
				}
				my $search = quotemeta "[".$cat."]";
				$cogdescr =~ s/$search//;
				push(@{$h_nbOrfByCog{$curFosId}{$cat}{"COGs"}}, $cogdescr);
				push(@{$h_nbOrfByCog{$curFosId}{$cat}{"ORFs"}}, $curOrfId);
				$firstHit = 0;
			}
			elsif($line =~ /No hits found/) {
				push(@{$h_nbOrfByCog{$curFosId}{"No hits found"}{"ORFs"}}, $curOrfId);
				push(@{$h_nbOrfByCog{$curFosId}{"No hits found"}{"COGs"}}, "");
			}
		}
		close COGOUT;
		#use Data::Dumper qw(Dumper);
		#print Dumper \%h_nbOrfByCog;

		print LOG "done\n";
	}
	
	if($diamond) {
		print LOG "## Run diamond against your own fasta file...";
		`diamond makedb --in $diamond --db $outputDir/private_tmp_ >& $outputDir/private.diamond_makedb_tmp_.stderr`;
		`diamond blastx --db $outputDir/private_tmp_ --query $prokkaTFA --threads $threads --outfmt 6 qseqid sseqid pident nident length mismatch gaps gapopen qstart qend sstart send evalue bitscore stitle qcovhsp  -k 1000 --query-cover $diamond_queryCover -e $diamond_evalue --out $outputDir/private.diamond.tsv >& $outputDir/private.diamond_tmp_.stderr`;
		# Valid output ?
		my $error = `grep Error $outputDir/private.diamond_tmp_.stderr`;
		if($error) {
			chomp $error;
			print LOG colored(['bold red'], $error, "\n");
			#print LOG "$error\n";
		}
		else {
			print LOG "done\n";	
		}
	}

	
	#####
	## HTML Report
	#####
	print LOG "## Create HTML report";
	
	# Create fa and fai files
	copy($inputSeq, "$outputDir/fosmids.fa") or die "Unabled to copy $inputSeq input $outputDir: $!";
	`samtools faidx $outputDir/fosmids.fa`;
	
	# bgzip tabix for prokka.gff
	`bgzip -fc $prokkaGFF > $prokkaGFF.gz; tabix -p gff $prokkaGFF.gz`;
	
	my $dianrRTrack;     #Diamond VS NR 
	my $diauniRTrack;    #Diamond VS Uniprot 
	my $dianrorfRTrack;  #Diamond ORF VS NR
	my $diauniorfRTrack; #Diamond ORF VS Uniprot
	if($funAndTaxo) {
		# Tracks : Diamond VS NR (alignment and annotation tracks)
		# Alignment track  => Build sort.bam and bai from Diamond fosmid against NR tsv
		# Annotation track => for the "resume" track build $dianrRTrack
		# sam header from fai
		my $samheader = "";
		open(FAI,"$outputDir/fosmids.fa.fai") || die "Error: Unabled to open $outputDir/fosmids.fa.fai\n";
		while(my $line=<FAI>) {
			my @a_line = split(/\t/, $line);
			$samheader .= "\@SQ\tSN:$a_line[0]\tLN:$a_line[1]\n";
		}
		close(FAI);
		# sam body from contigs-nr.diamond.tsv
		open(SAM,">$outputDir/contigs-nr.diamond_tmp_.sam") || die "Error: Unabled to create $outputDir/contigs-nr.diamond_tmp_.sam\n";
		print SAM $samheader;
		open(DIANR, "$outputDir/contigs-nr.diamond.tsv") || die "Error: Unabled to open $outputDir/contigs-nr.diamond.tsv\n";
		my $cmp = 0;
		while(my $line=<DIANR>) {
			chomp $line;
			my ($qseqid,$sseqid,$pident,$nident,$length,$mismatch,$gaps,$gapopen,$qstart,$qend,$sstart,$send,$evalue,$bitscore,$stitle,$qcovhsp) = split(/\t/, $line);
			my $strand = 0;
			if($qstart>$qend) {
				my $tmp = $qstart;
				$qstart = $qend;
				$qend   = $tmp;
				$strand = 16;
			}
			print SAM "$sseqid\t$strand\t$qseqid\t$qstart\t30\t".($qend-$qstart+1)."M\t*\t0\t0\t*\t*\t";
			print SAM "XI:Z:$pident\tXJ:Z:$nident\tXL:i:$length\tXM:i:$mismatch\tXH:i:$gaps\tXG:i:$gapopen\tXC:Z:$qcovhsp\tXE:Z:$evalue\tXB:Z:$bitscore\tXA:Z:$stitle\n";
			my $link = $sseqid;
			$link =~s/^gi\|\d+\|[A-Za-z]*\|//;
			$link =~s/\|$//;
			$dianrRTrack .= "," if($cmp);
			$dianrRTrack .= "{chr:\"$qseqid\",Name:\"<a href=\\\"https://www.ncbi.nlm.nih.gov/protein/$link\\\" target=\\\"_blank\\\">$sseqid</a>\",start:$qstart,end:$qend,evalue:$evalue,color:\"rgba(64,140,193,0.2)\",description:\"$stitle\"}";
			$cmp++;
		}
		close(DIANR);
		close(SAM);
		`samtools view -bS $outputDir/contigs-nr.diamond_tmp_.sam | samtools sort - -o $outputDir/contigs-nr.diamond.bam ; samtools index $outputDir/contigs-nr.diamond.bam`;
		
		# Tracks : Diamond VS Uniprot (alignment and annotation tracks)
		# Alignmnet track  => Build sort.bam and bai from Diamond fosmid against uniprot tsv
		# Annotation track => for the "resume" track build $diauniRTrack
		# sam body from contigs-uniprot.diamond.tsv
		open(SAM,">$outputDir/contigs-uniprot.diamond_tmp_.sam") || die "Error: Unabled to create $outputDir/contigs-uniprot.diamond_tmp_.sam\n";
		print SAM $samheader;
		open(DIAUNI, "$outputDir/contigs-uniprot.diamond.tsv") || die "Error: Unabled to open $outputDir/contigs-uniprot.diamond.tsv\n";
		$cmp = 0;
		while(my $line=<DIAUNI>) {
			chomp $line;
			my ($qseqid,$sseqid,$pident,$nident,$length,$mismatch,$gaps,$gapopen,$qstart,$qend,$sstart,$send,$evalue,$bitscore,$stitle,$qcovhsp) = split(/\t/, $line);
			my $strand = 0;
			if($qstart>$qend) {
				my $tmp = $qstart;
				$qstart = $qend;
				$qend   = $tmp;
				$strand = 16;
			}
			print SAM "$sseqid\t$strand\t$qseqid\t$qstart\t30\t".($qend-$qstart+1)."M\t*\t0\t0\t*\t*\t";
			print SAM "XI:Z:$pident\tXJ:Z:$nident\tXL:i:$length\tXM:i:$mismatch\tXH:i:$gaps\tXG:i:$gapopen\tXC:Z:$qcovhsp\tXE:Z:$evalue\tXB:Z:$bitscore\tXA:Z:$stitle\n";
			my @a_link = split(/\|/,$sseqid);
			$diauniRTrack .= "," if($cmp);
			$diauniRTrack .= "{chr:\"$qseqid\",Name:\"<a href=\\\"https://www.ncbi.nlm.nih.gov/protein/$a_link[2]\\\" target=\\\"_blank\\\">$sseqid</a>\",start:$qstart,end:$qend,evalue:$evalue,color:\"rgba(88,162,88,0.2)\",description:\"$stitle\"}";
			$cmp++;
		}
		close(DIAUNI);
		close(SAM);
		`samtools view -bS $outputDir/contigs-uniprot.diamond_tmp_.sam | samtools sort - -o $outputDir/contigs-uniprot.diamond.bam ; samtools index $outputDir/contigs-uniprot.diamond.bam`;
			
		# Tracks : Diamond ORF VS NR (alignment and annotation tracks)
		# Alignmnet track  => Build sort.bam and bai from Diamond ORF against NR tsv
		# Annotation track => for the "resume" track build $dianrorfRTrack
		# Coord must be shift: orf to fosmid
		open(SAM,">$outputDir/prokka-nr.diamond_tmp_.sam") || die "Error: Unabled to create $outputDir/prokka-nr.diamond_tmp_.sam\n";
		print SAM $samheader;
		open(DIANR, "$outputDir/prokka-nr.diamond.tsv") || die "Error: Unabled to open $outputDir/prokka-nr.diamond.tsv\n";
		$cmp = 0;
		while(my $line=<DIANR>) {
			chomp $line;
			my ($qseqid,$sseqid,$pident,$nident,$length,$mismatch,$gaps,$gapopen,$qstart,$qend,$sstart,$send,$evalue,$bitscore,$stitle,$qcovhsp) = split(/\t/, $line);
			my $strand = 0;
			if($qstart>$qend) {
				my $tmp = $qstart;
				$qstart = $qend;
				$qend   = $tmp;
				$strand = 16;
			}
			my $orfname = "";
			if($qseqid =~ /^(\S*)#(ORF\d+)$/) {
				$qseqid  = $1;
				$orfname = $2;
			}
			my $len = $qend - $qstart + 1;
			my ($shiftstart, $shiftend) = (0, 0);
			if($h_orf{$qseqid}{$orfname}{"Strand"} == 0) {
				$shiftstart = $qstart + $h_orf{$qseqid}{$orfname}{"Start"};
				$shiftend   = $shiftstart + $len;
			}
			else {
				$shiftstart = $h_orf{$qseqid}{$orfname}{"Length"} - $qend + $h_orf{$qseqid}{$orfname}{"Start"};
				$shiftend   = $shiftstart + $len;
			}
			print SAM "$sseqid\t".$h_orf{$qseqid}{$orfname}{"Strand"}."\t$qseqid\t$shiftstart\t30\t$len"."M\t*\t0\t0\t*\t*\t";
			print SAM "XI:Z:$pident\tXJ:Z:$nident\tXL:i:$length\tXM:i:$mismatch\tXH:i:$gaps\tXG:i:$gapopen\tXC:Z:$qcovhsp\tXE:Z:$evalue\tXB:Z:$bitscore\tXA:Z:$stitle\n";
			my $link = $sseqid;
			$link =~s/^gi\|\d+\|[A-Za-z]*\|//;
			$link =~s/\|$//;
			$dianrorfRTrack .= "," if($cmp);
			$dianrorfRTrack .= "{chr:\"$qseqid\",Name:\"<a href=\\\"https://www.ncbi.nlm.nih.gov/protein/$link\\\" target=\\\"_blank\\\">$sseqid</a>\",start:$shiftstart,end:$shiftend,evalue:$evalue,color:\"rgba(193,64,64,0.2)\",description:\"$stitle\"}";
			
			# Save orf hsp, cov... info
			$h_orf{$qseqid}{$orfname}{"nrNbHSP"}++;
			$h_orf{$qseqid}{"nrORF"}++ if($h_orf{$qseqid}{$orfname}{"nrNbHSP"} == 1);
			$h_orf{$qseqid}{$orfname}{"nrMaxCov"} = $qcovhsp if($qcovhsp > $h_orf{$qseqid}{$orfname}{"nrMaxCov"});
			if($qcovhsp*$pident > $h_orf{$qseqid}{$orfname}{"nrCov"}*$h_orf{$qseqid}{$orfname}{"nrIden"}) {
				$h_orf{$qseqid}{$orfname}{"nrCov"}  = $qcovhsp;
				$h_orf{$qseqid}{$orfname}{"nrIden"} = $pident;
				if($stitle =~ /\[(.*?)\]/) {
					$h_orf{$qseqid}{$orfname}{"nrSpecies"} = $1;
				}
			}
			$cmp++;
		}
		close(DIANR);
		close(SAM);	
		`samtools view -bS $outputDir/prokka-nr.diamond_tmp_.sam | samtools sort - -o $outputDir/prokka-nr.diamond.bam ; samtools index $outputDir/prokka-nr.diamond.bam`;
		
		# Tracks : Diamond ORF VS Uniprot (alignment and annotation tracks)
		# Alignmnet track  => Build sort.bam and bai from Diamond ORF against uniprot tsv
		# Annotation track => for the "resume" track build $diauniorfRTrack
		# Coord must be shift: orf to fosmid
		open(SAM,">$outputDir/prokka-uniprot.diamond_tmp_.sam") || die "Error: Unabled to create $outputDir/prokka-uniprot.diamond_tmp_.sam\n";
		print SAM $samheader;
		open(DIAUNI, "$outputDir/prokka-uniprot.diamond.tsv") || die "Error: Unabled to open $outputDir/prokka-uniprot.diamond.tsv\n";
		$cmp = 0;
		while(my $line=<DIAUNI>) {
			chomp $line;
			my ($qseqid,$sseqid,$pident,$nident,$length,$mismatch,$gaps,$gapopen,$qstart,$qend,$sstart,$send,$evalue,$bitscore,$stitle,$qcovhsp) = split(/\t/, $line);
			my $strand = 0;
			if($qstart>$qend) {
				my $tmp = $qstart;
				$qstart = $qend;
				$qend   = $tmp;
				$strand = 16;
			}
			my $orfname = "";
			if($qseqid =~ /^(\S*)#(ORF\d+)$/) {
				$qseqid  = $1;
				$orfname = $2;
			}
			my $len = $qend - $qstart + 1;
			my ($shiftstart, $shiftend) = (0, 0);
			if($h_orf{$qseqid}{$orfname}{"Strand"} == 0) {
				$shiftstart = $qstart + $h_orf{$qseqid}{$orfname}{"Start"};
				$shiftend   = $shiftstart + $len;
			}
			else {
				$shiftstart = $h_orf{$qseqid}{$orfname}{"Length"} - $qend + $h_orf{$qseqid}{$orfname}{"Start"};
				$shiftend   = $shiftstart + $len;
			}
			print SAM "$sseqid\t".$h_orf{$qseqid}{$orfname}{"Strand"}."\t$qseqid\t$shiftstart\t30\t$len"."M\t*\t0\t0\t*\t*\t";
			print SAM "XI:Z:$pident\tXJ:Z:$nident\tXL:i:$length\tXM:i:$mismatch\tXH:i:$gaps\tXG:i:$gapopen\tXC:Z:$qcovhsp\tXE:Z:$evalue\tXB:Z:$bitscore\tXA:Z:$stitle\n";
			my @a_link = split(/\|/,$sseqid);
			$diauniorfRTrack .= "," if($cmp);
			$diauniorfRTrack .= "{chr:\"$qseqid\",Name:\"<a href=\\\"https://www.ncbi.nlm.nih.gov/protein/$a_link[2]\\\" target=\\\"_blank\\\">$sseqid</a>\",start:$shiftstart,end:$shiftend,evalue:$evalue,color:\"rgba(224,136,42,0.2)\",description:\"$stitle\"}";
			
			# Save orf hsp, cov... info
			$h_orf{$qseqid}{$orfname}{"uniNbHSP"}++;
			$h_orf{$qseqid}{"uniORF"}++ if($h_orf{$qseqid}{$orfname}{"uniNbHSP"} == 1);
			$h_orf{$qseqid}{$orfname}{"uniMaxCov"} = $qcovhsp if($qcovhsp > $h_orf{$qseqid}{$orfname}{"uniMaxCov"});
			if($qcovhsp*$pident > $h_orf{$qseqid}{$orfname}{"uniCov"}*$h_orf{$qseqid}{$orfname}{"uniIden"}) {
				$h_orf{$qseqid}{$orfname}{"uniCov"}  = $qcovhsp;
				$h_orf{$qseqid}{$orfname}{"uniIden"} = $pident;
				if($stitle =~ /OS=(\w+\s+\w+)/) {
					$h_orf{$qseqid}{$orfname}{"uniSpecies"} = $1;
				}
			}
			$cmp++;
		}
		close(DIAUNI);
		close(SAM);	
		`samtools view -bS $outputDir/prokka-uniprot.diamond_tmp_.sam | samtools sort - -o $outputDir/prokka-uniprot.diamond.bam ; samtools index $outputDir/prokka-uniprot.diamond.bam`;
	}
	
	open (HTML, ">$outputDir/$outputHtml") || die "Error: Unabled to create $outputDir/$outputHtml\n";
	my $html = $HTML_HEADER;
	$html =~ s/###TITLE###/Mintia annotation report/g;
		
	my $menu = '
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
				<li class="nav-item">
					<a class="nav-link" href="#annotation">
						<span class="oi oi-eye" aria-hidden="true"></span>
						Annotation
					</a>
				</li>
					<li class="nav-item" style="padding-left:10px">
						<a class="nav-link" href="#annotation-table">
							<span class="oi oi-spreadsheet" aria-hidden="true"></span>
							Result table
						</a>
					</li>
					<li class="nav-item" style="padding-left:10px">
						<a class="nav-link" href="#annotation-tablefos">
							<span class="oi oi-spreadsheet" aria-hidden="true"></span>
							Result table by fosmid
						</a>
					</li>
					<li class="nav-item" style="padding-left:10px">
						<a class="nav-link" href="#annotation-browser">
							<span class="oi oi-monitor" aria-hidden="true"></span>
							Fosmid browser
						</a>
					</li>
            </ul>
          </div>';
	
	$html =~ s/###MENU###/$menu/;
	print HTML $html;
	
	print HTML '
		<div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center pb-2 border-bottom">
			<h1 class="h4">Inputs and parameters</h1>
			<span class="anchor" id="inputs-parameters"></span>
		</div>';
	print HTML '
		<div class="d-flex">
			<div class="mt-4 mr-4 pl-0 col-md-4">
				<h5>Parameters</h5>
				<span class="anchor" id="parameters"></span>
				<ul class="list-group">';
	my $tmp = "No";
	my $badge = "badge-danger";
	if($funAndTaxo) { $tmp = "Yes"; $badge = "badge-success"; }
	print HTML "<li class=\"list-group-item d-flex justify-content-between align-items-center\">Run functional and taxonomic annotations <span class=\"badge $badge badge-pill ml-4\">$tmp</span></li>";
	$tmp = "No";
	$badge = "badge-danger";
	if($megan) { $tmp = "Yes"; $badge = "badge-success"; }
	print HTML "<li class=\"list-group-item d-flex justify-content-between align-items-center\">Run MEGAN <span class=\"badge $badge badge-pill ml-4\">$tmp</span></li>";
	$tmp = "No";
	$badge = "badge-danger";
	if($cog) { $tmp = "Yes"; $badge = "badge-success"; }
	print HTML "<li class=\"list-group-item d-flex justify-content-between align-items-center\">Run annotations with COGs <span class=\"badge $badge badge-pill ml-4\">$tmp</span></li>";
	$tmp = "No";
	$badge = "badge-danger";
	if($ncbi) { $tmp = "Yes"; $badge = "badge-success"; }
	print HTML "<li class=\"list-group-item d-flex justify-content-between align-items-center\">Build submission files <span class=\"badge $badge badge-pill ml-4\">$tmp</span></li>";
	$tmp = "No";
	$badge = "badge-danger";
	if($diamond) { $tmp = "Yes"; $badge = "badge-success"; }
	print HTML "<li class=\"list-group-item d-flex justify-content-between align-items-center\">Run blast against your own fasta file <span class=\"badge $badge badge-pill ml-4\">$tmp</span></li>";
	print HTML	'
				</ul>
			</div>
		</div>

		<div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center mt-5 pb-2 border-bottom">
			<h1 class="h4">Annotation</h1>
			<span class="anchor" id="annotation"></span>
		</div>
		
		<h5 class="mt-4">Result table</h5>
		<span class="anchor" id="annotation-table"></span>
		<table class="table table-striped table-bordered mb-0" style="width:100%;">
			<thead>
				<tr>
					<th style="text-align:center">Prokka Fasta files</th>
					<th style="text-align:center">Prokka annotation</th>
	';
	if($funAndTaxo) {
		print HTML '
					<th style="text-align:center">Diamond fosmid against NR</th>
					<th style="text-align:center">Diamond fosmid against Uniprot</th>
					<th style="text-align:center">Diamond ORFs against NR</th>
					<th style="text-align:center">Diamond ORFs against Uniprot</th>
		';
	}
	if($cog) {
		print HTML '<th style="text-align:center">COGs annotation</th>';
	}
	if($diamond) {
		print HTML '<th style="text-align:center">Diamond ORFs against private DB</th>';
	}
	print HTML '
				</tr>
			</thead>
			<tbody>
				<tr>
					<td nowrap style="text-align:center">
						<a href="./prokka.fasta"      target="_blank" class="btn btn-sm btn-outline-secondary">DNA</a> -
						<a href="./prokka-prot.fasta" target="_blank" class="btn btn-sm btn-outline-secondary">PROT</a></td>
					<td nowrap style="text-align:center"><a href="./prokka.gff" target="_blank" class="btn btn-sm btn-outline-secondary">GFF</a></td>
	';
	if($funAndTaxo) {
		print HTML '
					<td nowrap style="text-align:center"><a href="./contigs-nr.diamond.tsv" target="_blank" class="btn btn-sm btn-outline-secondary">TSV</a></td>
					<td nowrap style="text-align:center"><a href="./contigs-uniprot.diamond.tsv" target="_blank" class="btn btn-sm btn-outline-secondary">TSV</a></td>
					<td nowrap style="text-align:center"><a href="./prokka-nr.diamond.tsv" target="_blank" class="btn btn-sm btn-outline-secondary">TSV</a></td>
					<td nowrap style="text-align:center"><a href="./prokka-uniprot.diamond.tsv" target="_blank" class="btn btn-sm btn-outline-secondary">TSV</a></td>
		';
	}
	if($cog) {
		print HTML '<td nowrap style="text-align:center"><a href="./cogs_1e-07.rpsblast" target="_blank" class="btn btn-sm btn-outline-secondary">RPSBLAST</a></td>';
	}
	if($diamond) {
		print HTML '<td nowrap style="text-align:center"><a href="./private.diamond.tsv" target="_blank" class="btn btn-sm btn-outline-secondary">TSV</a></td>';
	}
	print HTML ' 
				</tr>
			</tbody>
		</table>
		
		<h5 class="mt-4">Result table by fosmid</h5>
		<span class="anchor" id="annotation-tablefos"></span>
		<table class="table table-striped table-bordered mb-0" style="width:100%;">
			<thead>
			<tr>
				<th nowrap class="valn" style="text-align:center">Sample</th>
				<th nowrap class="valn" style="text-align:center">Length</th>
				<th nowrap class="valn" style="text-align:center">GC%</th>
				<th nowrap class="valn" style="text-align:center">A</th>
				<th nowrap class="valn" style="text-align:center">T</th>
				<th nowrap class="valn" style="text-align:center">G</th>
				<th nowrap class="valn" style="text-align:center">C</th>
				<th nowrap class="valn" style="text-align:center">Others</th>
				<th nowrap class="valn" style="text-align:center">Prokka</th>
	';
	if($ncbi) {
		print HTML '
				<th nowrap class="valn" style="text-align:center">GenBank</th>
		';
	}
	if($funAndTaxo) {
		print HTML '
				<th nowrap class="valn" style="text-align:center" colspan="2">Annotated ORFs<br/>against Uniprot</th>
				<th nowrap class="valn" style="text-align:center" colspan="2">Annotated ORFs<br/>against NR</th>
		';
	}
	print HTML '<th nowrap class="valn" style="text-align:center">MEGAN</th>' if($megan);
	print HTML '<th nowrap class="valn" style="text-align:center" colspan="2">COGs</th>' if($cog);
	print HTML '
			</tr>
			</thead>
			<tbody>
	';
	
	foreach my $i (@a_fosmid) {
		my $nbseq = `grep -c ">" $outputDir/$i/fosmid.fasta`;
		chomp $nbseq;
		print HTML "<tr class='text-right'><td class='valn text-left' rowspan=$nbseq>$i</td>\n";
		open(FOSMID, "$outputDir/$i/fosmid.fasta") || die "Error: Unabled to open $outputDir/$i/fosmid.fasta\n";
		my ($id,$len,$gc,$A,$T,$C,$G, $o) = ("", 0, 0, 0, 0, 0, 0, 0);
		while(my $line=<FOSMID>) {
			chomp $line;
			if($line =~/^>(\S+)/) { 
				if($id ne "") {
					$gc = sprintf("%.2f", ($G+$C)/$len*100);
					$o  = $len-($A+$T+$C+$G);					
					print HTML "<td class='valn'>$len</td> <td class='valn'>$gc</td>"
							.	"<td class='valn'>$A</td>  <td class='valn'>$T</td>"
							.	"<td class='valn'>$G</td>  <td class='valn'>$C</td>"
							.	"<td class='valn'>$o</td>"
							.   "<td class='valn text-center' rowspan=$nbseq>"
							.		"<a href='./$i/prokka.fasta' class='btn btn-sm btn-outline-secondary' target='_blank'>Fasta</a></td>";
					if($ncbi) {
						print HTML "<td class='valn text-center' rowspan=$nbseq>"
							.		"<a href='./$i/prokka.gbk' class='btn btn-sm btn-outline-secondary' target='_blank'>GBK</a></td>";
					}
					if($funAndTaxo) {
						print HTML "<td class='valn'>".$h_orf{$id}{'uniORF'}." / ".$h_orf{$id}{'ORF'}."</td><td style='width:59px'>";
						if($h_orf{$id}{"uniORF"} > 0) {
							print HTML "<button type='button' data-toggle='modal' data-target='#uniAnnot-$id' class='btn btn-sm btn-outline-secondary'>"
									  ."<span class='oi oi-zoom-in' aria-hidden='true'></span></button>";
						}
						print HTML "</td>"
								.  "<td class='valn'>".$h_orf{$id}{'nrORF'}." / ".$h_orf{$id}{'ORF'}."</td><td style='width:59px'>";
						if($h_orf{$id}{"nrORF"} > 0) {
							print HTML "<button type='button' data-toggle='modal' data-target='#nrAnnot-$id' class='btn btn-sm btn-outline-secondary'>"
								.		"<span class='oi oi-zoom-in' aria-hidden='true'></span></button>";
						}
						print HTML "</td>";
					}
					if($megan) {
						print HTML "<td class='valn text-center' rowspan=$nbseq>"
								.		"<button type='button' data-toggle='modal' data-target='#meganview' class='pop btn btn-sm btn-outline-secondary' style='padding:0 4px 0 4px'>"
								.		"<img style='width:auto;height:24px' data-title='".$id."' src='./$i/megan.png'></button>"
								.	"</td>";
					}
					if($cog) {
						print HTML "<td class='valn textoverflow' rowspan=$nbseq>";
						my $tmpCat = "";
						my $tmpMax = 0;
						foreach my $k (keys %{$h_nbOrfByCog{$id}}) {
							if($k ne "No hits found" && scalar(@{$h_nbOrfByCog{$id}{$k}{"ORFs"}}) > $tmpMax) {
								$tmpCat = $k;
								$tmpMax = scalar(@{$h_nbOrfByCog{$id}{$k}{"ORFs"}});
							}
						}
						print HTML $tmpCat;
						print HTML "</td><td style='width:59px'>"
								."<button type='button' data-toggle='modal' data-target='#cogAnnot-$id' class='btn btn-sm btn-outline-secondary'>"
								."<span class='oi oi-pie-chart' aria-hidden='true'></span></button>";
						print HTML "</td>";
					}
					print HTML "</tr>";
					($len,$gc,$A,$T,$C,$G,$o) = (0, 0, 0, 0, 0, 0, 0);
					print HTML "</tr><tr class='text-right'>";
				}
				$id = $1;
			}
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
		print HTML "<td class='valn'>$len</td> <td class='valn'>$gc</td>"
				.	"<td class='valn'>$A</td>  <td class='valn'>$T</td>"
				.	"<td class='valn'>$G</td>  <td class='valn'>$C</td>"
				.	"<td class='valn'>$o</td>";
		if($nbseq==1) {
			print HTML "<td class='valn text-center' rowspan=$nbseq>"
					.		"<a href='./$i/prokka.fasta' class='btn btn-sm btn-outline-secondary' target='_blank'>Fasta</a></td>";
		}
		if($ncbi && $nbseq==1) {
			print HTML "<td class='valn text-center' rowspan=$nbseq>"
					.       "<a href='./$i/prokka.gbk' class='btn btn-sm btn-outline-secondary' target='_blank'>GBK</a></td>";
		}
		if($funAndTaxo) {
			print HTML "<td class='valn'>".$h_orf{$id}{'uniORF'}." / ".$h_orf{$id}{'ORF'}."</td><td style='width:59px'>";
			if($h_orf{$id}{"uniORF"} > 0) {
				print HTML "<button type='button' data-toggle='modal' data-target='#uniAnnot-$id' class='btn btn-sm btn-outline-secondary'>"
						  ."<span class='oi oi-zoom-in' aria-hidden='true'></span></button>";
			}
			print HTML "</td>"
					.  "<td class='valn'>".$h_orf{$id}{'nrORF'}." / ".$h_orf{$id}{'ORF'}."</td><td style='width:59px'>";
			if($h_orf{$id}{"nrORF"} > 0) {
				print HTML "<button type='button' data-toggle='modal' data-target='#nrAnnot-$id' class='btn btn-sm btn-outline-secondary'>"
						  ."<span class='oi oi-zoom-in' aria-hidden='true'></span></button>";
			}
			print HTML "</td>";
		}
		if($megan && $nbseq==1) {
			print HTML "<td class='valn text-center' rowspan=$nbseq>"
					.		"<button type='button' data-toggle='modal' data-target='#meganview' class='pop btn btn-sm btn-outline-secondary' style='padding:0 4px 0 4px'>"
					.		"<img style='width:auto;height:24px' data-title='".$id."' src='./$i/megan.png'></button>"
					.	"</td>";
		}
		if($cog) {
			print HTML "<td class='valn textoverflow' rowspan=$nbseq>";
			my $tmpCat = "";
			my $tmpMax = 0;
			foreach my $k (keys %{$h_nbOrfByCog{$id}}) {
				if($k ne "No hits found" && scalar(@{$h_nbOrfByCog{$id}{$k}{"ORFs"}}) > $tmpMax) {
					$tmpCat = $k;
					$tmpMax = scalar(@{$h_nbOrfByCog{$id}{$k}{"ORFs"}});
				}
			}
			print HTML $tmpCat;
			print HTML "</td><td style='width:59px'>"
					."<button type='button' data-toggle='modal' data-target='#cogAnnot-$id' class='btn btn-sm btn-outline-secondary'>"
					."<span class='oi oi-pie-chart' aria-hidden='true'></span></button>";
			print HTML "</td>";
		}
		print HTML "</tr>";
		close(FOSMID);
	}
	print HTML '
			</tbody>
		</table>
	';
	
	# Modal megan
	print HTML '
		<div class="modal fade" id="meganview" tabindex="-1" role="dialog" aria-labelledby="myModalLabel" aria-hidden="true">
			<div class="modal-dialog" data-dismiss="modal" style="max-width:1232px;">
				<div class="modal-content">
					<div class="modal-header">
						<h5 class="modal-title">
							<span id="megantitle"></span>
							<br/>
							<small class="text-muted">Produced by MEGAN5</small>
						</h5>
						<button type="button" class="close" data-dismiss="modal" aria-label="Close">
						  <span aria-hidden="true">&times;</span>
						</button>
					</div>
					<div class="modal-body">
						<img src="" id="meganimg" style="width: 100%;" >
					</div>
				</div>          
			</div>
		</div>
	';	
	
	# Modal Uniprot annotated ORFs
	foreach my $id (keys (%h_orf)) {
		print HTML '
		<div class="modal fade" id="nrAnnot-'.$id.'" tabindex="-1" role="dialog" aria-labelledby="myModalLabel" aria-hidden="true">
			<div class="modal-dialog modal-lg" role="document">
				<div class="modal-content">
					<div class="modal-header">
						<h5 class="modal-title">'.$id.'
							<br/>
							<small class="text-muted">Based on Diamond results of each ORF against NR, details for each annotated ORFs ('.$h_orf{$id}{"nrORF"}.'/'.$h_orf{$id}{"ORF"}.')</small>
						</h5>
						<button type="button" class="close" data-dismiss="modal" aria-label="Close">
						  <span aria-hidden="true">&times;</span>
						</button>
					</div>
					<div class="modal-body">
				
						<table class="table table-striped table-bordered mb-0" style="width:100%;">
							<thead>
								<tr>
									<th rowspan=2 class="valn">ORF</th>
									<th rowspan=2 class="valn" nowrap style="text-align:center">HSP</th>
									<th rowspan=2 class="valn" nowrap style="text-align:center;width:175px">Max coverage</th>
									<th colspan=3 class="valn" nowrap style="text-align:center;border-bottom-width:1px">Best HSP (max. of %iden * coverage)</th>
								</tr>
								<tr>
									<th nowrap class="valn" style="text-align:center">%iden.</th>
									<th nowrap class="valn" style="text-align:center;width:175px">Coverage</th>
									<th class="valn textoverflow" style="text-align:center">Species</th>
								</tr>
							</thead>
							<tbody>';
			
		for(my $c=1;$c<$h_orf{$id}{"ORF"};$c++) {
			next if($h_orf{$id}{"ORF".$c}{"nrMaxCov"} == 0);
			print HTML '
								<tr>
									<td >ORF'.$c.'</td>
									<td style="text-align:right">'.$h_orf{$id}{"ORF".$c}{"nrNbHSP"}.'</td>
									<td> 
										<div class="progress" style="margin-top:3px">
											<div class="progress-bar" role="progressbar" style="width: '.
												$h_orf{$id}{"ORF".$c}{"nrMaxCov"}.'%;background-color:rgba(193,64,64,0.8) !important;" aria-valuenow="'.
												$h_orf{$id}{"ORF".$c}{"nrMaxCov"}.'" aria-valuemin="0" aria-valuemax="100">'.
												$h_orf{$id}{"ORF".$c}{"nrMaxCov"}.'%</div>
											</div>
										</div>
									</td>
									<td style="text-align:right">'.$h_orf{$id}{"ORF".$c}{"nrIden"}.'</td>
									<td>
										<div class="progress" style="margin-top:3px">
											<div class="progress-bar" role="progressbar" style="width: '.
												$h_orf{$id}{"ORF".$c}{"nrCov"}.'%;background-color:rgba(193,64,64,0.8) !important;" aria-valuenow="'.
												$h_orf{$id}{"ORF".$c}{"nrCov"}.'" aria-valuemin="0" aria-valuemax="100">'.
												$h_orf{$id}{"ORF".$c}{"nrCov"}.'%</div>
											</div>
										</div>
									</td>
									<td class="textoverflow">'.$h_orf{$id}{"ORF".$c}{"nrSpecies"}.'</td>
								</tr>
							';
		}
		print HTML '
							<tbody>
						</table>
					</div>
				</div>
			</div>
		</div>
		';
	}
	
	# Modal NR annotated ORFs
	foreach my $id (keys (%h_orf)) {
		print HTML '
		<div class="modal fade" id="uniAnnot-'.$id.'" tabindex="-1" role="dialog" aria-labelledby="myModalLabel" aria-hidden="true">
			<div class="modal-dialog modal-lg" role="document">
				<div class="modal-content">
					<div class="modal-header">
						<h5 class="modal-title">'.$id.'
							<br/>
							<small class="text-muted">Based on Diamond results of each ORF against NR, details for each annotated ORFs ('.$h_orf{$id}{"uniORF"}.'/'.$h_orf{$id}{"ORF"}.')</small>
						</h5>
						<button type="button" class="close" data-dismiss="modal" aria-label="Close">
						  <span aria-hidden="true">&times;</span>
						</button>
					</div>
					<div class="modal-body">
				
						<table class="table table-striped table-bordered mb-0" style="width:100%;">
							<thead>
								<tr>
									<th rowspan=2 class="valn">ORF</th>
									<th rowspan=2 class="valn" nowrap style="text-align:center">HSP</th>
									<th rowspan=2 class="valn" nowrap style="text-align:center;width:175px">Max coverage</th>
									<th colspan=3 class="valn" nowrap style="text-align:center;border-bottom-width:1px">Best HSP (max. of %iden * coverage)</th>
								</tr>
								<tr>
									<th nowrap class="valn" style="text-align:center">%iden.</th>
									<th nowrap class="valn" style="text-align:center;width:175px">Coverage</th>
									<th class="valn textoverflow" style="text-align:center">Species</th>
								</tr>
							</thead>
							<tbody>';
			
		for(my $c=1;$c<$h_orf{$id}{"ORF"};$c++) {
			next if($h_orf{$id}{"ORF".$c}{"uniMaxCov"} == 0);
			print HTML '
								<tr>
									<td >ORF'.$c.'</td>
									<td style="text-align:right">'.$h_orf{$id}{"ORF".$c}{"uniNbHSP"}.'</td>
									<td> 
										<div class="progress" style="margin-top:3px">
											<div class="progress-bar" role="progressbar" style="width: '.
												$h_orf{$id}{"ORF".$c}{"uniMaxCov"}.'%;background-color:rgba(224,136,42,0.8) !important;" aria-valuenow="'.
												$h_orf{$id}{"ORF".$c}{"uniMaxCov"}.'" aria-valuemin="0" aria-valuemax="100">'.
												$h_orf{$id}{"ORF".$c}{"uniMaxCov"}.'%</div>
											</div>
										</div>
									</td>
									<td style="text-align:right">'.$h_orf{$id}{"ORF".$c}{"uniIden"}.'</td>
									<td>
										<div class="progress" style="margin-top:3px">
											<div class="progress-bar" role="progressbar" style="width: '.
												$h_orf{$id}{"ORF".$c}{"uniCov"}.'%;background-color:rgba(224,136,42,0.8) !important;" aria-valuenow="'.
												$h_orf{$id}{"ORF".$c}{"uniCov"}.'" aria-valuemin="0" aria-valuemax="100">'.
												$h_orf{$id}{"ORF".$c}{"uniCov"}.'%</div>
											</div>
										</div>
									</td>
									<td class="textoverflow">'.$h_orf{$id}{"ORF".$c}{"uniSpecies"}.'</td>
								</tr>
							';
		}
		print HTML '
							<tbody>
						</table>
					</div>
				</div>
			</div>
		</div>
		';
	}
	
	# Modal COGs annotated
	foreach my $id (keys (%h_nbOrfByCog)) {
		print HTML '
		<div class="modal fade" id="cogAnnot-'.$id.'" tabindex="-1" role="dialog" aria-labelledby="myModalLabel" aria-hidden="true">
			<div class="modal-dialog modal-lg" role="document">
				<div class="modal-content">
					<div class="modal-header">
						<h5 class="modal-title">'.$id.'
							<br/>
							<small class="text-muted">Based on Rps-Blast results against COGs database, details for each COGs categories found.</small>
						</h5>
						<button type="button" class="close" data-dismiss="modal" aria-label="Close">
						  <span aria-hidden="true">&times;</span>
						</button>
					</div>
					<div class="modal-body">
						<div id="cog-pie-'.$id.'" style="min-width: 310px; height: 400px; max-width: 600px; margin: 0 auto"></div>
						<script type="text/javascript" language="javascript">						
							var pieColors = (function () {
								var colors = [],
									base = Highcharts.getOptions().colors[0],
									i;

								for (i = 0; i < 10; i += 0.75) {
									// Start out with a darkened base color (negative brighten), and end
									// up with a much brighter color
									colors.push(Highcharts.Color(base).brighten((i - 3) / 7).get());
								}
								return colors;
							}());
							// Build the chart
							Highcharts.chart("cog-pie-'.$id.'", {
								chart: {
									plotBackgroundColor: null,
									plotBorderWidth: null,
									plotShadow: false,
									type: "pie"
								},
								title: {
									text: "Number of ORF per COG categories - '.$id.'"
								},
								tooltip: {
									headerFormat: "",
									pointFormat: "<b>{point.y}</b> ORF(s)</b>"
									
								},
								subtitle: {
									text: "Produced by '.$MINTIA_VERSION.'"
								},
								credits: false,	
								exporting: {
									buttons: {
										contextButton: {
											menuItems: ["printChart",
														"downloadPNG",
														"downloadJPEG",
														"downloadPDF",
														"downloadSVG"]
										}
									}
								},
								plotOptions: {
									pie: {
										allowPointSelect: true,
										cursor: "pointer",
										colors: pieColors,
										dataLabels: {
											enabled: true,
											format: "<b>{point.name}</b>: {point.percentage:.1f} %",
											connectorColor: "silver"
										}
									}
								},
								series: [{
									name: "Number of ORF",
									data: [';
								my $nofirst = 0;
								foreach my $k ( sort ({ scalar(@{$h_nbOrfByCog{$id}{$b}{"COGs"}}) <=> scalar(@{$h_nbOrfByCog{$id}{$a}{"COGs"}}) } keys %{$h_nbOrfByCog{$id}})){
									#if($k ne "No hits found") {
										my $len = scalar(@{$h_nbOrfByCog{$id}{$k}{"COGs"}});
										if($nofirst) { print HTML ","; }
										print HTML '{name: "'.$k.'", y:'.$len;
										print HTML ', sliced: true, selected: true' if($nofirst == 1);
										print HTML '}';
										$nofirst++;
									#}
								}
								print HTML '
									]
								}]
							});
						</script>						
								
						<table class="table table-sm table-striped table-bordered mb-0" style="width:100%;font-size:small">
							<thead>
								<tr>
									<th class="valn" style="text-align:center">COGs category</th>
									<th class="valn" nowrap style="text-align:center">ORFs</th>
									<th class="valn" nowrap style="text-align:center">COGs (Best hit)</th>
								</tr>
							</thead>
							<tbody>';
		foreach my $k ( sort ({ scalar(@{$h_nbOrfByCog{$id}{$b}{"COGs"}}) <=> scalar(@{$h_nbOrfByCog{$id}{$a}{"COGs"}}) } keys %{$h_nbOrfByCog{$id}})){
			if($k ne "No hits found") {
				my $len = scalar(@{$h_nbOrfByCog{$id}{$k}{"COGs"}});
				print HTML '<tr><td rowspan='.$len.' class="textoverflow valn" style="width:30%">'.$k.'</td>';
				for(my $i=0;$i<$len;$i++) {
					print HTML '<td style="text-align:center">'.$h_nbOrfByCog{$id}{$k}{"ORFs"}[$i].'</td>';
					print HTML '<td class="textoverflow" style="width:60%">'.$h_nbOrfByCog{$id}{$k}{"COGs"}[$i].'</td>';
					print HTML '</tr>';
				}
			}
		}
		print HTML '
							<tbody>
						</table>
					</div>
				</div>
			</div>
		</div>
		';
	}
	
	print HTML '<h5 class="mt-4">Fosmid browser</h5>
		<span class="anchor" id="annotation-browser"></span>
		<div id="igv-div" class="mb-3 pt-2 pb-2" style="border:1px solid lightgray; width:100%"></div>
	';
	
	print HTML '
	<script type="text/javascript">
		document.addEventListener("DOMContentLoaded", function () {
			if($(location).attr("href").match("^file:")) {
				var options = {
					"minimumBases": 1500,
					reference:{
						"id": "Mintia",
						"name": "Mintia",
						"fastaURL": "./fosmids.fa",
						"indexURL": "./fosmids.fa.fai",
						"tracks": [
						{
							type: "annotation",
							format: "gff",
							url: "./prokka.gff.gz",
							indexURL: "./prokka.gff.gz.tbi",
							displayMode: "EXPANDED",
							name: "Prokka",
							visibilityWindow: 1000000,
							color:"rgba(104,64,193,0.8)"
						}';
	if($funAndTaxo) {
		print HTML ',
						{
							name: "Diamond VS NR",
							type: "annotation",
							visibilityWindow: 50000,
							height: 40,
							features: ['.$dianrRTrack.']
						},
						{
							name: "Diamond VS Uniprot",
							type: "annotation",
							visibilityWindow: 50000,
							height: 40,
							features: ['.$diauniRTrack.']
						},
						{
							name: "Diamond ORF VS NR",
							type: "annotation",
							visibilityWindow: 50000,
							height: 40,
							features: ['.$dianrorfRTrack.']
						},
						{
							name: "Diamond ORF VS Uniprot",
							type: "annotation",
							visibilityWindow: 50000,
							height: 40,
							features: ['.$diauniorfRTrack.']
						}';
	}
	print HTML '
				]}};
			}
			else {
				var options = {
					reference:{
						"id": "Mintia",
						"name": "Mintia",
						"fastaURL": "./fosmids.fa",
						"indexURL": "./fosmids.fa.fai",
						"tracks": [
						{
							type: "annotation",
							format: "gff",
							url: "./prokka.gff.gz",
							indexURL: "./prokka.gff.gz.tbi",
							displayMode: "EXPANDED",
							name: "Prokka",
							visibilityWindow: 1000000,
							color:"rgba(104,64,193,0.8)"
						}';
	if($funAndTaxo) {
		print HTML ',
						{
							name: "Diamond VS NR",
							type: "alignment",
							format: "bam",
							url: "./contigs-nr.diamond.bam",
							indexURL: "./contigs-nr.diamond.bam.bai",
							height: 125,
							visibilityWindow: 50000,
							color:"rgba(64,140,193,0.8)"
						},
						{
							name: "Diamond VS NR",
							type: "annotation",
							visibilityWindow: 50000,
							height: 40,
							features: ['.$dianrRTrack.']
						},
						{
							name: "Diamond VS Uniprot",
							type: "alignment",
							format: "bam",
							url: "./contigs-uniprot.diamond.bam",
							indexURL: "./contigs-uniprot.diamond.bam.bai",
							height: 125,
							visibilityWindow: 50000,
							color:"rgba(88,162,88,0.8)"
						},
						{
							name: "Diamond VS Uniprot",
							type: "annotation",
							visibilityWindow: 50000,
							height: 40,
							features: ['.$diauniRTrack.']
						},
						{
							name: "Diamond ORF VS NR",
							type: "alignment",
							format: "bam",
							url: "./prokka-nr.diamond.bam",
							indexURL: "./prokka-nr.diamond.bam.bai",
							height: 125,
							visibilityWindow: 50000,
							color:"rgba(193,64,64,0.8)"
						},
						{
							name: "Diamond ORF VS NR",
							type: "annotation",
							visibilityWindow: 50000,
							height: 40,
							features: ['.$dianrorfRTrack.']
						},
						{
							name: "Diamond ORF VS Uniprot",
							type: "alignment",
							format: "bam",
							url: "./prokka-uniprot.diamond.bam",
							indexURL: "./prokka-uniprot.diamond.bam.bai",
							height: 125,
							visibilityWindow: 50000,
							color:"rgba(224,136,42,0.8)"
						},
						{
							name: "Diamond ORF VS Uniprot",
							type: "annotation",
							visibilityWindow: 50000,
							height: 40,
							features: ['.$diauniorfRTrack.']
						}';
	}
	print HTML '
				]}};
			}
			
			const igvDiv = document.getElementById("igv-div");			
			igv.createBrowser(igvDiv, options).then(function (browser) {
                    browser.on(\'trackclick\', function (track, popoverData) {

                        let markup = "<table class=\"igv-popover-table\">";

                        // Don\'t show a pop-over when there\'s no data.
                        if (!popoverData || !popoverData.length) {
                            return false;
                        }
						
						firstheader = 0;
                        popoverData.forEach(function (nameValue) {
							if (nameValue.name) {
								let name  = nameValue.name;
								let value = nameValue.value;
								if(track.type == "annotation") {
									if(name == "id" || name == "chr") { return; }
								}
								else if(track.type == "alignment") {
									if( name == "Cigar"              || name == "Mapped"           ||
										name == "Secondary"          || name == "Supplementary"    ||
										name == "Duplicate"          || name == "Failed QC"        ||
										name == "Mapping Quality"    || name.startsWith("Genomic Location") ||
										name.startsWith("Read Base") || name.startsWith("Base Quality") ||
										name == "A" || name == "T" || name == "C" || name == "G" || name == "N")
										{ return; }
									if( name == "Alignment Start" ) {
										name = "Start-End";
										value =  value.replace(",","");
										let end = parseInt(popoverData[4].value.slice(0, -1)) + parseInt(value) - 1;
										value = value+"-"+end;
									}
									if( name == "Read Strand" )     { name = "Strand" }
									if( name == "XI" )              { name = "% identity" }
									if( name == "XJ" )              { name = "Nb. identity" }
									if( name == "XL" )              { name = "Alignment length" }
									if( name == "XM" )              { name = "Nb. mismatch" }
									if( name == "XH" )              { name = "Nb. gaps" }
									if( name == "XG" )              { name = "Nb. gap openings" }
									if( name == "XC" )              { name = "Coverage Per HSP" }
									if( name == "XE" )              { name = "Expect value" }
									if( name == "XB" )              { name = "Bit score" }
									if( name == "XA" )              { name = "Description" }
									if( name == "Read Name" ) {
										name = "Name";
										let a_value = value.split("|");
										let url = "<a href=\"https://www.ncbi.nlm.nih.gov/protein/";
										if(a_value[a_value.length-1] != "") {
											url += a_value[a_value.length-1];
										}
										else {
											url += a_value[a_value.length-2];
										}
										url += "\" target=\"_blank\">" + value + "</a>";
										value = url;
									}
								}
								
								let header = "";
								if((name == "ID" || name == "Name") && firstheader==0) {
										header = " style=\"background-color:#eee;border-bottom:1px solid darkgrey\"";
										if(name == "ID") { firstheader = 1; }
								}
								
								markup += "<tr" + header + "><td class=\"igv-popover-item\">" + name + "</td>"
										+ "<td class=\"igv-popover-value\">" + value + "</td>"
										+ "</td></tr>";
							}
                            else {
                                // not a name/value pair
                                markup += "<tr><td>" + nameValue.toString() + "</td></tr>";
                            }
                        });

                        markup += "</table>";

                        // By returning a string from the trackclick handler we\'re asking IGV to use our custom HTML in its pop-over.
                        return markup;
                    });
            });
		});
	</script>
	';
	
	print HTML $HTML_FOOTER;
	close(HTML);
	print LOG ".......done\n";
		
	if(!$keep) {
		print LOG "## Remove temporary files";
		unlink glob "$outputDir/*_tmp_*";
		unlink glob "$outputDir/*/*_tmp_*";
		print LOG "...done\n";
	}
	close(LOG);
}


############################################################################
# MAIN 
############################################################################
MAIN:
{
	my @ARGV_SAVE = @ARGV;
	my $argKO = 0;
	
	# No arg
	if($#ARGV == -1) { $argKO = 1; }
	# No command
	elsif($ARGV[0] ne "check" && $ARGV[0] ne "assemble" && $ARGV[0] ne "annotate") {
		$argKO = 1;
		# Global help
		if ( join("#", @ARGV) =~ /\-{1,2}h[elp]{0,3}/ ) {
			pod2usage(
			-verbose => 99,
			-sections => "NAME|COMMANDS|DESCRIPTION|MAIN OPTIONS|VERSION|AUTHORS|COPYRIGHT")
		}
		# Version
		elsif ( join("#", @ARGV) =~ /\-{1,2}v[ersion]{0,6}/ ) {
			pod2usage(
			-verbose => 99,
			-sections => "VERSION")
		}
		print "[main] unrecognized command: $ARGV[0]\n\n";
	}
	pod2usage(
		-verbose => 99,
		-sections => "NAME|COMMANDS|DESCRIPTION|MAIN OPTIONS|VERSION|AUTHORS|COPYRIGHT"
	) if($argKO);
	
	## Check dependencies
	if (defined($ARGV[0])  &&  $ARGV[0] eq "check") {
		shift @ARGV_SAVE;
		@ARGV = @ARGV_SAVE;
		check();
	}
	## Module1: Assemble
	elsif (defined($ARGV[0])  &&  $ARGV[0] eq "assemble") {
		shift @ARGV_SAVE;
		@ARGV = @ARGV_SAVE;
		assemble();
	}
	## Module2 : Annotate
	elsif (defined($ARGV[0])  &&  $ARGV[0] eq "annotate") {
		shift @ARGV_SAVE;
		@ARGV = @ARGV_SAVE;
		annotate();
	}
}

use strict;
use warnings;
use Getopt::Long;

use Net::FTP;
use Bio::SeqIO;

require File::Temp;
use File::Temp ();
use File::Temp qw/ :seekable /;

use File::Path qw(make_path remove_tree);
use File::Basename;
use File::Copy;

use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use Term::UI;

our $debug;

our @SERVERS = qw/ftp.ensembl.org
		  ftp.ensemblgenomes.org/;

our @SUFS = qw/.dna.toplevel.fa.gz
		   .dna.toplevel.fa
		   .fa
		   .fa.gz
		   .gtf
		   .gtf.gz/;

### This could be expended to include new aligners
our %INDEX =  (
	       bowtie  => 'bowtie-build',
	       bowtie2 => 'bowtie2-build'
	      );

### Little dataframe to keep track of important info
our %genomeInfo = map{$_,'--'} qw|fullName
				  shortName
				  group
				  base
				  id
				  releaseDate|;

our $TERM = Term::ReadLine->new('brand');

MAIN:{
  init();
  my $files = download_Genome();
  create_idx($_,$files) for keys %INDEX;
  create_tophat_idx($files);
  forgeBSgenome($files);
  ##create_IGV_genome($files->{GTF}); ## TODO fix this
  exit 0;
}
  
sub download_Genome{
  my $server = $TERM->get_reply(
				prompt => 'Pick a server',
				choices => \@SERVERS,
				default =>  $SERVERS[0],
			       );
  
  print "Downloading the list of available genomes from the current release at $server...\n";
  my $ftp = Net::FTP->new($server,
			  Passive => 1)
    or die "Can't connect: $@\n";
  $ftp->login
    or die "Couldn't login\n";
  $ftp->binary();

  my @releases = sort{($a=~/-(\d+)/)[0] <=> ($b=~/-(\d+)/)[0]} map{basename $_} grep{/release-/} $ftp->ls('pub');
  
  my %releases = map{/(.+)-(.+)/;"$1 $2" => $_} @releases[0..$#releases-1];
  $releases{current} = $releases[-1];
  my @rel_names = ((map{/(.+)-(.+)/;"$1 $2"} @releases[0..$#releases-1]),"current");
  
  my $release = $TERM->get_reply(
				 prompt => 'Pick a release',
				 choices => \@rel_names,
				 default => 'current',
				);
  if ($server eq $SERVERS[0]){
    my @genome_name = my @avail_sp = map{basename($_)} $ftp->ls("/pub/$releases{$release}/gtf");
    map{s/(^\w)/\U$1/g;s/_/ /g} @avail_sp;
    my %genomes_id = map{$avail_sp[$_],$_} 0..$#avail_sp;
    
    my $genome = $TERM->get_reply(
				  prompt => 'Pick a genome: ',
				  choices => \@avail_sp,
				 );
    $genomeInfo{release}   = $releases{$release};
    $genomeInfo{fullName}  = $genome;
    $genomeInfo{shortName} = $genome_name[$genomes_id{$genome}];
  
  } elsif ($server eq $SERVERS[1]){
    my %availGenome;
    my ($in, $species_file);
    open($in, '>', \$species_file);
    
    $ftp->get("/pub/$releases{$release}/species.txt", $in)
      or die "get failed ", $ftp->message;
    
    for (split "\n", $species_file){
    chomp;
      my @line = split /\t/;
      my ($gr) = ($line[2] =~/Ensembl(.+)/);
      next unless $gr;
      push @{$availGenome{$gr}},\@line;
    }
    
    my $kingdom = $TERM->get_reply(
				   prompt => 'Pick a kingdom: ',
				   choices => [keys %availGenome],
				  );
    
    my $genome = $TERM->get_reply(
				  prompt => 'Pick a genome: ',
				  choices => [map{$_->[0]} @{$availGenome{$kingdom}}]
				 );

    my ($idx) = grep{${$availGenome{$kingdom}->[$_]}[0] eq $genome} 0 .. $#{$availGenome{$kingdom}};
    
    my @genome = @{$availGenome{$kingdom}[$idx]};
    
    $genomeInfo{release}     = $releases{$release};
    $genomeInfo{fullName}    = $genome[0];
    $genomeInfo{shortName}   = $genome[1];
    $genomeInfo{group}       = $kingdom;
    $genomeInfo{id}          = $genome[3];
    $genomeInfo{releaseDate} = $genome[5];
  }
  
  my $files = download_files($ftp,$server);
  filterChromosomes($files);
  return $files;
}

sub download_files {
  my $ftp = shift;
  my $server = shift;
  
  my %files;
  
  print "Connecting to $server for ftp download\n";
  my (@fastas,@gtfs);
  if ($server eq $SERVERS[0]){
    @fastas = $ftp->ls("pub/$genomeInfo{release}/fasta/$genomeInfo{shortName}/dna");
    @gtfs = $ftp->ls("pub/$genomeInfo{release}/gtf/$genomeInfo{shortName}");
  } 
  elsif($server eq $SERVERS[1]){
    my $group = lc $genomeInfo{group};
    @fastas = $ftp->ls("pub/$genomeInfo{release}/$group/fasta/$genomeInfo{shortName}/dna");
    @gtfs =  $ftp->ls("pub/$genomeInfo{release}/$group/gtf/$genomeInfo{shortName}");
  }
  ($files{fasta}->{remote}) = grep{/\.dna\.toplevel/} @fastas;
  ($files{gtf}->{remote}) = grep{/\.gtf/} @gtfs;
  die "Could not find the genome files to download\n" unless $files{fasta} && $files{gtf};
  ($genomeInfo{base}) = (basename($files{fasta}->{remote}) =~ /(.+)\.dna/);
  
  my ($repos) = ($server =~ /\.(.+)\.org/);
  $files{curr_dir} = "current_genome/$repos/$genomeInfo{release}/$genomeInfo{shortName}";
  
  for my $file_type (qw|gtf fasta|){
    $files{$file_type}->{URL}   = "ftp://${server}/$files{$file_type}->{remote}";
    $files{$file_type}->{dir}   = "$files{curr_dir}/$file_type";
    $files{$file_type}->{local} = $files{$file_type}->{dir} . "/" . basename($files{$file_type}->{remote});
    make_path($files{$file_type}->{dir});
    print "Downloading $files{$file_type}->{remote}\n";
    
    ## This FTP will restart a previously failed download attemp or 
    ## skip if the file exists and have same size has remote.
    my $ftpsize = $ftp->size($files{$file_type}->{remote});
    my $localsize= stat($files{$file_type}->{local}) ?(stat(_))[7]:(0);
    next if $ftpsize == $localsize;
    $ftp->get($files{$file_type}->{remote},$files{$file_type}->{local},$localsize) 
      or die "Couldn't dowload $_\n",$ftp->message;
  }
  print "Done downloading the files\n";
  return \%files;
}

sub filterChromosomes {
  my $files = shift;
  
  my $g = readFasta($files->{fasta}->{local});
  my @chrs = selectChrs(keys %{$g});
  $files->{fasta}->{minimal} = "$files->{fasta}->{dir}/$genomeInfo{base}.min.fa";
  open OUT,">$files->{fasta}->{minimal}";
  # Print a minimal fasta with only the selected chromosomes
  print OUT (join("\n",">$_",@{$g->{$_}}),"\n") for @chrs;
  close OUT;
}

sub readFasta {
  my $genome = shift;
  my $old_IRS = $/;
  $/='>';
  print "Deciding what are the chromosomes versus unassembled contigs...\n";
  if ($genome =~ /\.gz$/){
    open IN, "gunzip -c $genome |";
  } else {
    open IN, $genome;
  }
  my @g = <IN>;
  #remove the first empty def line
  shift @g; 
  # Remove the > at the end of line (currently set as EOL character); 
  # Remove extra info from def line
  map{chomp;s/(^.+?)(\s.+)\n/$1\n/} @g; 
  my %g = map{my @t = split /\s+/;$t[0] => [@t[1..$#t]]} @g;
  $/=$old_IRS;
  return \%g;
}

sub selectChrs {
  my @chrs = map {/(.+?)(\s|$)/;$1} @_;
  
  my %selected = map{$_=>1} grep{/(^chr|^)(\d+|Mito|mito|MT|Z|W|Y|[IVX]+)([LR]|(|Het|het|$))/} @chrs;
  my @not_selected = grep{!exists $selected{$_}} @chrs;

  print("These are the suggested chromosomes that would be used to build the alignment indexes:\n");
  print_in_n_cols(5,keys %selected);
  print("\n");
  
  if (@not_selected){
    print("And these will not be included:\n");
    my $i = scalar @not_selected > 100?100:$#not_selected;
    print_in_n_cols(5,@not_selected[0..$i]);
    print("... (".scalar @not_selected." more)\n") if scalar @not_selected > 100;
  }
  my $bool = $TERM->ask_yn( prompt   => 'Continue with that selection?',
			    default  => 'y',
			  );
  print("All entries in the fasta files will then be used to build the indexes\n") unless $bool;
  return($bool?keys %selected: @chrs);
}

sub print_in_n_cols {
    my ($count, @items) = @_;
    while (@items){
        for (1 .. $count){
            print "\t", shift @items;
            last unless @items;
        }
        print "\n";
    }
}

sub create_idx {  
  my $idx = shift;
  my $files = shift;
    
  my $idx_dir = "$files->{curr_dir}/${idx}_indexes";
  make_path($idx_dir);
  
  my $fasta = $files->{fasta}->{minimal};
  print "indexing ",basename($fasta), " for $idx\n";
  
  my ($name,$path,$suffix) = fileparse($fasta,@SUFS);
  $files->{$idx} = "$idx_dir/$name";
  
  system("$INDEX{$idx} $fasta $files->{$idx} >$idx_dir/${idx}.log 2>&1") unless $debug;
}

sub create_tophat_idx{
  my $files = shift;
  my $gtf = $files->{gtf}->{local};
  
  print "Creating the transcript index for TopHat splice site alignments\n"; 
  my $tx_idx_d = "$files->{curr_dir}/tophat_tx_idx";
  make_path($tx_idx_d);
  
  my $fh_gtf = "$tx_idx_d/".join("",($genomeInfo{fullName} =~ /^(.).+\s(.+)/)).".gtf";
  
  if ($gtf =~ /\.gz$/) {
      gunzip $gtf => $fh_gtf
	or die "gunzip failed: $GunzipError\n";
      $gtf=$fh_gtf;
    }
  $files->{GTF} = $fh_gtf;
    
  ### create a little temp fastq file just to fire up tophat
  my $fh = File::Temp->new(SUFFIX => ".fastq");
  my $fname = $fh->filename;
  print $fh <DATA>;
  close $fh;
  my $name = basename($files->{fasta}->{minimal},'.fa');
  system("tophat -G $gtf --transcriptome-index=$tx_idx_d/$name.tc  $files->{bowtie2} $fname >$tx_idx_d/tx.log 2>&1");
  remove_tree 'tophat_out';
}

sub create_IGV_genome{
  my $files = shift;
  #TO FIX
  # my $fasta = $files->{fasta}->{minimal};    
  # $files->{igv}->{dir} = "$curr_dir/for_IGV";
  # make_path($files->{igv}->{dir});
  
  # copy $fasta,"$files->{igv}->{dir}/".basename( $fasta) 
  #   or die "Could not move $fasta to IGV directory\n";
  
  # copy $file,"$files->{igv}->{dir}/".basename($file)
  #   or die "Could not move $file to IGV directory\n"; 

  # my $igv_d = $files->{igv}->{dir};
  
  # system("tar -cvf ${igv_d}.tar.gz $igv_d >$igv_d/compression.log 2>&1");
}

sub forgeBSgenome{
  my $files = shift;
  my %fa = %{$files->{fasta}};
  
  my ($name,$path,$suffix) = fileparse($fa{local},@SUFS);
  my ($orgID,$assembly,$release) = ($name =~ /^(.*?)\.(.*)\.(.*)/);
  my ($genus,$species) = split /_/,$orgID;
  my $shortName = substr($genus, 0, 1) . $species;
  
  print "Forging BSgenome $shortName\n";

  my $package = join('.','BSgenome',$shortName,$assembly,$release);
  $package =~ s/_//g;
  $package =~ s/-/\./g;
  my $title = "$genus $species ($genomeInfo{group}) full genome (Ensembl assembly $assembly version $release)";
  
  my $bsg_base_dir="$files->{curr_dir}/BSgenomeForge";
  my $bsg_dir = "$bsg_base_dir/srcdata/$package";
  my $seqs_srcdir  = "$bsg_dir/seqs";
  my $masks_srcdir = "$bsg_dir/masks";
  make_path($seqs_srcdir,$masks_srcdir);
  
  my (@chr_name,@contig_name);
  my $in = Bio::SeqIO->new(-file => "gunzip -c $fa{local} |");
  my $out_contig = Bio::SeqIO->new(-file   => ">${seqs_srcdir}/contigs.fa",
				   -format => 'fasta');
  while (my $seq = $in->next_seq){
    if ($seq->display_id =~ /contig/i || $seq->desc =~ /contig/i){
      ## Remove the description if exists!
      $seq->desc('');
      $out_contig->write_seq($seq);
      push @contig_name, $seq->display_id;
    } else {
      ## Remove the description if exists!
      $seq->desc('');
      push @chr_name, $seq->display_id;
      my $out_name = $seqs_srcdir."/".$seq->display_id.".fa";
      my $out = Bio::SeqIO->new(-file   => ">$out_name",
				-format => 'fasta');
      $out->write_seq($seq);
    }
  }
  
  my $mseqnames = @contig_name?"c(\"contigs\")":"NULL";
  my $seqnames = @chr_name?"c(".join(",",map{"\"$_\""} @chr_name).")":"NULL";

  my $seed_file = "$ENV{PWD}/$bsg_base_dir/BSgenomeSeed.txt";
  open OUT, ">",$seed_file;
  
  print OUT <<SEED
Package: $package
Title: $title
Description: $title and stored in Biostrings objects.
Version: 1.0.0
organism: $genus $species
species: $genomeInfo{group}
provider: Ensembl
provider_version: $release
release_date: --
release_name: $assembly
source_url: $fa{URL}
organism_biocview: $genomeInfo{shortName}
BSgenomeObjname: $shortName
seqnames: $seqnames
mseqnames: $mseqnames
seqs_srcdir: $ENV{PWD}/$seqs_srcdir
masks_srcdir: $ENV{PWD}/$masks_srcdir
SEED
    ;
  close OUT;  

  ### Now to the forge
  open FORGE,"| R --vanilla --slave >$bsg_base_dir/forging.log 2>&1" or die "Can't open R to forge BSgenome";
  print FORGE <<TOFORGE
setwd(\"$bsg_base_dir\")
library(BSgenome)
forgeBSgenomeDataPkg(\"$seed_file\")
q()
TOFORGE
    ;
  close FORGE;
  
  my $oldPWD = $ENV{PWD};
  chdir "$ENV{PWD}/$bsg_base_dir" or die "Could not change directory to build the BSgenome\n";
  system "R CMD build $package >RBuild.log 2>&1";
  
  open OUT, ">README.txt" || die "Can't create README.txt file\n";
  print OUT ("To install the forged genome, run the following line in the terminal\n",
	     "R CMD install $ENV{PWD}/$package\n",
	     "\n",
	     "To use in R, use the following line:\n",
	     "library($package)\n",
	    ); 
  close OUT; 
  chdir $oldPWD;
}

sub init{
  GetOptions('debug'   => \$debug,
	    );
  print "#### RUNNING IN DEBUGING MODE ####\n" if $debug;
}


__DATA__
@ILLUMINA-C4D679_0026_FC:3:1:11884:1072#0/1
CTCGGGTATACATCAAGTGATGGATTATGCAAATTTTGGG
+ILLUMINA-C4D679_0026_FC:3:1:11884:1072#0/1
f[dfdd`\Jcfffff]cOdc^ad[IcdaJdfadfffff]a
@ILLUMINA-C4D679_0026_FC:3:1:16648:1081#0/1
GGACATCTACAAGCCAGTTGACAAAGTCGAACCCGGTACC
+ILLUMINA-C4D679_0026_FC:3:1:16648:1081#0/1
acfbfffc^fYc_ffdf]cfffcfc[bSKcaaadKS^ccc
@ILLUMINA-C4D679_0026_FC:3:1:9324:1106#0/1
CAGCAGTACCAATAACTCCAATAATATGAAGGATCAAGTC
+ILLUMINA-C4D679_0026_FC:3:1:9324:1106#0/1
f_fggfggYggcgagfdggdgdcggaf^[addfffcffRd
@ILLUMINA-C4D679_0026_FC:3:1:10004:1101#0/1
TCGAAAGATGAAAAGAACTTTGAAAAGAGAGTTAAATAGT
+ILLUMINA-C4D679_0026_FC:3:1:10004:1101#0/1
gfcagg_d^ggg_g^f`cda\dffcaa_S\afaead\YYW

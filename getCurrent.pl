use strict;
use warnings;
use Net::FTP;
use Bio::SeqIO;

require File::Temp;
use File::Temp ();
use File::Temp qw/ :seekable /;

use File::Path qw(make_path remove_tree);
use File::Basename;
use File::Copy;

use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use Archive::Tar;

use Term::UI;
use Term::ReadLine;
use Getopt::Long;

our $debug;

our @servers = qw/ftp.ensembl.org
		  ftp.ensemblgenomes.org/;

our @toDownload = qw/fasta gtf/;

our @suffixes = qw/.dna.toplevel.fa.gz
		   .dna.toplevel.fa
		   .fa
		   .fa.gz
		   .gtf
		   .gtf.gz/;


our %genomeInfo = map{$_,'--'} qw|fullName
				  shortName
				  group
				  id
				  releaseDate|;

our (%files,$curr_dir);


MAIN:{
  GetOptions('debug' => \$debug);
  print "#### RUNNING IN DEBUGING MODE ####\n" if $debug;

  my $term = Term::ReadLine->new('brand');
  my $server = $term->get_reply(
				prompt => 'Pick a server',
				choices => \@servers,
				default => $servers[0],
			       );

  openServer($server,$term);
  download_files($server);
  forgeBSgenome(%{$files{fasta}});
  create_bowtie_idx(@{$files{fasta}->{files}});
  create_bowtie2_idx(@{$files{fasta}->{files}});
  create_tophat_idx(@{$files{gtf}->{files}});
  create_IGV_genome($files{GTF});
  exit 1;
}

sub openServer{
  my $server = shift;
  my $term = shift;
  
  print "Downloading the list of available genomes from $server...\n";
  my $ftp = Net::FTP->new($server)
    or die "Can't connect: $@\n";
  $ftp->login
    or die "Couldn't login\n";
  $ftp->binary();
  
  if ($server eq $servers[0]){

    my @genome_name = my @avail_sp = $ftp->ls('/pub/current_gtf');
    map{s/(^\w)/\U$1/g;s/_/ /g} @avail_sp;
    my %genomes_id = map{$avail_sp[$_],$_} 0..$#avail_sp;
    
    my $genome = $term->get_reply(
				  prompt => 'Pick a genome: ',
				  choices => \@avail_sp,
				 );
    
    $genomeInfo{fullName} = $genome;
    $genomeInfo{shortName} = $genome_name[$genomes_id{$genome}];
    
    
    print("here \n");
  } elsif ($server eq $servers[1]){
    my %availGenome;
    open OUT, $ftp->get("/pub/current/species.txt");
    while (<OUT>){
      chomp;
	my @line = split /\t/;
      my ($gr) = ($line[2] =~/Ensembl(.+)/);
      next unless $gr;
      push @{$availGenome{$gr}},\@line;
    }
    
    my $kingdom = $term->get_reply(
				   prompt => 'Pick a kingdom: ',
				   choices => [keys %availGenome],
				  );
    
    my $genome = $term->get_reply(
				  prompt => 'Pick a genome: ',
				  choices => [map{$_->[0]} @{$availGenome{$kingdom}}]
				 );

    my ($idx) = grep{${$availGenome{$kingdom}->[$_]}[0] eq $genome} 0 .. $#{$availGenome{$kingdom}};
    
    my @genome = @{$availGenome{$kingdom}[$idx]};
    
    %genomeInfo = (fullName    => $genome[0],
		   shortName   => $genome[1],
		   group       => $kingdom,
		   id          => $genome[3],
		   releaseDate => $genome[5]);
  }
}

sub download_files {
  my $server = shift;
  print "Connecting to $server for ftp download\n";
  
  my $ftp = Net::FTP->new($server)
    or die "Can't connect: $@\n";
  $ftp->login
    or die "Couldn't login\n";
  $ftp->binary();
  
  for my $file_type (@toDownload) {
    
    my $group = lc $genomeInfo{group};

    my $remote_dir;
    if ($server eq $servers[0]){
      my ($subtype,$subdir) = $file_type eq 'fasta'?('current_fasta','dna'):('current_gtf','');
      $remote_dir="/pub/$subtype/$genomeInfo{shortName}/$subdir";
      
    }elsif($server eq $servers[1]){
      my $subtype = $file_type eq 'fasta'?'dna':'';
      $remote_dir="/pub/current/$group/$file_type/$genomeInfo{shortName}/$subtype";
    }
    
    $ftp->cwd($remote_dir)
      or die "Couldn't change directory\n";
    
    unless (exists $genomeInfo{release}){
      #($genomeInfo{release}) = ($ftp->pwd() =~ /pub\/(.+?)\//);
      ($genomeInfo{release}) = ($ftp->pwd() =~ /(release.+?)\//);
    }
    
    my ($repos) = ($server =~ /\.(.+)\.org/);
    $curr_dir = "current_genome/$repos/$genomeInfo{release}/$genomeInfo{shortName}";
    my $local_dir="$curr_dir/$file_type";
    make_path($local_dir);
    
    $files{$file_type}->{dir}=$local_dir;
    
    my (@files) = $ftp->ls();
    @files = grep{/\.dna\.toplevel/} @files if $file_type eq 'fasta';
    @files = grep{/\.gtf/} @files if $file_type eq 'gtf';
    die "Coulnd not find a file to download\n" unless @files;
  
    for ( @files) {
      my $local_file = $local_dir."/".$_;
      print "Downloading $_\n";
      $files{$file_type}->{URL} = "ftp://${server}${remote_dir}/$_";
      ($ftp->get($_,$local_file)
       or die "Couldn't dowload $_\n") unless $debug;
      push @{$files{$file_type}->{files}},$local_file
    }
  }
}



sub create_bowtie_idx {  
  my @files = @_;
  
  print "indexing $fasta for BowTie\n";
  my $bwt_idx_dir = "$curr_dir/bowtie_indexes";
  make_path($bwt_idx_dir);
  
  for my $fasta (@files) {
    my $file = $fasta;
    
    my ($name,$path,$suffix) = fileparse($file,@suffixes);
    $files{bwt_index} = "$bwt_idx_dir/$name";
  
    if ($fasta =~/\.gz/) {
      $file = "$bwt_idx_dir/$name.fa";
      system("gunzip -c $fasta > $file") == 0 || die "Can't unzip $fasta\n";
      $files{unzipped_fasta} = "$ENV{PWD}/$file";
    }
    system("bowtie-build $file $files{bwt_index} >$bwt_idx_dir/bwt_idx.log 2>&1") unless $debug;
  }
}

sub create_bowtie2_idx {  
  my @files = @_;
  
  print "indexing $fasta for BowTie2\n";
  my $bwt2_idx_dir = "$curr_dir/bowtie2_indexes";
  make_path($bwt2_idx_dir);
  
  for my $fasta (@files) {
    my $file = $fasta;
    
    my ($name,$path,$suffix) = fileparse($file,@suffixes);
    $files{bwt2_index} = "$bwt2_idx_dir/$name";
  
    if ($fasta =~/\.gz/) {
      $file = "$bwt2_idx_dir/$name.fa";
      system("gunzip -c $fasta > $file") == 0 || die "Can't unzip $fasta\n";
      $files{unzipped_fasta} = "$ENV{PWD}/$file";
    }
    system("bowtie2-build $file $files{bwt2_index} >$bwt2_idx_dir/bwt2_idx.log 2>&1") unless $debug;
  }
}

sub create_tophat_idx{
  my @files = @_;

  print "Creating the transcript index for TopHat splice site alignments\n"; 
  my $tx_idx_d = "$curr_dir/tophat_tx_idx";
  make_path($tx_idx_d);
    
  for my $gtf (@files) {
    next unless $gtf =~ /\.gtf/;
    my $file = $gtf;
    
    my $fh_gtf = "$tx_idx_d/".join("",($genomeInfo{fullName} =~ /^(.).+\s(.+)/)).".gtf";

    if ($gtf =~ /\.gz$/) {
      gunzip $gtf => $fh_gtf
	or die "gunzip failed: $GunzipError\n";
      #$gtf = $fh_gtf->filename;
      $gtf=$fh_gtf;
    }
    $files{GTF} = $fh_gtf;
    
    ### create a little temp fastq file just to fire up tophat
    my $fh = File::Temp->new(SUFFIX => ".fastq");
    my $fname = $fh->filename;
    print $fh <DATA>;
    close $fh;

    system("tophat -p12 -G $gtf --transcriptome-index=$tx_idx_d/known_tc  $files{bwt2_index} $fname >$tx_idx_d/tx.log 2>&1")
      unless $debug;
    remove_tree 'tophat_out';
  }
}

sub create_IGV_genome{
  my $file = shift;
      
  $files{igv}->{dir} = "$curr_dir/for_IGV";
  make_path($files{igv}->{dir});
  
  copy $files{unzipped_fasta},"$files{igv}->{dir}/".basename( $files{unzipped_fasta}) 
    or die "Could not move $files{unzipped_fasta} to IGV directory\n";
  
  copy $file,"$files{igv}->{dir}/".basename($file)
    or die "Could not move $file to IGV directory\n"; 

  my $igv_d = $files{igv}->{dir};
  
  system("tar -cvf ${igv_d}.tar.gz $igv_d >$igv_d/compression.log 2>&1");
}

sub forgeBSgenome{
  my %fa = @_;
  
  print "Forging BSgenome $shortName\n";

  my ($name,$path,$suffix) = fileparse($fa{files}[0],@suffixes);
  my ($orgID,$assembly,$release) = ($name =~ /^(.*?)\.(.*)\.(.*)/);
  my ($genus,$species) = split /_/,$orgID;
  my $shortName = substr($genus, 0, 1) . $species;

  my $package = join('.','BSgenome',$shortName,$assembly,$release);
  $package =~ s/_//g;
  $package =~ s/-/\./g;
  my $title = "$genus $species ($genomeInfo{group}) full genome (Ensembl assembly $assembly version $release)";
  
  my $bsg_base_dir="$curr_dir/BSgenomeForge";
  my $bsg_dir = "$bsg_base_dir/srcdata/$package";
  my $seqs_srcdir  = "$bsg_dir/seqs";
  my $masks_srcdir = "$bsg_dir/masks";
  make_path($seqs_srcdir,$masks_srcdir);
  
  my (@chr_name,@contig_name);
  my $in = Bio::SeqIO->new(-file => "gunzip -c $fa{files}[0] |");
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
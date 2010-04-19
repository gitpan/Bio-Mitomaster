#!perl -T

use Test::More tests => 19;

BEGIN {
	use_ok( 'Bio::Mitomaster::SeqIO' );
	use_ok( 'Bio::Mitomaster::FastaIO' );
	use_ok( 'Bio::Mitomaster::GenbankIO' );
	use_ok( 'Bio::Mitomaster::VariantsIO' );
}

diag( "Testing Bio::Mitomaster::SeqIO $Bio::Mitomaster::SeqIO::VERSION, Perl $], $^X" );

ok ( my $io = Bio::Mitomaster::SeqIO->new() );
$io = Bio::Mitomaster::SeqIO->new();
ok( !defined( eval {$io->read_seqs('/tmp/non_existent_file.fasta')} ) );  # should throw an error


# FastaIO
my $fasta1 = ">seq1|haplogroup|A1\nGACTTAAACTA";
my $fasta2 = ">seq2|haplogroup|A2\nGACTTGAACCTA";
my $both_fasta = $fasta1 . "\n" . $fasta2;

# When the test flag is set, the file contents can be passed directly
# to read_seqs as a string.  
ok ( my $fasta_io = Bio::Mitomaster::FastaIO->new(test=>1) );
my $fasta_seqs = $fasta_io->read_seqs($both_fasta);
my $fasta_info1 = shift (@{$fasta_seqs});
ok ( $fasta_info1->{name} eq 'seq1' );
ok ( $fasta_info1->{haplogroup} eq 'A1' );
my $fasta_info2 = shift (@{$fasta_seqs});
ok ( $fasta_info2->{name} eq 'seq2' );
ok ( $fasta_info2->{haplogroup} eq 'A2' );

# Use this to test an actual fasta file
#my $fasta_io = Bio::Mitomaster::FastaIO->new();
#my $seqs = $fasta_io->read_seqs('/tmp/seqs2.fasta');
#for (@{$seqs}) {
#    diag ("Seq: ", $_->{name}, ' ', $_->{hap}, "\n",  $_->{seq}, "\n");
#}



# GenbankIO
ok ( my $genbank_io = Bio::Mitomaster::GenbankIO->new() );

# Use this to test an actual genbank file
#my $gen_io = Bio::Mitomaster::GenbankIO->new();
#my $seq = $gen_io->read_seqs('/tmp/sequences.gb');
#diag (
#    "locus: ", $seq->{locus}, "\n", 
#    "def: ", $seq->{definition}, "\n",  
#    "accession: ", $seq->{accession}, "\n",  
#    "version: ", $seq->{version}, "\n",  
#    "dblink: ", $seq->{dblink}, "\n",  
#    "keywords: ", $seq->{keywords}, "\n",  
#    "source: ", $seq->{source}, "\n",  
#    "organism: ", $seq->{organism}, "\n");
#my $refs = $seq->{references};
#for (@{$refs}) {
#    diag ( 
#        "authors: ", $_->{AUTHORS}, "\n",
#        "title: ", $_->{TITLE}, "\n",
#        "journal: ", $_->{JOURNAL}, "\n",
#        "medline: ", $_->{MEDLINE}, "\n",
#        "pubmed: ", $_->{PUBMED}, "\n"
#    );
#}
#diag ($seq->{seq}, "\n");



# VariantsIO
my $variants1 = ">seq1\n100A\n200G\n400.001iCC\n500-\n600---\n700-7";
my $variants2 = ">seq2\n800a\n900g\n1000.002icc\n1100-\n1200---\n1300-7";
my $both_variants = $variants1 . "\n" . $variants2;

# When the test flag is set, the file contents can be passed directly
# to read_seqs as a string.  
ok ( my $variant_io = Bio::Mitomaster::VariantsIO->new(test=>1) );
$fasta_seqs = $variant_io->read_seqs($both_variants);
$fasta_info1 = shift (@{$fasta_seqs});
ok ( $fasta_info1->{name} eq 'seq1' );
my $fasta_var1 = $fasta_info1->{variants};
ok ( $fasta_var1->{'200'} eq 'G' );
ok ( $fasta_var1->{'400.001'} eq 'CC' );
$fasta_info2 = shift (@{$fasta_seqs});
ok ( $fasta_info2->{name} eq 'seq2' );
my $fasta_var2 = $fasta_info2->{variants};
ok ( $fasta_var2->{'900'} eq 'G' );
ok ( $fasta_var2->{'1200'} eq '---' );


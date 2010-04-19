#!perl -T

use Test::More tests => 15;

BEGIN {
	use_ok( 'Bio::Mitomaster::RNASeq' );
}


diag( "Testing Bio::Mitomaster::MitoSeq $Bio::Mitomaster::MitoSeq::VERSION, Perl $], $^X" );

use Bio::Mitomaster::Seq;
use Bio::Mitomaster::SpeciesRef;
my $ref = Bio::Mitomaster::SpeciesRef->new(species=>'human', reference=>'rCRS');


my $rna_seq = Bio::Mitomaster::RNASeq->new(species_ref=>$ref, locus_id=>33, variants=>{});

# Instantiation, attributes, and delegation
ok( !$rna_seq->wrapping );  # RNASeq objects are non-wrapping
$rna_seq->wrapping(1);
ok( $rna_seq->wrapping );  # RNASeq objects are non-wrapping
$rna_seq->wrapping(0);
#diag('length: ', length($rna_seq->seq), "\n");
#diag('seq: ', $rna_seq->seq, "\n");


# Seq Manipulation

# Simple insertion
$rna_seq = Bio::Mitomaster::RNASeq->new(species_ref=>$ref, locus_id=>33, variants=>{2=>'C', 4.01=>'GG'});
# This rna_seq is the same as this dna seq transcribed.
# my $seq = Bio::Mitomaster::Seq->new(species_ref=>$ref, variants=>{14669.01=>'CC', 14672=>'G'});
# my $transcript = $seq->transcribe(33);
ok ( $rna_seq->seq(1,5) eq 'ACGAGGU' );
#diag('length: ', length($rna_seq->seq), "\n");
#diag('seq: ', $rna_seq->seq, "\n");


# Simple deletion
$rna_seq = Bio::Mitomaster::RNASeq->new(species_ref=>$ref, locus_id=>33, variants=>{4=>'--'});
ok ( $rna_seq->seq(4) eq '-' );
my $raw_seq = $rna_seq->seq;
$raw_seq =~ /-/g;
ok ( length($raw_seq) == 525 );


# Boundary checking
ok ( !defined( eval{$rna_seq = Bio::Mitomaster::RNASeq->new(species_ref=>$ref, locus_id=>33, variants=>{}, start=>5, end=>5000)}));


# Partial Seq
ok ( $rna_seq = Bio::Mitomaster::RNASeq->new(species_ref=>$ref, locus_id=>33, variants=>{}, start=>5, end=>15));
ok ( $rna_seq->seq eq 'UGUAUGCUUUG');
ok ( !defined( eval{$rna_seq->seq(1,5)}) );  # boundary checking of seq method




# Translation
$rna_seq = Bio::Mitomaster::RNASeq->new(species_ref=>$ref, locus_id=>37, variants=>{});
ok( !defined(eval {$rna_seq->translate}) );  # Translating a non-coding region should throw an error.
$rna_seq = Bio::Mitomaster::RNASeq->new(species_ref=>$ref, locus_id=>33, variants=>{});
ok( my $aa_seq = $rna_seq->translate );
ok( $aa_seq->ref_seq(1,5) eq 'MMYAL' );


$rna_seq = Bio::Mitomaster::RNASeq->new(species_ref=>$ref, locus_id=>33, variants=>{2=>'C'});
ok($aa_seq = $rna_seq->translate);
ok( $aa_seq->seq(1,5) eq 'TMYAL' );

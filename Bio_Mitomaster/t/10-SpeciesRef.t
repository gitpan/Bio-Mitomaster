#!perl -T

use Test::More tests => 19;

BEGIN {
	use_ok( 'Bio::Mitomaster::SpeciesRef' );
}

diag( "Testing Bio::Mitomaster::SpeciesRef $Bio::Mitomaster::SpeciesRef::VERSION, Perl $], $^X" );

my $ref = Bio::Mitomaster::SpeciesRef->new(species=>'human', reference=>'rCRS');
ok( $ref->ref_seq );
ok( $ref->dna_wrapping() == 1 );
ok( $ref->wrapping() == 1 );  # alias for dna_wrapping
ok( $ref->rna_wrapping() == 0 );
ok( $ref->aa_wrapping() == 0 );
ok( $ref->start() == 1 );
ok( $ref->end() == 16569 );
ok( $ref->codon_code );
ok( $ref->codon_code('TTT') eq 'F' );
ok( $ref->codon_code('UUU') eq 'F' );
ok( $ref->locus(16)->{'start'} == 5904 );
ok( $ref->protein(16)->{'predicted_weight'} == 57000 );
#diag( $ref->ref_seq(100,105) );
ok( $ref->ref_seq(100,105) eq 'GGAGCC' );
ok( $ref->transcript(16) );
ok( $ref->transcript(16,1,3) eq 'AUG' );
ok( $ref->translation(16) );
ok( $ref->translation(16,1,3) eq 'MFA' );
ok( $ref->translation(16,3,3) eq 'A' );

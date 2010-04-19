#!perl -T

use Test::More tests => 17;

BEGIN {
	use_ok( 'Bio::Mitomaster::AASeq' );
}


diag( "Testing Bio::Mitomaster::AASeq $Bio::Mitomaster::MitoSeq::VERSION, Perl $], $^X" );

use Bio::Mitomaster::SpeciesRef;
my $ref = Bio::Mitomaster::SpeciesRef->new(species=>'human', reference=>'rCRS');



my $aa_seq = Bio::Mitomaster::AASeq->new(species_ref=>$ref, locus_id=>33, variants=>{2=>'UUG', 3=>'AUG +1', 4=>'CUU +1', 5=>'UGG'});
ok( !$aa_seq->wrapping );  # AASeq objects are non-wrapping
ok ( !$aa_seq->show_codons );  # show_codons should be false by default
ok ( !$aa_seq->show_frames );  # show_frames should be false by default

my $aav = $aa_seq->variants;
ok ( $aav->{2} eq 'L' );
ok ( $aav->{3} eq 'M' );
ok ( $aav->{4} eq 'L' );
ok ( $aav->{5} eq 'W' );

$aa_seq->show_codons(1);
$aav = $aa_seq->variants;
ok ( $aav->{2} eq 'UUG' );
ok ( $aav->{3} eq 'AUG' );
ok ( $aav->{4} eq 'CUU' );
ok ( $aav->{5} eq 'UGG' );

$aa_seq->show_frames(1);
$aav = $aa_seq->variants;
ok ( $aav->{2} eq 'UUG' );
ok ( $aav->{3} eq 'AUG +1' );
ok ( $aav->{4} eq 'CUU +1' );
ok ( $aav->{5} eq 'UGG' );

ok ( $aa_seq->seq(2,5) eq 'LMLW' );
#diag ( $aa_seq->ref_seq(1,5) );
#diag ( $aa_seq->seq(1,5) );

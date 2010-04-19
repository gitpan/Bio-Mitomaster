#!perl -T

use Test::More tests => 22;

BEGIN {
	use_ok( 'Bio::Mitomaster::Seq' );
}

diag( "Testing Bio::Mitomaster::Seq $Bio::Mitomaster::Seq::VERSION, Perl $], $^X" );


use Bio::Mitomaster::SpeciesRef;
my $ref = Bio::Mitomaster::SpeciesRef->new(species=>'human', reference=>'rCRS');


# Seq objects do very little besides delegate to their SpeciesRef object or inherit
# methods from MitoSeq.


# Instantiation, attributes, and delegation
my $seq = Bio::Mitomaster::Seq->new(species_ref=>$ref, variants=>{2=>'-', 4=>'G', 6.01=>'CC'});
ok( $seq->species eq 'human' );
ok( $seq->ref_seq(14660, 14673) eq 'AAAGCATACATCAT' );
ok( $seq->codon_code('TTT') eq 'F' );
#diag( $seq->seq(1,10));
ok( $seq->seq(1,10) eq 'G-TGACCCAGGT');
ok( $seq->wrapping );  # The wrapping flag should have been automatically set to 1
ok( $seq->start == 1);
ok( $seq->end == 16569);



# Seq manipulation

# Basic seq with deletion and insertions
$seq = Bio::Mitomaster::Seq->new(species_ref=>$ref, variants=>{2=>'-', 4.01=>'G', 6.01=>'CC', 8=>'T'});
#diag ('seq: ', $seq->seq(3,10), "\n" );
#diag ('ref: ', $seq->ref_seq(3,10), "\n" );

ok( $seq->seq(5,15) eq 'ACCCATGTCTATC');
#diag ('seq: ', $seq->seq(5,15), "\n" );
#diag ('ref: ', $seq->ref_seq(5,15), "\n" );

# Deletions and insertions and crosses the origin
$seq = Bio::Mitomaster::Seq->new(species_ref=>$ref, variants=>{16567=>'G', 2=>'-', 4.01=>'CC'});
ok( $seq->seq(16565,6) eq 'CGGTGG-TCCCAC');
#diag ('seq: ', $seq->seq(16565,6), "\n" );
#diag ('ref: ', $seq->ref_seq(16565,6), "\n" );

# Deletions, insertions, and crosses the origin
$seq = Bio::Mitomaster::Seq->new(species_ref=>$ref, variants=>{100.01=>'A', 16567.01=>'CC', 2=>'-', 4=>'G'});
ok( $seq->seq(16565,6) eq 'CGACCTGG-TGAC');
#diag ('seq: ', $seq->seq(16565,6), "\n" );
#diag ('ref: ', $seq->ref_seq(16565,6), "\n" );


# Variant access 
$seq = Bio::Mitomaster::Seq->new(species_ref=>$ref, variants=>{100.01=>'A', 16567.01=>'CC', 2=>'-', 3=>'-',4=>'G'});
my $v = $seq->variants;
ok( $v->{100.01} eq 'A' );
ok( $v->{16567.01} eq 'CC' );
ok( $v->{2} eq '-' );
ok( $v->{3} eq '-' );
ok( $v->{4} eq 'G' );


# Transcription
ok( !defined(eval {$seq->transcribe(72)}) );  # transcribing a non-coding region should throw an error.
ok( $transcript = $seq->transcribe(33) );
ok( ref $transcript and $transcript =~ /RNASeq/ );
$seq = Bio::Mitomaster::Seq->new(species_ref=>$ref, variants=>{14669.01=>'CC', 14672=>'G'});
my $transcript = $seq->transcribe(33);
ok ( $transcript->seq(1,5) eq 'ACGAGGU' );

$seq = Bio::Mitomaster::Seq->new(species_ref=>$ref, variants=>{14669=>'--', 14672=>'G'});
ok( $transcript = $seq->transcribe(33) );
ok ( $transcript->seq(1,5) eq 'ACG--' );


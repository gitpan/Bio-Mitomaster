#!perl -T


use Test::More tests => 15;

BEGIN {
	use_ok( 'Bio::Mitomaster::MitoSeq' );
}

diag( "Testing Bio::Mitomaster::MitoSeq $Bio::Mitomaster::MitoSeq::VERSION, Perl $], $^X" );

use Bio::Mitomaster::SpeciesRef;
my $ref = Bio::Mitomaster::SpeciesRef->new(species=>'human', reference=>'rCRS');

# Check instantiation, attributes, and delegation
my $ms = Bio::Mitomaster::MitoSeq->new(species_ref=>$ref, variants=>{2=>'A', 3=>'A'}, info=>{name=>'seq1', donor=>'Inuit'});
ok( $ms->species eq 'human' );
ok( $ms->start == 1 );
ok( $ms->end == 16569 );
ok( $ms->info('name') eq 'seq1' );
$ms->info('name','seq2');
ok( $ms->info('name') eq 'seq2' );
ok( $ms->wrapping == 1 );
ok( $ms->ref_seq(16567,3) eq 'ATGGAT' );
ok( $ms->seq(1,4) eq 'GAAC' );
#diag( $ms->seq(1,4) );
ok( $ms->seq(16567,3) eq 'ATGGAA' );
#diag( 'sub_seq: ', $ms->seq(16567,3), "\n" );



##### Seq Manipulation #####

# Crossing the origin with an insertion
$ms = Bio::Mitomaster::MitoSeq->new(species_ref=>$ref, variants=>{2=>'A', 3=>'A', 
    16569.01=>'C'}, info=>{name=>'seq1', donor=>'Inuit'});
ok( $ms->seq(16567,3) eq 'ATGCGAA' );
#diag( 'sub_seq: ', $ms->seq(16567,3) );


# Partial sequence and out-of-bounds checking
$ms = Bio::Mitomaster::MitoSeq->new(species_ref=>$ref, start=>1, end=>10, 
    variants=>{3=>'A'}, info=>{name=>'seq1', donor=>'Inuit'});
ok( $ms->seq eq 'GAACACAGGT' );
ok ( !defined(eval { $ms->seq(1,11) }));


# Instantiation boundary checking
ok ( !defined( eval { $ms = Bio::Mitomaster::MitoSeq->new(species_ref=>$ref, 
    start=>1, end=>16570, variants=>{3=>'A'}, info=>{name=>'seq1', donor=>'Inuit'}) } ) );


# Variant boundary checking
ok ( !defined(eval{$ms = Bio::Mitomaster::MitoSeq->new(species_ref=>$ref, start=>1, end=>10, 
        variants=>{11=>'A'}, info=>{name=>'seq1', donor=>'Inuit'})}) );


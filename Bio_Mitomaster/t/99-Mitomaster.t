#!perl -T

use Test::More tests => 14;

BEGIN {
	use_ok( 'Bio::Mitomaster' );
}

diag( "Testing Bio::Mitomaster $Bio::Mitomaster::VERSION, Perl $], $^X" );

my $mm = Bio::Mitomaster->new(species=>'human', reference=>'rCRS');
ok( $mm->reference eq 'rCRS' );
ok( $mm->species eq 'human' );

# Test delegation to the SpeciesRef object
ok( $mm->transcript(16,1,5) eq 'AUGUU' );
ok( $mm->transcript(16,5) eq 'U' );
ok( $mm->transcript(16,5,5) eq 'U' );
ok( $mm->translation(16,1,5) eq 'MFADR' );
ok( $mm->translation(16,5) eq 'R' );
ok( $mm->translation(16,5,5) eq 'R' );
#diag( $mm->translation(16,1,5) );


# Test ref_seq retrieval
ok( $mm->ref_seq());
ok( $mm->ref_seq(16568) eq 'T');
ok( $mm->ref_seq(3, 5) eq 'TCA');

# Test instantiation of Seq objects
my $seq;
ok($seq = $mm->seq(start=>1, end=>5, variants=>{1=>'T',4=>'T'}));
ok($seq->seq eq 'TATTA');

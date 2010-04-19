package Bio::Mitomaster::AASeq;

use Moose;
use Moose::Util::TypeConstraints;
use Bio::Mitomaster::Types;
use Carp;
extends 'Bio::Mitomaster::MitoSeq';



=head1 NAME

Bio::Mitomaster::AASeq - An mtDNA translation.

=cut

our $VERSION = '0.10';


=head1 SYNOPSIS

Generally, these objects should be instantiated by calling the translate method on a Bio::Mitomaster::RNASeq object.  If for some reason you want to create an AASeq object directly, you must specify a SpeciesRef object, locus ID, and a hash reference of codon variants with frameshift information appended to the end of each variant (separated by a space).  If the sequence is partial, then you must also specify the start and end values.  When created using the translate method, an AASeq object will inherit the SpeciesRef object and locus_id from the RNASeq object that created it.  


    use Bio::Mitomaster::SpeciesRef;
    use Bio::Mitomaster::AASeq;

    my $ref = Bio::Mitomaster::SpeciesRef->new(species=>'human', reference=>'rCRS');
    my $aa_seq1 = Bio::Mitomaster::AASeq->new(species_ref=>$ref, locus_id=>33, variants=>{2=>'UUG', 3=>'AUG +1', 4=>'CUU +1', 5=>'UGG'});
    my $aa_seq2 = Bio::Mitomaster::AASeq->new(species_ref=>$ref, locus_id=>33, start=>1, end=>15, variants=>{2=>'UUG', 3=>'AUG +1', 4=>'CUU +1', 5=>'UGG'});  # same thing except that we only have the first 15 amino acids

    my $s = $aa_seq1->seq(1,5);  # $s is 'MMYAL', the first 5 amino acids
    my $r = $aa_seq1->ref_seq(1,5);  # $s is 'MLMLW', the first 5 amino acids of its ref translation



=head1 ATTRIBUTES

See Bio::Mitomaster::MitoSeq for inherited attributes.

=cut


# The species_ref object is inherited from MitoSeq, but we need to make a few adjustments.
has 'species_ref' => (#FOLDBEG
    handles =>{
        wrapping => 'aa_wrapping',
    },
);#FOLDEND


=head2 end

The end position of the sequence.  Defaults to the length of the reference sequence translation.  Set this in order to analyze partial sequences.  It will also cause the ref_seq value to end at this position.

=cut

has 'end' => (#FOLDBEG
    is => 'ro',
    isa => 'PositiveInt',
    required => 1,
    lazy_build => 1,
);

sub _build_end {
    # If no end value is supplied we use the length of the reference sequence.

    my $self = shift;
    return length($self->species_ref->translation($self->locus_id));

}#FOLDEND



=head2 locus_id 

The locus_id for this AA sequence.  Use a Mitomaster object to see a list of the valid locus_id values that can be used (see the Bio::Mitomaster documentation).  This is a read-only value and must be set during object construction.

=cut

has 'locus_id' => (#FOLDBEG
    is => 'ro',
    isa => 'PositiveInt',
    required => 1,
);#FOLDEND


=head2 variants

A hash reference to a list of variants, indexed by translation position.  Note that AASeq stores variants as codons, along with frame information, and converts the values to amino acids using the codon code in the SpeciesRef object.  This behavior can be altered by setting the show_codons and show_frames attributes.


With this flag set, a single base DNA insertion would cause a -1 frame for the encoded amino acid, while a single base deletion would result in a +1 frame.  Variants that are in frame do not have any frame information recorded.  This flag is often used in combination with the show_codons flag.   

=cut


=head2 show_codons

A flag that alters the variants method.  When set to 1, the variants will be returned as codons, rather than single letter amino acid values.  See 'variants'.

=cut

has 'show_codons' => (#FOLDBEG
    is => 'rw',
    isa => 'Bool',
    default => 0,
);#FOLDEND


=head2 show_frames

A flag that alters the variants method.  When set to 1, the variants will have frame information appended to the end (separated by a space).  See 'variants'. 

=cut

has 'show_frames' => (#FOLDBEG
    is => 'rw',
    isa => 'Bool',
    default => 0,
);#FOLDEND


# The species_ref object is inherited from MitoSeq, but we need to make a few adjustments.
has 'species_ref' => (#FOLDBEG
    handles =>{
        wrapping => 'aa_wrapping',
    },
);#FOLDEND



=head2 start

The start position of the sequence.  Defaults to 1.  Set this in order to analyze partial sequences.  It will also alter the behavior of the ref_seq method.

=cut


has 'start' => (#FOLDBEG
    is => 'ro',
    isa => 'PositiveInt',
    required => 1,
    default => 1,
);#FOLDEND


=head1 METHODS

See Bio::Mitomaster::MitoSeq for inherited methods, particularly the seq method, which is probably the most important one provided.

=cut



around 'seq' => \&_get_translation_seq;#FOLDBEG

sub _get_translation_seq {

    my ($orig, $self, $s, $e) = @_;

    # Store the current values of show_codons and show_frames
    my $stored_show_codons = $self->show_codons;
    my $stored_show_frames = $self->show_frames;

    # Now, turn them off so that the seq will contain AA values
    $self->show_codons(0);
    $self->show_frames(0);

    my $translation = $self->$orig($s, $e);

    # Restore the stored settings
    $self->show_codons($stored_show_codons);
    $self->show_frames($stored_show_frames);


    return $translation;

}#FOLDEND





around 'variants' => \&_variants_mod;#FOLDBEG

sub _variants_mod {
    my ($orig, $self) = @_;

    my $var_ref = $self->$orig;

    my %variants;
    for (keys %{$var_ref}) {
        my($c, $f) = split /\s/, $var_ref->{$_};

        my $v;
        if ($self->show_codons) {
           $v = $c;
        }
        else {
           $v = $self->codon_code($c); 
        }

        if ($self->show_frames) {
            $v = $v . ' ' . $f if $f;
        }

        $variants{$_} = $v;
    }

    return \%variants;

}#FOLDEND



sub BUILD {#FOLDBEG
    # Does some basic sanity checks when the object is constructed.

    my $self = shift;

    # retrieve locus info
    my $locus_info = $self->locus($self->locus_id);


    # check that the locus requested is an mRNA coding locus
    unless ( ($locus_info->{'type'} eq 'm') ) {
        confess "Locus ID: ", $self->locus_id, " does not refer to a coding locus";
        return undef;
    }



    return 1;

}#FOLDEND



=head1 AUTHOR

Marty Brandon, C<< <mbrandon at uci.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-bio-mitomaster-seq at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Bio-Mitomaster-AASeq>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Bio::Mitomaster::AASeq


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Bio-Mitomaster-AASeq>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Bio-Mitomaster-AASeq>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/Bio-Mitomaster-AASeq>

=item * Search CPAN

L<http://search.cpan.org/dist/Bio-Mitomaster-AASeq/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 COPYRIGHT & LICENSE

Copyright 2009 Marty Brandon, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.


=cut


#These two line are a recommended speedup.
no Moose;
__PACKAGE__->meta->make_immutable();
1; # End of Bio::Mitomaster::AASeq

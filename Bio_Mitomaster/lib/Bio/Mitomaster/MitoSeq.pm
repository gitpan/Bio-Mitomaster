package Bio::Mitomaster::MitoSeq;
use Moose;
use Moose::Util::TypeConstraints;
use Bio::Mitomaster::Types;
with 'Bio::Mitomaster::SeqManipRole';
use Carp;


=head1 NAME

Bio::Mitomaster::MitoSeq - Superclass for Mitomaster seq objects.

=cut

our $VERSION = '0.10';



=head1 SYNOPSIS

Don't instantiate this directly.  It's meant to provide common functionality to those SEQ classes that inherit from it, and those SEQ classes should always be instantiated via a Bio::Mitomaster object.


=head1 ATTRIBUTES




=head2 end

The end position of the sequence.  Defaults to the length of the reference sequence.  Set this in order to analyze partial sequences.  It will also cause the ref_seq value to end at this position.

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
    return length($self->species_ref->ref_seq);

}#FOLDEND


=head2 info

A hash reference for associating a sequence object to any value desired.  Values are made available through the info method.

 $ms->info('name', 'AJX002');    # set the seq name
 $seq_name = $ms->info('name');  # get the seq name

=cut

has 'info' => (#FOLDBEG
    is => 'rw',
    isa => 'HashRef',
    required => 0,
);

around 'info' => sub {
    my $orig = shift;
    my $self = shift;
    my $key = shift;

    my $info_ref = $self->$orig;
    if (@_) {
        # a setter
        $info_ref->{$key} = shift;        
    }
    else {
        # a getter
        return $info_ref->{$key};
    }

};#FOLDEND



=head2 start

The start position of the sequence.  Defaults to 1.  Set this in order to analyze partial sequences.  It will also cause the ref_seq value to begin at this position.

=cut


has 'start' => (#FOLDBEG
    is => 'ro',
    isa => 'PositiveInt',
    required => 1,
    default => 1,
);#FOLDEND


=head2 species_ref

A Bio::Mitomaster::SpeciesRef object.  Most of the methods available in those objects are delegated straight-through in this one.  The ref_seq method is enhanced with the Bio::Mitomaster::SeqManipRole (see below).

=cut


has 'species_ref' => (#FOLDBEG
    is => 'ro',
    isa => 'Bio::Mitomaster::SpeciesRef',
    required => 1,
    handles =>{
        species => 'species', 
        reference => 'reference',
        codon_code => 'codon_code',
        locus => 'locus',
        protein => 'protein',
        transcript => 'transcript',
        translation => 'translation',
        wrapping => 'dna_wrapping', 
    },
);#FOLDEND


=head2 variants

A hash reference holding the list of variants found in the sequence.  Key values are the positions at which the variant occurs.

=cut

has 'variants' => (#FOLDBEG
    is => 'ro',
    isa => 'HashRef',
    required => 1,
);#FOLDEND



=head1 METHODS

=cut


sub BUILD {#FOLDBEG

    # Does some basic sanity checks when the object is constructed.

    my $self = shift;
    my $wrapping = $self->wrapping;
    my $start = $self->start;
    my $end = $self->end;

    my $ref_start = 1;
    my $ref_end = $start + (length($self->ref_seq) - 1);


    # Check the start and end values to make sure they are in the same range
    # as the reference sequence.
    confess "Start position: $start is out of bounds"
        if $start < $ref_start or $start > $ref_end;

    confess "end position: $end is out of bounds  ref_end: $ref_end"
        if $end < $ref_start or $end > $ref_end;


    if (!$wrapping) {
        # A linear molecule.
        confess "Start position: $start cannot be greater than end: $end for a linear molecule"
            if $start >= $end;
    }


    # Check the position value of each variant to make sure they are within the
    # same range as the start and end for the sequence.  Make sure to allow for
    # insertions at the very last position (use $_ >= ($end+1)).
    my $variants = $self->variants;
    for (keys %{$variants}){
        if ( ($_ < $start) or ($_ >= ($end + 1)) ) {
            confess "Variant position $_ is beyond seq boundaries: $start - $end";
        }
    }

    return 1;

}#FOLDEND



=head2 ref_seq($s, $e, $f)

Retrieves all or part of the SEQ objects ref sequence.  It analyzes the caller to determine whether to return DNA, RNA, or AA.  Provide start and end idices to retrieve a partial ref seq.

=cut

sub ref_seq {#FOLDBEG

    my ($self, $s, $e) = @_;

    my ($start, $end) = $self->get_seq_indices($s, $e);
    $self->check_seq_indices($start, $end, $self->wrapping);

    my $seq_string;
    if ( $self =~ /RNA/ ) {
        return $self->species_ref->transcript($self->locus_id, $start, $end);
    }
    elsif ( $self =~ /AA/ ) {
        return $self->species_ref->translation($self->locus_id, $start, $end);
    }
    else {
        return $self->species_ref->ref_seq($start, $end);
    }

}#FOLDEND



=head2 seq
=head2 seq($position)
=head2 seq($start, $end)

Returns a string value representing the sequence.  Works by merging the list of variants with the reference sequence.  Part of a sequence can be returned by specifying start and end values.  Note that these values are absolute position numbers (relative to the reference sequence used).  So, if you have a sequence that begins at position 577, you would call seq(577,581) to get the first 5 moieties.  Calling the method with one argument will retrieve the single moiety at that position.  When called with indicies outside the declared start and end values of the sequence, an error is returned.

 $s->seq;  # the entire seq string
 $s->seq(11);  # the moiety at position 11 
 $s->seq(11,11);  # same thing
 $s->seq(11,20);  # 10 moieties from positions 11 to 20

When the wrapping attribute is set to a true value (usually set automatically from species meta data during SEQ object construction) then we can also do this:

 $s->seq(16567,3);  # moieties taken from the end and concatenated with ones from the beginning 

=cut



sub seq {#FOLDBEG

    my ($self, $s, $e) = @_;

    # We construct a complete sequence according to the ref_seq provided by
    # our SpeciesRef object and the list of variants, then we trim the seq
    # according to the indices submitted using the sub_seq method.


    # Get validated indices
    my ($start, $end) = $self->get_seq_indices($s, $e);
    $self->check_seq_indices($start, $end, $self->wrapping);


    # Retrieve the variants hash ref and if there are no variants defined
    # for this sequence then just return the corresponding reference sequence.
    my $variants = $self->variants;
    return $self->ref_seq($start, $end) unless keys %{$variants};


    # Define another set of indices that will be used when the sequence is
    # built and then needs to be trimmed.
    my ($trim_start, $trim_end) = ($start, $end);


    my $seq_string;
    # seq_string is made by looping through all the variants.  For each variant, 
    # add any preceeding ref sequence, then add the variant, and finish with any 
    # final ref segment.  The final step is to trim the sequence down to what was 
    # actually requested.  Variant positions are validated during object 
    # construction, so they should be trustworthy.  

    # Building an entire sequence can be a little expensive if we are only using
    # partial seqs, but it made it easier to implement wrapping and mostly we 
    # retrieve full seqs.

    my $i = 1;
    my $ref_end;
    if ($start > $end) {
        # wrapping
        $ref_end = length($self->species_ref->ref_genome);
    }
    else {
        # non-wrapping seq
        $ref_end = length($self->ref_seq);
    }
    

    for (sort {$a <=> $b} keys %{$variants}){

        # add the ref_seq portion
        if ($i < $_){
            $seq_string .= int($_) == $_ ? 
                $self->ref_seq($i, $_ - 1) : $self->ref_seq($i, int($_));
        }    
        # add the variant
        $seq_string .= $variants->{$_};
        

        # Adjust $trim_start and $trim_end for insertions.  Deletions are not a 
        # problem since we are using '-' to represent them.
        if ( int($_) < $_ ) {
            # An insertion

            if ($_ < $start) {
                $trim_start = $trim_start + length($variants->{$_});
            }
            if ($_ < $end) {
                $trim_end = $trim_end + length($variants->{$_});
            }
        }


        # increment the current position
        if ($variants->{$_} =~ /^-+$/) {
            $i = int($_) + length($variants->{$_});
        }
        else {
            $i = int($_) + 1;
        }

    }

	# add the final ref_seq segment unless the last position was variable
    $seq_string .= $self->ref_seq($i, $ref_end) 
        unless $i > $ref_end;


    return $self->sub_seq($seq_string, $trim_start, $trim_end);

}#FOLDEND



=head1 AUTHOR

Marty Brandon, C<< <mbrandon at uci.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-bio-mitomaster-mitoseq at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Bio-Mitomaster-MitoSeq>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Bio::Mitomaster::MitoSeq


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Bio-Mitomaster-MitoSeq>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Bio-Mitomaster-MitoSeq>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/Bio-Mitomaster-MitoSeq>

=item * Search CPAN

L<http://search.cpan.org/dist/Bio-Mitomaster-MitoSeq/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 COPYRIGHT & LICENSE

Copyright 2009 Marty Brandon, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.


=cut

no Moose;
__PACKAGE__->meta->make_immutable();
1; # End of Bio::Mitomaster::MitoSeq

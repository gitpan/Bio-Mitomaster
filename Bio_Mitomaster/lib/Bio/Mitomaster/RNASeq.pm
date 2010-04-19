package Bio::Mitomaster::RNASeq;

use Moose;
use Moose::Util::TypeConstraints;
use Bio::Mitomaster::Types;
extends 'Bio::Mitomaster::MitoSeq';



=head1 NAME

Bio::Mitomaster::RNASeq - An mtDNA transcript.

=cut

our $VERSION = '0.10';


=head1 SYNOPSIS

Generally, these objects are instantiated by calling the transcribe method on a Bio::Mitomaster::Seq object.  If for some reason you want to create an RNASeq object directly, you must specify a SpeciesRef object, locus ID, and a hash reference of variants.  If the sequence is partial, then you must also specify the start and end values.  When created using the transcribe method, an RNASeq object will inherit the SpeciesRef object from the Seq object that creates it.  


    use Bio::Mitomaster::SpeciesRef;
    use Bio::Mitomaster::RNASeq;

    my $ref = Bio::Mitomaster::SpeciesRef->new(species=>'human', reference=>'rCRS');
    my $rna_seq1 = Bio::Mitomaster::RNASeq->new(species_ref=>$ref, locus_id=>5, variants=>{1=>'T',2=>'T'});
    my $rna_seq2 = Bio::Mitomaster::RNASeq->new(species_ref=>$ref, locus_id=>5, variants=>{1=>'T',2=>'T'}, start=>1, end=>15); #same thing except that we only have the first 15 nucleotides

    my $s = $rna_seq1->seq(1,5);  #$s is 'TTTCA', the first 5 RNA nucleotides of the transcript from locus 5;



=head1 ATTRIBUTES

See Bio::Mitomaster::MitoSeq for inherited attributes.


=head2 end

The end position of the sequence.  Defaults to the length of the reference sequence transcript.  Set this in order to analyze partial sequences.  It will also cause the ref_seq value to end at this position.

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
    return length($self->species_ref->transcript($self->locus_id));

}#FOLDEND



=head2 locus_id 

The locus_id for this RNA sequence.  Use a Mitomaster object to see a list of the valid locus_id values that can be used (see the Bio::Mitomaster documentation).

=cut

has 'locus_id' => (
    is => 'ro',
    isa => 'PositiveInt',
    required => 1,
);



# The species_ref object is inherited from MitoSeq, but we need to make a few adjustments.
has 'species_ref' => (#FOLDBEG
    handles =>{
        wrapping => 'rna_wrapping',
    },
);


=head2 start

The start position of the sequence.  Defaults to 1.  Set this in order to analyze partial sequences.  It will also alter the behavior of the ref_seq method.

=cut


has 'start' => (
    is => 'ro',
    isa => 'PositiveInt',
    required => 1,
    default => 1,
);


=head1 METHODS

See Bio::Mitomaster::MitoSeq for inherited methods, particularly the seq method, which is probably the most important one provided.

=cut 


around 'seq' => \&_get_transcript_seq;#FOLDBEG

sub _get_transcript_seq {

    my ($orig, $self, $s, $e) = @_;

    my $transcript = $self->$orig($s, $e);

    my $ref_end = length($self->species_ref->transcript($self->locus_id));

    if (!$s and !$e and $self->start == 1 and $self->end == $ref_end) {
        # poly-adenylate if this is a full transcript and no start or end value was requested
        $transcript = $self->_poly_adenylate($transcript);
    }

    return $transcript;

}#FOLDEND


sub _poly_adenylate {#FOLDBEG

    my ($self, $raw_seq) = @_;

    # Poly-adenylation (PA).  How long is the actual transcript in nature?  PA completes 
    # a stop codon, but presumably there are additional unused adenosines.  If the frame 
    # is disrupted by indels, then we would naturally expect to have the potential for a 
    # runon poly-peptide.  But let's not go there for now.  We needn't model reality to that
    # level for this purpose.  Here we'll assume that the upper bound for the length of any
    # transcript is determined by the corresponding reference transcript plus enough 
    # adenosines to complete a codon.  This means that our level of poly-adenylation will change
    # according to the presence of indels -- deletions will cause additional adensines, 
    # while insertions will cause fewer, perhaps even truncation of the regular sequence.
    # The other affect is that our in silico transcript is not guaranteed to end with a 
    # stop codon.  This should not be too problematic.  If during an analysis we should
    # have a transcript that does not end in a stop codon, then we can infer what is
    # occurring.


    # Get the species_ref first, and use it to retrieve the ref transcript so that we 
    # can avoid truncating the ref seq according to the start and end values.
    my $sr = $self->species_ref;
    my $ref_length = length($sr->transcript($self->locus_id));
    my $target_length = length($ref_length) % 3 ?
        $ref_length + (3 - ( length($ref_length) % 3 ) ) :
        $ref_length;
    
    # Remove any deletion placeholders from the position to get an accurate length
    $raw_seq =~ s/-//g;

    if (length($raw_seq) < $target_length) {
        # poly-adenylate
        $raw_seq .= 'A' x ( $target_length - length($raw_seq) );
    }
    elsif (length($raw_seq) > $target_length) {
        # truncate
        $raw_seq = substr($raw_seq, 0, $target_length );
    }

    return $raw_seq;

}#FOLDEND

 

sub BUILD {#FOLDBEG

    # Does some basic sanity checks when the object is constructed.

    my $self = shift;
    my $start = $self->start;
    my $end = $self->end;
    my $ref_start = 1;
    my $ref_end = length($self->species_ref->transcript($self->locus_id));


    # Check the start and end values to make sure they are in the same range
    # as the reference sequence.
    confess "Start position: $start is out of bounds"
        if $start < $ref_start or $start > $ref_end;

    confess "end position: $end is out of bounds  ref_end: $ref_end"
        if $end < $ref_start or $end > $ref_end;

    return 1;

}#FOLDEND


=head2 translate

Returns a Bio::Mitomaster::AASeq object representing the translation that would be produced by translating the transcript sequence of this object.  This method cannot be called on partial transcript sequences.  Normally, you would begin with a Bio::Mitomaster object and use its seq method to produce a Seq object.  Then use the transcribe method of the Seq object to produce and RNASeq object, and finally use the translation method of the RNASeq object to produce an AASeq object.

=cut

sub translate {

    my ($self) = @_;


    # Compare the locus boundaries with our start and end values to make sure the 
    # requested locus is covered.  Fortunately no coding loci span the origin and
    # "beginning" values are always with respect to the lowest numbered nucleotide
    # on the heavy strand.
    if ($self->start > 1) {
        confess "Start value: ", $self->start, " should begin at 1 for translate method";
        return undef;
    }
    if ($self->end < length($self->ref_seq) ) {
        confess "End value: ", $self->end, " must be greater than the ref transcript length ",
            "for translate method";
        return undef;
    }


    # Traverse the transcript, building codons, and recording variant codons.
    # Essentially, we keep three queues one for insertion nucleotides, one for
    # variant positions, and one for ref positions (deletions only adjust the
    # frame value, so they don't need a queue).  Each codon is built by taking
    # the nucleotide from the queue with the highest precedence.  The insert
    # queue has highest precedence, followed by polymorphic positions, and 
    # lastly reference nucleotides.
    my $rna_variants = $self->variants;
    my @r_positions = (1 .. length($self->species_ref->transcript($self->locus_id)));
    my @v_positions = sort{$a <=> $b} keys %{$rna_variants};
    my %aa_variants;
    my $insert_queue;
    my $aa_pos = 1;
    my $frame = 0;
    my $codon;
    while ( @v_positions ) {

        if ($insert_queue) {
            # Preferentially take a base from the insert queue
            $codon .= substr($insert_queue, 0, 1);
            $insert_queue = substr($insert_queue, 1);
            $frame--;
        }
        elsif ($v_positions[0] <= $r_positions[0]) {
            # If there is no insert_queue, then check to see if there
            # is a variant at that position (could also be an insertion
            # between the positions).  Also give special treatment to 
            # deletions.

            if ($rna_variants->{$v_positions[0]} =~ /-/) {
                # A deletion
                $frame += length($rna_variants->{$v_positions[0]});
                shift @v_positions;
                next;
            }


            if (int($v_positions[0]) < $v_positions[0]) {
                # An insertion
                $frame -= length($rna_variants->{$v_positions[0]});
                my $insert = $rna_variants->{$v_positions[0]};
                $codon .= substr($insert, 0, 1);
                $insert_queue .= substr($insert, 1) if $insert;
            }
            else {
                # A polymorphism
                $codon .= $rna_variants->{$v_positions[0]};
            }

            shift @v_positions;

        }
        else {
            # No insert_queue or other variant at this position so grab
            # a reference nucleotide.
            $codon .= $self->ref_seq($r_positions[0]);
            shift @r_positions;
        }


        if (length($codon) == 3) {
            # Check to see if we've completed a codon
            $aa_variants{$aa_pos} = "$codon";
            $aa_variants{$aa_pos} .= " $frame" if $frame;
            $codon = '';
            $aa_pos++;
        }

        unless (@v_positions) {
            # When there are no more variants, complete the codon and exit.

            if ($codon) {
                # complete the partial codon
                $codon .= $self->ref_seq($r_positions[0]);
                if (length($codon) < 3) {
                    shift @r_positions;
                    $codon .= $self->ref_seq($r_positions[0]);
                }

                $aa_variants{$aa_pos} = "$codon";
                $aa_variants{$aa_pos} .= " $frame" if $frame;
                $codon = '';
            }
            last;
        }
    }

    eval {require Bio::Mitomaster::AASeq;};
    my $aa_seq = Bio::Mitomaster::AASeq->new(species_ref=>$self->species_ref, 
        locus_id=>$self->locus_id, variants=>\%aa_variants);

    return $aa_seq;

}



=head1 AUTHOR

Marty Brandon, C<< <mbrandon at uci.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-bio-mitomaster-seq at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Bio-Mitomaster-RNASeq>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Bio::Mitomaster::RNASeq


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Bio-Mitomaster-RNASeq>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Bio-Mitomaster-RNASeq>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/Bio-Mitomaster-RNASeq>

=item * Search CPAN

L<http://search.cpan.org/dist/Bio-Mitomaster-RNASeq/>

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
1; # End of Bio::Mitomaster::RNASeq

package Bio::Mitomaster::Seq;

use Moose;
extends 'Bio::Mitomaster::MitoSeq';


=head1 NAME

Bio::Mitomaster::Seq - An mtDNA sequence.

=cut

our $VERSION = '0.10';


=head1 SYNOPSIS

These objects should almost always be instantiated by calling the seq method on a Bio::Mitomaster object.  Each Seq object has an associated SpeciesRef object (inherited attribute from MitoSeq).  The SpeciesRef is an internal reference for meta-data, such as locus positions and the reference sequence.  A Seq is specific for a particular species and reference sequence, and once instantiated the context cannot be changed.  To create a Seq object you must provide the SpeciesRef, a hash of variants, and the start and end values of the sequence.  Bio::Mitomaster encapsulates the logic of creating the right SpeciesRef object, reading seqs from files, and performing the needed alignment.  


    use Bio::Mitomaster::SpeciesRef;
    use Bio::Mitomaster::Seq;

    my $ref = Bio::Mitomaster::SpeciesRef->new(species=>'human', reference=>'rCRS');
    my $seq1 = Bio::Mitomaster::Seq->new(species_ref=>$ref, variants=>{1=>'T',2=>'T'});
    my $seq2 = Bio::Mitomaster::Seq->new(species_ref=>$ref, variants=>{1=>'T',2=>'T'}, start=>1, end=>577); #same thing except that we only have the first 577 nucleotides

    my $s = $seq1->seq(1,5);  #$s is 'TTTCA';
    my $r_seq = $seq1->transcribe(1);  #$r_seq is an RNASeq object containing the transcript for the locus with ID 1.  



=head1 ATTRIBUTES

See Bio::Mitomaster::MitoSeq for inherited attributes.

=cut 

# The species_ref object is inherited from MitoSeq, but we need to make a few adjustments.
has 'species_ref' => (#FOLDBEG
    handles =>{
        wrapping => 'dna_wrapping',
    },
);#FOLDEND




=head1 METHODS

See Bio::Mitomaster::MitoSeq for inherited methods, particularly the seq method, which is probably the most important one provided.

=cut


sub BUILD {#FOLDBEG
    # Does some basic sanity checks when the object is constructed.

    my $self = shift;
    return 1;

}#FOLDEND


=head2 transcribe($locus_id)

Returns a Bio::Mitomaster::RNASeq object representing the transcript that would be produced by transcribing the locus identified by $locus_id.  The boundaries of the locus to transcribe are retrieved from the species_ref attribute (a Bio::Mitomaster::SpeciesRef object).  If the locus boundaries fall outside the sequence start and end values an error is returned.  The RNA sequence returned within the RNASeq object is adjusted for all expected natural processes (poly-adenylation, uracil replaces thymine, etc.). 

=cut

sub transcribe {#FOLDBEG
    my ($self, $locus_id) = @_;

    # Note that poly-adenlation concerns are handled by a modifier for the seq
    # method in RNASeq.

    # retrieve locus info
    my $locus_info = $self->locus($locus_id);


    # check that the locus requested is a coding locus
    unless ( ($locus_info->{'type'} eq 'm') || ($locus_info->{'type'} eq 'r')
            || ($locus_info->{'type'} eq 't') ) {

        confess "Locus ID: $locus_id does not refer to a coding locus";
        return undef;
    }


    # Compare the locus boundaries with our start and end values to make sure the 
    # requested locus is covered.  Fortunately no coding loci span the origin and
    # "beginning" values are always with respect to the lowest numbered nucleotide
    # on the heavy strand.
    if ($self->start > $locus_info->{'start'}) {
        confess "Start value: ", $self->start, " is greater than the locus start: ", $locus_info->{'start'}, "\n";
        return undef;
    }
    if ($self->end < $locus_info->{'end'}) {
        confess "End value: ", $self->end, " is less than the locus end: ", $locus_info->{'end'}, "\n";
        return undef;
    }


    # Get some locus info
    my $strand = $locus_info->{strand};
    my $start = $locus_info->{start};
    my $end = $locus_info->{end};
    my $dna = $self->ref_seq($start, $end);


    # Loop through the DNA variants converting them to RNA variants
    my $dna_variants = $self->variants;
    my %rna_variants;

    for my $p (keys %{$dna_variants}) {

        # Check to see if the variant resides within this locus
        if ($p < $start or $p > $end) {
            next;
        }

        my $v = $dna_variants->{$p};

        if ($strand eq 'L') {
            # make adjustments for light strand loci
            $v =~ tr/ACGT/TGCA/;

            if ( ($p - int($p)) > 0 ) {
                # an insertion
                # This formatting will get rid of trailing zeros due to rounding the
                # subtraction of the position values.  It also limits the max number 
                # of digits (before and after the decimal) to 12.
                $p = sprintf("%.12g", ( ($end - int($p)) + ($p - int($p)) ) );
            }
            else {
                # a polymorphism or deletion
                $p = ( $end - $p ) + 1;

                if ($v =~ /--/) {
                    # multiple deletions on the light strand cause the positions
                    # to be shifted during transcription
                    $p = $p - length($v) + 1; 
                }
            }
        }
        elsif ($strand eq 'H') {
            # make adjustments for heavy strand loci
            $p = ($p - $start) + 1;
        }
        else {
            # Strand is not 'H' or 'L' - an error.
            confess "Invalid strand value: $strand";
        }
        
        # Replace thymine with uracil and add the variant to our collection of RNA variants
        $v =~ s/T/U/g;
        $rna_variants{$p} = $v;
    }


    eval {require Bio::Mitomaster::RNASeq;};
    my $rna_seq = Bio::Mitomaster::RNASeq->new(species_ref=>$self->species_ref, 
        locus_id=>$locus_id, variants=>\%rna_variants);

    return $rna_seq;

}#FOLDEND



=head1 AUTHOR

Marty Brandon, C<< <mbrandon at uci.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-bio-mitomaster-seq at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Bio-Mitomaster-Seq>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Bio::Mitomaster::Seq


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Bio-Mitomaster-Seq>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Bio-Mitomaster-Seq>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/Bio-Mitomaster-Seq>

=item * Search CPAN

L<http://search.cpan.org/dist/Bio-Mitomaster-Seq/>

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
1; # End of Bio::Mitomaster::Seq

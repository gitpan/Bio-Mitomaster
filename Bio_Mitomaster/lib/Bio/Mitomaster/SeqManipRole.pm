package Bio::Mitomaster::SeqManipRole;
use Moose::Role;


=head1 NAME

Bio::Mitomaster::SeqManipRole - A role for objects that manipulate indexed mitochondrial sequences.

=cut


our $VERSION = '0.10';


=head1 SYNOPSIS

This role provides methods for checking position indices and retrieving partial sequences using absolute position values.  Mostly, these are utility methods intended to be used by the seq method of an object class (look at code of MitoSeq.pm)

 $seq = MySeqClass->new;  # MySeqClass consumes this role and provides 'start' and 'end' methods.
 ($start, $end) = $seq->get_seq_indices($s,$e);
 $seq->check_seq_indices($start,$end,$seq->wrapping); 
 $seq->sub_seq($start, $end);

=cut



=head1 ATTRIBUTES



=head1 METHODS

=head2 check_seq_indices

A utility method for doing sanity checks on sequence indices.  Usually, called after get_seq_indices to check the indices it returns.

=cut

sub check_seq_indices {#FOLDBEG
    # A utility method for doing sanity checks on sequence indices.  Usually,
    # called after get_seq_indices to check the indices it returns.
    

    my ($self, $s, $e, $wrapping) = @_;

    confess "check_seq_indices requires two integer arguments" 
        unless $s and $e and $s =~ /\d+/ and $e =~ /\d+/;

    if (!$wrapping) {
        # A linear molecule.  Either we have a molecule whose positions wrap 
        # or one in which the positions do not wrap.
        if ($self->start <= $self->end) {
            # Linear molecule without wrapping positions numbers
            confess "Start index: $s is greater than ending index: $e"
                if ($s > $e);

            confess "Start index: $s is before the start position"
                if ($s < $self->start);

            confess "End index: $e is beyond the end position"
                if ($e > $self->end);
        }
        else {
            # Linear molecule with wrapping position numbers, i.e. go to 
            # the end and startover at 1.  The region greater than
            # end and less than start is out of bounds.

            confess "End index: $e is not part of the sequence"
                if ($e > $self->start and $e < $self->end);
            
            confess "Start index: $s is not part of the sequence"
                if ($s > $self->start and $s < $self->end);
        }
    }
    else {
        # A wraped molecule.  Check that both indices are within its range.

        confess "Start index: $s is not part of the sequence"
            if ($s < $self->start or $s > $self->end);

        confess "End index: $e is not part of the sequence"
            if ($e < $self->start or $e > $self->end);

    }

    return 1;

}#FOLDEND


=head2 get_seq_indices

A utility method giving polymorphic behavior to the seq method, allowing the user to call the method with one, two, or no arguments. 

=cut


sub get_seq_indices {#FOLDBEG
    # A utility method giving polymorphic behavior to the seq method, 
    # allowing the user to call the method with one, two, or no arguments.
    # Generally, this method is used only by the seq method to return a
    # pair of indices.  


    my ($self, $s, $e) = @_;
    my ($start, $end);

    if ( $s ) {

        confess "Beginning index: $s is not a positive integer." 
            unless $s =~ /\d+/; 
            
        confess "Beginning index: $s is less than zero." if $s <= 0;

        if ( $e ) {

            confess "Ending index: $e is not a positive integer." 
                unless $e =~ /\d+/;

            confess "Ending index: $e is less than zero." if $e <= 0;

            # Both positions were given
            $start = $s;
            $end = $e;
        }
        else {
            # A single position was given
            $start = $s;
            $end = $start;
        }
    }
    else {
        # No arguments were given.
        $start = $self->start;
        $end = $self->end;
    }

    return ($start, $end);

}#FOLDEND

    

=head2 sub_seq

Submit a string and indices representing position numbers (first letter is position 1) and it returns a sub-string.  It assumes wrapping if the beginning indice is greater than the ending indice.  No checking is done on the indices since it is assumed that you call this method following a calls to get_seq_indices and check_seq_indices (See the code of the seq method in MitoSeq.pm).

 $r = 'GATCCATAAT';
 $seq->sub_seq($r,1,3);  # returns 'GAT'
 $seq->sub_seq($r,10,2);  # returns 'TGA'

=cut

sub sub_seq {#FOLDBEG
    my ($self, $string, $start, $end) = @_;

    confess "A string and two position values are needed by sub_seq" 
        unless $string and $start and $end;

    my $len = length($string);
    my $sub_seq;

    if ($start > $end) {
        # We are wrapping at the end of the indices
        $sub_seq = substr($string, $start - ($len + 1));
        $sub_seq .= substr($string, 0, ( $end ));

    }
    else {
        $sub_seq = substr($string, $start - 1, ($end - $start + 1));
    }

    return $sub_seq; 


}#FOLDEND



=head1 AUTHOR

Marty Brandon, C<< <marty.brandon at gmail.com> >>



=head1 BUGS

Please report any bugs or feature requests to C<bug-bio-mm-refmitohuman at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Bio-Bio::Mitomaster-RefMitoHuman>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.



=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Bio::Mitomaster::SeqManipRole

=over 4

=item * RT: CPAN's request tracker

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Bio-Bio::Mitomaster-RefMitoHuman>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Bio-Bio::Mitomaster-RefMitoHuman>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/Bio-Bio::Mitomaster-RefMitoHuman>

=item * Search CPAN

L<http://search.cpan.org/dist/Bio-Bio::Mitomaster-RefMitoHuman/>

=back


=head1 COPYRIGHT & LICENSE

Copyright 2008 Marty Brandon, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.


=cut

1; # End of Bio::Mitomaster::SeqManipRole

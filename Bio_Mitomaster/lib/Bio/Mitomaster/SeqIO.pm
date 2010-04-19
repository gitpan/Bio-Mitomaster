package Bio::Mitomaster::SeqIO;
use Moose;



=head1 NAME

Bio::Mitomaster::SeqIO - Superclass for all file IO objects.

=cut

our $VERSION = '0.10';



=head1 SYNOPSIS

Don't instantiate this directly.  It's meant to provide common functionality to file IO classes that inherit from it.


=head1 ATTRIBUTES

=cut

# a flag used to give differential behavior during testing
has 'test' => (
    is => 'rw',
    isa => 'Bool',
    default => 0,
);

=head1 METHODS

These are all basically helper methods for reading/writing files.  Some are expected to be overridden or modified by the subclasses.

=cut



=head1 convert_to_unix

Takes a string and converts any DOS or Apple line endings to unix line endings.

=cut

sub convert_to_unix {#FOLDBEG
    # Takes a string and converts DOS and Apple line endings to unix line endings.

    my $self = shift;
    my $seq_text = shift;
    $seq_text =~ s/\r\n/\n/g;
    $seq_text =~ s/\r/\n/g;

    return $seq_text;
}#FOLDEND


=head2 read_check

Checks to see if a file location is readable.

=cut

sub read_check {#FOLDBEG
    # Basic checks for the readability of a seq file

    my $self = shift;
    my $file_name = shift;
        
    # Check that the file exists, is not empty, is readable, and is a text file.
 	-e $file_name || confess "File: $file_name does not exist";
	-s $file_name || confess "File: $file_name has no data"; 
	-r $file_name || confess "File: $file_name is not readable"; 
	-T $file_name || confess "File: $file_name is not a text file"; 
	
    return 1;
}#FOLDEND


=head1 read_seqs

Used to read a seq file.  It takes the file name and a subroutine to a seq parser (defined in the subclass) and returns a hashref whose indices are the seq identifiers, each of which points to a hash of information about a particular seq. 

=cut

sub read_seqs {#FOLDBEG
    my $self = shift;
    my $file_name = shift;
    my $seq_parser = shift;

    confess "No seq parser given to read_seqs.  Are you trying to call the superclass method?" 
        unless $seq_parser;


    my $seq_text;

    if ($self->test) {
        # To test this module we send in a string of text for $file_name
        $seq_text = $file_name;
    }
    else {
        # Do basic file checks
        $self->read_check($file_name);
    
        # Read the entire file into a string
        open my $INFILE, '<', $file_name  or
            confess "Can't open file: $file_name";
        my @lines = <$INFILE>;
        $seq_text = join '', @lines;
    }



    # Convert line endings 
    $seq_text = $self->convert_to_unix($seq_text);

    # Return a reference to an array of hashes each containing info for a sequence
    return &$seq_parser($seq_text);

}#FOLDEND


=head2 write_check

Checks to see if a file location is writable.

=cut


sub write_check {#FOLDBEG

    # Basic write checks for creating a new seq file

    my $file_name = shift;
        
    # Check that the new file location is writable.
 	-w $file_name || confess "File: $file_name is not writable";
	
    return 1;

}#FOLDEND


sub write_seqs {#FOLDBEG
    confess "You reached a superclass method that should have been overridden";
}#FOLDEND



=head1 AUTHOR

Marty Brandon, C<< <marty.brandon at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-bio-mitomaster-mitoseq at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Bio-Mitomaster-SeqIO>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.



=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Bio::Mitomaster::SeqIO


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Bio-Mitomaster-SeqIO>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Bio-Mitomaster-SeqIO>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/Bio-Mitomaster-SeqIO>

=item * Search CPAN

L<http://search.cpan.org/dist/Bio-Mitomaster-SeqIO/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 COPYRIGHT & LICENSE

Copyright 2009 Marty Brandon, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.


=cut

no Moose;
__PACKAGE__->meta->make_immutable();
1; # End of Bio::Mitomaster::SeqIO

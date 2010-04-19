package Bio::Mitomaster::FastaIO;
use Moose;
extends 'Bio::Mitomaster::SeqIO';
use Readonly;


=head1 NAME

Bio::Mitomaster::FastaIO - Read and write fasta files.

=cut

our $VERSION = '0.10';



=head1 SYNOPSIS

Provides IO functionality for fasta formatted sequence files.


=head1 ATTRIBUTES



=head1 METHODS

=cut


around 'read_seqs' => sub {#FOLDBEG

    my $orig = shift;
    my $self = shift;
    my $file_name = shift;

    $self->$orig($file_name, \&_parse_seq);

};#FOLDEND


# Each SeqIO subclass must define a _parse_seq subroutine, which is passed to read_seqs.
sub _parse_seq {#FOLDBEG

    # A callback used by the read_seqs method.

    my $seq_text = shift;
    my @seqs;
    Readonly my $ID => qr{[^\n]+};
    Readonly my $SEQ => qr{[^>]+};

    # Define an RE for validating a sequence.  By the time we match this the 
    # seq will have been cleansed of spaces, newlines, and placeholders, and
    # will have been made upper-case.
    Readonly my $VALIDSEQ => qr{^[ACGTUNRYWSMKHBVDK]+$};


    while ($seq_text =~ m{>($ID)\n($SEQ)}gxm) {

        my ($identifier, $sequence) = ($1, $2);

        # remove spaces, line breaks, and placeholders from the seq
        $sequence =~ s/\s//g;
        $sequence =~ s/X//g;  # X is sometimes used for the placeholder at 3107
        $sequence =~ s/://g;
        $sequence =~ s/\-//g;

        # make it uppercase
        $sequence =~ tr/acgtunrywsmkhbvdx/ACGTUNRYWSMKHBVDX/;

        # check that the sequence is valid
        unless ($sequence =~ /$VALIDSEQ/) {
            confess "Sequence: $sequence is not valid\n";
        }

        # Extract information in the identifier line which can contain any values 
        # (separated by a |).  The name should be the first value after which we
        # assume that the values occur in key-value pairs.
        my ($seq_name, %seq_info) = split('\|', $identifier);
        $seq_info{name} = $seq_name;
        $seq_info{seq} = $sequence;
        
        push @seqs, \%seq_info;
    }

    return \@seqs;

}#FOLDEND



=head1 AUTHOR

Marty Brandon, C<< <marty.brandon at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-bio-mitomaster-mitoseq at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Bio-Mitomaster-FastaIO>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Bio::Mitomaster::FastaIO


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Bio-Mitomaster-FastaIO>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Bio-Mitomaster-FastaIO>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/Bio-Mitomaster-FastaIO>

=item * Search CPAN

L<http://search.cpan.org/dist/Bio-Mitomaster-FastaIO/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 COPYRIGHT & LICENSE

Copyright 2009 Marty Brandon, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.


=cut

no Moose;
__PACKAGE__->meta->make_immutable();
1; # End of Bio::Mitomaster::FastaIO

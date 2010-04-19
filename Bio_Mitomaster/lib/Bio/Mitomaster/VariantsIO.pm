package Bio::Mitomaster::VariantsIO;
use Moose;
extends 'Bio::Mitomaster::SeqIO';
use Readonly;


=head1 NAME

Bio::Mitomaster::VariantsIO - Read and write fasta files.

=cut

our $VERSION = '0.10';



=head1 SYNOPSIS

Provides IO functionality for fasta formatted sequence files.


=head1 ATTRIBUTES



=head1 METHODS

These are all basically helper methods for reading/writing files.  It is expected that they would be overridden or modified by the subclasses.

=cut


around 'read_seqs' => sub {#FOLDBEG
    my $orig = shift;
    my $self = shift;
    my $file_name = shift;

    $self->$orig($file_name, \&_parse_seq);
};#FOLDEND


sub _parse_seq {#FOLDBEG

    # A callback used by the read_seqs method.

    my $seq_text = shift;
    my @seqs;
    Readonly my $ID => qr{[^\n]+};
    Readonly my $SEQ => qr{[^>]+};


    # Define REs for the various types of variants that might be found
	Readonly my $POSITION => qr/\d{1,6}/;
	Readonly my $DECPOSITION => qr/\d{1,6}\.\d{1,6}/;
	Readonly my $POLY => qr{[ACGTUNRYWSMKHBVDK]}i;
	Readonly my $INS => qr{i[ACGTUNRYWSMKHBVDK]+}i;
	Readonly my $DEL => qr{\-+\d\{0,6\}};
	Readonly my $QUANTIFIER => qr/\d{1,6}/;
	Readonly my $VARIANT => qr{($POSITION|$DECPOSITION) ($POLY|$INS|$DEL) $QUANTIFIER?}x;


    # Loop through the seqs
    while ($seq_text =~ m{>($ID)\n($SEQ)}gxm) {

        my ($identifier, $variants) = ($1, $2);

        # Extract information in the identifier line which can contain any values 
        # (separated by a |).  The name should be the first value after which we
        # assume that the values occur in key-value pairs.
        my ($seq_name, %seq_info) = split('\|', $identifier);
        $seq_info{name} = $seq_name;


        # Loop through the variants
        my %variants;
        while ($variants =~ m{($VARIANT)\n}gxm) {

            my $variant = $1;
            $variant =~ s/\s//g;
            
            if ($variant =~ m{($POSITION)($POLY)}) {
                # a polymorphism
                my ($p, $v) = ($1, $2);
                $v =~ tr/[a-z]/[A-Z]/;
                $variants{$p} = $v;
            }
            elsif ($variant =~ m{($DECPOSITION)($INS)}) {
                # an insertion
                my ($p, $v) = ($1, $2);
                $v = substr($v, 1);
                $v =~ tr/[a-z]/[A-Z]/;
                $variants{$p} = $v;
            }
            elsif ($variant =~ m{($POSITION)-($QUANTIFIER)}) {
                # a deletion expressed with a quantifier
                my ($p, $q) = ($1, $2);
                $variants{$p} = '-' x $q;
            }
            elsif ($variant =~ m{($POSITION)(-+)}) {
                # a deletion expressed literally
                my ($p, $d) = ($1, $2);
                $variants{$p} = $d;
            }
            else {
                warn "variant: $variant is not expressed correctly - skipping\n";
            }

        }

        $seq_info{variants} = \%variants;
        push @seqs, \%seq_info;
    }

    return \@seqs;

}#FOLDEND



=head1 AUTHOR

Marty Brandon, C<< <mbrandon at uci.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-bio-mitomaster-mitoseq at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Bio-Mitomaster-VariantsIO>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Bio::Mitomaster::VariantsIO


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Bio-Mitomaster-VariantsIO>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Bio-Mitomaster-VariantsIO>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/Bio-Mitomaster-VariantsIO>

=item * Search CPAN

L<http://search.cpan.org/dist/Bio-Mitomaster-VariantsIO/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 COPYRIGHT & LICENSE

Copyright 2009 Marty Brandon, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.


=cut

no Moose;
__PACKAGE__->meta->make_immutable();
1; # End of Bio::Mitomaster::VariantsIO

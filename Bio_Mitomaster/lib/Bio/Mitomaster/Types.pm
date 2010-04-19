package Bio::Mitomaster::Types;
use Moose::Util::TypeConstraints;


=head1 NAME

Bio::Mitomaster::Types - Just a package of type constraints.  Not useful by itself

=cut

our $VERSION = '0.10';



=head1 SYNOPSIS

Use this for custom type constraint definitions on the attributes in other packages.

 use Bio::Mitomaster::Types;

 has 'position' => (
    is => 'rw',
    isa => 'PositiveInt',
 );

=cut




=head1 TYPES


=head2 PositiveInt

Use this for positive integers.

=cut

subtype 'PositiveInt'
    => as 'Int',
    => where { $_ > 0 },
    => message { 'Seq start: $_ must be greater than zero' };




=head1 AUTHOR

Marty Brandon, C<< <mbrandon at uci.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-bio-mitomaster-mitoseq at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Bio-Mitomaster-Types>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Bio::Mitomaster::Types


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Bio-Mitomaster-Types>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Bio-Mitomaster-Types>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/Bio-Mitomaster-Types>

=item * Search CPAN

L<http://search.cpan.org/dist/Bio-Mitomaster-Types/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 COPYRIGHT & LICENSE

Copyright 2009 Marty Brandon, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.


=cut

1; # End of Bio::Mitomaster::Types

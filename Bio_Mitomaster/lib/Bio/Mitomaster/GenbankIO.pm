package Bio::Mitomaster::GenbankIO;
use Moose;
extends 'Bio::Mitomaster::SeqIO';
use Readonly;


=head1 NAME

Bio::Mitomaster::GenbankIO - Read and write fasta files.

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
sub _parse_seq{#FOLDBEG
    my $seq_text = shift;

    #declare some syntactical blocks to match
    Readonly my $LOCUS => qr{^\s*LOCUS.+?}xms;
    Readonly my $DEFINITION => qr{^\s*DEFINITION.+?}xms;
    Readonly my $ACCESSION => qr{^\s*ACCESSION.+?}xms;
    Readonly my $VERSION => qr{^\s*VERSION.+?}xms;
    Readonly my $DBLINK => qr{^\s*DBLINK.+?}xms;
    Readonly my $KEYWORDS => qr{^\s*KEYWORDS.+?}xms;
    Readonly my $SOURCE => qr{^\s*SOURCE.+?}xms;
    Readonly my $REFERENCE => qr{^\s*REFERENCE.+?}xms;
    Readonly my $COMMENT => qr{^\s*COMMENT.+?}xms;
    Readonly my $FEATURES => qr{^\s*FEATURES.+?}xms;

    # Genbank sequences have the full sequence located between ORIGIN and //.
    # We accept all of the IUPAC characters in upper or lower case, ':'
    # which is used by Sequencher for spaces in an alignment, '-' which
    # is also sometimes used for a space, and also 'X'
	# which is used as a placeholder for position 3107.    
    Readonly my $ORIGIN => qr{ORIGIN[ACGTUNRYWSMKHBVDKX:\-\s\d\n]+//}ixms;
    
    # validate the file
    $seq_text =~ m{
                ($LOCUS)
                ($DEFINITION)
                ($ACCESSION)
                ($VERSION)
                ($DBLINK)
                ($KEYWORDS)
                ($SOURCE)
                ($REFERENCE)
                ($COMMENT)?
                ($FEATURES)
                ($ORIGIN)}xs || return 0;
                    
    my ($locus,     $definition,  $accession, $version, $dblink,
        $keywords,  $source,      $reference, $comment, 
        $features,  $origin) = ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11);      

    # Extract and clean the values
    my %seq_info;

    # extract the locus
    $locus =~ m{LOCUS\s*(.+)}xms;
    chomp ($seq_info{locus} = $1);

    # extract the definition
    $definition =~ m{DEFINITION\s*(.+)}xms;
    chomp ($seq_info{definition} = $1);

    # extract the accession
    $accession =~ m{ACCESSION\s*(.+)}xms;
    chomp ($seq_info{accession} = $1);

    # extract the version and gi
    $version =~ m{VERSION\s*(.+)}xms;
    chomp ($seq_info{version} = $1);
    $seq_info{version} =~ m{(.+?)\s*GI:?\s*(.+)};
    $seq_info{gi} = $2;

    # extract the dblink
    $accession =~ m{DBLINK\s*(.+)}xms;
    chomp ($seq_info{dblink} = $1);

    # extract the keywords
    $keywords =~ m{KEYWORDS\s*(.+)}xms;
    chomp ($seq_info{keywords} = $1);

    # extract the source
    $source =~ s/\n//g;
    $source =~ m{SOURCE\s*(.+)}xms;
    $seq_info{source} = $1;
    $seq_info{source} =~ m{(.+?)\s*ORGANISM\s*(.+)};
    $seq_info{source} = $1 if $1;
    $seq_info{organism} = $2 if $2;
    $seq_info{organism} =~ s/\s+/ /g;

    # extract the sequence
    $origin =~ m{ORIGIN(.+)//}xms;
    $origin = $1;
    $origin =~ s/\d//g; #remove digits
    $origin =~ s/\s//g; #remove spaces
    $origin =~ s/X//g;  #get rid of placeholder
    $origin =~ s/://g;  #get rid of gaps
    $origin =~ s/\-//g; #get rid of gaps   
    $origin =~ tr/acgtunrywsmkhbvdx/ACGTUNRYWSMKHBVDX/;  #convert to uppercase
    $seq_info{seq} = $origin;

    # extract the references
    my @references;
    my @ref_strings = split(/\s*REFERENCE[^\n]+\n\s*/, $reference);
    shift @ref_strings;
    for my $ref_string (@ref_strings){
    
        # define the sub-fields used for genbank REFERENCE entries
        # these may change in the future, so just adjust the list
        Readonly my $REF_FIELDS => qr{AUTHORS|TITLE|JOURNAL|MEDLINE|PUBMED}x;
        my @field_entries = $ref_string   
           =~ (m{\s*((?:$REF_FIELDS) .+? (?=(?:$REF_FIELDS|\Z)))}gxms);
    
        # seperate the field from the value
        my %reflisting;
        map { 
            my ($field, $value) = (m{\s*(\S+)\s+(.+)}xms);
            $value =~ s/\n//g; 
            $reflisting{$field} = $value; } @field_entries;

        # references are stored as an array of hash references
        push (@references, \%reflisting);
    }
    $seq_info{references} = \@references if @references;

    return \%seq_info;

}#FOLDEND


=head1 AUTHOR

Marty Brandon, C<< <marty.brandon at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-bio-mitomaster-mitoseq at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Bio-Mitomaster-GenbankIO>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Bio::Mitomaster::GenbankIO


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Bio-Mitomaster-GenbankIO>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Bio-Mitomaster-GenbankIO>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/Bio-Mitomaster-GenbankIO>

=item * Search CPAN

L<http://search.cpan.org/dist/Bio-Mitomaster-GenbankIO/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 COPYRIGHT & LICENSE

Copyright 2009 Marty Brandon, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.


=cut

no Moose;
__PACKAGE__->meta->make_immutable();
1; # End of Bio::Mitomaster::GenbankIO

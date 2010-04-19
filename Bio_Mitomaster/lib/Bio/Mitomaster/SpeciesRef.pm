package Bio::Mitomaster::SpeciesRef;
use Moose;
with 'Bio::Mitomaster::SeqManipRole';

=head1 NAME

Bio::Mitomaster::SpeciesRef - Reference object for species and ref sequence meta-data.

=cut


our $VERSION = '0.10';


=head1 SYNOPSIS

  use Bio::Mitomaster::SpeciesRef;
  $ref = Bio::Mitomaster::SpeciesRef->new(species=>'human', reference=>'rCRS');
  $ref = Bio::Mitomaster::SpeciesRef->new();  #same thing

  $amino_acid = $ref->codon_code('TTC');
  # $amino_acid is now F
  $amino_acid = $ref->codon_code('UUC');  
  # Same thing but with Uracil

  print "The mtDNA reference seq is: ", $ref->ref_seq();
  $transcript = $ref->transcript($locus_id);  # Get the ND1 transcript
  $translation = $ref->translation($locus_id); # Get the ND1 translation
  $locus_info = $ref->locus(5);
  print $locus_info, ' begins at ', $locus_info->{start};  # Prints 'MTTL1 begins at 3230'

=cut




=head1 ATTRIBUTES

=cut

=head2 dna_wrapping
=head2 rna_wrapping
=head2 aa_wrapping

These attributes are used to indicate if position values wrap, i.e. the last moiety is a neighbor to the first.  Bio::Mitomaster was designed with support for circular molecules and setting these flag will cause the ending position to simply wrap back on the beginning.  Note that the molecule may still be linear.  Such would be the case if we had a partial human mtDNA sequences that began in the middle of the sequence and extended a little past the origin of replication.  That situation would warrant setting the wrapping value to 1, so that the seq would be treated as one continuous molecule.  Generally, this value is read from the species meta-data and we don't have to manipulate it directly.  The dna_wrapping attribute has an alias accessor called wrapping since that is the value that we normally want to manipulate.  The rna_wrapping and aa_wrapping were included for completeness and internal use by the software, but perhaps there is some unimagined scenario where a transcript or translation could be thought of as wrapping.

=cut

has 'dna_wrapping' => (#FOLDBEG
    is => 'rw',
    isa => 'Bool',
    default => '1',
);

has 'rna_wrapping' => (
    is => 'rw',
    isa => 'Bool',
    default => '0',
);

has 'aa_wrapping' => (
    is => 'rw',
    isa => 'Bool',
    default => '0',
);#FOLDEND


=head2 Constructor Attributes

The reference and species are required for instantiation, but if you're studying human seqs, then the default values are probably okay.

=cut



=head3 end

The last moiety position.  This value is contextual.  While it defaults to the last genome ref position, if the method is called by an RNASeq object (RNASeq delegates to this method), it will return the last position of the corresponding ref transcript.  Same goes for AASeq objects and translations. 

=cut

has 'end' => (#FOLDBEG
    is => 'ro',
    isa => 'Int',
    default => 16569,
);

around 'end' => \&_get_end;

sub _get_end {
    # Alters the behavior of the end accessor to return a value 
    # depending on the type of object that invoked the method.

    my ($orig, $self) = @_;
    if ($self =~ /RNA/) {
        return length($self->transcript($self->locus_id));
    }
    elsif ($self =~ /AA/) {
        return length($self->translation($self->locus_id));
    }
    else {
        return $self->$orig;
    }
}#FOLDEND



=head3 reference

The reference sequence for the species.  Defaults to 'rCRS' (the human revised Cambridge sequence).

=cut

has 'reference' => (#FOLDBEG
    is => 'ro',
    isa => 'Str',
    default => 'rCRS',
    required => 1,
);#FOLDEND




=head2 ref_seq
=head2 ref_seq($position)
=head2 ref_seq($start, $end)

A string for the the MITOMAP reference sequence (the revised Cambridge sequence).  Note that the revised Cambridge sequence (rCRS) has a placeholder at position 3107, which is denoted by an 'N'.  Also note that the placeholder base is represented by an 'N'. 


Returns a string representing the ref sequence value.  This is the same sequence provided by the contained SpeciesRef object, but enhanced with the Bio::Mitomaster::SeqManipRole and trimmed according to the start and end values of this particlar sequence object.  When called with indicies outside the declared start and end values of the sequence, an error is returned.

=cut

has 'ref_genome' => (#FOLDBEG
    is => 'ro',
    isa => 'Str',
    lazy_build => 1,
);

sub _build_ref_genome {
    my $self = shift;
    my $r = <<ref_genome;
GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTT
CGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTC
GCAGTATCTGTCTTTGATTCCTGCCTCATCCTATTATTTATCGCACCTACGTTCAATATT
ACAGGCGAACATACTTACTAAAGTGTGTTAATTAATTAATGCTTGTAGGACATAATAATA
ACAATTGAATGTCTGCACAGCCACTTTCCACACAGACATCATAACAAAAAATTTCCACCA
AACCCCCCCTCCCCCGCTTCTGGCCACAGCACTTAAACACATCTCTGCCAAACCCCAAAA
ACAAAGAACCCTAACACCAGCCTAACCAGATTTCAAATTTTATCTTTTGGCGGTATGCAC
TTTTAACAGTCACCCCCCAACTAACACATTATTTTCCCCTCCCACTCCCATACTACTAAT
CTCATCAATACAACCCCCGCCCATCCTACCCAGCACACACACACCGCTGCTAACCCCATA
CCCCGAACCAACCAAACCCCAAAGACACCCCCCACAGTTTATGTAGCTTACCTCCTCAAA
GCAATACACTGAAAATGTTTAGACGGGCTCACATCACCCCATAAACAAATAGGTTTGGTC
CTAGCCTTTCTATTAGCTCTTAGTAAGATTACACATGCAAGCATCCCCGTTCCAGTGAGT
TCACCCTCTAAATCACCACGATCAAAAGGAACAAGCATCAAGCACGCAGCAATGCAGCTC
AAAACGCTTAGCCTAGCCACACCCCCACGGGAAACAGCAGTGATTAACCTTTAGCAATAA
ACGAAAGTTTAACTAAGCTATACTAACCCCAGGGTTGGTCAATTTCGTGCCAGCCACCGC
GGTCACACGATTAACCCAAGTCAATAGAAGCCGGCGTAAAGAGTGTTTTAGATCACCCCC
TCCCCAATAAAGCTAAAACTCACCTGAGTTGTAAAAAACTCCAGTTGACACAAAATAGAC
TACGAAAGTGGCTTTAACATATCTGAACACACAATAGCTAAGACCCAAACTGGGATTAGA
TACCCCACTATGCTTAGCCCTAAACCTCAACAGTTAAATCAACAAAACTGCTCGCCAGAA
CACTACGAGCCACAGCTTAAAACTCAAAGGACCTGGCGGTGCTTCATATCCCTCTAGAGG
AGCCTGTTCTGTAATCGATAAACCCCGATCAACCTCACCACCTCTTGCTCAGCCTATATA
CCGCCATCTTCAGCAAACCCTGATGAAGGCTACAAAGTAAGCGCAAGTACCCACGTAAAG
ACGTTAGGTCAAGGTGTAGCCCATGAGGTGGCAAGAAATGGGCTACATTTTCTACCCCAG
AAAACTACGATAGCCCTTATGAAACTTAAGGGTCGAAGGTGGATTTAGCAGTAAACTAAG
AGTAGAGTGCTTAGTTGAACAGGGCCCTGAAGCGCGTACACACCGCCCGTCACCCTCCTC
AAGTATACTTCAAAGGACATTTAACTAAAACCCCTACGCATTTATATAGAGGAGACAAGT
CGTAACATGGTAAGTGTACTGGAAAGTGCACTTGGACGAACCAGAGTGTAGCTTAACACA
AAGCACCCAACTTACACTTAGGAGATTTCAACTTAACTTGACCGCTCTGAGCTAAACCTA
GCCCCAAACCCACTCCACCTTACTACCAGACAACCTTAGCCAAACCATTTACCCAAATAA
AGTATAGGCGATAGAAATTGAAACCTGGCGCAATAGATATAGTACCGCAAGGGAAAGATG
AAAAATTATAACCAAGCATAATATAGCAAGGACTAACCCCTATACCTTCTGCATAATGAA
TTAACTAGAAATAACTTTGCAAGGAGAGCCAAAGCTAAGACCCCCGAAACCAGACGAGCT
ACCTAAGAACAGCTAAAAGAGCACACCCGTCTATGTAGCAAAATAGTGGGAAGATTTATA
GGTAGAGGCGACAAACCTACCGAGCCTGGTGATAGCTGGTTGTCCAAGATAGAATCTTAG
TTCAACTTTAAATTTGCCCACAGAACCCTCTAAATCCCCTTGTAAATTTAACTGTTAGTC
CAAAGAGGAACAGCTCTTTGGACACTAGGAAAAAACCTTGTAGAGAGAGTAAAAAATTTA
ACACCCATAGTAGGCCTAAAAGCAGCCACCAATTAAGAAAGCGTTCAAGCTCAACACCCA
CTACCTAAAAAATCCCAAACATATAACTGAACTCCTCACACCCAATTGGACCAATCTATC
ACCCTATAGAAGAACTAATGTTAGTATAAGTAACATGAAAACATTCTCCTCCGCATAAGC
CTGCGTCAGATTAAAACACTGAACTGACAATTAACAGCCCAATATCTACAATCAACCAAC
AAGTCATTATTACCCTCACTGTCAACCCAACACAGGCATGCTCATAAGGAAAGGTTAAAA
AAAGTAAAAGGAACTCGGCAAATCTTACCCCGCCTGTTTACCAAAAACATCACCTCTAGC
ATCACCAGTATTAGAGGCACCGCCTGCCCAGTGACACATGTTTAACGGCCGCGGTACCCT
AACCGTGCAAAGGTAGCATAATCACTTGTTCCTTAAATAGGGACCTGTATGAATGGCTCC
ACGAGGGTTCAGCTGTCTCTTACTTTTAACCAGTGAAATTGACCTGCCCGTGAAGAGGCG
GGCATAACACAGCAAGACGAGAAGACCCTATGGAGCTTTAATTTATTAATGCAAACAGTA
CCTAACAAACCCACAGGTCCTAAACTACCAAACCTGCATTAAAAATTTCGGTTGGGGCGA
CCTCGGAGCAGAACCCAACCTCCGAGCAGTACATGCTAAGACTTCACCAGTCAAAGCGAA
CTACTATACTCAATTGATCCAATAACTTGACCAACGGAACAAGTTACCCTAGGGATAACA
GCGCAATCCTATTCTAGAGTCCATATCAACAATAGGGTTTACGACCTCGATGTTGGATCA
GGACATCCCGATGGTGCAGCCGCTATTAAAGGTTCGTTTGTTCAACGATTAAAGTCCTAC
GTGATCTGAGTTCAGACCGGAGTAATCCAGGTCGGTTTCTATCTACXTTCAAATTCCTCC
CTGTACGAAAGGACAAGAGAAATAAGGCCTACTTCACAAAGCGCCTTCCCCCGTAAATGA
TATCATCTCAACTTAGTATTATACCCACACCCACCCAAGAACAGGGTTTGTTAAGATGGC
AGAGCCCGGTAATCGCATAAAACTTAAAACTTTACAGTCAGAGGTTCAATTCCTCTTCTT
AACAACATACCCATGGCCAACCTCCTACTCCTCATTGTACCCATTCTAATCGCAATGGCA
TTCCTAATGCTTACCGAACGAAAAATTCTAGGCTATATACAACTACGCAAAGGCCCCAAC
GTTGTAGGCCCCTACGGGCTACTACAACCCTTCGCTGACGCCATAAAACTCTTCACCAAA
GAGCCCCTAAAACCCGCCACATCTACCATCACCCTCTACATCACCGCCCCGACCTTAGCT
CTCACCATCGCTCTTCTACTATGAACCCCCCTCCCCATACCCAACCCCCTGGTCAACCTC
AACCTAGGCCTCCTATTTATTCTAGCCACCTCTAGCCTAGCCGTTTACTCAATCCTCTGA
TCAGGGTGAGCATCAAACTCAAACTACGCCCTGATCGGCGCACTGCGAGCAGTAGCCCAA
ACAATCTCATATGAAGTCACCCTAGCCATCATTCTACTATCAACATTACTAATAAGTGGC
TCCTTTAACCTCTCCACCCTTATCACAACACAAGAACACCTCTGATTACTCCTGCCATCA
TGACCCTTGGCCATAATATGATTTATCTCCACACTAGCAGAGACCAACCGAACCCCCTTC
GACCTTGCCGAAGGGGAGTCCGAACTAGTCTCAGGCTTCAACATCGAATACGCCGCAGGC
CCCTTCGCCCTATTCTTCATAGCCGAATACACAAACATTATTATAATAAACACCCTCACC
ACTACAATCTTCCTAGGAACAACATATGACGCACTCTCCCCTGAACTCTACACAACATAT
TTTGTCACCAAGACCCTACTTCTAACCTCCCTGTTCTTATGAATTCGAACAGCATACCCC
CGATTCCGCTACGACCAACTCATACACCTCCTATGAAAAAACTTCCTACCACTCACCCTA
GCATTACTTATATGATATGTCTCCATACCCATTACAATCTCCAGCATTCCCCCTCAAACC
TAAGAAATATGTCTGATAAAAGAGTTACTTTGATAGAGTAAATAATAGGAGCTTAAACCC
CCTTATTTCTAGGACTATGAGAATCGAACCCATCCCTGAGAATCCAAAATTCTCCGTGCC
ACCTATCACACCCCATCCTAAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCG
AAAATGTTGGTTATACCCTTCCCGTACTAATTAATCCCCTGGCCCAACCCGTCATCTACT
CTACCATCTTTGCAGGCACACTCATCACAGCGCTAAGCTCGCACTGATTTTTTACCTGAG
TAGGCCTAGAAATAAACATGCTAGCTTTTATTCCAGTTCTAACCAAAAAAATAAACCCTC
GTTCCACAGAAGCTGCCATCAAGTATTTCCTCACGCAAGCAACCGCATCCATAATCCTTC
TAATAGCTATCCTCTTCAACAATATACTCTCCGGACAATGAACCATAACCAATACTACCA
ATCAATACTCATCATTAATAATCATAATAGCTATAGCAATAAAACTAGGAATAGCCCCCT
TTCACTTCTGAGTCCCAGAGGTTACCCAAGGCACCCCTCTGACATCCGGCCTGCTTCTTC
TCACATGACAAAAACTAGCCCCCATCTCAATCATATACCAAATCTCTCCCTCACTAAACG
TAAGCCTTCTCCTCACTCTCTCAATCTTATCCATCATAGCAGGCAGTTGAGGTGGATTAA
ACCAAACCCAGCTACGCAAAATCTTAGCATACTCCTCAATTACCCACATAGGATGAATAA
TAGCAGTTCTACCGTACAACCCTAACATAACCATTCTTAATTTAACTATTTATATTATCC
TAACTACTACCGCATTCCTACTACTCAACTTAAACTCCAGCACCACGACCCTACTACTAT
CTCGCACCTGAAACAAGCTAACATGACTAACACCCTTAATTCCATCCACCCTCCTCTCCC
TAGGAGGCCTGCCCCCGCTAACCGGCTTTTTGCCCAAATGGGCCATTATCGAAGAATTCA
CAAAAAACAATAGCCTCATCATCCCCACCATCATAGCCACCATCACCCTCCTTAACCTCT
ACTTCTACCTACGCCTAATCTACTCCACCTCAATCACACTACTCCCCATATCTAACAACG
TAAAAATAAAATGACAGTTTGAACATACAAAACCCACCCCATTCCTCCCCACACTCATCG
CCCTTACCACGCTACTCCTACCTATCTCCCCTTTTATACTAATAATCTTATAGAAATTTA
GGTTAAATACAGACCAAGAGCCTTCAAAGCCCTCAGTAAGTTGCAATACTTAATTTCTGT
AACAGCTAAGGACTGCAAAACCCCACTCTGCATCAACTGAACGCAAATCAGCCACTTTAA
TTAAGCTAAGCCCTTACTAGACCAATGGGACTTAAACCCACAAACACTTAGTTAACAGCT
AAGCACCCTAATCAACTGGCTTCAATCTACTTCTCCCGCCGCCGGGAAAAAAGGCGGGAG
AAGCCCCGGCAGGTTTGAAGCTGCTTCTTCGAATTTGCAATTCAATATGAAAATCACCTC
GGAGCTGGTAAAAAGAGGCCTAACCCCTGTCTTTAGATTTACAGTCCAATGCTTCACTCA
GCCATTTTACCTCACCCCCACTGATGTTCGCCGACCGTTGACTATTCTCTACAAACCACA
AAGACATTGGAACACTATACCTATTATTCGGCGCATGAGCTGGAGTCCTAGGCACAGCTC
TAAGCCTCCTTATTCGAGCCGAGCTGGGCCAGCCAGGCAACCTTCTAGGTAACGACCACA
TCTACAACGTTATCGTCACAGCCCATGCATTTGTAATAATCTTCTTCATAGTAATACCCA
TCATAATCGGAGGCTTTGGCAACTGACTAGTTCCCCTAATAATCGGTGCCCCCGATATGG
CGTTTCCCCGCATAAACAACATAAGCTTCTGACTCTTACCTCCCTCTCTCCTACTCCTGC
TCGCATCTGCTATAGTGGAGGCCGGAGCAGGAACAGGTTGAACAGTCTACCCTCCCTTAG
CAGGGAACTACTCCCACCCTGGAGCCTCCGTAGACCTAACCATCTTCTCCTTACACCTAG
CAGGTGTCTCCTCTATCTTAGGGGCCATCAATTTCATCACAACAATTATCAATATAAAAC
CCCCTGCCATAACCCAATACCAAACGCCCCTCTTCGTCTGATCCGTCCTAATCACAGCAG
TCCTACTTCTCCTATCTCTCCCAGTCCTAGCTGCTGGCATCACTATACTACTAACAGACC
GCAACCTCAACACCACCTTCTTCGACCCCGCCGGAGGAGGAGACCCCATTCTATACCAAC
ACCTATTCTGATTTTTCGGTCACCCTGAAGTTTATATTCTTATCCTACCAGGCTTCGGAA
TAATCTCCCATATTGTAACTTACTACTCCGGAAAAAAAGAACCATTTGGATACATAGGTA
TGGTCTGAGCTATGATATCAATTGGCTTCCTAGGGTTTATCGTGTGAGCACACCATATAT
TTACAGTAGGAATAGACGTAGACACACGAGCATATTTCACCTCCGCTACCATAATCATCG
CTATCCCCACCGGCGTCAAAGTATTTAGCTGACTCGCCACACTCCACGGAAGCAATATGA
AATGATCTGCTGCAGTGCTCTGAGCCCTAGGATTCATCTTTCTTTTCACCGTAGGTGGCC
TGACTGGCATTGTATTAGCAAACTCATCACTAGACATCGTACTACACGACACGTACTACG
TTGTAGCCCACTTCCACTATGTCCTATCAATAGGAGCTGTATTTGCCATCATAGGAGGCT
TCATTCACTGATTTCCCCTATTCTCAGGCTACACCCTAGACCAAACCTACGCCAAAATCC
ATTTCACTATCATATTCATCGGCGTAAATCTAACTTTCTTCCCACAACACTTTCTCGGCC
TATCCGGAATGCCCCGACGTTACTCGGACTACCCCGATGCATACACCACATGAAACATCC
TATCATCTGTAGGCTCATTCATTTCTCTAACAGCAGTAATATTAATAATTTTCATGATTT
GAGAAGCCTTCGCTTCGAAGCGAAAAGTCCTAATAGTAGAAGAACCCTCCATAAACCTGG
AGTGACTATATGGATGCCCCCCACCCTACCACACATTCGAAGAACCCGTATACATAAAAT
CTAGACAAAAAAGGAAGGAATCGAACCCCCCAAAGCTGGTTTCAAGCCAACCCCATGGCC
TCCATGACTTTTTCAAAAAGGTATTAGAAAAACCATTTCATAACTTTGTCAAAGTTAAAT
TATAGGCTAAATCCTATATATCTTAATGGCACATGCAGCGCAAGTAGGTCTACAAGACGC
TACTTCCCCTATCATAGAAGAGCTTATCACCTTTCATGATCACGCCCTCATAATCATTTT
CCTTATCTGCTTCCTAGTCCTGTATGCCCTTTTCCTAACACTCACAACAAAACTAACTAA
TACTAACATCTCAGACGCTCAGGAAATAGAAACCGTCTGAACTATCCTGCCCGCCATCAT
CCTAGTCCTCATCGCCCTCCCATCCCTACGCATCCTTTACATAACAGACGAGGTCAACGA
TCCCTCCCTTACCATCAAATCAATTGGCCACCAATGGTACTGAACCTACGAGTACACCGA
CTACGGCGGACTAATCTTCAACTCCTACATACTTCCCCCATTATTCCTAGAACCAGGCGA
CCTGCGACTCCTTGACGTTGACAATCGAGTAGTACTCCCGATTGAAGCCCCCATTCGTAT
AATAATTACATCACAAGACGTCTTGCACTCATGAGCTGTCCCCACATTAGGCTTAAAAAC
AGATGCAATTCCCGGACGTCTAAACCAAACCACTTTCACCGCTACACGACCGGGGGTATA
CTACGGTCAATGCTCTGAAATCTGTGGAGCAAACCACAGTTTCATGCCCATCGTCCTAGA
ATTAATTCCCCTAAAAATCTTTGAAATAGGGCCCGTATTTACCCTATAGCACCCCCTCTA
CCCCCTCTAGAGCCCACTGTAAAGCTAACTTAGCATTAACCTTTTAAGTTAAAGATTAAG
AGAACCAACACCTCTTTACAGTGAAATGCCCCAACTAAATACTACCGTATGGCCCACCAT
AATTACCCCCATACTCCTTACACTATTCCTCATCACCCAACTAAAAATATTAAACACAAA
CTACCACCTACCTCCCTCACCAAAGCCCATAAAAATAAAAAATTATAACAAACCCTGAGA
ACCAAAATGAACGAAAATCTGTTCGCTTCATTCATTGCCCCCACAATCCTAGGCCTACCC
GCCGCAGTACTGATCATTCTATTTCCCCCTCTATTGATCCCCACCTCCAAATATCTCATC
AACAACCGACTAATCACCACCCAACAATGACTAATCAAACTAACCTCAAAACAAATGATA
ACCATACACAACACTAAAGGACGAACCTGATCTCTTATACTAGTATCCTTAATCATTTTT
ATTGCCACAACTAACCTCCTCGGACTCCTGCCTCACTCATTTACACCAACCACCCAACTA
TCTATAAACCTAGCCATGGCCATCCCCTTATGAGCGGGCACAGTGATTATAGGCTTTCGC
TCTAAGATTAAAAATGCCCTAGCCCACTTCTTACCACAAGGCACACCTACACCCCTTATC
CCCATACTAGTTATTATCGAAACCATCAGCCTACTCATTCAACCAATAGCCCTGGCCGTA
CGCCTAACCGCTAACATTACTGCAGGCCACCTACTCATGCACCTAATTGGAAGCGCCACC
CTAGCAATATCAACCATTAACCTTCCCTCTACACTTATCATCTTCACAATTCTAATTCTA
CTGACTATCCTAGAAATCGCTGTCGCCTTAATCCAAGCCTACGTTTTCACACTTCTAGTA
AGCCTCTACCTGCACGACAACACATAATGACCCACCAATCACATGCCTATCATATAGTAA
AACCCAGCCCATGACCCCTAACAGGGGCCCTCTCAGCCCTCCTAATGACCTCCGGCCTAG
CCATGTGATTTCACTTCCACTCCATAACGCTCCTCATACTAGGCCTACTAACCAACACAC
TAACCATATACCAATGATGGCGCGATGTAACACGAGAAAGCACATACCAAGGCCACCACA
CACCACCTGTCCAAAAAGGCCTTCGATACGGGATAATCCTATTTATTACCTCAGAAGTTT
TTTTCTTCGCAGGATTTTTCTGAGCCTTTTACCACTCCAGCCTAGCCCCTACCCCCCAAT
TAGGAGGGCACTGGCCCCCAACAGGCATCACCCCGCTAAATCCCCTAGAAGTCCCACTCC
TAAACACATCCGTATTACTCGCATCAGGAGTATCAATCACCTGAGCTCACCATAGTCTAA
TAGAAAACAACCGAAACCAAATAATTCAAGCACTGCTTATTACAATTTTACTGGGTCTCT
ATTTTACCCTCCTACAAGCCTCAGAGTACTTCGAGTCTCCCTTCACCATTTCCGACGGCA
TCTACGGCTCAACATTTTTTGTAGCCACAGGCTTCCACGGACTTCACGTCATTATTGGCT
CAACTTTCCTCACTATCTGCTTCATCCGCCAACTAATATTTCACTTTACATCCAAACATC
ACTTTGGCTTCGAAGCCGCCGCCTGATACTGGCATTTTGTAGATGTGGTTTGACTATTTC
TGTATGTCTCCATCTATTGATGAGGGTCTTACTCTTTTAGTATAAATAGTACCGTTAACT
TCCAATTAACTAGTTTTGACAACATTCAAAAAAGAGTAATAAACTTCGCCTTAATTTTAA
TAATCAACACCCTCCTAGCCTTACTACTAATAATTATTACATTTTGACTACCACAACTCA
ACGGCTACATAGAAAAATCCACCCCTTACGAGTGCGGCTTCGACCCTATATCCCCCGCCC
GCGTCCCTTTCTCCATAAAATTCTTCTTAGTAGCTATTACCTTCTTATTATTTGATCTAG
AAATTGCCCTCCTTTTACCCCTACCATGAGCCCTACAAACAACTAACCTGCCACTAATAG
TTATGTCATCCCTCTTATTAATCATCATCCTAGCCCTAAGTCTGGCCTATGAGTGACTAC
AAAAAGGATTAGACTGAACCGAATTGGTATATAGTTTAAACAAAACGAATGATTTCGACT
CATTAAATTATGATAATCATATTTACCAAATGCCCCTCATTTACATAAATATTATACTAG
CATTTACCATCTCACTTCTAGGAATACTAGTATATCGCTCACACCTCATATCCTCCCTAC
TATGCCTAGAAGGAATAATACTATCGCTGTTCATTATAGCTACTCTCATAACCCTCAACA
CCCACTCCCTCTTAGCCAATATTGTGCCTATTGCCATACTAGTCTTTGCCGCCTGCGAAG
CAGCGGTGGGCCTAGCCCTACTAGTCTCAATCTCCAACACATATGGCCTAGACTACGTAC
ATAACCTAAACCTACTCCAATGCTAAAACTAATCGTCCCAACAATTATATTACTACCACT
GACATGACTTTCCAAAAAACACATAATTTGAATCAACACAACCACCCACAGCCTAATTAT
TAGCATCATCCCTCTACTATTTTTTAACCAAATCAACAACAACCTATTTAGCTGTTCCCC
AACCTTTTCCTCCGACCCCCTAACAACCCCCCTCCTAATACTAACTACCTGACTCCTACC
CCTCACAATCATGGCAAGCCAACGCCACTTATCCAGTGAACCACTATCACGAAAAAAACT
CTACCTCTCTATACTAATCTCCCTACAAATCTCCTTAATTATAACATTCACAGCCACAGA
ACTAATCATATTTTATATCTTCTTCGAAACCACACTTATCCCCACCTTGGCTATCATCAC
CCGATGAGGCAACCAGCCAGAACGCCTGAACGCAGGCACATACTTCCTATTCTACACCCT
AGTAGGCTCCCTTCCCCTACTCATCGCACTAATTTACACTCACAACACCCTAGGCTCACT
AAACATTCTACTACTCACTCTCACTGCCCAAGAACTATCAAACTCCTGAGCCAACAACTT
AATATGACTAGCTTACACAATAGCTTTTATAGTAAAGATACCTCTTTACGGACTCCACTT
ATGACTCCCTAAAGCCCATGTCGAAGCCCCCATCGCTGGGTCAATAGTACTTGCCGCAGT
ACTCTTAAAACTAGGCGGCTATGGTATAATACGCCTCACACTCATTCTCAACCCCCTGAC
AAAACACATAGCCTACCCCTTCCTTGTACTATCCCTATGAGGCATAATTATAACAAGCTC
CATCTGCCTACGACAAACAGACCTAAAATCGCTCATTGCATACTCTTCAATCAGCCACAT
AGCCCTCGTAGTAACAGCCATTCTCATCCAAACCCCCTGAAGCTTCACCGGCGCAGTCAT
TCTCATAATCGCCCACGGGCTTACATCCTCATTACTATTCTGCCTAGCAAACTCAAACTA
CGAACGCACTCACAGTCGCATCATAATCCTCTCTCAAGGACTTCAAACTCTACTCCCACT
AATAGCTTTTTGATGACTTCTAGCAAGCCTCGCTAACCTCGCCTTACCCCCCACTATTAA
CCTACTGGGAGAACTCTCTGTGCTAGTAACCACGTTCTCCTGATCAAATATCACTCTCCT
ACTTACAGGACTCAACATACTAGTCACAGCCCTATACTCCCTCTACATATTTACCACAAC
ACAATGGGGCTCACTCACCCACCACATTAACAACATAAAACCCTCATTCACACGAGAAAA
CACCCTCATGTTCATACACCTATCCCCCATTCTCCTCCTATCCCTCAACCCCGACATCAT
TACCGGGTTTTCCTCTTGTAAATATAGTTTAACCAAAACATCAGATTGTGAATCTGACAA
CAGAGGCTTACGACCCCTTATTTACCGAGAAAGCTCACAAGAACTGCTAACTCATGCCCC
CATGTCTAACAACATGGCTTTCTCAACTTTTAAAGGATAACAGCTATCCATTGGTCTTAG
GCCCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATAACCATGCACACTACTATAACC
ACCCTAACCCTGACTTCCCTAATTCCCCCCATCCTTACCACCCTCGTTAACCCTAACAAA
AAAAACTCATACCCCCATTATGTAAAATCCATTGTCGCATCCACCTTTATTATCAGTCTC
TTCCCCACAACAATATTCATGTGCCTAGACCAAGAAGTTATTATCTCGAACTGACACTGA
GCCACAACCCAAACAACCCAGCTCTCCCTAAGCTTCAAACTAGACTACTTCTCCATAATA
TTCATCCCTGTAGCATTGTTCGTTACATGGTCCATCATAGAATTCTCACTGTGATATATA
AACTCAGACCCAAACATTAATCAGTTCTTCAAATATCTACTCATCTTCCTAATTACCATA
CTAATCTTAGTTACCGCTAACAACCTATTCCAACTGTTCATCGGCTGAGAGGGCGTAGGA
ATTATATCCTTCTTGCTCATCAGTTGATGATACGCCCGAGCAGATGCCAACACAGCAGCC
ATTCAAGCAATCCTATACAACCGTATCGGCGATATCGGTTTCATCCTCGCCTTAGCATGA
TTTATCCTACACTCCAACTCATGAGACCCACAACAAATAGCCCTTCTAAACGCTAATCCA
AGCCTCACCCCACTACTAGGCCTCCTCCTAGCAGCAGCAGGCAAATCAGCCCAATTAGGT
CTCCACCCCTGACTCCCCTCAGCCATAGAAGGCCCCACCCCAGTCTCAGCCCTACTCCAC
TCAAGCACTATAGTTGTAGCAGGAATCTTCTTACTCATCCGCTTCCACCCCCTAGCAGAA
AATAGCCCACTAATCCAAACTCTAACACTATGCTTAGGCGCTATCACCACTCTGTTCGCA
GCAGTCTGCGCCCTTACACAAAATGACATCAAAAAAATCGTAGCCTTCTCCACTTCAAGT
CAACTAGGACTCATAATAGTTACAATCGGCATCAACCAACCACACCTAGCATTCCTGCAC
ATCTGTACCCACGCCTTCTTCAAAGCCATACTATTTATGTGCTCCGGGTCCATCATCCAC
AACCTTAACAATGAACAAGATATTCGAAAAATAGGAGGACTACTCAAAACCATACCTCTC
ACTTCAACCTCCCTCACCATTGGCAGCCTAGCATTAGCAGGAATACCTTTCCTCACAGGT
TTCTACTCCAAAGACCACATCATCGAAACCGCAAACATATCATACACAAACGCCTGAGCC
CTATCTATTACTCTCATCGCTACCTCCCTGACAAGCGCCTATAGCACTCGAATAATTCTT
CTCACCCTAACAGGTCAACCTCGCTTCCCCACCCTTACTAACATTAACGAAAATAACCCC
ACCCTACTAAACCCCATTAAACGCCTGGCAGCCGGAAGCCTATTCGCAGGATTTCTCATT
ACTAACAACATTTCCCCCGCATCCCCCTTCCAAACAACAATCCCCCTCTACCTAAAACTC
ACAGCCCTCGCTGTCACTTTCCTAGGACTTCTAACAGCCCTAGACCTCAACTACCTAACC
AACAAACTTAAAATAAAATCCCCACTATGCACATTTTATTTCTCCAACATACTCGGATTC
TACCCTAGCATCACACACCGCACAATCCCCTATCTAGGCCTTCTTACGAGCCAAAACCTG
CCCCTACTCCTCCTAGACCTAACCTGACTAGAAAAGCTATTACCTAAAACAATTTCACAG
CACCAAATCTCCACCTCCATCATCACCTCAACCCAAAAAGGCATAATTAAACTTTACTTC
CTCTCTTTCTTCTTCCCACTCATCCTAACCCTACTCCTAATCACATAACCTATTCCCCCG
AGCAATCTCAATTACAATATATACACCAACAAACAATGTTCAACCAGTAACTACTACTAA
TCAACGCCCATAATCATACAAAGCCCCCGCACCAATAGGATCCTCCCGAATCAACCCTGA
CCCCTCTCCTTCATAAATTATTCAGCTTCCTACACTATTAAAGTTTACCACAACCACCAC
CCCATCATACTCTTTCACCCACAGCACCAATCCTACCTCCATCGCTAACCCCACTAAAAC
ACTCACCAAGACCTCAACCCCTGACCCCCATGCCTCAGGATACTCCTCAATAGCCATCGC
TGTAGTATATCCAAAGACAACCATCATTCCCCCTAAATAAATTAAAAAAACTATTAAACC
CATATAACCTCCCCCAAAATTCAGAATAATAACACACCCGACCACACCGCTAACAATCAA
TACTAAACCCCCATAAATAGGAGAAGGCTTAGAAGAAAACCCCACAAACCCCATTACTAA
ACCCACACTCAACAGAAACAAAGCATACATCATTATTCTCGCACGGACTACAACCACGAC
CAATGATATGAAAAACCATCGTTGTATTTCAACTACAAGAACACCAATGACCCCAATACG
CAAAACTAACCCCCTAATAAAATTAATTAACCACTCATTCATCGACCTCCCCACCCCATC
CAACATCTCCGCATGATGAAACTTCGGCTCACTCCTTGGCGCCTGCCTGATCCTCCAAAT
CACCACAGGACTATTCCTAGCCATGCACTACTCACCAGACGCCTCAACCGCCTTTTCATC
AATCGCCCACATCACTCGAGACGTAAATTATGGCTGAATCATCCGCTACCTTCACGCCAA
TGGCGCCTCAATATTCTTTATCTGCCTCTTCCTACACATCGGGCGAGGCCTATATTACGG
ATCATTTCTCTACTCAGAAACCTGAAACATCGGCATTATCCTCCTGCTTGCAACTATAGC
AACAGCCTTCATAGGCTATGTCCTCCCGTGAGGCCAAATATCATTCTGAGGGGCCACAGT
AATTACAAACTTACTATCCGCCATCCCATACATTGGGACAGACCTAGTTCAATGAATCTG
AGGAGGCTACTCAGTAGACAGTCCCACCCTCACACGATTCTTTACCTTTCACTTCATCTT
GCCCTTCATTATTGCAGCCCTAGCAACACTCCACCTCCTATTCTTGCACGAAACGGGATC
AAACAACCCCCTAGGAATCACCTCCCATTCCGATAAAATCACCTTCCACCCTTACTACAC
AATCAAAGACGCCCTCGGCTTACTTCTCTTCCTTCTCTCCTTAATGACATTAACACTATT
CTCACCAGACCTCCTAGGCGACCCAGACAATTATACCCTAGCCAACCCCTTAAACACCCC
TCCCCACATCAAGCCCGAATGATATTTCCTATTCGCCTACACAATTCTCCGATCCGTCCC
TAACAAACTAGGAGGCGTCCTTGCCCTATTACTATCCATCCTCATCCTAGCAATAATCCC
CATCCTCCATATATCCAAACAACAAAGCATAATATTTCGCCCACTAAGCCAATCACTTTA
TTGACTCCTAGCCGCAGACCTCCTCATTCTAACCTGAATCGGAGGACAACCAGTAAGCTA
CCCTTTTACCATCATTGGACAAGTAGCATCCGTACTATACTTCACAACAATCCTAATCCT
AATACCAACTATCTCCCTAATTGAAAACAAAATACTCAAATGGGCCTGTCCTTGTAGTAT
AAACTAATACACCAGTCTTGTAAACCGGAGATGAAAACCTTTTTCCAAGGACAAATCAGA
GAAAAAGTCTTTAACTCCACCATTAGCACCCAAAGCTAAGATTCTAATTTAAACTATTCT
CTGTTCTTTCATGGGGAAGCAGATTTGGGTACCACCCAAGTATTGACTCACCCATCAACA
ACCGCTATGTATTTCGTACATTACTGCCAGCCACCATGAATATTGTACGGTACCATAAAT
ACTTGACCACCTGTAGTACATAAAAACCCAATCCACATCAAAACCCCCTCCCCATGCTTA
CAAGCAAGTACAGCAATCAACCCTCAACTATCACACATCAACTGCAACTCCAAAGCCACC
CCTCACCCACTAGGATACCAACAAACCTACCCACCCTTAACAGTACATAGTACATAAAGC
CATTTACCGTACATAGCACATTACAGTCAAATCCCTTCTCGTCCCCATGGATGACCCCCC
TCAGATAGGGGTCCCTTGACCACCATCCTCCGTGAAATCAATATCCCGCACAAGAGTGCT
ACTCTCCTCGCTCCGGGCCCATAACACTTGGGGGTAGCTAAAGTGAACTGTATCCGACAT
CTGGTTCCTACTTCAGGGTCATAAAGCCTAAATAGCCCACACGTTCCCCTTAAATAAGAC
ATCACGATG
ref_genome

    $r =~ s/\n//g;
    # remove the 3107 placeholder position unless we want the revised seq
    # $r = substr($r, 0, 3106) . substr($r, 3107);
    return $r;
}


sub ref_seq {

    my ($self, $s, $e) = @_;

    my ($start, $end) = $self->get_seq_indices($s, $e);
    $self->check_seq_indices($start, $end, 1);

    return $self->sub_seq($self->ref_genome, $start, $end);

}#FOLDEND



=head3 species

The species from which the sequence is derived.  Defaults to 'human'.

=cut

has 'species' => (#FOLDBEG
    is => 'ro',
    isa => 'Str',
    default => 'human',
    required => 1,
);#FOLDEND



=head3 start

The first nucleotide position.  Defaults to one.

=cut

has 'start' => (#FOLDBEG
    is => 'ro',
    isa => 'Int',
    default => 1,
);#FOLDEND


=head2 codon_code

The human mitochondrial translation code that maps codons to their amino acid.  Note that the mitochondrial code has some exceptions to the standard code.  There are key values for both thymine and uracil codes, so you can use either one.

 $r = Bio::Mitomaster::SpeciesRef->new();
 $c = $r->codon_code('TTT');  # $c is now 'F'
 $c = $r->codon_code('UUU');  # same thing
 $c = $r->codon_code();  # $c is a hash reference to a codon lookup table

=cut

has 'codon_code' => (#FOLDBEG
    is => 'ro',
    lazy_build => 1,
);

around 'codon_code' => \&_codon_lookup;

sub _build_codon_code {
    my $self = shift;
    my %code = (
	'TTT'=>'F','TTC'=>'F','TTA'=>'L','TTG'=>'L','CTT'=>'L','CTC'=>'L','CTA'=>'L',
	'CTG'=>'L','ATT'=>'I','ATC'=>'I','ATA'=>'M','ATG'=>'M','GTT'=>'V','GTC'=>'V',
	'GTA'=>'V','GTG'=>'V','TCT'=>'S','TCC'=>'S','TCA'=>'S','TCG'=>'S','CCT'=>'P',
	'CCC'=>'P','CCA'=>'P','CCG'=>'P','ACT'=>'T','ACC'=>'T','ACA'=>'T','ACG'=>'T',
	'GCT'=>'A','GCC'=>'A','GCA'=>'A','GCG'=>'A','TAT'=>'Y','TAC'=>'Y','TAA'=>'TERM',
	'TAG'=>'TERM','CAT'=>'H','CAC'=>'H','CAA'=>'Q','CAG'=>'Q','AAT'=>'N','AAC'=>'N',
	'AAA'=>'K','AAG'=>'K','GAT'=>'D','GAC'=>'D','GAA'=>'E','GAG'=>'E','TGT'=>'C',
	'TGC'=>'C','TGA'=>'W','TGG'=>'W','CGT'=>'R','CGC'=>'R','CGA'=>'R','CGG'=>'R',
	'AGT'=>'S','AGC'=>'S','AGA'=>'TERM','AGG'=>'TERM','GGT'=>'G','GGC'=>'G',
	'GGA'=>'G','GGG'=>'G',
	'UUU'=>'F','UUC'=>'F','UUA'=>'L','UUG'=>'L','CUU'=>'L','CUC'=>'L','CUA'=>'L',
	'CUG'=>'L','AUU'=>'I','AUC'=>'I','AUA'=>'M','AUG'=>'M','GUU'=>'V','GUC'=>'V',
	'GUA'=>'V','GUG'=>'V','UCU'=>'S','UCC'=>'S','UCA'=>'S','UCG'=>'S','CCU'=>'P',
	'CCC'=>'P','CCA'=>'P','CCG'=>'P','ACU'=>'T','ACC'=>'T','ACA'=>'T','ACG'=>'T',
	'GCU'=>'A','GCC'=>'A','GCA'=>'A','GCG'=>'A','UAU'=>'Y','UAC'=>'Y','UAA'=>'TERM',
	'UAG'=>'TERM','CAU'=>'H','CAC'=>'H','CAA'=>'Q','CAG'=>'Q','AAU'=>'N','AAC'=>'N',
	'AAA'=>'K','AAG'=>'K','GAU'=>'D','GAC'=>'D','GAA'=>'E','GAG'=>'E','UGU'=>'C',
	'UGC'=>'C','UGA'=>'W','UGG'=>'W','CGU'=>'R','CGC'=>'R','CGA'=>'R','CGG'=>'R',
	'AGU'=>'S','AGC'=>'S','AGA'=>'TERM','AGG'=>'TERM','GGU'=>'G','GGC'=>'G',
	'GGA'=>'G','GGG'=>'G'
    );

    return \%code;
}

sub _codon_lookup {
    my ($orig, $self, $code) = @_;

    my $c = $self->$orig;
    return $c unless defined $code;

    confess "$code is not a valid codon" unless $c->{$code};

    return $c->{$code};
}#FOLDEND



=head2 locus

The human mitochondrial locus regions.  Contains information for these regions.  Locus regions are indexed by ID# and name.

 $r = Bio::Mitomaster::SpeciesRef->new();
 $t = $r->locus();  # $l is a hashref to a hash that contains all the locus hashrefs indexed by ID#
 $l = $r->locus(16);  # $l is a hashref for information about the Cytochrome oxidase subunit I locus
 $s = $l->{start};  # $s is 5903, the start for position for MTCO1

Each locus hashref contains the keys: name, common_name, start, ending, strand, type, product.

=cut

has 'locus' => (#FOLDBEG
    is => 'ro',
    lazy_build => 1,
);

around 'locus' => \&_locus_lookup;


sub _build_locus {

    my $self = shift;

    my %locus_table = (
	38=>{name=>'MTHV2',common_name=>'HVS2',start=>'57',end=>'372',strand=>'N',type=>'n',product=>'Hypervariable segment 2'},
	39=>{name=>'MTOHR',common_name=>'OH',start=>'110',end=>'441',strand=>'H',type=>'n',product=>'H-strand origin'},
	40=>{name=>'MTCSB1',common_name=>'CSB1',start=>'213',end=>'235',strand=>'N',type=>'n',product=>'Conserved sequence block 1'},
	41=>{name=>'MTTFX',common_name=>'TFX',start=>'233',end=>'260',strand=>'N',type=>'n',product=>'mtTF1 binding site'},
	42=>{name=>'MTTFY',common_name=>'TFY',start=>'276',end=>'303',strand=>'N',type=>'n',product=>'mtTF1 binding site'},
	43=>{name=>'MTCSB2',common_name=>'CSB2',start=>'299',end=>'315',strand=>'N',type=>'n',product=>'Conserved sequence block 2'},
	44=>{name=>'MTHPR',common_name=>'HPR',start=>'317',end=>'321',strand=>'N',type=>'n',product=>'replication primer'},
	45=>{name=>'MTCSB3',common_name=>'CSB3',start=>'346',end=>'363',strand=>'N',type=>'n',product=>'Conserved sequence block 3'},
	46=>{name=>'MTMT4H',common_name=>'mt4H',start=>'371',end=>'379',strand=>'H',type=>'n',product=>'mt4 H-strand control element'},
	47=>{name=>'MTMT3H',common_name=>'mt3H',start=>'384',end=>'391',strand=>'H',type=>'n',product=>'mt3 H-strand control element'},
	48=>{name=>'MTLSP',common_name=>'PL',start=>'392',end=>'445',strand=>'L',type=>'n',product=>'L-strand promoter'},
	49=>{name=>'MTTFL',common_name=>'TFL',start=>'418',end=>'445',strand=>'N',type=>'n',product=>'mtTF1 binding site'},
	50=>{name=>'MTTFH',common_name=>'TFH',start=>'523',end=>'550',strand=>'N',type=>'n',product=>'mtTF1 binding site'},
	51=>{name=>'MTHSP1',common_name=>'PH1',start=>'545',end=>'567',strand=>'H',type=>'n',product=>'Major H-strand promoter'},
	1=>{name=>'MTTF',common_name=>'F',start=>'577',end=>'647',strand=>'H',type=>'t',product=>'tRNA phenylanine'},
	52=>{name=>'MTHSP2',common_name=>'PH2',start=>'645',end=>'645',strand=>'H',type=>'n',product=>'Minor H-strand promoter'},
	2=>{name=>'MTRNR1',common_name=>'12S',start=>'648',end=>'1601',strand=>'H',type=>'r',product=>'12S ribosomal RNA'},
	3=>{name=>'MTTV',common_name=>'V',start=>'1602',end=>'1670',strand=>'H',type=>'t',product=>'tRNA valine'},
	4=>{name=>'MTRNR2',common_name=>'16S',start=>'1671',end=>'3229',strand=>'H',type=>'r',product=>'16S ribosomal RNA'},
	53=>{name=>'MTRNR3',common_name=>'RNR3',start=>'3205',end=>'3229',strand=>'N',type=>'n',product=>'5S-like sequence'},
	54=>{name=>'MTTER',common_name=>'TER',start=>'3228',end=>'3255',strand=>'N',type=>'n',product=>'Transcription terminator'},
	5=>{name=>'MTTL1',common_name=>'L(UUA/G)',start=>'3229',end=>'3303',strand=>'H',type=>'t',product=>'tRNA leucine 1'},
	55=>{name=>'MTNC1',common_name=>'NC1',start=>'3304',end=>'3305',strand=>'N',type=>'n',product=>'non-coding nucleotides between MTTL1 and MTND1'},
	6=>{name=>'MTND1',common_name=>'ND1',start=>'3306',end=>'4261',strand=>'H',type=>'m',product=>'NADH Dehydrogenase subunit 1'},
	7=>{name=>'MTTI',common_name=>'I',start=>'4262',end=>'4330',strand=>'H',type=>'t',product=>'tRNA isoleucine'},
	8=>{name=>'MTTQ',common_name=>'Q',start=>'4328',end=>'4399',strand=>'L',type=>'t',product=>'tRNA glutamine'},
	56=>{name=>'MTNC2',common_name=>'NC2',start=>'4400',end=>'4400',strand=>'N',type=>'n',product=>'non-coding nucleotide between MTTQ and MTTM'},
	9=>{name=>'MTTM',common_name=>'M',start=>'4401',end=>'4468',strand=>'H',type=>'t',product=>'tRNA methionine'},
	10=>{name=>'MTND2',common_name=>'ND2',start=>'4469',end=>'5510',strand=>'H',type=>'m',product=>'NADH dehydrogenase subunit 2'},
	11=>{name=>'MTTW',common_name=>'W',start=>'5511',end=>'5578',strand=>'H',type=>'t',product=>'tRNA tryptophan'},
	57=>{name=>'MTNC3',common_name=>'NC3',start=>'5576',end=>'5585',strand=>'N',type=>'n',product=>'non-coding nucleotides between MTTW and MTTA'},
	12=>{name=>'MTTA',common_name=>'A',start=>'5586',end=>'5654',strand=>'L',type=>'t',product=>'tRNA alanine'},
	58=>{name=>'MTNC4',common_name=>'NC4',start=>'5655',end=>'5655',strand=>'N',type=>'n',product=>'non-coding nucleotide between MTTA and MTTN'},
	13=>{name=>'MTTN',common_name=>'N',start=>'5656',end=>'5728',strand=>'L',type=>'t',product=>'tRNA asparagine'},
	59=>{name=>'MTOLR',common_name=>'OL',start=>'5720',end=>'5797',strand=>'L',type=>'n',product=>'L-strand origin'},
	14=>{name=>'MTTC',common_name=>'C',start=>'5760',end=>'5825',strand=>'L',type=>'t',product=>'tRNA cysteine'},
	15=>{name=>'MTTY',common_name=>'Y',start=>'5825',end=>'5890',strand=>'L',type=>'t',product=>'tRNA tyrosine'},
	60=>{name=>'MTNC5',common_name=>'NC5',start=>'5891',end=>'5902',strand=>'N',type=>'n',product=>'non-coding nucleotides between MTTY and MTCO1'},
	16=>{name=>'MTCO1',common_name=>'COI',start=>'5903',end=>'7444',strand=>'H',type=>'m',product=>'Cytochrome c oxidase subunit I'},
	17=>{name=>'MTTS1',common_name=>'S(UCN)',start=>'7445',end=>'7515',strand=>'L',type=>'t',product=>'tRNA serine 1'},
	61=>{name=>'MTNC6',common_name=>'NC6',start=>'7516',end=>'7516',strand=>'N',type=>'n',product=>'non-coding nucleotide between MTTS1 and MTTD'},
	18=>{name=>'MTTD',common_name=>'D',start=>'7517',end=>'7584',strand=>'H',type=>'t',product=>'tRNA aspartic acid'},
	19=>{name=>'MTCO2',common_name=>'COII',start=>'7585',end=>'8268',strand=>'H',type=>'m',product=>'Cytochrome c oxidase subunit II'},
	62=>{name=>'MTNC7',common_name=>'NC7',start=>'8269',end=>'8293',strand=>'N',type=>'n',product=>'non-coding nucleotides between MTCO2 and MTTK'},
	20=>{name=>'MTTK',common_name=>'K',start=>'8294',end=>'8363',strand=>'H',type=>'t',product=>'tRNA lysine'},
	63=>{name=>'MTNC8',common_name=>'NC8',start=>'8364',end=>'8364',strand=>'N',type=>'n',product=>'non-coding nucleotide between MTTK and MTATP8'},
	21=>{name=>'MTATP8',common_name=>'ATPase8',start=>'8365',end=>'8571',strand=>'H',type=>'m',product=>'ATP synthase F0 subunit 8'},
	22=>{name=>'MTATP6',common_name=>'ATPase6',start=>'8526',end=>'9206',strand=>'H',type=>'m',product=>'ATP synthase F0 subunit 6'},
	23=>{name=>'MTCO3',common_name=>'COIII',start=>'9206',end=>'9989',strand=>'H',type=>'m',product=>'Cytochrome c oxidase subunit III'},
	24=>{name=>'MTTG',common_name=>'G',start=>'9990',end=>'10057',strand=>'H',type=>'t',product=>'tRNA glycine'},
	25=>{name=>'MTND3',common_name=>'ND3',start=>'10058',end=>'10403',strand=>'H',type=>'m',product=>'NADH dehydrogenase subunit 3'},
	26=>{name=>'MTTR',common_name=>'R',start=>'10404',end=>'10468',strand=>'H',type=>'t',product=>'tRNA arginine'},
	27=>{name=>'MTND4L',common_name=>'ND4L',start=>'10469',end=>'10765',strand=>'H',type=>'m',product=>'NADH dehydrogenase subunit 4L'},
	28=>{name=>'MTND4',common_name=>'ND4',start=>'10759',end=>'12136',strand=>'H',type=>'m',product=>'NADH dehydrogenase subunit 4'},
	29=>{name=>'MTTH',common_name=>'H',start=>'12137',end=>'12205',strand=>'H',type=>'t',product=>'tRNA histidine'},
	30=>{name=>'MTTS2',common_name=>'S(AGY)',start=>'12206',end=>'12264',strand=>'H',type=>'t',product=>'tRNA serine2'},
	31=>{name=>'MTTL2',common_name=>'L(CUN)',start=>'12265',end=>'12335',strand=>'H',type=>'t',product=>'tRNA leucine2'},
	32=>{name=>'MTND5',common_name=>'ND5',start=>'12336',end=>'14147',strand=>'H',type=>'m',product=>'NADH dehydrogenase subunit 5'},
	33=>{name=>'MTND6',common_name=>'ND6',start=>'14148',end=>'14672',strand=>'L',type=>'m',product=>'NADH dehydrogenase subunit 6'},
	34=>{name=>'MTTE',common_name=>'E',start=>'14673',end=>'14741',strand=>'L',type=>'t',product=>'tRNA glutamic acid'},
	64=>{name=>'MTNC9',common_name=>'NC9',start=>'14742',end=>'14745',strand=>'N',type=>'n',product=>'non-coding nucleotides between MTTE and MTCYB'},
	35=>{name=>'MTCYB',common_name=>'Cytb',start=>'14746',end=>'15886',strand=>'H',type=>'m',product=>'Cytochrome b'},
	36=>{name=>'MTTT',common_name=>'T',start=>'15887',end=>'15952',strand=>'H',type=>'t',product=>'tRNA threonine'},
	65=>{name=>'MTATT',common_name=>'ATT',start=>'15924',end=>'499',strand=>'N',type=>'n',product=>'membrane attachment site'},
	66=>{name=>'MTNC10',common_name=>'NC10',start=>'15953',end=>'15953',strand=>'N',type=>'n',product=>'non-coding nucleotide between MTTT and MTTP'},
	37=>{name=>'MTTP',common_name=>'P',start=>'15954',end=>'16022',strand=>'L',type=>'t',product=>'tRNA proline'},
	67=>{name=>'MTDLOOP',common_name=>'D-Loop',start=>'16023',end=>'576',strand=>'N',type=>'n',product=>'Displacement loop'},
	68=>{name=>'MTHV1',common_name=>'HVS1',start=>'16023',end=>'16382',strand=>'N',type=>'n',product=>'Hypervariable segment 1'},
	69=>{name=>'MT7SDNA',common_name=>'7S DNA',start=>'16105',end=>'191',strand=>'N',type=>'n',product=>'7S DNA'},
	70=>{name=>'MTTAS',common_name=>'TAS',start=>'16156',end=>'16171',strand=>'N',type=>'n',product=>'termination-associated sequence'},
	71=>{name=>'MTMT5',common_name=>'mt5',start=>'16193',end=>'16207',strand=>'N',type=>'n',product=>'control element'},
	72=>{name=>'MTMT3L',common_name=>'mt3L',start=>'16498',end=>'16505',strand=>'L',type=>'n',product=>'L-strand control element'},
    );

    # adjust position numbers for the rCRS numbering
    for (1..72) {
        $locus_table{$_}->{'start'}++ if $locus_table{$_}->{'start'} > 3106;
        $locus_table{$_}->{'end'}++ if $locus_table{$_}->{'end'} > 3106;
    }

    return \%locus_table;
}

sub _locus_lookup {

    my ($orig, $self, $index) = @_;

    my $l = $self->$orig;
    return $l unless defined $index;

    confess "locus index: $index is not a digit" unless $index =~ /^\d+$/;
    confess "locus index: $index is not valid" unless $l->{$index};

    return $l->{$index};

}#FOLDEND



=head2 protein

The human mitochondrially encoded proteins.  Contains information for these.  Proteins are indexed by locus ID# and name.
 $r = Bio::Mitomaster::SpeciesRef->new();
 $p = $r->protein(16);  # $p is a hashref for information about the Cytochrome oxidase subunit I
 $p = $r->protein('COI');  # same thing
 $w = $p->{predicted_weight};  # $w is 57000, the molecular weight predicted for Cytochrome oxidase subunit I
 $p = $r->protein();  # $l is a hashref to a hash that contains all the protein hashrefs indexed by locus ID#
Each protein hashref contains the keys:  id, predicted_weight, glycine_number, glycine_weight, urea_number, urea_weight, name, npid, and giid.

=cut


has 'protein' => (#FOLDBEG
    is => 'ro',
    lazy_build => 1,
);

around 'protein' => \&_protein_lookup;


sub _build_protein {

    my $self = shift;

    my %protein_table = (
	    32=>{id=>1, predicted_weight=>66600, glycine_number=>1, glycine_weight=>43.5, 
	        urea_number=>1, urea_weight=>51, name=>'NADH dehydrogenase subunit 5', 
	        npid=>'NP_536853.1', giid=>'GI:17981863'}, 
	    16=>{id=>2, predicted_weight=>57000, glycine_number=>2, glycine_weight=>39.5, 
	        urea_number=>23, urea_weight=>45, name=>'cytochrome c oxidase subunit I', 
	        npid=>'NP_536845.1', giid=>'GI:17981855'}, 
	    28=>{id=>3, predicted_weight=>51400, glycine_number=>3, glycine_weight=>36.5, 
	        urea_number=>45, urea_weight=>'39, 36', name=>'NADH dehydrogenase subunit 4', 
	        npid=>'NP_536852.1', giid=>'GI:17981862'}, 
	    35=>{id=>4, predicted_weight=>42700, glycine_number=>6, glycine_weight=>27.5, 
	        urea_number=>'7-10', urea_weight=>29, name=>'cytochrome b', npid=>'NP_536855.1', 
	        giid=>'GI:17981865'}, 
	    10=>{id=>5, predicted_weight=>38900, glycine_number=>4, glycine_weight=>31.5, 
	        urea_number=>11, urea_weight=>25, name=>'NADH dehydrogenase subunit 2', 
	        npid=>'NP_536844.1', giid=>'GI:17981854'}, 
	    6=>{id=>6, predicted_weight=>35600, glycine_number=>5, glycine_weight=>29.5, 
	        urea_number=>12, urea_weight=>24, 
	        name=>'NADH dehydrogenase subunit 1', npid=>'NP_536843.1', giid=>'GI:17981853'}, 
	    23=>{id=>7, predicted_weight=>30000, glycine_number=>8, glycine_weight=>22.5, 
	        urea_number=>'15, 16', urea_weight=>18, name=>'cytochrome c oxidase subunit 3', 
	        npid=>'NP_536849.1', giid=>'GI:17981859'}, 
	    19=>{id=>8, predicted_weight=>25500, glycine_number=>7, glycine_weight=>23.6, 
	        urea_number=>'13, 14', urea_weight=>20, name=>'cytochrome c oxidase subunit 2', 
	        npid=>'NP_536846.1', giid=>'GI:17981856'}, 
	    22=>{id=>9, predicted_weight=>24800, glycine_number=>9, glycine_weight=>21.6, 
	        urea_number=>'17, 18', urea_weight=>16, name=>'ATP synthase F0 subunit 6', 
	        npid=>'NP_536848.1', giid=>'GI:17981858'}, 
	    33=>{id=>10, predicted_weight=>18600, glycine_number=>10, glycine_weight=>16.7, 
	        urea_number=>, urea_weight=>, name=>'NADH dehydrogenase subunit 6', 
	        npid=>'NP_536854.1', giid=>'GI:17981864'}, 
	    25=>{id=>11, predicted_weight=>13200, glycine_number=>12, glycine_weight=>13.5, 
	        urea_number=>'23, 24', urea_weight=>6, name=>'NADH dehydrogenase subunit 3', 
	        npid=>'NP_536850.1', giid=>'GI:17981860'}, 
	    27=>{id=>12, predicted_weight=>10700, glycine_number=>11, glycine_weight=>14.8, 
	        urea_number=>26, urea_weight=>3.5, name=>'NADH dehydrogenase subunit 4L', 
	        npid=>'NP_536851.1', giid=>'GI:17981861'}, 
	    21=>{id=>13, predicted_weight=>7900, glycine_number=>'13+13a', glycine_weight=>9.8, 
	        urea_number=>25, urea_weight=>4.5, name=>'ATP synthase F0 subunit 8', 
	        npid=>'NP_536847.1', giid=>'GI:17981857'}
    );

    return \%protein_table;
}

sub _protein_lookup {
    my ($orig, $self, $index) = @_;

    my $p = $self->$orig;
    return $p unless defined $index;

    confess "locus index: $index is not a digit" unless $index =~ /^\d+$/;
    confess "locus index: $index is not valid" unless $p->{$index};

    return $p->{$index};

}#FOLDEND



=head2 transcript

The human mitochondrial transcript sequences indexed by their locus ID value.  Values returned are the RNA strings that would be seen following transcription: compensation has been made for the strand on which the locus is encoded, uracil replaces thymine, and any poly-adenlyation.

 $r = Bio::Mitomaster::SpeciesRef->new();
 $t = $r->transcript(16);  # $t contains the transcript string for the Cytochrome oxidase subunit I
 $t = $r->transcript(16,1,10);  # same thing but only the first 10 nucleotides
 $t = $r->transcript();  # $t is a hashref to a hash that contains all the mitochondrially encoded transcripts indexed by locus ID#

The options for retrieving sub string values of a transcript is documented in the section Seq String Attributes.

=cut


has 'transcript' => (#FOLDBEG
    is => 'ro',
    lazy_build => 1,
);

around 'transcript' => \&_get_transcript_seq;


sub _build_transcript {

    my $self = shift;

    my %transcript_table = (

        1=>'GUUUAUGUAGCUUACCUCCUCAAAGCAAUACACUGAAAAUGUUUAGACGGGCUCACAUCACCCCAUAAACA',
        2=>'AAUAGGUUUGGUCCUAGCCUUUCUAUUAGCUCUUAGUAAGAUUACACAUGCAAGCAUCCCCGUUCCAGUGAGUUCACCCUCUAAAUCACCACGAUCAAAAGGAACAAGCAUCAAGCACGCAGCAAUGCAGCUCAAAACGCUUAGCCUAGCCACACCCCCACGGGAAACAGCAGUGAUUAACCUUUAGCAAUAAACGAAAGUUUAACUAAGCUAUACUAACCCCAGGGUUGGUCAAUUUCGUGCCAGCCACCGCGGUCACACGAUUAACCCAAGUCAAUAGAAGCCGGCGUAAAGAGUGUUUUAGAUCACCCCCUCCCCAAUAAAGCUAAAACUCACCUGAGUUGUAAAAAACUCCAGUUGACACAAAAUAGACUACGAAAGUGGCUUUAACAUAUCUGAACACACAAUAGCUAAGACCCAAACUGGGAUUAGAUACCCCACUAUGCUUAGCCCUAAACCUCAACAGUUAAAUCAACAAAACUGCUCGCCAGAACACUACGAGCCACAGCUUAAAACUCAAAGGACCUGGCGGUGCUUCAUAUCCCUCUAGAGGAGCCUGUUCUGUAAUCGAUAAACCCCGAUCAACCUCACCACCUCUUGCUCAGCCUAUAUACCGCCAUCUUCAGCAAACCCUGAUGAAGGCUACAAAGUAAGCGCAAGUACCCACGUAAAGACGUUAGGUCAAGGUGUAGCCCAUGAGGUGGCAAGAAAUGGGCUACAUUUUCUACCCCAGAAAACUACGAUAGCCCUUAUGAAACUUAAGGGUCGAAGGUGGAUUUAGCAGUAAACUAAGAGUAGAGUGCUUAGUUGAACAGGGCCCUGAAGCGCGUACACACCGCCCGUCACCCUCCUCAAGUAUACUUCAAAGGACAUUUAACUAAAACCCCUACGCAUUUAUAUAGAGGAGACAAGUCGUAACAUGGUAAGUGUACUGGAAAGUGCACUUGGACGAAC',

        3=>'CAGAGUGUAGCUUAACACAAAGCACCCAACUUACACUUAGGAGAUUUCAACUUAACUUGACCGCUCUGA',

        4=>'GCUAAACCUAGCCCCAAACCCACUCCACCUUACUACCAGACAACCUUAGCCAAACCAUUUACCCAAAUAAAGUAUAGGCGAUAGAAAUUGAAACCUGGCGCAAUAGAUAUAGUACCGCAAGGGAAAGAUGAAAAAUUAUAACCAAGCAUAAUAUAGCAAGGACUAACCCCUAUACCUUCUGCAUAAUGAAUUAACUAGAAAUAACUUUGCAAGGAGAGCCAAAGCUAAGACCCCCGAAACCAGACGAGCUACCUAAGAACAGCUAAAAGAGCACACCCGUCUAUGUAGCAAAAUAGUGGGAAGAUUUAUAGGUAGAGGCGACAAACCUACCGAGCCUGGUGAUAGCUGGUUGUCCAAGAUAGAAUCUUAGUUCAACUUUAAAUUUGCCCACAGAACCCUCUAAAUCCCCUUGUAAAUUUAACUGUUAGUCCAAAGAGGAACAGCUCUUUGGACACUAGGAAAAAACCUUGUAGAGAGAGUAAAAAAUUUAACACCCAUAGUAGGCCUAAAAGCAGCCACCAAUUAAGAAAGCGUUCAAGCUCAACACCCACUACCUAAAAAAUCCCAAACAUAUAACUGAACUCCUCACACCCAAUUGGACCAAUCUAUCACCCUAUAGAAGAACUAAUGUUAGUAUAAGUAACAUGAAAACAUUCUCCUCCGCAUAAGCCUGCGUCAGAUUAAAACACUGAACUGACAAUUAACAGCCCAAUAUCUACAAUCAACCAACAAGUCAUUAUUACCCUCACUGUCAACCCAACACAGGCAUGCUCAUAAGGAAAGGUUAAAAAAAGUAAAAGGAACUCGGCAAAUCUUACCCCGCCUGUUUACCAAAAACAUCACCUCUAGCAUCACCAGUAUUAGAGGCACCGCCUGCCCAGUGACACAUGUUUAACGGCCGCGGUACCCUAACCGUGCAAAGGUAGCAUAAUCACUUGUUCCUUAAAUAGGGACCUGUAUGAAUGGCUCCACGAGGGUUCAGCUGUCUCUUACUUUUAACCAGUGAAAUUGACCUGCCCGUGAAGAGGCGGGCAUAACACAGCAAGACGAGAAGACCCUAUGGAGCUUUAAUUUAUUAAUGCAAACAGUACCUAACAAACCCACAGGUCCUAAACUACCAAACCUGCAUUAAAAAUUUCGGUUGGGGCGACCUCGGAGCAGAACCCAACCUCCGAGCAGUACAUGCUAAGACUUCACCAGUCAAAGCGAACUACUAUACUCAAUUGAUCCAAUAACUUGACCAACGGAACAAGUUACCCUAGGGAUAACAGCGCAAUCCUAUUCUAGAGUCCAUAUCAACAAUAGGGUUUACGACCUCGAUGUUGGAUCAGGACAUCCCGAUGGUGCAGCCGCUAUUAAAGGUUCGUUUGUUCAACGAUUAAAGUCCUACGUGAUCUGAGUUCAGACCGGAGUAAUCCAGGUCGGUUUCUAUCUACUUCAAAUUCCUCCCUGUACGAAAGGACAAGAGAAAUAAGGCCUACUUCACAAAGCGCCUUCCCCCGUAAAUGAUAUCAUCUCAACUUAGUAUUAUACCCACACCCACCCAAGAACAGGGUUUG',

        5=>'GUUAAGAUGGCAGAGCCCGGUAAUCGCAUAAAACUUAAAACUUUACAGUCAGAGGUUCAAUUCCUCUUCUUAACA',
	6=>'AUACCCAUGGCCAACCUCCUACUCCUCAUUGUACCCAUUCUAAUCGCAAUGGCAUUCCUAAUGCUUACCGAACGAAAAAUUCUAGGCUAUAUACAACUACGCAAAGGCCCCAACGUUGUAGGCCCCUACGGGCUACUACAACCCUUCGCUGACGCCAUAAAACUCUUCACCAAAGAGCCCCUAAAACCCGCCACAUCUACCAUCACCCUCUACAUCACCGCCCCGACCUUAGCUCUCACCAUCGCUCUUCUACUAUGAACCCCCCUCCCCAUACCCAACCCCCUGGUCAACCUCAACCUAGGCCUCCUAUUUAUUCUAGCCACCUCUAGCCUAGCCGUUUACUCAAUCCUCUGAUCAGGGUGAGCAUCAAACUCAAACUACGCCCUGAUCGGCGCACUGCGAGCAGUAGCCCAAACAAUCUCAUAUGAAGUCACCCUAGCCAUCAUUCUACUAUCAACAUUACUAAUAAGUGGCUCCUUUAACCUCUCCACCCUUAUCACAACACAAGAACACCUCUGAUUACUCCUGCCAUCAUGACCCUUGGCCAUAAUAUGAUUUAUCUCCACACUAGCAGAGACCAACCGAACCCCCUUCGACCUUGCCGAAGGGGAGUCCGAACUAGUCUCAGGCUUCAACAUCGAAUACGCCGCAGGCCCCUUCGCCCUAUUCUUCAUAGCCGAAUACACAAACAUUAUUAUAAUAAACACCCUCACCACUACAAUCUUCCUAGGAACAACAUAUGACGCACUCUCCCCUGAACUCUACACAACAUAUUUUGUCACCAAGACCCUACUUCUAACCUCCCUGUUCUUAUGAAUUCGAACAGCAUACCCCCGAUUCCGCUACGACCAACUCAUACACCUCCUAUGAAAAAACUUCCUACCACUCACCCUAGCAUUACUUAUAUGAUAUGUCUCCAUACCCAUUACAAUCUCCAGCAUUCCCCCUCAAACCUAA',

        7=>'AGAAAUAUGUCUGAUAAAAGAGUUACUUUGAUAGAGUAAAUAAUAGGAGCUUAAACCCCCUUAUUUCUA',

        8=>'UAGGAUGGGGUGUGAUAGGUGGCACGGAGAAUUUUGGAUUCUCAGGGAUGGGUUCGAUUCUCAUAGUCCUAG',

        9=>'AGUAAGGUCAGCUAAAUAAGCUAUCGGGCCCAUACCCCGAAAAUGUUGGUUAUACCCUUCCCGUACUA',
	10=>'AUUAAUCCCCUGGCCCAACCCGUCAUCUACUCUACCAUCUUUGCAGGCACACUCAUCACAGCGCUAAGCUCGCACUGAUUUUUUACCUGAGUAGGCCUAGAAAUAAACAUGCUAGCUUUUAUUCCAGUUCUAACCAAAAAAAUAAACCCUCGUUCCACAGAAGCUGCCAUCAAGUAUUUCCUCACGCAAGCAACCGCAUCCAUAAUCCUUCUAAUAGCUAUCCUCUUCAACAAUAUACUCUCCGGACAAUGAACCAUAACCAAUACUACCAAUCAAUACUCAUCAUUAAUAAUCAUAAUAGCUAUAGCAAUAAAACUAGGAAUAGCCCCCUUUCACUUCUGAGUCCCAGAGGUUACCCAAGGCACCCCUCUGACAUCCGGCCUGCUUCUUCUCACAUGACAAAAACUAGCCCCCAUCUCAAUCAUAUACCAAAUCUCUCCCUCACUAAACGUAAGCCUUCUCCUCACUCUCUCAAUCUUAUCCAUCAUAGCAGGCAGUUGAGGUGGAUUAAACCAAACCCAGCUACGCAAAAUCUUAGCAUACUCCUCAAUUACCCACAUAGGAUGAAUAAUAGCAGUUCUACCGUACAACCCUAACAUAACCAUUCUUAAUUUAACUAUUUAUAUUAUCCUAACUACUACCGCAUUCCUACUACUCAACUUAAACUCCAGCACCACGACCCUACUACUAUCUCGCACCUGAAACAAGCUAACAUGACUAACACCCUUAAUUCCAUCCACCCUCCUCUCCCUAGGAGGCCUGCCCCCGCUAACCGGCUUUUUGCCCAAAUGGGCCAUUAUCGAAGAAUUCACAAAAAACAAUAGCCUCAUCAUCCCCACCAUCAUAGCCACCAUCACCCUCCUUAACCUCUACUUCUACCUACGCCUAAUCUACUCCACCUCAAUCACACUACUCCCCAUAUCUAACAACGUAAAAAUAAAAUGACAGUUUGAACAUACAAAACCCACCCCAUUCCUCCCCACACUCAUCGCCCUUACCACGCUACUCCUACCUAUCUCCCCUUUUAUACUAAUAAUCUUAUAA',

        11=>'AGAAAUUUAGGUUAAAUACAGACCAAGAGCCUUCAAAGCCCUCAGUAAGUUGCAAUACUUAAUUUCUG',

        12=>'AAGGGCUUAGCUUAAUUAAAGUGGCUGAUUUGCGUUCAGUUGAUGCAGAGUGGGGUUUUGCAGUCCUUA',

        13=>'UAGAUUGAAGCCAGUUGAUUAGGGUGCUUAGCUGUUAACUAAGUGUUUGUGGGUUUAAGUCCCAUUGGUCUAG',

        14=>'AGCUCCGAGGUGAUUUUCAUAUUGAAUUGCAAAUUCGAAGAAGCAGCUUCAAACCUGCCGGGGCUU',

        15=>'GGUAAAAUGGCUGAGUGAAGCAUUGGACUGUAAAUCUAAAGACAGGGGUUAGGCCUCUUUUUACCA',
	16=>'AUGUUCGCCGACCGUUGACUAUUCUCUACAAACCACAAAGACAUUGGAACACUAUACCUAUUAUUCGGCGCAUGAGCUGGAGUCCUAGGCACAGCUCUAAGCCUCCUUAUUCGAGCCGAGCUGGGCCAGCCAGGCAACCUUCUAGGUAACGACCACAUCUACAACGUUAUCGUCACAGCCCAUGCAUUUGUAAUAAUCUUCUUCAUAGUAAUACCCAUCAUAAUCGGAGGCUUUGGCAACUGACUAGUUCCCCUAAUAAUCGGUGCCCCCGAUAUGGCGUUUCCCCGCAUAAACAACAUAAGCUUCUGACUCUUACCUCCCUCUCUCCUACUCCUGCUCGCAUCUGCUAUAGUGGAGGCCGGAGCAGGAACAGGUUGAACAGUCUACCCUCCCUUAGCAGGGAACUACUCCCACCCUGGAGCCUCCGUAGACCUAACCAUCUUCUCCUUACACCUAGCAGGUGUCUCCUCUAUCUUAGGGGCCAUCAAUUUCAUCACAACAAUUAUCAAUAUAAAACCCCCUGCCAUAACCCAAUACCAAACGCCCCUCUUCGUCUGAUCCGUCCUAAUCACAGCAGUCCUACUUCUCCUAUCUCUCCCAGUCCUAGCUGCUGGCAUCACUAUACUACUAACAGACCGCAACCUCAACACCACCUUCUUCGACCCCGCCGGAGGAGGAGACCCCAUUCUAUACCAACACCUAUUCUGAUUUUUCGGUCACCCUGAAGUUUAUAUUCUUAUCCUACCAGGCUUCGGAAUAAUCUCCCAUAUUGUAACUUACUACUCCGGAAAAAAAGAACCAUUUGGAUACAUAGGUAUGGUCUGAGCUAUGAUAUCAAUUGGCUUCCUAGGGUUUAUCGUGUGAGCACACCAUAUAUUUACAGUAGGAAUAGACGUAGACACACGAGCAUAUUUCACCUCCGCUACCAUAAUCAUCGCUAUCCCCACCGGCGUCAAAGUAUUUAGCUGACUCGCCACACUCCACGGAAGCAAUAUGAAAUGAUCUGCUGCAGUGCUCUGAGCCCUAGGAUUCAUCUUUCUUUUCACCGUAGGUGGCCUGACUGGCAUUGUAUUAGCAAACUCAUCACUAGACAUCGUACUACACGACACGUACUACGUUGUAGCCCACUUCCACUAUGUCCUAUCAAUAGGAGCUGUAUUUGCCAUCAUAGGAGGCUUCAUUCACUGAUUUCCCCUAUUCUCAGGCUACACCCUAGACCAAACCUACGCCAAAAUCCAUUUCACUAUCAUAUUCAUCGGCGUAAAUCUAACUUUCUUCCCACAACACUUUCUCGGCCUAUCCGGAAUGCCCCGACGUUACUCGGACUACCCCGAUGCAUACACCACAUGAAACAUCCUAUCAUCUGUAGGCUCAUUCAUUUCUCUAACAGCAGUAAUAUUAAUAAUUUUCAUGAUUUGAGAAGCCUUCGCUUCGAAGCGAAAAGUCCUAAUAGUAGAAGAACCCUCCAUAAACCUGGAGUGACUAUAUGGAUGCCCCCCACCCUACCACACAUUCGAAGAACCCGUAUACAUAAAAUCUAGA',

        17=>'UUGAAAAAGUCAUGGAGGCCAUGGGGUUGGCUUGAAACCAGCUUUGGGGGGUUCGAUUCCUUCCUUUUUUG',

        18=>'AAGGUAUUAGAAAAACCAUUUCAUAACUUUGUCAAAGUUAAAUUAUAGGCUAAAUCCUAUAUAUCUUA',
	19=>'AUGGCACAUGCAGCGCAAGUAGGUCUACAAGACGCUACUUCCCCUAUCAUAGAAGAGCUUAUCACCUUUCAUGAUCACGCCCUCAUAAUCAUUUUCCUUAUCUGCUUCCUAGUCCUGUAUGCCCUUUUCCUAACACUCACAACAAAACUAACUAAUACUAACAUCUCAGACGCUCAGGAAAUAGAAACCGUCUGAACUAUCCUGCCCGCCAUCAUCCUAGUCCUCAUCGCCCUCCCAUCCCUACGCAUCCUUUACAUAACAGACGAGGUCAACGAUCCCUCCCUUACCAUCAAAUCAAUUGGCCACCAAUGGUACUGAACCUACGAGUACACCGACUACGGCGGACUAAUCUUCAACUCCUACAUACUUCCCCCAUUAUUCCUAGAACCAGGCGACCUGCGACUCCUUGACGUUGACAAUCGAGUAGUACUCCCGAUUGAAGCCCCCAUUCGUAUAAUAAUUACAUCACAAGACGUCUUGCACUCAUGAGCUGUCCCCACAUUAGGCUUAAAAACAGAUGCAAUUCCCGGACGUCUAAACCAAACCACUUUCACCGCUACACGACCGGGGGUAUACUACGGUCAAUGCUCUGAAAUCUGUGGAGCAAACCACAGUUUCAUGCCCAUCGUCCUAGAAUUAAUUCCCCUAAAAAUCUUUGAAAUAGGGCCCGUAUUUACCCUAUAG',

        20=>'CACUGUAAAGCUAACUUAGCAUUAACCUUUUAAGUUAAAGAUUAAGAGAACCAACACCUCUUUACAGUGA',
	    21=>'AUGCCCCAACUAAAUACUACCGUAUGGCCCACCAUAAUUACCCCCAUACUCCUUACACUAUUCCUCAUCACCCAACUAAAAAUAUUAAACACAAACUACCACCUACCUCCCUCACCAAAGCCCAUAAAAAUAAAAAAUUAUAACAAACCCUGAGAACCAAAAUGAACGAAAAUCUGUUCGCUUCAUUCAUUGCCCCCACAAUCCUAG',
	22=>'AUGAACGAAAAUCUGUUCGCUUCAUUCAUUGCCCCCACAAUCCUAGGCCUACCCGCCGCAGUACUGAUCAUUCUAUUUCCCCCUCUAUUGAUCCCCACCUCCAAAUAUCUCAUCAACAACCGACUAAUCACCACCCAACAAUGACUAAUCAAACUAACCUCAAAACAAAUGAUAACCAUACACAACACUAAAGGACGAACCUGAUCUCUUAUACUAGUAUCCUUAAUCAUUUUUAUUGCCACAACUAACCUCCUCGGACUCCUGCCUCACUCAUUUACACCAACCACCCAACUAUCUAUAAACCUAGCCAUGGCCAUCCCCUUAUGAGCGGGCACAGUGAUUAUAGGCUUUCGCUCUAAGAUUAAAAAUGCCCUAGCCCACUUCUUACCACAAGGCACACCUACACCCCUUAUCCCCAUACUAGUUAUUAUCGAAACCAUCAGCCUACUCAUUCAACCAAUAGCCCUGGCCGUACGCCUAACCGCUAACAUUACUGCAGGCCACCUACUCAUGCACCUAAUUGGAAGCGCCACCCUAGCAAUAUCAACCAUUAACCUUCCCUCUACACUUAUCAUCUUCACAAUUCUAAUUCUACUGACUAUCCUAGAAAUCGCUGUCGCCUUAAUCCAAGCCUACGUUUUCACACUUCUAGUAAGCCUCUACCUGCACGACAACACAUAA',
	23=>'AUGACCCACCAAUCACAUGCCUAUCAUAUAGUAAAACCCAGCCCAUGACCCCUAACAGGGGCCCUCUCAGCCCUCCUAAUGACCUCCGGCCUAGCCAUGUGAUUUCACUUCCACUCCAUAACGCUCCUCAUACUAGGCCUACUAACCAACACACUAACCAUAUACCAAUGAUGGCGCGAUGUAACACGAGAAAGCACAUACCAAGGCCACCACACACCACCUGUCCAAAAAGGCCUUCGAUACGGGAUAAUCCUAUUUAUUACCUCAGAAGUUUUUUUCUUCGCAGGAUUUUUCUGAGCCUUUUACCACUCCAGCCUAGCCCCUACCCCCCAAUUAGGAGGGCACUGGCCCCCAACAGGCAUCACCCCGCUAAAUCCCCUAGAAGUCCCACUCCUAAACACAUCCGUAUUACUCGCAUCAGGAGUAUCAAUCACCUGAGCUCACCAUAGUCUAAUAGAAAACAACCGAAACCAAAUAAUUCAAGCACUGCUUAUUACAAUUUUACUGGGUCUCUAUUUUACCCUCCUACAAGCCUCAGAGUACUUCGAGUCUCCCUUCACCAUUUCCGACGGCAUCUACGGCUCAACAUUUUUUGUAGCCACAGGCUUCCACGGACUUCACGUCAUUAUUGGCUCAACUUUCCUCACUAUCUGCUUCAUCCGCCAACUAAUAUUUCACUUUACAUCCAAACAUCACUUUGGCUUCGAAGCCGCCGCCUGAUACUGGCAUUUUGUAGAUGUGGUUUGACUAUUUCUGUAUGUCUCCAUCUAUUGAUGAGGGUCUUAA',

        24=>'ACUCUUUUAGUAUAAAUAGUACCGUUAACUUCCAAUUAACUAGUUUUGACAACAUUCAAAAAAGAGUA',
	25=>'AUAAACUUCGCCUUAAUUUUAAUAAUCAACACCCUCCUAGCCUUACUACUAAUAAUUAUUACAUUUUGACUACCACAACUCAACGGCUACAUAGAAAAAUCCACCCCUUACGAGUGCGGCUUCGACCCUAUAUCCCCCGCCCGCGUCCCUUUCUCCAUAAAAUUCUUCUUAGUAGCUAUUACCUUCUUAUUAUUUGAUCUAGAAAUUGCCCUCCUUUUACCCCUACCAUGAGCCCUACAAACAACUAACCUGCCACUAAUAGUUAUGUCAUCCCUCUUAUUAAUCAUCAUCCUAGCCCUAAGUCUGGCCUAUGAGUGACUACAAAAAGGAUUAGACUGAACCGAAUAA',

        26=>'UGGUAUAUAGUUUAAACAAAACGAAUGAUUUCGACUCAUUAAAUUAUGAUAAUCAUAUUUACCAA',

        27=>'AUGCCCCUCAUUUACAUAAAUAUUAUACUAGCAUUUACCAUCUCACUUCUAGGAAUACUAGUAUAUCGCUCACACCUCAUAUCCUCCCUACUAUGCCUAGAAGGAAUAAUACUAUCGCUGUUCAUUAUAGCUACUCUCAUAACCCUCAACACCCACUCCCUCUUAGCCAAUAUUGUGCCUAUUGCCAUACUAGUCUUUGCCGCCUGCGAAGCAGCGGUGGGCCUAGCCCUACUAGUCUCAAUCUCCAACACAUAUGGCCUAGACUACGUACAUAACCUAAACCUACUCCAAUGCUAA',
	28=>'AUGCUAAAACUAAUCGUCCCAACAAUUAUAUUACUACCACUGACAUGACUUUCCAAAAAACACAUAAUUUGAAUCAACACAACCACCCACAGCCUAAUUAUUAGCAUCAUCCCUCUACUAUUUUUUAACCAAAUCAACAACAACCUAUUUAGCUGUUCCCCAACCUUUUCCUCCGACCCCCUAACAACCCCCCUCCUAAUACUAACUACCUGACUCCUACCCCUCACAAUCAUGGCAAGCCAACGCCACUUAUCCAGUGAACCACUAUCACGAAAAAAACUCUACCUCUCUAUACUAAUCUCCCUACAAAUCUCCUUAAUUAUAACAUUCACAGCCACAGAACUAAUCAUAUUUUAUAUCUUCUUCGAAACCACACUUAUCCCCACCUUGGCUAUCAUCACCCGAUGAGGCAACCAGCCAGAACGCCUGAACGCAGGCACAUACUUCCUAUUCUACACCCUAGUAGGCUCCCUUCCCCUACUCAUCGCACUAAUUUACACUCACAACACCCUAGGCUCACUAAACAUUCUACUACUCACUCUCACUGCCCAAGAACUAUCAAACUCCUGAGCCAACAACUUAAUAUGACUAGCUUACACAAUAGCUUUUAUAGUAAAGAUACCUCUUUACGGACUCCACUUAUGACUCCCUAAAGCCCAUGUCGAAGCCCCCAUCGCUGGGUCAAUAGUACUUGCCGCAGUACUCUUAAAACUAGGCGGCUAUGGUAUAAUACGCCUCACACUCAUUCUCAACCCCCUGACAAAACACAUAGCCUACCCCUUCCUUGUACUAUCCCUAUGAGGCAUAAUUAUAACAAGCUCCAUCUGCCUACGACAAACAGACCUAAAAUCGCUCAUUGCAUACUCUUCAAUCAGCCACAUAGCCCUCGUAGUAACAGCCAUUCUCAUCCAAACCCCCUGAAGCUUCACCGGCGCAGUCAUUCUCAUAAUCGCCCACGGGCUUACAUCCUCAUUACUAUUCUGCCUAGCAAACUCAAACUACGAACGCACUCACAGUCGCAUCAUAAUCCUCUCUCAAGGACUUCAAACUCUACUCCCACUAAUAGCUUUUUGAUGACUUCUAGCAAGCCUCGCUAACCUCGCCUUACCCCCCACUAUUAACCUACUGGGAGAACUCUCUGUGCUAGUAACCACGUUCUCCUGAUCAAAUAUCACUCUCCUACUUACAGGACUCAACAUACUAGUCACAGCCCUAUACUCCCUCUACAUAUUUACCACAACACAAUGGGGCUCACUCACCCACCACAUUAACAACAUAAAACCCUCAUUCACACGAGAAAACACCCUCAUGUUCAUACACCUAUCCCCCAUUCUCCUCCUAUCCCUCAACCCCGACAUCAUUACCGGGUUUUCCUCUUAA',

        29=>'GUAAAUAUAGUUUAACCAAAACAUCAGAUUGUGAAUCUGACAACAGAGGCUUACGACCCCUUAUUUACC',

        30=>'GAGAAAGCUCACAAGAACUGCUAACUCAUGCCCCCAUGUCUAACAACAUGGCUUUCUCA',

        31=>'ACUUUUAAAGGAUAACAGCUAUCCAUUGGUCUUAGGCCCCAAAAAUUUUGGUGCAACUCCAAAUAAAAGUA',
	32=>'AUAACCAUGCACACUACUAUAACCACCCUAACCCUGACUUCCCUAAUUCCCCCCAUCCUUACCACCCUCGUUAACCCUAACAAAAAAAACUCAUACCCCCAUUAUGUAAAAUCCAUUGUCGCAUCCACCUUUAUUAUCAGUCUCUUCCCCACAACAAUAUUCAUGUGCCUAGACCAAGAAGUUAUUAUCUCGAACUGACACUGAGCCACAACCCAAACAACCCAGCUCUCCCUAAGCUUCAAACUAGACUACUUCUCCAUAAUAUUCAUCCCUGUAGCAUUGUUCGUUACAUGGUCCAUCAUAGAAUUCUCACUGUGAUAUAUAAACUCAGACCCAAACAUUAAUCAGUUCUUCAAAUAUCUACUCAUCUUCCUAAUUACCAUACUAAUCUUAGUUACCGCUAACAACCUAUUCCAACUGUUCAUCGGCUGAGAGGGCGUAGGAAUUAUAUCCUUCUUGCUCAUCAGUUGAUGAUACGCCCGAGCAGAUGCCAACACAGCAGCCAUUCAAGCAAUCCUAUACAACCGUAUCGGCGAUAUCGGUUUCAUCCUCGCCUUAGCAUGAUUUAUCCUACACUCCAACUCAUGAGACCCACAACAAAUAGCCCUUCUAAACGCUAAUCCAAGCCUCACCCCACUACUAGGCCUCCUCCUAGCAGCAGCAGGCAAAUCAGCCCAAUUAGGUCUCCACCCCUGACUCCCCUCAGCCAUAGAAGGCCCCACCCCAGUCUCAGCCCUACUCCACUCAAGCACUAUAGUUGUAGCAGGAAUCUUCUUACUCAUCCGCUUCCACCCCCUAGCAGAAAAUAGCCCACUAAUCCAAACUCUAACACUAUGCUUAGGCGCUAUCACCACUCUGUUCGCAGCAGUCUGCGCCCUUACACAAAAUGACAUCAAAAAAAUCGUAGCCUUCUCCACUUCAAGUCAACUAGGACUCAUAAUAGUUACAAUCGGCAUCAACCAACCACACCUAGCAUUCCUGCACAUCUGUACCCACGCCUUCUUCAAAGCCAUACUAUUUAUGUGCUCCGGGUCCAUCAUCCACAACCUUAACAAUGAACAAGAUAUUCGAAAAAUAGGAGGACUACUCAAAACCAUACCUCUCACUUCAACCUCCCUCACCAUUGGCAGCCUAGCAUUAGCAGGAAUACCUUUCCUCACAGGUUUCUACUCCAAAGACCACAUCAUCGAAACCGCAAACAUAUCAUACACAAACGCCUGAGCCCUAUCUAUUACUCUCAUCGCUACCUCCCUGACAAGCGCCUAUAGCACUCGAAUAAUUCUUCUCACCCUAACAGGUCAACCUCGCUUCCCCACCCUUACUAACAUUAACGAAAAUAACCCCACCCUACUAAACCCCAUUAAACGCCUGGCAGCCGGAAGCCUAUUCGCAGGAUUUCUCAUUACUAACAACAUUUCCCCCGCAUCCCCCUUCCAAACAACAAUCCCCCUCUACCUAAAACUCACAGCCCUCGCUGUCACUUUCCUAGGACUUCUAACAGCCCUAGACCUCAACUACCUAACCAACAAACUUAAAAUAAAAUCCCCACUAUGCACAUUUUAUUUCUCCAACAUACUCGGAUUCUACCCUAGCAUCACACACCGCACAAUCCCCUAUCUAGGCCUUCUUACGAGCCAAAACCUGCCCCUACUCCUCCUAGACCUAACCUGACUAGAAAAGCUAUUACCUAAAACAAUUUCACAGCACCAAAUCUCCACCUCCAUCAUCACCUCAACCCAAAAAGGCAUAAUUAAACUUUACUUCCUCUCUUUCUUCUUCCCACUCAUCCUAACCCUACUCCUAAUCACAUAA',
	33=>'AUGAUGUAUGCUUUGUUUCUGUUGAGUGUGGGUUUAGUAAUGGGGUUUGUGGGGUUUUCUUCUAAGCCUUCUCCUAUUUAUGGGGGUUUAGUAUUGAUUGUUAGCGGUGUGGUCGGGUGUGUUAUUAUUCUGAAUUUUGGGGGAGGUUAUAUGGGUUUAAUAGUUUUUUUAAUUUAUUUAGGGGGAAUGAUGGUUGUCUUUGGAUAUACUACAGCGAUGGCUAUUGAGGAGUAUCCUGAGGCAUGGGGGUCAGGGGUUGAGGUCUUGGUGAGUGUUUUAGUGGGGUUAGCGAUGGAGGUAGGAUUGGUGCUGUGGGUGAAAGAGUAUGAUGGGGUGGUGGUUGUGGUAAACUUUAAUAGUGUAGGAAGCUGAAUAAUUUAUGAAGGAGAGGGGUCAGGGUUGAUUCGGGAGGAUCCUAUUGGUGCGGGGGCUUUGUAUGAUUAUGGGCGUUGAUUAGUAGUAGUUACUGGUUGAACAUUGUUUGUUGGUGUAUAUAUUGUAAUUGAGAUUGCUCGGGGGAAUAGG',

        34=>'GUUCUUGUAGUUGAAAUACAACGAUGGUUUUUCAUAUCAUUGGUCGUGGUUGUAGUCCGUGCGAGAAUA',
	35=>'AUGACCCCAAUACGCAAAACUAACCCCCUAAUAAAAUUAAUUAACCACUCAUUCAUCGACCUCCCCACCCCAUCCAACAUCUCCGCAUGAUGAAACUUCGGCUCACUCCUUGGCGCCUGCCUGAUCCUCCAAAUCACCACAGGACUAUUCCUAGCCAUGCACUACUCACCAGACGCCUCAACCGCCUUUUCAUCAAUCGCCCACAUCACUCGAGACGUAAAUUAUGGCUGAAUCAUCCGCUACCUUCACGCCAAUGGCGCCUCAAUAUUCUUUAUCUGCCUCUUCCUACACAUCGGGCGAGGCCUAUAUUACGGAUCAUUUCUCUACUCAGAAACCUGAAACAUCGGCAUUAUCCUCCUGCUUGCAACUAUAGCAACAGCCUUCAUAGGCUAUGUCCUCCCGUGAGGCCAAAUAUCAUUCUGAGGGGCCACAGUAAUUACAAACUUACUAUCCGCCAUCCCAUACAUUGGGACAGACCUAGUUCAAUGAAUCUGAGGAGGCUACUCAGUAGACAGUCCCACCCUCACACGAUUCUUUACCUUUCACUUCAUCUUGCCCUUCAUUAUUGCAGCCCUAGCAACACUCCACCUCCUAUUCUUGCACGAAACGGGAUCAAACAACCCCCUAGGAAUCACCUCCCAUUCCGAUAAAAUCACCUUCCACCCUUACUACACAAUCAAAGACGCCCUCGGCUUACUUCUCUUCCUUCUCUCCUUAAUGACAUUAACACUAUUCUCACCAGACCUCCUAGGCGACCCAGACAAUUAUACCCUAGCCAACCCCUUAAACACCCCUCCCCACAUCAAGCCCGAAUGAUAUUUCCUAUUCGCCUACACAAUUCUCCGAUCCGUCCCUAACAAACUAGGAGGCGUCCUUGCCCUAUUACUAUCCAUCCUCAUCCUAGCAAUAAUCCCCAUCCUCCAUAUAUCCAAACAACAAAGCAUAAUAUUUCGCCCACUAAGCCAAUCACUUUAUUGACUCCUAGCCGCAGACCUCCUCAUUCUAACCUGAAUCGGAGGACAACCAGUAAGCUACCCUUUUACCAUCAUUGGACAAGUAGCAUCCGUACUAUACUUCACAACAAUCCUAAUCCUAAUACCAACUAUCUCCCUAAUUGAAAACAAAAUACUCAAAUGGGCCUAA',

        36=>'GUCCUUGUAGUAUAAACUAAUACACCAGUCUUGUAAACCGGAGAUGAAAACCUUUUUCCAAGGACA',

        37=>'CAGAGAAUAGUUUAAAUUAGAAUCUUAGCUUUGGGUGCUAAUGGUGGAGUUAAAGACUUUUUCUCUGAU',
    );

    return \%transcript_table;
}

sub _get_transcript_seq {

    # Method modifier to make the transcript method have the polymorphic and 
    # sub_seq behavior defined in SeqManipRole.pm

    my ($orig, $self, $index, $s, $e) = @_;

    confess "locus id: $index is not a digit" unless $index =~ /^\d+$/;

    my ($start, $end) = $self->get_seq_indices($s, $e);
    $self->check_seq_indices($start, $end, 0);

    my $t = $self->$orig;
    confess "locus id: $index is not a valid transcript index" unless $t->{$index};

    return $self->sub_seq($t->{$index}, $start, $end);

}#FOLDEND



=head2 translation

The human mitochondrial translation sequences.  Values are the amino acid strings (IUPAC single-letter symbols) that would be seen following translation using the transcript strings returned by the transcript method.

    $r = Bio::Mitomaster::SpeciesRef->new();
    $t = $r->translation(16);  # $t contains the translation string for the Cytochrome oxidase subunit I
    $t = $r->translation(16,1,10);  # same thing but only the first 10 amino acids
    $t = $r->translation();  # $t is a hashref to a hash that contains all the mitochondrially encoded translations indexed by locus ID#

The options for retrieving sub string values of a translation is documented in the section Seq String Attributes.

=cut


has 'translation' => (#FOLDBEG
    is => 'ro',
    lazy_build => 1,
);

around 'translation' => \&_get_translation_seq;


sub _build_translation {

    my $self = shift;

    my %translation_table = (
	32=>'MTMHTTMTTLTLTSLIPPILTTLVNPNKKNSYPHYVKSIVASTFIISLFPTTMFMCLDQEVIISNWHWATTQTTQLSLSFKLDYFSMMFIPVALFVTWSIMEFSLWYMNSDPNINQFFKYLLIFLITMLILVTANNLFQLFIGWEGVGIMSFLLISWWYARADANTAAIQAILYNRIGDIGFILALAWFILHSNSWDPQQMALLNANPSLTPLLGLLLAAAGKSAQLGLHPWLPSAMEGPTPVSALLHSSTMVVAGIFLLIRFHPLAENSPLIQTLTLCLGAITTLFAAVCALTQNDIKKIVAFSTSSQLGLMMVTIGINQPHLAFLHICTHAFFKAMLFMCSGSIIHNLNNEQDIRKMGGLLKTMPLTSTSLTIGSLALAGMPFLTGFYSKDHIIETANMSYTNAWALSITLIATSLTSAYSTRMILLTLTGQPRFPTLTNINENNPTLLNPIKRLAAGSLFAGFLITNNISPASPFQTTIPLYLKLTALAVTFLGLLTALDLNYLTNKLKMKSPLCTFYFSNMLGFYPSITHRTIPYLGLLTSQNLPLLLLDLTWLEKLLPKTISQHQISTSIITSTQKGMIKLYFLSFFFPLILTLLLIT',
	16=>'MFADRWLFSTNHKDIGTLYLLFGAWAGVLGTALSLLIRAELGQPGNLLGNDHIYNVIVTAHAFVMIFFMVMPIMIGGFGNWLVPLMIGAPDMAFPRMNNMSFWLLPPSLLLLLASAMVEAGAGTGWTVYPPLAGNYSHPGASVDLTIFSLHLAGVSSILGAINFITTIINMKPPAMTQYQTPLFVWSVLITAVLLLLSLPVLAAGITMLLTDRNLNTTFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIVTYYSGKKEPFGYMGMVWAMMSIGFLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAIPTGVKVFSWLATLHGSNMKWSAAVLWALGFIFLFTVGGLTGIVLANSSLDIVLHDTYYVVAHFHYVLSMGAVFAIMGGFIHWFPLFSGYTLDQTYAKIHFTIMFIGVNLTFFPQHFLGLSGMPRRYSDYPDAYTTWNILSSVGSFISLTAVMLMIFMIWEAFASKRKVLMVEEPSMNLEWLYGCPPPYHTFEEPVYMKS',
	28=>'MLKLIVPTIMLLPLTWLSKKHMIWINTTTHSLIISIIPLLFFNQINNNLFSCSPTFSSDPLTTPLLMLTTWLLPLTIMASQRHLSSEPLSRKKLYLSMLISLQISLIMTFTATELIMFYIFFETTLIPTLAIITRWGNQPERLNAGTYFLFYTLVGSLPLLIALIYTHNTLGSLNILLLTLTAQELSNSWANNLMWLAYTMAFMVKMPLYGLHLWLPKAHVEAPIAGSMVLAAVLLKLGGYGMMRLTLILNPLTKHMAYPFLVLSLWGMIMTSSICLRQTDLKSLIAYSSISHMALVVTAILIQTPWSFTGAVILMIAHGLTSSLLFCLANSNYERTHSRIMILSQGLQTLLPLMAFWWLLASLANLALPPTINLLGELSVLVTTFSWSNITLLLTGLNMLVTALYSLYMFTTTQWGSLTHHINNMKPSFTRENTLMFMHLSPILLLSLNPDIITGFSS',
	35=>'MTPMRKTNPLMKLINHSFIDLPTPSNISAWWNFGSLLGACLILQITTGLFLAMHYSPDASTAFSSIAHITRDVNYGWIIRYLHANGASMFFICLFLHIGRGLYYGSFLYSETWNIGIILLLATMATAFMGYVLPWGQMSFWGATVITNLLSAIPYIGTDLVQWIWGGYSVDSPTLTRFFTFHFILPFIIAALATLHLLFLHETGSNNPLGITSHSDKITFHPYYTIKDALGLLLFLLSLMTLTLFSPDLLGDPDNYTLANPLNTPPHIKPEWYFLFAYTILRSVPNKLGGVLALLLSILILAMIPILHMSKQQSMMFRPLSQSLYWLLAADLLILTWIGGQPVSYPFTIIGQVASVLYFTTILILMPTISLIENKMLKWA',
	10=>'INPLAQPVIYSTIFAGTLITALSSHWFFTWVGLEMNMLAFIPVLTKKMNPRSTEAAIKYFLTQATASMILLMAILFNNMLSGQWTMTNTTNQYSSLMIMMAMAMKLGMAPFHFWVPEVTQGTPLTSGLLLLTWQKLAPISIMYQISPSLNVSLLLTLSILSIMAGSWGGLNQTQLRKILAYSSITHMGWMMAVLPYNPNMTILNLTIYIILTTTAFLLLNLNSSTTTLLLSRTWNKLTWLTPLIPSTLLSLGGLPPLTGFLPKWAIIEEFTKNNSLIIPTIMATITLLNLYFYLRLIYSTSITLLPMSNNVKMKWQFEHTKPTPFLPTLIALTTLLLPISPFMLMIL',
	6=>'MPMANLLLLIVPILIAMAFLMLTERKILGYMQLRKGPNVVGPYGLLQPFADAMKLFTKEPLKPATSTITLYITAPTLALTIALLLWTPLPMPNPLVNLNLGLLFILATSSLAVYSILWSGWASNSNYALIGALRAVAQTISYEVTLAIILLSTLLMSGSFNLSTLITTQEHLWLLLPSWPLAMMWFISTLAETNRTPFDLAEGESELVSGFNIEYAAGPFALFFMAEYTNIIMMNTLTTTIFLGTTYDALSPELYTTYFVTKTLLLTSLFLWIRTAYPRFRYDQLMHLLWKNFLPLTLALLMWYVSMPITISSIPPQT',
	23=>'MTHQSHAYHMVKPSPWPLTGALSALLMTSGLAMWFHFHSMTLLMLGLLTNTLTMYQWWRDVTRESTYQGHHTPPVQKGLRYGMILFITSEVFFFAGFFWAFYHSSLAPTPQLGGHWPPTGITPLNPLEVPLLNTSVLLASGVSITWAHHSLMENNRNQMIQALLITILLGLYFTLLQASEYFESPFTISDGIYGSTFFVATGFHGLHVIIGSTFLTICFIRQLMFHFTSKHHFGFEAAAWYWHFVDVVWLFLYVSIYWWGS',
	19=>'MAHAAQVGLQDATSPIMEELITFHDHALMIIFLICFLVLYALFLTLTTKLTNTNISDAQEMETVWTILPAIILVLIALPSLRILYMTDEVNDPSLTIKSIGHQWYWTYEYTDYGGLIFNSYMLPPLFLEPGDLRLLDVDNRVVLPIEAPIRMMITSQDVLHSWAVPTLGLKTDAIPGRLNQTTFTATRPGVYYGQCSEICGANHSFMPIVLELIPLKIFEMGPVFTL',
	22=>'MNENLFASFIAPTILGLPAAVLIILFPPLLIPTSKYLINNRLITTQQWLIKLTSKQMMTMHNTKGRTWSLMLVSLIIFIATTNLLGLLPHSFTPTTQLSMNLAMAIPLWAGTVIMGFRSKIKNALAHFLPQGTPTPLIPMLVIIETISLLIQPMALAVRLTANITAGHLLMHLIGSATLAMSTINLPSTLIIFTILILLTILEIAVALIQAYVFTLLVSLYLHDNT',
	33=>'MMYALFLLSVGLVMGFVGFSSKPSPIYGGLVLIVSGVVGCVIILNFGGGYMGLMVFLIYLGGMMVVFGYTTAMAIEEYPEAWGSGVEVLVSVLVGLAMEVGLVLWVKEYDGVVVVVNFNSVGSWMIYEGEGSGLIREDPIGAGALYDYGRWLVVVTGWTLFVGVYIVIEIARGN',
	    25=>'MNFALILMINTLLALLLMIITFWLPQLNGYMEKSTPYECGFDPMSPARVPFSMKFFLVAITFLLFDLEIALLLPLPWALQTTNLPLMVMSSLLLIIILALSLAYEWLQKGLDWTE',
	    27=>'MPLIYMNIMLAFTISLLGMLVYRSHLMSSLLCLEGMMLSLFIMATLMTLNTHSLLANIVPIAMLVFAACEAAVGLALLVSISNTYGLDYVHNLNLLQC',
	    21=>'MPQLNTTVWPTMITPMLLTLFLITQLKMLNTNYHLPPSPKPMKMKNYNKPWEPKWTKICSLHSLPPQS'
	);   

    return \%translation_table;
}


sub _get_translation_seq {

    # A utility method to give the translation method the polymorphic and sub_seq 
    # functionality defined in SeqManipRole.pm

    my ($orig, $self, $index, $s, $e) = @_;

    confess "index: $index is not a valid translation index" unless $index =~ /^\d+$/;

    my ($start, $end) = $self->get_seq_indices($s, $e);
    $self->check_seq_indices($start, $end, 0);

    my $t = $self->$orig;
    confess "index: $index is not a valid translation index" unless $t->{$index};

    return $self->sub_seq($t->{$index}, $start, $end);

}#FOLDEND


sub _seq_string {#FOLDBEG
    my ($self, $s, $b, $e) = @_;

    confess "No string value supplied." unless defined $s;
    my $l = length($s);


    my $seq_string;

    if (defined $b) {
        confess "Beginning value must be an integer." 
            unless $b =~ /-?\d+/; 
            
        confess "Beginning value: $b not in the range -$l to $l"
            unless ($b <= $l or $b >= -$l);

        confess "Beginning value cannot be zero." if $b == 0;


        if (defined $e) {

            confess "Ending value must be an integer." 
                unless $e =~ /-?\d+/;

            confess "Ending value: $e not in the range -$l to $l"
                unless ($e <= $l or $e >= -$l);

            confess "Ending value cannot be zero." 
                if $e == 0;


            # A beginning and end were both defined and not zero

            if ( $b > 0 ) {
                # positive beginning
            
                if ( $e > 0 ) {
                    # positive ending
                    $seq_string = substr($s, $b - 1, ( $e - $b + 1));
                }
                else {
                    # negative ending
                    $e++;
                    $seq_string = $e ? substr($s, $b - 1, $e) : substr($s, $b - 1);
                }
            }
            else {
                # negative beginning

                if ( $e > 0 ) {
                    # positive ending
                    confess "Seq_string indices: $b and $e have no overlap"
                        unless ($e + abs($b)) >= $l;

                    $seq_string = substr($s, $l + $b, $e - ($l + $b));
                }
                else {
                    # negative ending

                    confess "Seq_string indices: $b and $e have no overlap"
                        unless $b <= $e;

                    $seq_string = substr($s, $l + $b, ($e - $b) + 1);
                }

            }

        }

        else {

            # A beginning with no end was given
            if ( $b > 0 ) {
                $seq_string = substr($s, $b - 1, 1);
            }
            else {
                $seq_string = substr($s, $b, 1);
            }

        }
    }
    elsif (defined $e) {
        # Getting here means that undef was passed for the beginning, but we
        # cannot have an end without defining a beginning, so this is an
        # error.
        confess "Ending value given without beginning value"
            if defined $e and not defined $b;
    }
    else {
        # Neither a beginning or end was supplied - a wasted call, so just
        # return the unmodified string.
        $seq_string = $s;
    }

    return $seq_string; 


}#FOLDEND



=head1 METHODS

=head2 wrapping

Alias accessor/setter for the dna_wrapping attribute.  We define an alias since dna_wrapping is normally the value we want for wrapping.

=cut


sub wrapping {#FOLDBEG
    # This method is an alias for accessing/setting the dna_wrapping attribute.  
    my $self = shift;

    if (@_) {
        $self->dna_wrapping(shift);
        return 1;
    }
    else {
        return $self->dna_wrapping;
    }

}#FOLDEND



=head1 AUTHOR

Marty Brandon, C<< <marty.brandon at gmail.com> >>



=head1 BUGS

Please report any bugs or feature requests to C<bug-bio-mm-refmitohuman at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Bio-Bio::Mitomaster-RefMitoHuman>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.



=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Bio::Mitomaster::SpeciesRef

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

no Moose;
__PACKAGE__->meta->make_immutable;
1; # End of Bio::Mitomaster::SpeciesRef

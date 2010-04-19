use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'Bio::Mitomaster',
    license             => 'perl',
    dist_author         => 'Marty Brandon <mbrandon@uci.edu>',
    dist_version_from   => 'lib/Bio/Mitomaster.pm',
    build_requires => {
        'Test::More' => 0,
    },
    add_to_cleanup      => [ 'Bio-Mitomaster-*' ],
    create_makefile_pl => 'traditional',
);

$builder->create_build_script();

#!/home/ivan/bin/perl -s

use blib;
use Chemistry::File::SMILES;
use Chemistry::File::MDLMol;

#my $mol = Chemistry::Mol->parse($ARGV[0] || 'CCC', format => 'smiles',
    #kekulize => $kekulize);
my $mol = Chemistry::Mol->read($ARGV[0] || die "pls give a .mol\n");

printf "%s (%s)\n", $mol->name, $mol->formula;
$_->printf("%s H%h V%v\n") for $mol->atoms;

my $smiles = $mol->print(format => 'smiles', 
    unique => $unique, aromatic => $aromatic);

print "$smiles\n";



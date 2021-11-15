package Chemistry::File;
$VERSION = '0.35';

=head1 NAME

Chemistry::File - Molecule file I/O base class

=head1 SYNOPSIS

    # As a convenient interface for several mol readers:
    use Chemistry::File qw(PDB MDLMol); # load PDB and MDL modules
    
    # or try to use every file I/O module installed in the system:
    use Chemistry::File ':auto';

    my $mol1 = Chemistry::Mol->read("file.pdb");
    my $mol2 = Chemistry::Mol->read("file.mol");


    # as a base for a mol reader:

    package Chemistry::File::Myfile;
    use base qw(Chemistry::File);
    Chemistry::Mol->register_type("myfile", __PACKAGE__);

    # override the read_mol method
    sub read_mol {
        my ($self, $fh, %opts) = shift;
        my $mol_class = $opts{mol_class} || "Chemistry::Mol";
        my $mol = $mol_class->new;
        # ... do some stuff with $fh and $mol ...
        return $mol;
    }

    # override the write_mol method
    sub write_mol {
        my ($self, $fh, $mol, %opts) = shift;
        print $fh $mol->name, "\n"; 
        # ... do some stuff with $fh and $mol ...
    }

=head1 DESCRIPTION

The main use of this module is as a base class for other molecule file I/O
modules (for example, Chemistry::File::PDB). Such modules should override and
extend the Chemistry::File methods as needed. You only need to care about the
methods here if if you are writing a file I/O module or if you want a finer
degree of control than what is offered by the simple read and write methods
in the Chemistry::Mol class.

From the user's point of view, this module can also be used as shorthand
for using several Chemistry::File modules at the same time.

    use Chemistry::File qw(PDB MDLMol);

is exactly equivalent to

    use Chemistry::File::PDB;
    use Chemistry::File::MDLMol;

If you use the :auto keyword, Chemistry::File will autodetect and load
all the Chemistry::File::* modules installed in your system.

    use Chemistry::File ':auto';

=head1 FILE I/O MODEL

Before version 0.30, file I/O modules typically used only parse_string,
write_string, parse_file, and write_file, and they were generally used as class
methods. A file could contain one or more molecules and only be read or written
whole; reading it would return every molecule on the file. This was problematic
when dealing with large multi-molecule files (such as SDF files), because all
the molecules would have to be loaded into memory at the same time.

While version 0.30 retains backward compatibility with that simple model, it
also allows a more flexible interface that allows reading one molecule at a
time, skipping molecules, and reading and writing file-level information that
is not associated with specific molecules. The following diagram shows the
global structure of a file according to the new model:

    +-----------+
    | header    |
    +-----------+
    | molecule  |
    +-----------+
    | molecule  |
    +-----------+
    | ...       |
    +-----------+
    | footer    |
    +-----------+

In cases where the header and the footer are empty, the model reduces to the
pre-0.30 version. The low-level steps to read a file are the following:

    $file = Chemistry::File::MyFormat->new(file => 'xyz.mol');
    $file->open('<');
    $file->read_header;
    while (my $mol = $self->read_mol($file->fh, %opts)) {
        # do something with $mol...
    }
    $self->read_footer;

The C<read> method does all the above automatically, and it stores all the
molecules read in the mols property.

=head1 STANDARD OPTIONS

All the methods below include a list of options %opts at the end of the
parameter list. Each class implementing this interface may have its own
particular options. However, the following options should be recognized by all
classes:

=over

=item mol_class

A class or object with a C<new> method that constructs a molecule. This is 
needed when the user want to specify a molecule subclass different from the
default. When this option is not defined, the module may use Chemistry::Mol 
or whichever class is appropriate for that file format.

=item format

The name of the file format being used, as registered by
Chemistry::Mol->register_format.

=item fatal

If true, parsing errors should throw an exception; if false, they should just
try to recover if possible. True by default.

=back

=head1 CLASS METHODS

The class methods in this class (or rather, its derived classes) are usually
not called directly. Instead, use Chemistry::Mol->read, write, print, parse,
and file. These methods also work if called as instance methods.

=over


=cut

use strict;
use warnings;
no warnings qw(uninitialized);
use Carp;
use FileHandle;
use base qw(Chemistry::Obj);
# don't blame our problems in the Chemistry::Mol module ;-)
our @CARP_NOT = qw(Chemistry::Mol);

# This subroutine implements the :auto functionality
sub import {
    my $pack = shift;
    for my $param (@_){
        if ($param eq ':auto') {
            for my $pmfile (map {glob "$_/Chemistry/File/*.pm"} @INC) {
                my ($pm) = $pmfile =~ m|(Chemistry/File/.*\.pm)$|;
                #warn "requiring $pm\n";
                eval { require $pm }; 
                die "Error in Chemistry::File: '$@'; pmfile='$pmfile'; pm='$pm'\n" if $@;
            }
        } else {
            eval "use ${pack}::$param";
            die "$@" if $@;
        }
    }
}

=item $class->parse_string($s, %options)

Parse a string $s and return one or mole molecule objects. This is an abstract
method, so it should be provided by all derived classes.

=cut

sub parse_string {
    my ($self, $s, %opts) = @_;
    if ($opts{_must_override}) {
        my $class = ref $self || $self;
        croak "parse_string() is not implemented for $class";
    }
    $self->new(file => \$s, opts => \%opts)->read;
}


=item $class->write_string($mol, %options)

Convert a molecule to a string. This is an abstract method, so it should be
provided by all derived classes.

=cut

sub write_string {
    my ($self, $mol, %opts) = @_;
    if ($opts{_must_override}) {
        my $class = ref $self || $self;
        croak "write_string() is not implemented for $class";
    }
    my $s;
    $self->new(file => \$s, mols => [$mol], opts => \%opts)->write;
    $s;
}

=item $class->parse_file($file, %options)

Reads the file $file and returns one or more molecules. The default method
slurps the whole file and then calls parse_string, but derived classes may
choose to override it. $file can be a filehandle, a filename, or a scalar
reference. See C<new> for details.

=cut

sub parse_file {
    my ($self, $file, %opts) = @_;
    $self->new(file => $file, opts => \%opts)->read;
}

=item $class->write_file($mol, $file, %options)

Writes a file $file containing the molecule $mol. The default method calls
write_string first and then saves the string to a file, but derived classes
may choose to override it. $file can be either a filehandle or a filename.

=cut

sub write_file {
    my ($self, $mol, $file, %opts) = @_;

    $self->new(file => $file, mols => [$mol], opts => \%opts)->write;
}

=item $class->name_is($fname, %options)

Returns true if a filename is of the format corresponding to the class.
It should look at the filename only, because it may be called with
non-existent files. It is used to determine with which format to save a file.
For example, the Chemistry::File::PDB returns true if the file ends in .pdb.

=cut

sub name_is {
    0;
}

=item $class->string_is($s, %options)

Examines the string $s and returns true if it has the format of the class.

=cut

sub string_is {
    0;
}

=item $class->file_is($file, %options)

Examines the file $file and returns true if it has the format of the class.
The default method slurps the whole file and then calls string_is, but derived
classes may choose to override it.

=cut

sub file_is {
    my ($self, $file, %opts) = @_;
    
    my $s = eval {
        $self->open('<');
        $self->slurp;
    };
    if ($s) {
        $self->string_is($s, %opts);
    } elsif (! ref $file) {
        $self->name_is($file, %opts);
    }
}

=item $class->slurp($file %opts)

Reads a file into a scalar. Automatic decompression of gzipped files is
supported if the Compress::Zlib module is installed. Files ending in .gz are
assumed to be compressed; otherwise it is possible to force decompression by
passing the gzip => 1 option (or no decompression with gzip => 0).

=cut

# slurp a file into a scalar, with transparent decompression
sub slurp {
    my ($self) = @_;

    my $fh = $self->fh;
    local $/; 
    <$fh>;
}

=item $class->new(file => $file, opts => \%opts)

Create a new file object. This method is usually called indirectly via
the Chemistry::Mol->file method. $file may be a scalar with a filename, an
open filehandle, or a reference to a scalar. If a reference to a scalar is 
used, the string contained in the scalar is used as an in-memory file.

=cut

sub new {
    my $self = shift->SUPER::new(@_);
    $self->{opts}{fatal} = 1 unless exists $self->{opts}{fatal};
    $self;
}

Chemistry::Obj::accessor(qw(file fh opts mols mode));

=back

=head1 INSTANCE METHODS

=head2 Accessors

Chemistry::File objects are derived from Chemistry::Obj and have the same
properties (name, id, and type), as well as the following ones:

=over

=item file

The "file" as described above under C<new>.

=item fh

The filehandle used for reading and writing molecules. It is opened by C<open>.

=item opts

A hashref containing the options that are passed through to the old-style class
methods. They are also passed to the instance method to keep a similar
interface, but they could access them via $self->opts anyway.

=item mode

'>' if the file is open for writing, '<' for reading, and false if not open.

=item mols

C<read> stores all the molecules that were read in this property as an array
reference. C<write> gets the molecules to write from here.

=back

=head2 Abstract methods

These methods should be overridden, because they don't really do much by
default.

=over

=item $file->read_header

Read whatever information is available in the file before the first molecule.
Does nothing by default.

=cut

sub read_header { }

=item $file->read_footer

Read whatever information is available in the file after the last molecule.
Does nothing by default.

=cut

sub read_footer { }

=item $self->slurp_mol($fh)

Reads from the input string until the end of the current molecule and returns
the "slurped" string. It does not parse the string. It returns undefined if
there are no more molecules in the file. This method should be overridden if
needed; by default, it slurps until the end of the file.

=cut

sub slurp_mol {
    my ($self, $fh) = @_;
    local $/; <$fh>;
}

=item $self->skip_mol($fh)

Similar to slurp_mol, but it doesn't need to return anything except true or 
false. It should also be overridden if needed; by default, it just calls 
slurp_mol.

=cut

sub skip_mol { shift->slurp_mol(@_) }

=item $file->read_mol($fh, %opts)

Read the next molecule in the input stream. It returns false if there are no
more molecules in the file. This method should be overridden by derived
classes; otherwise it will call slurp_mol and parse_string (for backwards
compatibility; it is recommended to override read_mol directly in new modules).

Note: some old file I/O modules (written before the 0.30 interface) may return
more than one molecule anyway, so it is recommended to call read_mol in list
context to be safe:

    ($mol) = $file->read_mol($fh, %opts);

=cut

sub read_mol {
    my ($self, $fh, %opts) = @_;
    my $s = $self->slurp_mol($fh);
    return unless defined $s and length $s;
    $self->parse_string($s, %opts, _must_override => 1);
}
=item $file->write_header

Write whatever information is needed before the first molecule.
Does nothing by default.

=cut

sub write_header { }

=item $file->write_footer

Write whatever information is needed after the last molecule.
Does nothing by default.

=cut

sub write_footer { }

=item $self->write_mol($fh, $mol, %opts)

Write one molecule to $fh. By default and for backward compatibility, it just
calls C<write_string> and prints its return value to $self->fh. New classes
should override it.

=cut

sub write_mol {
    my ($self, $fh, $mol, %opts) = @_;
    print $fh $self->write_string($mol, %opts, _must_override => 1);
}

########################## OTHER ##################################

=back

=head2 Other methods

=over 

=item $self->open($mode) 

Opens the file (held in $self->file) for reading by default, or for writing if
$mode eq '>'. This method sets $self->fh transparently regardless of whether
$self->file is a filename (compressed or not), a scalar reference, or a
filehandle.

=cut

sub open {
    my ($self, $mode) = @_;
    my $fh;
    my $s;
    $mode ||= '<';
    $self->mode($mode);
    my $file = $self->file;
    croak "Chemistry::File::open: no file supplied" unless defined $file;
    if (ref $file eq 'SCALAR') {
        croak "decompression only supported for files" if $self->{opts}{gzip};
        if ($] >= 5.008) {
            open $fh, $mode, $file;
        } else { 
            require IO::String;
            $fh = IO::String->new($$file);
        }
    } elsif (ref $file) {
        croak "decompression only supported for files" if $self->{opts}{gzip};
        $fh = $file;
    } elsif ($self->{opts}{gzip} 
        or !defined $self->{opts}{gzip} and $file =~ /.gz$/) 
    {
        eval { require Compress::Zlib } # Carp
            or croak "Compress::Zlib not installed!";
        require File::Temp;

        $fh = File::Temp::tempfile();
        $self->{opts}{gzip} ||= 1;
        unless ($mode eq '>') { 
            my $gz = Compress::Zlib::gzopen($file, "rb") 
                 or croak "Cannot open compressed $file: "
                     . "$Compress::Zlib::gzerrno\n";

            my $buffer;
            print $fh $buffer while $gz->gzread($buffer) > 0;
        
            if ($Compress::Zlib::gzerrno != Compress::Zlib::Z_STREAM_END()) {
                croak "Error reading from $file: $Compress::Zlib::gzerrno"
                    . ($Compress::Zlib::gzerrno+0) . "\n";
            }
            $gz->gzclose();
            seek $fh, 0, 0;
        }
    } else {
        $fh = FileHandle->new("$mode$file") 
            or croak "Could not open file $file: $!";
    }
    $self->fh($fh);
    $self;
}

=item $self->close

Close the file. For regular files this just closes the filehandle, but for
gzipped files it does some additional postprocessing. This method is called
automatically on object destruction, so it is not mandatory to call it
explicitly.

=cut

sub close {
    my ($self) = @_;
    my $fh = $self->fh;
    if ($fh and $self->mode eq '>' and $self->{opts}{gzip}) {
        my $level = $self->{opts}{gzip} || 6;
        $level = 6 if $level == 1;
        my $file = $self->file;
        if (ref $file) { 
            croak "compression only supported for files";
        } else {
            seek $fh, 0, 0;
            my $gz = Compress::Zlib::gzopen($file, "wb$level")
                or croak "Cannot open $file $Compress::Zlib::gzerrno\n";
            local $_;
            while (<$fh>) {
                $gz->gzwrite($_) 
                    or croak "error writing: $Compress::Zlib::gzerrno\n";
            }
            $gz->gzclose;
        }
    }
    if ($self->mode) {
        if ($fh) { $fh->close or croak "$!" };
        $self->mode('');
    }
}

sub DESTROY { shift->close }

=item $file->read

Read the whole file. This calls open, read_header, read_mol until there are no
more molecules left, read_footer, and close. Returns a list of molecules if
called in list context, or the first molecule in scalar context.

=cut

sub read {
    my ($self) = @_;
    $self->open('<');
    $self->read_header;
    my @all_mols;
    while (my @mols = $self->read_mol($self->fh, %{$self->{opts}})) {
        push @all_mols, @mols;
    }
    $self->read_footer;
    $self->close;
    $self->mols(\@all_mols);
    wantarray ? @all_mols : $all_mols[0];
}

=item $self->write

Write all the molecules in $self->mols. It just calls open, write_header, 
write_mol (per each molecule), write_footer, and close.

=cut

sub write {
    my ($self) = @_;
    $self->open('>');
    $self->write_header;
    for my $mol (@{$self->mols}) {
        $self->write_mol($self->fh, $mol, %{$self->{opts}});
    }
    $self->write_footer;
    $self->close;
}

1;

=back

=head1 CAVEATS

The :auto feature may not be entirely portable, but it is known to work under
Unix and Windows (either Cygwin or Activestate).

=head1 VERSION

0.35

=head1 SEE ALSO

L<Chemistry::Mol>

The PerlMol website L<http://www.perlmol.org/>

=head1 AUTHOR

Ivan Tubert-Brohman-Brohman <itub@cpan.org>

=head1 COPYRIGHT

Copyright (c) 2005 Ivan Tubert-Brohman. All rights reserved. This program is
free software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut


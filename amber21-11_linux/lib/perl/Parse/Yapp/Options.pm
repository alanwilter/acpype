#
# Module Parse::Yapp::Options
#
# (c) Copyright 1999-2001 Francois Desarmenien, all rights reserved.
# (see the pod text in Parse::Yapp module for use and distribution rights)
#
package Parse::Yapp::Options;

use strict;
use Carp;

############################################################################
#Definitions of options
#
# %known_options    allowed options
#
# %default_options  default
#
# %actions          sub refs to execute if option is set with ($self,$value)
#                   as parameters
############################################################################
#
#A value of '' means any value can do
#
my(%known_options)= (
    language    =>  {
        perl    => "Ouput parser for Perl language",
# for future use...
#       'c++'   =>  "Output parser for C++ language",
#       c       =>  "Output parser for C language"
    },
    linenumbers =>  {
        0       =>  "Don't embbed line numbers in parser",
        1       =>  "Embbed source line numbers in parser"
    },
    inputfile   =>  {
        ''      =>  "Input file name: will automagically fills input"
    },
    classname   =>  {
        ''      =>  "Class name of parser object (Perl and C++)"
    },
    standalone  =>  {
        0       =>  "Don't create a standalone parser (Perl and C++)",
        1       =>  "Create a standalone parser"
    },
    input       =>  {
        ''      =>  "Input text of grammar"
    },
    template    => {
        ''      =>  "Template text for generating grammar file"
    },
);

my(%default_options)= (
    language => 'perl',
    linenumbers => 1,
    inputfile => undef,
    classname   => 'Parser',
    standalone => 0,
    input => undef,
    template => undef,
    shebang => undef,
);

my(%actions)= (
    inputfile => \&__LoadFile
);

#############################################################################
#
# Actions
#
# These are NOT a method, although they look like...
#
# They are super-private routines (that's why I prepend __ to their names)
#
#############################################################################
sub __LoadFile {
    my($self,$filename)=@_;

        open(IN,"<$filename")
    or  croak "Cannot open input file '$filename' for reading";
    $self->{OPTIONS}{input}=join('',<IN>);
    close(IN);
}

#############################################################################
#
# Private methods
#
#############################################################################

sub _SetOption {
    my($self)=shift;
    my($key,$value)=@_;

    $key=lc($key);

        @_ == 2
    or  croak "Invalid number of arguments";

        exists($known_options{$key})
    or  croak "Unknown option: '$key'";

    if(exists($known_options{$key}{lc($value)})) {
        $value=lc($value);
    }
    elsif(not exists($known_options{$key}{''})) {
        croak "Invalid value '$value' for option '$key'";
    }

        exists($actions{$key})
    and &{$actions{$key}}($self,$value);

    $self->{OPTIONS}{$key}=$value;
}

sub _GetOption {
    my($self)=shift;
    my($key)=map { lc($_) } @_;

        @_ == 1
    or  croak "Invalid number of arguments";

        exists($known_options{$key})
    or  croak "Unknown option: '$key'";

    $self->{OPTIONS}{$key};
}

#############################################################################
#
# Public methods
#
#############################################################################

#
# Constructor
#
sub new {
    my($class)=shift;
    my($self)={ OPTIONS => { %default_options } };

        ref($class)
    and $class=ref($class);
    
    bless($self,$class);

    $self->Options(@_);

    $self;
}

#
# Specify one or more options to set
#
sub Options {
    my($self)=shift;
    my($key,$value);

        @_ % 2 == 0
    or  croak "Invalid number of arguments";

    while(($key,$value)=splice(@_,0,2)) {
        $self->_SetOption($key,$value);
    }
}

#
# Set (2 parameters) or Get (1 parameter) values for one option
#
sub Option {
    my($self)=shift;
    my($key,$value)=@_;

        @_ == 1
    and return $self->_GetOption($key);

        @_ == 2
    and return $self->_SetOption($key,$value);

    croak "Invalid number of arguments";

}

1;

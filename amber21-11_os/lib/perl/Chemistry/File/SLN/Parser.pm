####################################################################
#
#    This file was generated using Parse::Yapp version 1.05.
#
#        Don't edit this file, use source file instead.
#
#             ANY CHANGE MADE HERE WILL BE LOST !
#
####################################################################
package Chemistry::File::SLN::Parser;
use vars qw ( @ISA );
use strict;

@ISA= qw ( Parse::Yapp::Driver );
use Parse::Yapp::Driver;

#line 3 "Parser.yp"
 no warnings 'uninitialized'; 

sub new {
        my($class)=shift;
        ref($class)
    and $class=ref($class);

    my($self)=$class->SUPER::new( yyversion => '1.05',
                                  yystates =>
[
	{#State 0
		ACTIONS => {
			'H' => 1,
			"<" => 4,
			'UC_NON_H' => 6
		},
		DEFAULT => -34,
		GOTOS => {
			'symbol' => 3,
			'UC' => 2,
			'atom' => 5,
			'ctab' => 8,
			'ctab_bracket' => 7
		}
	},
	{#State 1
		DEFAULT => -21
	},
	{#State 2
		DEFAULT => -18,
		GOTOS => {
			'lc_str' => 9
		}
	},
	{#State 3
		ACTIONS => {
			"[" => 11
		},
		DEFAULT => -28,
		GOTOS => {
			'atom_bracket' => 10
		}
	},
	{#State 4
		ACTIONS => {
			'H' => 1,
			"-" => 12,
			'LC' => 15,
			"+" => 16,
			'UC_NON_H' => 6
		},
		DEFAULT => -39,
		GOTOS => {
			'UC' => 13,
			'charge' => 20,
			'key_val' => 14,
			'attr_list' => 17,
			'ALPHA' => 21,
			'attr' => 18,
			'key' => 19
		}
	},
	{#State 5
		ACTIONS => {
			'' => -3,
			"-" => 23,
			":" => 22,
			"<" => -3,
			"=" => 24,
			"(" => 26,
			"#" => 29,
			"." => 30
		},
		DEFAULT => -8,
		GOTOS => {
			'bond' => 25,
			'chain' => 28,
			'bondatom' => 27,
			'bond_symbol' => 31
		}
	},
	{#State 6
		DEFAULT => -20
	},
	{#State 7
		DEFAULT => -1
	},
	{#State 8
		ACTIONS => {
			'' => 32
		}
	},
	{#State 9
		ACTIONS => {
			'LC' => 33
		},
		DEFAULT => -17
	},
	{#State 10
		ACTIONS => {
			'H' => 34
		},
		DEFAULT => -36,
		GOTOS => {
			'hcount' => 35
		}
	},
	{#State 11
		ACTIONS => {
			'NUM' => 37,
			'H' => 1,
			"-" => 12,
			'LC' => 15,
			"+" => 16,
			'UC_NON_H' => 6
		},
		DEFAULT => -39,
		GOTOS => {
			'UC' => 13,
			'charge' => 20,
			'key_val' => 14,
			'attr_list' => 17,
			'ALPHA' => 21,
			'attr' => 36,
			'key' => 19
		}
	},
	{#State 12
		ACTIONS => {
			'NUM' => 38
		},
		DEFAULT => -54
	},
	{#State 13
		DEFAULT => -22
	},
	{#State 14
		DEFAULT => -41
	},
	{#State 15
		DEFAULT => -23
	},
	{#State 16
		ACTIONS => {
			'NUM' => 39
		},
		DEFAULT => -53
	},
	{#State 17
		ACTIONS => {
			";" => 40
		},
		DEFAULT => -40
	},
	{#State 18
		ACTIONS => {
			">" => 41
		}
	},
	{#State 19
		ACTIONS => {
			"=" => 42
		},
		DEFAULT => -44
	},
	{#State 20
		DEFAULT => -45
	},
	{#State 21
		DEFAULT => -49,
		GOTOS => {
			'key_tail' => 43
		}
	},
	{#State 22
		DEFAULT => -15
	},
	{#State 23
		DEFAULT => -11
	},
	{#State 24
		DEFAULT => -12
	},
	{#State 25
		ACTIONS => {
			'H' => 1,
			"\@" => 44,
			'UC_NON_H' => 6
		},
		GOTOS => {
			'symbol' => 3,
			'UC' => 2,
			'atom' => 46,
			'closure' => 45
		}
	},
	{#State 26
		ACTIONS => {
			"-" => 23,
			":" => 22,
			"=" => 24,
			"#" => 29,
			"." => 30
		},
		DEFAULT => -8,
		GOTOS => {
			'bond' => 25,
			'bondatom' => 47,
			'bond_symbol' => 31
		}
	},
	{#State 27
		ACTIONS => {
			"-" => 23,
			":" => 22,
			"\@" => -8,
			'UC_NON_H' => -8,
			"=" => 24,
			'H' => -8,
			"(" => 26,
			"#" => 29,
			"." => 30
		},
		DEFAULT => -3,
		GOTOS => {
			'bond' => 25,
			'chain' => 48,
			'bondatom' => 27,
			'bond_symbol' => 31
		}
	},
	{#State 28
		ACTIONS => {
			"<" => 4
		},
		DEFAULT => -34,
		GOTOS => {
			'ctab_bracket' => 49
		}
	},
	{#State 29
		DEFAULT => -13
	},
	{#State 30
		DEFAULT => -14
	},
	{#State 31
		ACTIONS => {
			"[" => 51
		},
		DEFAULT => -32,
		GOTOS => {
			'bond_bracket' => 50
		}
	},
	{#State 32
		DEFAULT => 0
	},
	{#State 33
		DEFAULT => -19
	},
	{#State 34
		ACTIONS => {
			'NUM' => 52
		},
		DEFAULT => -38
	},
	{#State 35
		DEFAULT => -16
	},
	{#State 36
		ACTIONS => {
			"]" => 53
		}
	},
	{#State 37
		ACTIONS => {
			":" => 54,
			"]" => 55
		}
	},
	{#State 38
		DEFAULT => -52
	},
	{#State 39
		DEFAULT => -51
	},
	{#State 40
		ACTIONS => {
			'H' => 1,
			"-" => 12,
			'LC' => 15,
			"+" => 16,
			'UC_NON_H' => 6
		},
		GOTOS => {
			'UC' => 13,
			'charge' => 20,
			'key_val' => 56,
			'ALPHA' => 21,
			'key' => 19
		}
	},
	{#State 41
		DEFAULT => -35
	},
	{#State 42
		ACTIONS => {
			'STRING' => 59
		},
		DEFAULT => -55,
		GOTOS => {
			'value' => 57,
			'string' => 58
		}
	},
	{#State 43
		ACTIONS => {
			'LC' => 15,
			'UC_NON_H' => 6,
			"_" => 62,
			'NUM' => 63,
			'H' => 1
		},
		DEFAULT => -48,
		GOTOS => {
			'WORD' => 61,
			'ALNUM' => 60,
			'UC' => 13,
			'ALPHA' => 64
		}
	},
	{#State 44
		ACTIONS => {
			'NUM' => 65
		}
	},
	{#State 45
		DEFAULT => -7
	},
	{#State 46
		DEFAULT => -6
	},
	{#State 47
		ACTIONS => {
			"-" => 23,
			":" => 22,
			"=" => 24,
			"(" => 26,
			"#" => 29,
			"." => 30,
			")" => -3
		},
		DEFAULT => -8,
		GOTOS => {
			'bond' => 25,
			'chain' => 66,
			'bondatom' => 27,
			'bond_symbol' => 31
		}
	},
	{#State 48
		DEFAULT => -4
	},
	{#State 49
		DEFAULT => -2
	},
	{#State 50
		DEFAULT => -9
	},
	{#State 51
		ACTIONS => {
			'H' => 1,
			"-" => 12,
			'LC' => 15,
			"+" => 16,
			'UC_NON_H' => 6
		},
		DEFAULT => -39,
		GOTOS => {
			'UC' => 13,
			'charge' => 20,
			'key_val' => 14,
			'attr_list' => 17,
			'ALPHA' => 21,
			'attr' => 67,
			'key' => 19
		}
	},
	{#State 52
		DEFAULT => -37
	},
	{#State 53
		DEFAULT => -30
	},
	{#State 54
		ACTIONS => {
			'H' => 1,
			"-" => 12,
			'LC' => 15,
			"+" => 16,
			'UC_NON_H' => 6
		},
		DEFAULT => -39,
		GOTOS => {
			'UC' => 13,
			'charge' => 20,
			'key_val' => 14,
			'attr_list' => 17,
			'ALPHA' => 21,
			'attr' => 68,
			'key' => 19
		}
	},
	{#State 55
		DEFAULT => -31
	},
	{#State 56
		DEFAULT => -42
	},
	{#State 57
		DEFAULT => -43
	},
	{#State 58
		ACTIONS => {
			"-" => 69,
			'LC' => 15,
			"+" => 71,
			'UC_NON_H' => 6,
			"_" => 62,
			'NUM' => 63,
			'H' => 1,
			'OTHER_CHAR' => 74,
			"." => 73
		},
		DEFAULT => -47,
		GOTOS => {
			'ALNUM' => 60,
			'WORD' => 70,
			'UC' => 13,
			'ALPHA' => 64,
			'string_char' => 72
		}
	},
	{#State 59
		DEFAULT => -46
	},
	{#State 60
		DEFAULT => -26
	},
	{#State 61
		DEFAULT => -50
	},
	{#State 62
		DEFAULT => -27
	},
	{#State 63
		DEFAULT => -25
	},
	{#State 64
		DEFAULT => -24
	},
	{#State 65
		DEFAULT => -10
	},
	{#State 66
		ACTIONS => {
			")" => 75
		}
	},
	{#State 67
		ACTIONS => {
			"]" => 76
		}
	},
	{#State 68
		ACTIONS => {
			"]" => 77
		}
	},
	{#State 69
		DEFAULT => -60
	},
	{#State 70
		DEFAULT => -57
	},
	{#State 71
		DEFAULT => -59
	},
	{#State 72
		DEFAULT => -56
	},
	{#State 73
		DEFAULT => -61
	},
	{#State 74
		DEFAULT => -58
	},
	{#State 75
		ACTIONS => {
			"-" => 23,
			":" => 22,
			"\@" => -8,
			'UC_NON_H' => -8,
			"=" => 24,
			'H' => -8,
			"(" => 26,
			"#" => 29,
			"." => 30
		},
		DEFAULT => -3,
		GOTOS => {
			'bond' => 25,
			'chain' => 78,
			'bondatom' => 27,
			'bond_symbol' => 31
		}
	},
	{#State 76
		DEFAULT => -33
	},
	{#State 77
		DEFAULT => -29
	},
	{#State 78
		DEFAULT => -5
	}
],
                                  yyrules  =>
[
	[#Rule 0
		 '$start', 2, undef
	],
	[#Rule 1
		 'ctab', 1,
sub
#line 10 "Parser.yp"
{ +{ chain => [], attr => $_[1] } }
	],
	[#Rule 2
		 'ctab', 3,
sub
#line 11 "Parser.yp"
{ +{ chain => [ $_[1], @{$_[2]} ],
                                    attr => $_[3] } }
	],
	[#Rule 3
		 'chain', 0,
sub
#line 14 "Parser.yp"
{ [] }
	],
	[#Rule 4
		 'chain', 2,
sub
#line 15 "Parser.yp"
{ [ @{$_[1]}, @{$_[2]}, ] }
	],
	[#Rule 5
		 'chain', 5,
sub
#line 16 "Parser.yp"
{ 
            [ $_[1], @{$_[2]}, @{$_[3]}, $_[4], @{$_[5]} ] }
	],
	[#Rule 6
		 'bondatom', 2,
sub
#line 20 "Parser.yp"
{ [@_[1,2]]}
	],
	[#Rule 7
		 'bondatom', 2,
sub
#line 21 "Parser.yp"
{ [@_[1,2]]}
	],
	[#Rule 8
		 'bond', 0,
sub
#line 24 "Parser.yp"
{ +{type => '-'} }
	],
	[#Rule 9
		 'bond', 2,
sub
#line 25 "Parser.yp"
{ +{
                                        type => $_[1],
                                        attr => $_[2]
                                        } 
                                    }
	],
	[#Rule 10
		 'closure', 2,
sub
#line 32 "Parser.yp"
{ +{closure => $_[2]} }
	],
	[#Rule 11
		 'bond_symbol', 1, undef
	],
	[#Rule 12
		 'bond_symbol', 1, undef
	],
	[#Rule 13
		 'bond_symbol', 1, undef
	],
	[#Rule 14
		 'bond_symbol', 1, undef
	],
	[#Rule 15
		 'bond_symbol', 1, undef
	],
	[#Rule 16
		 'atom', 3,
sub
#line 37 "Parser.yp"
{ +{
                                        symbol  => $_[1],
                                        id      => $_[2][0],         
                                        hcount  => $_[3],
                                        attr   => $_[2][1],
                                    } 
                                }
	],
	[#Rule 17
		 'symbol', 2,
sub
#line 46 "Parser.yp"
{ $_[1] . $_[2] }
	],
	[#Rule 18
		 'lc_str', 0, undef
	],
	[#Rule 19
		 'lc_str', 2,
sub
#line 50 "Parser.yp"
{ $_[1] . $_[2] }
	],
	[#Rule 20
		 'UC', 1, undef
	],
	[#Rule 21
		 'UC', 1, undef
	],
	[#Rule 22
		 'ALPHA', 1, undef
	],
	[#Rule 23
		 'ALPHA', 1, undef
	],
	[#Rule 24
		 'ALNUM', 1, undef
	],
	[#Rule 25
		 'ALNUM', 1, undef
	],
	[#Rule 26
		 'WORD', 1, undef
	],
	[#Rule 27
		 'WORD', 1, undef
	],
	[#Rule 28
		 'atom_bracket', 0, undef
	],
	[#Rule 29
		 'atom_bracket', 5,
sub
#line 59 "Parser.yp"
{ [@_[2,4]] }
	],
	[#Rule 30
		 'atom_bracket', 3,
sub
#line 60 "Parser.yp"
{ [undef, $_[2]] }
	],
	[#Rule 31
		 'atom_bracket', 3,
sub
#line 61 "Parser.yp"
{ [$_[2], undef] }
	],
	[#Rule 32
		 'bond_bracket', 0, undef
	],
	[#Rule 33
		 'bond_bracket', 3,
sub
#line 65 "Parser.yp"
{ $_[2] }
	],
	[#Rule 34
		 'ctab_bracket', 0, undef
	],
	[#Rule 35
		 'ctab_bracket', 3,
sub
#line 69 "Parser.yp"
{ $_[2] }
	],
	[#Rule 36
		 'hcount', 0, undef
	],
	[#Rule 37
		 'hcount', 2,
sub
#line 73 "Parser.yp"
{ $_[2] }
	],
	[#Rule 38
		 'hcount', 1,
sub
#line 74 "Parser.yp"
{ 1 }
	],
	[#Rule 39
		 'attr', 0, undef
	],
	[#Rule 40
		 'attr', 1, undef
	],
	[#Rule 41
		 'attr_list', 1, undef
	],
	[#Rule 42
		 'attr_list', 3,
sub
#line 81 "Parser.yp"
{ +{%{$_[1]}, %{$_[3]}} }
	],
	[#Rule 43
		 'key_val', 3,
sub
#line 84 "Parser.yp"
{ +{lc($_[1]) => $_[3]} }
	],
	[#Rule 44
		 'key_val', 1,
sub
#line 85 "Parser.yp"
{ +{lc($_[1]) => 'TRUE' } }
	],
	[#Rule 45
		 'key_val', 1,
sub
#line 86 "Parser.yp"
{ +{charge => $_[1] } }
	],
	[#Rule 46
		 'value', 1, undef
	],
	[#Rule 47
		 'value', 1, undef
	],
	[#Rule 48
		 'key', 2,
sub
#line 91 "Parser.yp"
{ $_[1] . $_[2] }
	],
	[#Rule 49
		 'key_tail', 0, undef
	],
	[#Rule 50
		 'key_tail', 2,
sub
#line 94 "Parser.yp"
{ $_[1] . $_[2] }
	],
	[#Rule 51
		 'charge', 2,
sub
#line 97 "Parser.yp"
{ $_[2] }
	],
	[#Rule 52
		 'charge', 2,
sub
#line 98 "Parser.yp"
{ -$_[2] }
	],
	[#Rule 53
		 'charge', 1,
sub
#line 99 "Parser.yp"
{ 1 }
	],
	[#Rule 54
		 'charge', 1,
sub
#line 100 "Parser.yp"
{ -1 }
	],
	[#Rule 55
		 'string', 0, undef
	],
	[#Rule 56
		 'string', 2,
sub
#line 104 "Parser.yp"
{ $_[1] . $_[2] }
	],
	[#Rule 57
		 'string_char', 1, undef
	],
	[#Rule 58
		 'string_char', 1, undef
	],
	[#Rule 59
		 'string_char', 1, undef
	],
	[#Rule 60
		 'string_char', 1, undef
	],
	[#Rule 61
		 'string_char', 1, undef
	]
],
                                  @_);
    bless($self,$class);
}

#line 112 "Parser.yp"


sub _Error {
        exists $_[0]->YYData->{ERRMSG}
    and do {
        warn $_[0]->YYData->{ERRMSG};
        delete $_[0]->YYData->{ERRMSG};
        return;
    };
    warn "Syntax error.\n";
}

sub _Lexer {
    my($parser)=shift;

        $parser->YYData->{INPUT}
    #or  $parser->YYData->{INPUT} = <STDIN>
    or  return('',undef);

    $parser->YYData->{INPUT}=~s/^[ \t]//;

    for ($parser->YYData->{INPUT}) {
        s/^([0-9]+(?:\.[0-9]+)?)//
                and return('NUM',$1);
        s/^(H)//
                and return('H',$1);
        s/^([A-Z])//
                and return('UC_NON_H',$1);
        s/^([a-z])//
                and return('LC',$1);
        s/^"(.*?)"//s
                and return('STRING',$1);
        s/^([<>[\];=_()\@#.:+-])//s 
                and return($1,$1);   # "special" character
        s/^(.)//s
                and return('OTHER_CHAR',$1);
    }
}

sub run {
    my($self)=shift;
    $self->YYData->{INPUT} = shift;
    $self->YYParse( yylex => \&_Lexer, yyerror => \&_Error );
}


1;

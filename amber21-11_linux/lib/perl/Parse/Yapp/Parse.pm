####################################################################
#
#    This file was generated using Parse::Yapp version 1.05.
#
#        Don't edit this file, use source file instead.
#
#             ANY CHANGE MADE HERE WILL BE LOST !
#
####################################################################
package Parse::Yapp::Parse;
use vars qw ( @ISA );
use strict;

@ISA= qw ( Parse::Yapp::Driver );
use Parse::Yapp::Driver;

#line 1 "YappParse.yp"
# (c) Copyright Francois Desarmenien 1998-2001, all rights reserved.
# (see COPYRIGHT in Parse::Yapp.pm pod section for use and distribution rights)
#
# Parse/Yapp/Parser.yp: Parse::Yapp::Parser.pm source file
#
# Use: yapp -m 'Parse::Yapp::Parse' -o Parse/Yapp/Parse.pm YappParse.yp
#
# to generate the Parser module.
# 
#line 12 "YappParse.yp"

require 5.004;

use Carp;

my($input,$lexlevel,@lineno,$nberr,$prec,$labelno);
my($syms,$head,$tail,$token,$term,$nterm,$rules,$precterm,$start,$nullable);
my($expect);



sub new {
        my($class)=shift;
        ref($class)
    and $class=ref($class);

    my($self)=$class->SUPER::new( yyversion => '1.05',
                                  yystates =>
[
	{#State 0
		ACTIONS => {
			"%%" => -6,
			'HEADCODE' => 3,
			'UNION' => 2,
			'TOKEN' => 5,
			'ASSOC' => 7,
			'START' => 6,
			'error' => 9,
			'TYPE' => 10,
			"\n" => 11,
			'EXPECT' => 13
		},
		GOTOS => {
			'head' => 1,
			'decls' => 12,
			'yapp' => 4,
			'decl' => 14,
			'headsec' => 8
		}
	},
	{#State 1
		ACTIONS => {
			'error' => 19,
			"%%" => 16,
			'IDENT' => 18
		},
		GOTOS => {
			'rules' => 15,
			'rulesec' => 20,
			'body' => 17
		}
	},
	{#State 2
		ACTIONS => {
			'CODE' => 21
		}
	},
	{#State 3
		ACTIONS => {
			"\n" => 22
		}
	},
	{#State 4
		ACTIONS => {
			'' => 23
		}
	},
	{#State 5
		ACTIONS => {
			"<" => 25
		},
		DEFAULT => -19,
		GOTOS => {
			'typedecl' => 24
		}
	},
	{#State 6
		ACTIONS => {
			'IDENT' => 26
		},
		GOTOS => {
			'ident' => 27
		}
	},
	{#State 7
		ACTIONS => {
			"<" => 25
		},
		DEFAULT => -19,
		GOTOS => {
			'typedecl' => 28
		}
	},
	{#State 8
		ACTIONS => {
			"%%" => 29
		}
	},
	{#State 9
		ACTIONS => {
			"\n" => 30
		}
	},
	{#State 10
		ACTIONS => {
			"<" => 25
		},
		DEFAULT => -19,
		GOTOS => {
			'typedecl' => 31
		}
	},
	{#State 11
		DEFAULT => -10
	},
	{#State 12
		ACTIONS => {
			"%%" => -7,
			'HEADCODE' => 3,
			'UNION' => 2,
			'TOKEN' => 5,
			'ASSOC' => 7,
			'START' => 6,
			'error' => 9,
			'TYPE' => 10,
			"\n" => 11,
			'EXPECT' => 13
		},
		GOTOS => {
			'decl' => 32
		}
	},
	{#State 13
		ACTIONS => {
			'NUMBER' => 33
		}
	},
	{#State 14
		DEFAULT => -9
	},
	{#State 15
		DEFAULT => -28
	},
	{#State 16
		DEFAULT => -26
	},
	{#State 17
		ACTIONS => {
			'TAILCODE' => 34
		},
		DEFAULT => -45,
		GOTOS => {
			'tail' => 35
		}
	},
	{#State 18
		ACTIONS => {
			":" => 36
		}
	},
	{#State 19
		ACTIONS => {
			";" => 37
		}
	},
	{#State 20
		ACTIONS => {
			'error' => 19,
			"%%" => 39,
			'IDENT' => 18
		},
		GOTOS => {
			'rules' => 38
		}
	},
	{#State 21
		ACTIONS => {
			"\n" => 40
		}
	},
	{#State 22
		DEFAULT => -14
	},
	{#State 23
		DEFAULT => -0
	},
	{#State 24
		ACTIONS => {
			'LITERAL' => 41,
			'IDENT' => 26
		},
		GOTOS => {
			'symlist' => 43,
			'ident' => 44,
			'symbol' => 42
		}
	},
	{#State 25
		ACTIONS => {
			'IDENT' => 45
		}
	},
	{#State 26
		DEFAULT => -4
	},
	{#State 27
		ACTIONS => {
			"\n" => 46
		}
	},
	{#State 28
		ACTIONS => {
			'LITERAL' => 41,
			'IDENT' => 26
		},
		GOTOS => {
			'symlist' => 47,
			'ident' => 44,
			'symbol' => 42
		}
	},
	{#State 29
		DEFAULT => -5
	},
	{#State 30
		DEFAULT => -18
	},
	{#State 31
		ACTIONS => {
			'IDENT' => 26
		},
		GOTOS => {
			'ident' => 48,
			'identlist' => 49
		}
	},
	{#State 32
		DEFAULT => -8
	},
	{#State 33
		ACTIONS => {
			"\n" => 50
		}
	},
	{#State 34
		DEFAULT => -46
	},
	{#State 35
		DEFAULT => -1
	},
	{#State 36
		ACTIONS => {
			'CODE' => 57,
			'LITERAL' => 41,
			'IDENT' => 26
		},
		DEFAULT => -35,
		GOTOS => {
			'rhselts' => 56,
			'rule' => 51,
			'code' => 52,
			'rhs' => 53,
			'ident' => 44,
			'rhselt' => 58,
			'rhss' => 55,
			'symbol' => 54
		}
	},
	{#State 37
		DEFAULT => -30
	},
	{#State 38
		DEFAULT => -27
	},
	{#State 39
		DEFAULT => -25
	},
	{#State 40
		DEFAULT => -15
	},
	{#State 41
		DEFAULT => -2
	},
	{#State 42
		DEFAULT => -22
	},
	{#State 43
		ACTIONS => {
			"\n" => 60,
			'LITERAL' => 41,
			'IDENT' => 26
		},
		GOTOS => {
			'ident' => 44,
			'symbol' => 59
		}
	},
	{#State 44
		DEFAULT => -3
	},
	{#State 45
		ACTIONS => {
			">" => 61
		}
	},
	{#State 46
		DEFAULT => -13
	},
	{#State 47
		ACTIONS => {
			"\n" => 62,
			'LITERAL' => 41,
			'IDENT' => 26
		},
		GOTOS => {
			'ident' => 44,
			'symbol' => 59
		}
	},
	{#State 48
		DEFAULT => -24
	},
	{#State 49
		ACTIONS => {
			"\n" => 63,
			'IDENT' => 26
		},
		GOTOS => {
			'ident' => 64
		}
	},
	{#State 50
		DEFAULT => -17
	},
	{#State 51
		DEFAULT => -32
	},
	{#State 52
		DEFAULT => -40
	},
	{#State 53
		ACTIONS => {
			'PREC' => 66
		},
		DEFAULT => -34,
		GOTOS => {
			'prec' => 65
		}
	},
	{#State 54
		DEFAULT => -39
	},
	{#State 55
		ACTIONS => {
			"|" => 68,
			";" => 67
		}
	},
	{#State 56
		ACTIONS => {
			'CODE' => 57,
			'LITERAL' => 41,
			'IDENT' => 26
		},
		DEFAULT => -36,
		GOTOS => {
			'code' => 52,
			'ident' => 44,
			'rhselt' => 69,
			'symbol' => 54
		}
	},
	{#State 57
		DEFAULT => -44
	},
	{#State 58
		DEFAULT => -38
	},
	{#State 59
		DEFAULT => -21
	},
	{#State 60
		DEFAULT => -11
	},
	{#State 61
		DEFAULT => -20
	},
	{#State 62
		DEFAULT => -12
	},
	{#State 63
		DEFAULT => -16
	},
	{#State 64
		DEFAULT => -23
	},
	{#State 65
		ACTIONS => {
			'CODE' => 57
		},
		DEFAULT => -42,
		GOTOS => {
			'code' => 70,
			'epscode' => 71
		}
	},
	{#State 66
		ACTIONS => {
			'LITERAL' => 41,
			'IDENT' => 26
		},
		GOTOS => {
			'ident' => 44,
			'symbol' => 72
		}
	},
	{#State 67
		DEFAULT => -29
	},
	{#State 68
		ACTIONS => {
			'CODE' => 57,
			'LITERAL' => 41,
			'IDENT' => 26
		},
		DEFAULT => -35,
		GOTOS => {
			'rhselts' => 56,
			'rule' => 73,
			'code' => 52,
			'rhs' => 53,
			'ident' => 44,
			'rhselt' => 58,
			'symbol' => 54
		}
	},
	{#State 69
		DEFAULT => -37
	},
	{#State 70
		DEFAULT => -43
	},
	{#State 71
		DEFAULT => -33
	},
	{#State 72
		DEFAULT => -41
	},
	{#State 73
		DEFAULT => -31
	}
],
                                  yyrules  =>
[
	[#Rule 0
		 '$start', 2, undef
	],
	[#Rule 1
		 'yapp', 3, undef
	],
	[#Rule 2
		 'symbol', 1,
sub
#line 30 "YappParse.yp"
{
                        exists($$syms{$_[1][0]})
                    or  do {
                        $$syms{$_[1][0]} = $_[1][1];
                        $$term{$_[1][0]} = undef;
                    };
                    $_[1]
                }
	],
	[#Rule 3
		 'symbol', 1, undef
	],
	[#Rule 4
		 'ident', 1,
sub
#line 41 "YappParse.yp"
{
                        exists($$syms{$_[1][0]})
                    or  do {
                        $$syms{$_[1][0]} = $_[1][1];
                        $$term{$_[1][0]} = undef;
                    };
                    $_[1]
                }
	],
	[#Rule 5
		 'head', 2, undef
	],
	[#Rule 6
		 'headsec', 0, undef
	],
	[#Rule 7
		 'headsec', 1, undef
	],
	[#Rule 8
		 'decls', 2, undef
	],
	[#Rule 9
		 'decls', 1, undef
	],
	[#Rule 10
		 'decl', 1, undef
	],
	[#Rule 11
		 'decl', 4,
sub
#line 66 "YappParse.yp"
{
                for (@{$_[3]}) {
                    my($symbol,$lineno)=@$_;

                        exists($$token{$symbol})
                    and do {
                        _SyntaxError(0,
                                "Token $symbol redefined: ".
                                "Previously defined line $$syms{$symbol}",
                                $lineno);
                        next;
                    };
                    $$token{$symbol}=$lineno;
                    $$term{$symbol} = [ ];
                }
                undef
            }
	],
	[#Rule 12
		 'decl', 4,
sub
#line 84 "YappParse.yp"
{
                for (@{$_[3]}) {
                    my($symbol,$lineno)=@$_;

                        defined($$term{$symbol}[0])
                    and do {
                        _SyntaxError(1,
                            "Precedence for symbol $symbol redefined: ".
                            "Previously defined line $$syms{$symbol}",
                            $lineno);
                        next;
                    };
                    $$token{$symbol}=$lineno;
                    $$term{$symbol} = [ $_[1][0], $prec ];
                }
                ++$prec;
                undef
            }
	],
	[#Rule 13
		 'decl', 3,
sub
#line 102 "YappParse.yp"
{ $start=$_[2][0]; undef }
	],
	[#Rule 14
		 'decl', 2,
sub
#line 103 "YappParse.yp"
{ push(@$head,$_[1]); undef }
	],
	[#Rule 15
		 'decl', 3,
sub
#line 104 "YappParse.yp"
{ undef }
	],
	[#Rule 16
		 'decl', 4,
sub
#line 106 "YappParse.yp"
{
                for ( @{$_[3]} ) {
                    my($symbol,$lineno)=@$_;

                        exists($$nterm{$symbol})
                    and do {
                        _SyntaxError(0,
                                "Non-terminal $symbol redefined: ".
                                "Previously defined line $$syms{$symbol}",
                                $lineno);
                        next;
                    };
                    delete($$term{$symbol});   #not a terminal
                    $$nterm{$symbol}=undef;    #is a non-terminal
                }
            }
	],
	[#Rule 17
		 'decl', 3,
sub
#line 122 "YappParse.yp"
{ $expect=$_[2][0]; undef }
	],
	[#Rule 18
		 'decl', 2,
sub
#line 123 "YappParse.yp"
{ $_[0]->YYErrok }
	],
	[#Rule 19
		 'typedecl', 0, undef
	],
	[#Rule 20
		 'typedecl', 3, undef
	],
	[#Rule 21
		 'symlist', 2,
sub
#line 130 "YappParse.yp"
{ push(@{$_[1]},$_[2]); $_[1] }
	],
	[#Rule 22
		 'symlist', 1,
sub
#line 131 "YappParse.yp"
{ [ $_[1] ] }
	],
	[#Rule 23
		 'identlist', 2,
sub
#line 134 "YappParse.yp"
{ push(@{$_[1]},$_[2]); $_[1] }
	],
	[#Rule 24
		 'identlist', 1,
sub
#line 135 "YappParse.yp"
{ [ $_[1] ] }
	],
	[#Rule 25
		 'body', 2,
sub
#line 140 "YappParse.yp"
{
                    $start
                or  $start=$$rules[1][0];

                    ref($$nterm{$start})
                or  _SyntaxError(2,"Start symbol $start not found ".
                                   "in rules section",$_[2][1]);

                $$rules[0]=[ '$start', [ $start, chr(0) ], undef, undef ];
            }
	],
	[#Rule 26
		 'body', 1,
sub
#line 150 "YappParse.yp"
{ _SyntaxError(2,"No rules in input grammar",$_[1][1]); }
	],
	[#Rule 27
		 'rulesec', 2, undef
	],
	[#Rule 28
		 'rulesec', 1, undef
	],
	[#Rule 29
		 'rules', 4,
sub
#line 157 "YappParse.yp"
{ _AddRules($_[1],$_[3]); undef }
	],
	[#Rule 30
		 'rules', 2,
sub
#line 158 "YappParse.yp"
{ $_[0]->YYErrok }
	],
	[#Rule 31
		 'rhss', 3,
sub
#line 161 "YappParse.yp"
{ push(@{$_[1]},$_[3]); $_[1] }
	],
	[#Rule 32
		 'rhss', 1,
sub
#line 162 "YappParse.yp"
{ [ $_[1] ] }
	],
	[#Rule 33
		 'rule', 3,
sub
#line 165 "YappParse.yp"
{ push(@{$_[1]}, $_[2], $_[3]); $_[1] }
	],
	[#Rule 34
		 'rule', 1,
sub
#line 166 "YappParse.yp"
{
                                my($code)=undef;

                                    defined($_[1])
                                and $_[1][-1][0] eq 'CODE'
                                and $code = ${pop(@{$_[1]})}[1];

                                push(@{$_[1]}, undef, $code);

                                $_[1]
                            }
	],
	[#Rule 35
		 'rhs', 0, undef
	],
	[#Rule 36
		 'rhs', 1, undef
	],
	[#Rule 37
		 'rhselts', 2,
sub
#line 183 "YappParse.yp"
{ push(@{$_[1]},$_[2]); $_[1] }
	],
	[#Rule 38
		 'rhselts', 1,
sub
#line 184 "YappParse.yp"
{ [ $_[1] ] }
	],
	[#Rule 39
		 'rhselt', 1,
sub
#line 187 "YappParse.yp"
{ [ 'SYMB', $_[1] ] }
	],
	[#Rule 40
		 'rhselt', 1,
sub
#line 188 "YappParse.yp"
{ [ 'CODE', $_[1] ] }
	],
	[#Rule 41
		 'prec', 2,
sub
#line 192 "YappParse.yp"
{
                       	defined($$term{$_[2][0]})
                    or  do {
                        _SyntaxError(1,"No precedence for symbol $_[2][0]",
                                         $_[2][1]);
                        return undef;
                    };

                    ++$$precterm{$_[2][0]};
                    $$term{$_[2][0]}[1];
				}
	],
	[#Rule 42
		 'epscode', 0,
sub
#line 205 "YappParse.yp"
{ undef }
	],
	[#Rule 43
		 'epscode', 1,
sub
#line 206 "YappParse.yp"
{ $_[1] }
	],
	[#Rule 44
		 'code', 1,
sub
#line 209 "YappParse.yp"
{ $_[1] }
	],
	[#Rule 45
		 'tail', 0, undef
	],
	[#Rule 46
		 'tail', 1,
sub
#line 215 "YappParse.yp"
{ $tail=$_[1] }
	]
],
                                  @_);
    bless($self,$class);
}

#line 218 "YappParse.yp"

sub _Error {
    my($value)=$_[0]->YYCurval;

    my($what)= $token ? "input: '$$value[0]'" : "end of input";

    _SyntaxError(1,"Unexpected $what",$$value[1]);
}

sub _Lexer {
 
    #At EOF
        pos($$input) >= length($$input)
    and return('',[ undef, -1 ]);

    #In TAIL section
        $lexlevel > 1
    and do {
        my($pos)=pos($$input);

        $lineno[0]=$lineno[1];
        $lineno[1]=-1;
        pos($$input)=length($$input);
        return('TAILCODE',[ substr($$input,$pos), $lineno[0] ]);
    };

    #Skip blanks
            $lexlevel == 0
        ?   $$input=~m{\G((?:
                                [\t\ ]+    # Any white space char but \n
                            |   \#[^\n]*  # Perl like comments
                            |   /\*.*?\*/ # C like comments
                            )+)}xsgc
        :   $$input=~m{\G((?:
                                \s+       # any white space char
                            |   \#[^\n]*  # Perl like comments
                            |   /\*.*?\*/ # C like comments
                            )+)}xsgc
    and do {
        my($blanks)=$1;

        #Maybe At EOF
            pos($$input) >= length($$input)
        and return('',[ undef, -1 ]);

        $lineno[1]+= $blanks=~tr/\n//;
    };

    $lineno[0]=$lineno[1];

        $$input=~/\G([A-Za-z_][A-Za-z0-9_]*)/gc
    and return('IDENT',[ $1, $lineno[0] ]);

        $$input=~/\G('(?:[^'\\]|\\\\|\\'|\\)+?')/gc
    and do {
            $1 eq "'error'"
        and do {
            _SyntaxError(0,"Literal 'error' ".
                           "will be treated as error token",$lineno[0]);
            return('IDENT',[ 'error', $lineno[0] ]);
        };
        return('LITERAL',[ $1, $lineno[0] ]);
    };

        $$input=~/\G(%%)/gc
    and do {
        ++$lexlevel;
        return($1, [ $1, $lineno[0] ]);
    };

        $$input=~/\G{/gc
    and do {
        my($level,$from,$code);

        $from=pos($$input);

        $level=1;
        while($$input=~/([{}])/gc) {
                substr($$input,pos($$input)-1,1) eq '\\' #Quoted
            and next;
                $level += ($1 eq '{' ? 1 : -1)
            or last;
        }
            $level
        and  _SyntaxError(2,"Unmatched { opened line $lineno[0]",-1);
        $code = substr($$input,$from,pos($$input)-$from-1);
        $lineno[1]+= $code=~tr/\n//;
        return('CODE',[ $code, $lineno[0] ]);
    };

    if($lexlevel == 0) {# In head section
            $$input=~/\G%(left|right|nonassoc)/gc
        and return('ASSOC',[ uc($1), $lineno[0] ]);
            $$input=~/\G%(start)/gc
        and return('START',[ undef, $lineno[0] ]);
            $$input=~/\G%(expect)/gc
        and return('EXPECT',[ undef, $lineno[0] ]);
            $$input=~/\G%{/gc
        and do {
            my($code);

                $$input=~/\G(.*?)%}/sgc
            or  _SyntaxError(2,"Unmatched %{ opened line $lineno[0]",-1);

            $code=$1;
            $lineno[1]+= $code=~tr/\n//;
            return('HEADCODE',[ $code, $lineno[0] ]);
        };
            $$input=~/\G%(token)/gc
        and return('TOKEN',[ undef, $lineno[0] ]);
            $$input=~/\G%(type)/gc
        and return('TYPE',[ undef, $lineno[0] ]);
            $$input=~/\G%(union)/gc
        and return('UNION',[ undef, $lineno[0] ]);
            $$input=~/\G([0-9]+)/gc
        and return('NUMBER',[ $1, $lineno[0] ]);

    }
    else {# In rule section
            $$input=~/\G%(prec)/gc
        and return('PREC',[ undef, $lineno[0] ]);
    }

    #Always return something
        $$input=~/\G(.)/sg
    or  die "Parse::Yapp::Grammar::Parse: Match (.) failed: report as a BUG";

        $1 eq "\n"
    and ++$lineno[1];

    ( $1 ,[ $1, $lineno[0] ]);

}

sub _SyntaxError {
    my($level,$message,$lineno)=@_;

    $message= "*".
              [ 'Warning', 'Error', 'Fatal' ]->[$level].
              "* $message, at ".
              ($lineno < 0 ? "eof" : "line $lineno").
              ".\n";

        $level > 1
    and die $message;

    warn $message;

        $level > 0
    and ++$nberr;

        $nberr == 20 
    and die "*Fatal* Too many errors detected.\n"
}

sub _AddRules {
    my($lhs,$lineno)=@{$_[0]};
    my($rhss)=$_[1];

        ref($$nterm{$lhs})
    and do {
        _SyntaxError(1,"Non-terminal $lhs redefined: ".
                       "Previously declared line $$syms{$lhs}",$lineno);
        return;
    };

        ref($$term{$lhs})
    and do {
        my($where) = exists($$token{$lhs}) ? $$token{$lhs} : $$syms{$lhs};
        _SyntaxError(1,"Non-terminal $lhs previously ".
                       "declared as token line $where",$lineno);
        return;
    };

        ref($$nterm{$lhs})      #declared through %type
    or  do {
            $$syms{$lhs}=$lineno;   #Say it's declared here
            delete($$term{$lhs});   #No more a terminal
    };
    $$nterm{$lhs}=[];       #It's a non-terminal now

    my($epsrules)=0;        #To issue a warning if more than one epsilon rule

    for my $rhs (@$rhss) {
        my($tmprule)=[ $lhs, [ ], splice(@$rhs,-2) ]; #Init rule

            @$rhs
        or  do {
            ++$$nullable{$lhs};
            ++$epsrules;
        };

        for (0..$#$rhs) {
            my($what,$value)=@{$$rhs[$_]};

                $what eq 'CODE'
            and do {
                my($name)='@'.++$labelno."-$_";
                push(@$rules,[ $name, [], undef, $value ]);
                push(@{$$tmprule[1]},$name);
                next;
            };
            push(@{$$tmprule[1]},$$value[0]);
        }
        push(@$rules,$tmprule);
        push(@{$$nterm{$lhs}},$#$rules);
    }

        $epsrules > 1
    and _SyntaxError(0,"More than one empty rule for symbol $lhs",$lineno);
}

sub Parse {
    my($self)=shift;

        @_ > 0
    or  croak("No input grammar\n");

    my($parsed)={};

    $input=\$_[0];

    $lexlevel=0;
    @lineno=(1,1);
    $nberr=0;
    $prec=0;
    $labelno=0;

    $head=();
    $tail="";

    $syms={};
    $token={};
    $term={};
    $nterm={};
    $rules=[ undef ];   #reserve slot 0 for start rule
    $precterm={};

    $start="";
    $nullable={};
    $expect=0;

    pos($$input)=0;


    $self->YYParse(yylex => \&_Lexer, yyerror => \&_Error);

        $nberr
    and _SyntaxError(2,"Errors detected: No output",-1);

    @$parsed{ 'HEAD', 'TAIL', 'RULES', 'NTERM', 'TERM',
              'NULL', 'PREC', 'SYMS',  'START', 'EXPECT' }
    =       (  $head,  $tail,  $rules,  $nterm,  $term,
               $nullable, $precterm, $syms, $start, $expect);

    undef($input);
    undef($lexlevel);
    undef(@lineno);
    undef($nberr);
    undef($prec);
    undef($labelno);

    undef($head);
    undef($tail);

    undef($syms);
    undef($token);
    undef($term);
    undef($nterm);
    undef($rules);
    undef($precterm);

    undef($start);
    undef($nullable);
    undef($expect);

    $parsed
}


1;

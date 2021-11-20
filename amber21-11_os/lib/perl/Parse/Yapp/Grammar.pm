#
# Module Parse::Yapp::Grammar
#
# (c) Copyright 1998-2001 Francois Desarmenien, all rights reserved.
# (see the pod text in Parse::Yapp module for use and distribution rights)
#
package Parse::Yapp::Grammar;
@ISA=qw( Parse::Yapp::Options );

require 5.004;

use Carp;
use strict;
use Parse::Yapp::Options;
use Parse::Yapp::Parse;

###############
# Constructor #
###############
sub new {
    my($class)=shift;
    my($values);

    my($self)=$class->SUPER::new(@_);

    my($parser)=new Parse::Yapp::Parse;

        defined($self->Option('input'))
    or  croak "No input grammar";

    $values = $parser->Parse($self->Option('input'));

    undef($parser);

    $$self{GRAMMAR}=_ReduceGrammar($values);

        ref($class)
    and $class=ref($class);

    bless($self, $class);
}

###########
# Methods #
###########
##########################
# Method To View Grammar #
##########################
sub ShowRules {
    my($self)=shift;
    my($rules)=$$self{GRAMMAR}{RULES};
    my($ruleno)=-1;
    my($text);

    for (@$rules) {
        my($lhs,$rhs)=@$_;

        $text.=++$ruleno.":\t".$lhs." -> ";
        if(@$rhs) {
            $text.=join(' ',map { $_ eq chr(0) ? '$end' : $_ } @$rhs);
        }
        else {
            $text.="/* empty */";
        }
        $text.="\n";
    }
    $text;
}

###########################
# Method To View Warnings #
###########################
sub Warnings {
    my($self)=shift;
    my($text);
    my($grammar)=$$self{GRAMMAR};

        exists($$grammar{UUTERM})
    and    do {
            $text="Unused terminals:\n\n";
            for (@{$$grammar{UUTERM}}) {
                $text.="\t$$_[0], declared line $$_[1]\n";    
            }
        $text.="\n";
        };
        exists($$grammar{UUNTERM})
    and    do {
            $text.="Useless non-terminals:\n\n";
            for (@{$$grammar{UUNTERM}}) {
                $text.="\t$$_[0], declared line $$_[1]\n";    
            }
        $text.="\n";
        };
        exists($$grammar{UURULES})
    and    do {
            $text.="Useless rules:\n\n";
            for (@{$$grammar{UURULES}}) {
                $text.="\t$$_[0] -> ".join(' ',@{$$_[1]})."\n";    
            }
        $text.="\n";
        };
    $text;
}

######################################
# Method to get summary about parser #
######################################
sub Summary {
    my($self)=shift;
    my($text);

    $text ="Number of rules         : ".
            scalar(@{$$self{GRAMMAR}{RULES}})."\n";
    $text.="Number of terminals     : ".
            scalar(keys(%{$$self{GRAMMAR}{TERM}}))."\n";
    $text.="Number of non-terminals : ".
            scalar(keys(%{$$self{GRAMMAR}{NTERM}}))."\n";
    $text;
}

###############################
# Method to Ouput rules table #
###############################
sub RulesTable {
    my($self)=shift;
    my($inputfile)=$self->Option('inputfile');
    my($linenums)=$self->Option('linenumbers');
    my($rules)=$$self{GRAMMAR}{RULES};
    my($ruleno);
    my($text);

        defined($inputfile)
    or  $inputfile = 'unkown';

    $text="[\n\t";

    $text.=join(",\n\t",
                map {
                    my($lhs,$rhs,$code)=@$_[0,1,3];
                    my($len)=scalar(@$rhs);
                    my($text);

                    $text.="[#Rule ".$ruleno++."\n\t\t '$lhs', $len,";
                    if($code) {
                        $text.= "\nsub".
                                (  $linenums
                                 ? qq(\n#line $$code[1] "$inputfile"\n)
                                 : " ").
                                "{$$code[0]}";
                    }
                    else {
                        $text.=' undef';
                    }
                    $text.="\n\t]";

                    $text;
                } @$rules);

    $text.="\n]";

    $text;
}

################################
# Methods to get HEAD and TAIL #
################################
sub Head {
    my($self)=shift;
    my($inputfile)=$self->Option('inputfile');
    my($linenums)=$self->Option('linenumbers');
    my($text);

        $$self{GRAMMAR}{HEAD}[0]
    or  return '';

        defined($inputfile)
    or  $inputfile = 'unkown';

    for (@{$$self{GRAMMAR}{HEAD}}) {
            $linenums
        and $text.=qq(#line $$_[1] "$inputfile"\n);
        $text.=$$_[0];
    }
    $text
}

sub Tail {
    my($self)=shift;
    my($inputfile)=$self->Option('inputfile');
    my($linenums)=$self->Option('linenumbers');
    my($text);

        $$self{GRAMMAR}{TAIL}[0]
    or  return '';

        defined($inputfile)
    or  $inputfile = 'unkown';

        $linenums
    and $text=qq(#line $$self{GRAMMAR}{TAIL}[1] "$inputfile"\n);
    $text.=$$self{GRAMMAR}{TAIL}[0];

    $text
}


#################
# Private Stuff #
#################

sub _UsefulRules {
    my($rules,$nterm) = @_;
    my($ufrules,$ufnterm);
    my($done);

    $ufrules=pack('b'.@$rules);
    $ufnterm={};

    vec($ufrules,0,1)=1;    #start rules IS always useful

    RULE:
    for (1..$#$rules) { # Ignore start rule
        for my $sym (@{$$rules[$_][1]}) {
                exists($$nterm{$sym})
            and next RULE;
        }
        vec($ufrules,$_,1)=1;
        ++$$ufnterm{$$rules[$_][0]};
    }

    do {
        $done=1;

        RULE:
        for (grep { vec($ufrules,$_,1) == 0 } 1..$#$rules) {
            for my $sym (@{$$rules[$_][1]}) {
                    exists($$nterm{$sym})
                and not exists($$ufnterm{$sym})
                and next RULE;
            }
            vec($ufrules,$_,1)=1;
                exists($$ufnterm{$$rules[$_][0]})
            or  do {
                $done=0;
                ++$$ufnterm{$$rules[$_][0]};
            };
        }

    }until($done);

    ($ufrules,$ufnterm)

}#_UsefulRules

sub _Reachable {
    my($rules,$nterm,$term,$ufrules,$ufnterm)=@_;
    my($reachable);
    my(@fifo)=( 0 );

    $reachable={ '$start' => 1 }; #$start is always reachable

    while(@fifo) {
        my($ruleno)=shift(@fifo);

        for my $sym (@{$$rules[$ruleno][1]}) {

                exists($$term{$sym})
            and do {
                ++$$reachable{$sym};
                next;
            };

                (   not exists($$ufnterm{$sym})
                 or exists($$reachable{$sym}) )
            and next;

            ++$$reachable{$sym};
            push(@fifo, grep { vec($ufrules,$_,1) } @{$$nterm{$sym}});
        }
    }

    $reachable

}#_Reachable

sub _SetNullable {
    my($rules,$term,$nullable) = @_;
    my(@nrules);
    my($done);

    RULE:
    for (@$rules) {
        my($lhs,$rhs)=@$_;

            exists($$nullable{$lhs})
        and next;

        for (@$rhs) {
                exists($$term{$_})
            and next RULE;
        }
        push(@nrules,[$lhs,$rhs]);
    }

    do {
        $done=1;

        RULE:
        for (@nrules) {
            my($lhs,$rhs)=@$_;

                    exists($$nullable{$lhs})
                and next;

                for (@$rhs) {
                        exists($$nullable{$_})
                    or  next RULE;
                }
            $done=0;
            ++$$nullable{$lhs};
        }

    }until($done);
}

sub _ReduceGrammar {
    my($values)=@_;
    my($ufrules,$ufnterm,$reachable);
    my($grammar)={ HEAD => $values->{HEAD},
                   TAIL => $values->{TAIL},
                   EXPECT => $values->{EXPECT} };
    my($rules,$nterm,$term) =  @$values {'RULES', 'NTERM', 'TERM'};

    ($ufrules,$ufnterm) = _UsefulRules($rules,$nterm);

        exists($$ufnterm{$values->{START}})
    or  die "*Fatal* Start symbol $values->{START} derives nothing, at eof\n";

    $reachable = _Reachable($rules,$nterm,$term,$ufrules,$ufnterm);

    $$grammar{TERM}{chr(0)}=undef;
    for my $sym (keys %$term) {
            (   exists($$reachable{$sym})
             or exists($values->{PREC}{$sym}) )
        and do {
            $$grammar{TERM}{$sym}
                = defined($$term{$sym}[0]) ? $$term{$sym} : undef;
            next;
        };
        push(@{$$grammar{UUTERM}},[ $sym, $values->{SYMS}{$sym} ]);
    }

    $$grammar{NTERM}{'$start'}=[];
    for my $sym (keys %$nterm) {
            exists($$reachable{$sym})
        and do {
                exists($values->{NULL}{$sym})
            and ++$$grammar{NULLABLE}{$sym};
            $$grammar{NTERM}{$sym}=[];
            next;
        };
        push(@{$$grammar{UUNTERM}},[ $sym, $values->{SYMS}{$sym} ]);
    }

    for my $ruleno (0..$#$rules) {
            vec($ufrules,$ruleno,1)
        and exists($$grammar{NTERM}{$$rules[$ruleno][0]})
        and do {
            push(@{$$grammar{RULES}},$$rules[$ruleno]);
            push(@{$$grammar{NTERM}{$$rules[$ruleno][0]}},$#{$$grammar{RULES}});
            next;
        };
        push(@{$$grammar{UURULES}},[ @{$$rules[$ruleno]}[0,1] ]);
    }

    _SetNullable(@$grammar{'RULES', 'TERM', 'NULLABLE'});

    $grammar;
}#_ReduceGrammar

1;

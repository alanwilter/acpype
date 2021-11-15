#
# Module Parse::Yapp::Lalr
#
# (c) Copyright 1998-2001 Francois Desarmenien, all rights reserved.
# (see the pod text in Parse::Yapp module for use and distribution rights)
#
package Parse::Yapp::Lalr;
@ISA=qw( Parse::Yapp::Grammar );

require 5.004;

use Parse::Yapp::Grammar;

=for nobody

Parse::Yapp::Compile Object Structure:
--------------------------------------
{
   GRAMMAR =>    Parse::Yapp::Grammar,
   STATES  =>    [ { CORE    => [ items... ],
                     ACTIONS  => { term => action }
                     GOTOS   => { nterm => stateno }
                   }... ]
   CONFLICTS=>{ SOLVED => { stateno  => [ ruleno, token, solved ] }, 
                FORCED => { TOTAL => [ nbsr, nbrr ],
                            DETAIL => { stateno => { TOTAL => [ nbsr, nbrr ] }
                                                     LIST => [ ruleno, token ]
                                                   }
                                      }
                          } 
}

'items' are of form: [ ruleno, dotpos ]
'term' in ACTIONS is '' means default action
'action' may be:
    undef:  explicit error (nonassociativity)
    0    :  accept
    >0   :  shift and go to state 'action'
    <0   :  reduce using rule -'action'
'solved' may have values of:
         'shift'  if solved as Shift
         'reduce' if solved as Reduce
         'error'  if solved by discarding both Shift and Reduce (nonassoc)

SOLVED is a set of states containing Solved conflicts
FORCED are forced conflict resolutions

nbsr and nbrr are number of shift/reduce and reduce/reduce conflicts

TOTAL is the total number of SR/RR conflicts for the parser

DETAIL is the detail of conflicts for each state
TOTAL is the total number of SR/RR conflicts for a state
LIST is the list of discarded reductions (for display purpose only)


=cut

use strict;

use Carp;

###############
# Constructor #
###############
sub new {
    my($class)=shift;

		ref($class)
	and	$class=ref($class);

	my($self)=$class->SUPER::new(@_);
    $self->_Compile();
    bless($self,$class);
}
###########
# Methods #
###########

###########################
# Method To View Warnings #
###########################
sub Warnings {
    my($self)=shift;
    my($text);
    my($nbsr,$nbrr)=@{$$self{CONFLICTS}{FORCED}{TOTAL}};

	$text=$self->SUPER::Warnings();

        $nbsr != $$self{GRAMMAR}{EXPECT}
    and $text.="$nbsr shift/reduce conflict".($nbsr > 1 ? "s" : "");

        $nbrr
    and do {
            $nbsr
        and $text.=" and ";
        $text.="$nbrr reduce/reduce conflict".($nbrr > 1 ? "s" : "");
    };

       (    $nbsr != $$self{GRAMMAR}{EXPECT}
        or  $nbrr)
    and $text.="\n";

    $text;
}
#############################
# Method To View DFA States #
#############################
sub ShowDfa {
    my($self)=shift;
    my($text);
    my($grammar,$states)=($$self{GRAMMAR}, $$self{STATES});

    for my $stateno (0..$#$states) {
        my(@shifts,@reduces,@errors,$default);

        $text.="State $stateno:\n\n";

        #Dump Kernel Items
        for (sort {     $$a[0] <=> $$b[0]
                    or  $$a[1] <=> $$b[1] } @{$$states[$stateno]{'CORE'}}) {
            my($ruleno,$pos)=@$_;
            my($lhs,$rhs)=@{$$grammar{RULES}[$ruleno]}[0,1];
            my(@rhscopy)=@$rhs;
        
                $ruleno
            or  $rhscopy[-1] = '$end';

            splice(@rhscopy,$pos,0,'.');
            $text.= "\t$lhs -> ".join(' ',@rhscopy)."\t(Rule $ruleno)\n";
        }

        #Prepare Actions
        for (keys(%{$$states[$stateno]{ACTIONS}})) {
            my($term,$action)=($_,$$states[$stateno]{ACTIONS}{$_});

                $term eq chr(0)
            and $term = '$end';

                not defined($action)
            and do {
                push(@errors,$term);
                next;
            };

                $action > 0
            and do {
                push(@shifts,[ $term, $action ]);
                next;
            };

            $action = -$action;

                $term
            or  do {
                $default= [ '$default', $action ];
                next;
            };

            push(@reduces,[ $term, $action ]);
        }

            #Dump shifts
            @shifts
        and do {
            $text.="\n";
            for (sort { $$a[0] cmp $$b[0] } @shifts) {
                my($term,$shift)=@$_;

                $text.="\t$term\tshift, and go to state $shift\n";
            }
        };

            #Dump errors
            @errors
        and do {
            $text.="\n";
            for my $term (sort { $a cmp $b } @errors) {
                $text.="\t$term\terror (nonassociative)\n";
            }
        };

        #Prepare reduces
            exists($$self{CONFLICTS}{FORCED}{DETAIL}{$stateno})
        and push(@reduces,@{$$self{CONFLICTS}{FORCED}{DETAIL}{$stateno}{LIST}});

        @reduces=sort { $$a[0] cmp $$b[0] or $$a[1] <=> $$b[1] } @reduces;

            defined($default)
        and push(@reduces,$default);

        #Dump reduces
            @reduces
        and do {
            $text.="\n";
            for (@reduces) {
                my($term,$ruleno)=@$_;
                my($discard);

                    $ruleno < 0
                and do {
                    ++$discard;
                    $ruleno = -$ruleno;
                };

                $text.= "\t$term\t".($discard  ? "[" : "");
                if($ruleno) {
                    $text.= "reduce using rule $ruleno ".
                            "($$grammar{RULES}[$ruleno][0])";
                }
                else {
                    $text.='accept';
                }
                $text.=($discard  ? "]" : "")."\n";
            }
        };

            #Dump gotos
            exists($$states[$stateno]{GOTOS})
        and    do {
                $text.= "\n";
                for (keys(%{$$states[$stateno]{GOTOS}})) {
                    $text.= "\t$_\tgo to state $$states[$stateno]{GOTOS}{$_}\n";
                }
            };

        $text.="\n";
    }
    $text;
}

######################################
# Method to get summary about parser #
######################################
sub Summary {
    my($self)=shift;
    my($text);

	$text=$self->SUPER::Summary();
    $text.="Number of states        : ".
            scalar(@{$$self{STATES}})."\n";
    $text;
}

#######################################
# Method To Get Infos about conflicts #
#######################################
sub Conflicts {
    my($self)=shift;
    my($states)=$$self{STATES};
    my($conflicts)=$$self{CONFLICTS};
    my($text);

    for my $stateno ( sort { $a <=> $b } keys(%{$$conflicts{SOLVED}})) {

        for (@{$$conflicts{SOLVED}{$stateno}}) {
            my($ruleno,$token,$how)=@$_;

                $token eq chr(0)
            and $token = '$end';

            $text.="Conflict in state $stateno between rule ".
                   "$ruleno and token $token resolved as $how.\n"; 
        }
    };

    for my $stateno ( sort { $a <=> $b } keys(%{$$conflicts{FORCED}{DETAIL}})) {
        my($nbsr,$nbrr)=@{$$conflicts{FORCED}{DETAIL}{$stateno}{TOTAL}};

        $text.="State $stateno contains ";

            $nbsr
        and $text.="$nbsr shift/reduce conflict".
                   ($nbsr > 1 ? "s" : "");

            $nbrr
        and do {
                $nbsr
            and $text.=" and ";

            $text.="$nbrr reduce/reduce conflict".
                   ($nbrr > 1 ? "s" : "");
        };
        $text.="\n";
    };

    $text;
}

#################################
# Method to dump parsing tables #
#################################
sub DfaTable {
    my($self)=shift;
    my($states)=$$self{STATES};
    my($stateno);
    my($text);

    $text="[\n\t{";

    $text.=join("\n\t},\n\t{",
                map {
                    my($state)=$_;
                    my($text);

                    $text="#State ".$stateno++."\n\t\t";

                       (    not exists($$state{ACTIONS}{''})
                        or  keys(%{$$state{ACTIONS}}) > 1)
                    and do {

                        $text.="ACTIONS => {\n\t\t\t";

                        $text.=join(",\n\t\t\t",
                                map {
                                    my($term,$action)=($_,$$state{ACTIONS}{$_});
                                    my($text);

                                    if(substr($term,0,1) eq "'") {
									    $term=~s/([\@\$\"])/\\$1/g;
                                        $term=~s/^'|'$/"/g;
                                    }
                                    else {
                                        $term=      $term eq chr(0)
                                                ?   "''" 
                                                :   "'$term'";
                                    }

                                    if(defined($action)) {
                                        $action=int($action);
                                    }
                                    else {
                                        $action='undef';
                                    }

                                    "$term => $action";
                                
                                } grep { $_ } keys(%{$$state{ACTIONS}}));

                        $text.="\n\t\t}";
                    };

                        exists($$state{ACTIONS}{''})
                    and do {
                            keys(%{$$state{ACTIONS}}) > 1
                        and $text.=",\n\t\t";

                        $text.="DEFAULT => $$state{ACTIONS}{''}";
                    };

                        exists($$state{GOTOS})
                    and do {
                        $text.=",\n\t\tGOTOS => {\n\t\t\t";
                        $text.=join(",\n\t\t\t",
                                map {
                                    my($nterm,$stateno)=($_,$$state{GOTOS}{$_});
                                    my($text);

                                    "'$nterm' => $stateno";
                                
                                } keys(%{$$state{GOTOS}}));
                        $text.="\n\t\t}";
                    };

                    $text;

                }@$states);

    $text.="\n\t}\n]";

    $text;

}


####################################
# Method to build Dfa from Grammar #
####################################
sub _Compile {
	my($self)=shift;
	my($grammar,$states);

	$grammar=$self->{GRAMMAR};

    $states = _LR0($grammar);

    $self->{CONFLICTS} = _LALR($grammar,$states);

    $self->{STATES}=$states;
}

#########################
# LR0 States Generation #
#########################
#
###########################
# General digraph routine #
###########################
sub _Digraph {
    my($rel,$F)=@_;
    my(%N,@S);
    my($infinity)=(~(1<<31));
    my($Traverse);

    $Traverse = sub {
        my($x,$d)=@_;
        my($y);

        push(@S,$x);
        $N{$x}=$d;

            exists($$rel{$x})
        and do {
            for $y (keys(%{$$rel{$x}})) {
                    exists($N{$y})
                or  &$Traverse($y,$d+1);

                    $N{$y} < $N{$x}
                and $N{$x} = $N{$y};

                $$F{$x}|=$$F{$y};
            }
        };

            $N{$x} == $d
        and do {
            for(;;) {
                $y=pop(@S);
                $N{$y}=$infinity;
                    $y eq $x
                and last;
                $$F{$y}=$$F{$x};
            }
        };
    };

    for (keys(%$rel)) {
            exists($N{$_})
        or  &$Traverse($_,1);
    }
}
#######################
# Generate LR0 states # 
#######################
=for nobody
Formula used for closures:

    CLOSE(A) = DCLOSE(A) u U (CLOSE(B) | A close B)

where:

    DCLOSE(A) = { [ A -> alpha ] in P }

    A close B iff [ A -> B gamma ] in P

=cut
sub _SetClosures {
	my($grammar)=@_;
    my($rel,$closures);

    for my $symbol (keys(%{$$grammar{NTERM}})) {
        $closures->{$symbol}=pack('b'.@{$$grammar{RULES}});

        for my $ruleno (@{$$grammar{NTERM}{$symbol}}) {
            my($rhs)=$$grammar{RULES}[$ruleno][1];

            vec($closures->{$symbol},$ruleno,1)=1;

                @$rhs > 0
            and exists($$grammar{NTERM}{$$rhs[0]})
            and ++$rel->{$symbol}{$$rhs[0]};
        }
    }
    _Digraph($rel,$closures);

	$closures
}

sub _Closures {
    my($grammar,$core,$closures)=@_;
    my($ruleset)=pack('b'.@{$$grammar{RULES}});

    for (@$core) {
        my($ruleno,$pos)=@$_;
        my($rhs)=$$grammar{RULES}[$ruleno][1];

            $pos < @$rhs
        and exists($closures->{$$rhs[$pos]})
        and $ruleset|=$closures->{$$rhs[$pos]};
    }
    [ @$core, map  { [ $_, 0 ] }
              grep { vec($ruleset,$_,1) }
              0..$#{$$grammar{RULES}} ];
}

sub _Transitions {
    my($grammar,$cores,$closures,$states,$stateno)=@_;
    my($core)=$$states[$stateno]{'CORE'};
    my(%transitions);

    for (@{_Closures($grammar,$core,$closures)}) {
        my($ruleno,$pos)=@$_;
        my($rhs)=$$grammar{RULES}[$ruleno][1];

            $pos == @$rhs
        and do {
            push(@{$$states[$stateno]{ACTIONS}{''}},$ruleno);
            next;
        };
        push(@{$transitions{$$rhs[$pos]}},[ $ruleno, $pos+1 ]);
    }

    for (keys(%transitions)) {
        my($symbol,$core)=($_,$transitions{$_});
        my($corekey)=join(',',map  { join('.',@$_) }
                              sort {    $$a[0] <=> $$b[0]
                                    or  $$a[1] <=> $$b[1] }
                              @$core);
        my($tostateno);

            exists($cores->{$corekey})
        or  do {
            push(@$states,{ 'CORE' => $core });
            $cores->{$corekey}=$#$states;
        };

        $tostateno=$cores->{$corekey};
        push(@{$$states[$tostateno]{FROM}},$stateno);

			exists($$grammar{TERM}{$_})
		and	do {
            $$states[$stateno]{ACTIONS}{$_} = [ $tostateno ];
			next;
		};
        $$states[$stateno]{GOTOS}{$_} = $tostateno;
    }
}

sub _LR0 {
	my($grammar)=@_;
	my($states) = [];
    my($stateno);
    my($closures);  #$closures={ nterm => ruleset,... }
	my($cores)={};  # { "itemlist" => stateno, ... }
					# where "itemlist" has the form:
					# "ruleno.pos,ruleno.pos" ordered by ruleno,pos

    $closures = _SetClosures($grammar);
    push(@$states,{ 'CORE' => [ [ 0, 0 ] ] });
    for($stateno=0;$stateno<@$states;++$stateno) {
        _Transitions($grammar,$cores,$closures,$states,$stateno);
    }

	$states
}

#########################################################
# Add Lookahead tokens where needed to make LALR states #
#########################################################
=for nobody
    Compute First sets for non-terminal using the following formula:

    FIRST(A) =      { a in T u { epsilon } | A l a }
                u
                U   { FIRST(B) | B in V and A l B }

    where:

    A l x iff [ A -> X1 X2 .. Xn x alpha ] in P and Xi =>* epsilon, 1 <= i <= n
=cut
sub _SetFirst {
	my($grammar,$termlst,$terminx)=@_;
    my($rel,$first)=( {}, {} );

    for my $symbol (keys(%{$$grammar{NTERM}})) {
        $first->{$symbol}=pack('b'.@$termlst);

        RULE:
        for my $ruleno (@{$$grammar{NTERM}{$symbol}}) {
            my($rhs)=$$grammar{RULES}[$ruleno][1];

            for (@$rhs) {
                    exists($terminx->{$_})
                and do {
                    vec($first->{$symbol},$terminx->{$_},1)=1;
                    next RULE;
                };
                ++$rel->{$symbol}{$_};
                    exists($$grammar{NULLABLE}{$_})
                or  next RULE;
            }
            vec($first->{$symbol},0,1)=1;
        }
    }
    _Digraph($rel,$first);

	$first
}

sub _Preds {
    my($states,$stateno,$len)=@_;
    my($queue, $preds);

        $len
    or  return [ $stateno ];

    $queue=[ [ $stateno, $len ] ];
    while(@$queue) {
        my($pred) = shift(@$queue);
        my($stateno, $len) = @$pred;

            $len == 1
        and do {
			push(@$preds,@{$states->[$stateno]{FROM}});
            next;
        };

        push(@$queue, map { [ $_, $len - 1 ] }
					  @{$states->[$stateno]{FROM}});
    }

    # Pass @$preds through a hash to ensure unicity
    [ keys( %{ +{ map { ($_,1) } @$preds } } ) ];
}

sub _FirstSfx {
    my($grammar,$firstset,$termlst,$terminx,$ruleno,$pos,$key)=@_;
    my($first)=pack('b'.@$termlst);
    my($rhs)=$$grammar{RULES}[$ruleno][1];

    for (;$pos < @$rhs;++$pos) {
            exists($terminx->{$$rhs[$pos]})
        and do {
            vec($first,$terminx->{$$rhs[$pos]},1)=1;
            return($first);
        };
        $first|=$firstset->{$$rhs[$pos]};

            vec($first,0,1)
        and vec($first,0,1)=0;

            exists($$grammar{NULLABLE}{$$rhs[$pos]})
        or  return($first);

    }
    vec($first,0,1)=1;
    $first;
}

=for noboby
    Compute Follow sets using following formula:

    FOLLOW(p,A) =       READ(p,A)
                    u
                    U   { FOLLOW(q,B) | (p,A) include (q,B)

    where:
 
    READ(p,A) = U { FIRST(beta) | [ A -> alpha A . beta ] in KERNEL(GOTO(p,A))
                  } - { epsilon }

    (p,a) include (q,B) iff [ B -> alpha A . beta ] in KERNEL(GOTO(p,A),
                            epsilon in FIRST(beta) and
                            q in PRED(p,alpha)
=cut
sub _ComputeFollows {
	my($grammar,$states,$termlst)=@_;
	my($firstset,$terminx);
	my($inconsistent, $rel, $follows, $sfx)= ( {}, {}, {}, {} );

    %$terminx= map { ($termlst->[$_],$_) } 0..$#$termlst;

    $firstset=_SetFirst($grammar,$termlst,$terminx);

    for my $stateno (0..$#$states) {
		my($state)=$$states[$stateno];

           	exists($$state{ACTIONS}{''})
        and (   @{$$state{ACTIONS}{''}} > 1
             or keys(%{$$state{ACTIONS}}) > 1 )
		and do {
			++$inconsistent->{$stateno};

			for my $ruleno (@{$$state{ACTIONS}{''}}) {
				my($lhs,$rhs)=@{$$grammar{RULES}[$ruleno]}[0,1];

                for my $predno (@{_Preds($states,$stateno,scalar(@$rhs))}) {
                    ++$rel->{"$stateno.$ruleno"}{"$predno.$lhs"};
                }
			}
		};

    		exists($$state{GOTOS})
		or	next;

        for my $symbol (keys(%{$$state{GOTOS}})) {
            my($tostate)=$$states[$$state{GOTOS}{$symbol}];
            my($goto)="$stateno.$symbol";

            $follows->{$goto}=pack('b'.@$termlst);

            for my $item (@{$$tostate{'CORE'}}) {
                my($ruleno,$pos)=@$item;
				my($key)="$ruleno.$pos";

					exists($sfx->{$key})
				or	$sfx->{$key} = _FirstSfx($grammar,$firstset,
											 $termlst,$terminx,
											 $ruleno,$pos,$key);

                $follows->{$goto}|=$sfx->{$key};

                    vec($follows->{$goto},0,1)
                and do {
                    my($lhs)=$$grammar{RULES}[$ruleno][0];

                    vec($follows->{$goto},0,1)=0;

                    for my $predno (@{_Preds($states,$stateno,$pos-1)}) {
                        ++$rel->{$goto}{"$predno.$lhs"};
                    }
                };
            }
        }
    }
    _Digraph($rel,$follows);

	($follows,$inconsistent)
}

sub _ComputeLA {
	my($grammar,$states)=@_;
	my($termlst)= [ '',keys(%{$$grammar{TERM}}) ];

    my($follows,$inconsistent) = _ComputeFollows($grammar,$states,$termlst);

    for my $stateno ( keys(%$inconsistent ) ) {
        my($state)=$$states[$stateno];
        my($conflict);

        #NB the sort is VERY important for conflicts resolution order
        for my $ruleno (sort { $a <=> $b }
                        @{$$state{ACTIONS}{''}}) {
            for my $term ( map { $termlst->[$_] } grep {
                           vec($follows->{"$stateno.$ruleno"},$_,1) }
                           0..$#$termlst) {
                    exists($$state{ACTIONS}{$term})
                and ++$conflict;
                push(@{$$state{ACTIONS}{$term}},-$ruleno);
            }
        }
        delete($$state{ACTIONS}{''});
            $conflict
        or  delete($inconsistent->{$stateno});
    }

	$inconsistent
}

#############################
# Solve remaining conflicts #
#############################

sub _SolveConflicts {
	my($grammar,$states,$inconsistent)=@_;
    my(%rulesprec,$RulePrec);
    my($conflicts)={    SOLVED  =>  {},
                    	FORCED  =>  {   TOTAL   =>  [ 0, 0 ],
                                    	DETAIL  =>  {}
                                 	}
                };

    $RulePrec = sub {
        my($ruleno)=@_;
        my($rhs,$rprec)=@{$$grammar{RULES}[$ruleno]}[1,2];
        my($lastterm);

            defined($rprec)
        and return($rprec);

            exists($rulesprec{$ruleno})
        and return($rulesprec{$ruleno});

        $lastterm=(grep { exists($$grammar{TERM}{$_}) } @$rhs)[-1];

            defined($lastterm)
        and ref($$grammar{TERM}{$lastterm})
        and do {
            $rulesprec{$ruleno}=$$grammar{TERM}{$lastterm}[1];
            return($rulesprec{$ruleno});
        };

        undef;
    };

    for my $stateno (keys(%$inconsistent)) {
        my($state)=$$states[$stateno];
        my($actions)=$$state{ACTIONS};
        my($nbsr,$nbrr);

        for my $term ( keys(%$actions) ) {
            my($act)=$$actions{$term};

                @$act > 1
            or  next;

                $$act[0] > 0
            and ref($$grammar{TERM}{$term})
            and do {
                my($assoc,$tprec)=@{$$grammar{TERM}{$term}};
                my($k,$error);

                for ($k=1;$k<@$act;++$k) {
                    my($ruleno)=-$$act[$k];
                    my($rprec)=&$RulePrec($ruleno);

                        defined($rprec)
                    or  next;

                        (     $tprec > $rprec
                         or ( $tprec == $rprec and $assoc eq 'RIGHT'))
                    and do {
                        push(@{$$conflicts{SOLVED}{$stateno}},
                             [ $ruleno, $term, 'shift' ]);
                        splice(@$act,$k--,1);
                        next;
                    };
                        (   $tprec < $rprec
                         or $assoc eq 'LEFT')
                    and do {
                        push(@{$$conflicts{SOLVED}{$stateno}},
                             [ $ruleno, $term, 'reduce' ]);
                            $$act[0] > 0
                        and do {
                            splice(@$act,0,1);
                            --$k;
                        };
                        next;
                    };
                    push(@{$$conflicts{SOLVED}{$stateno}},
                         [ $ruleno, $term, 'error' ]);
                    splice(@$act,$k--,1);
                        $$act[0] > 0
                    and do {
                        splice(@$act,0,1);
                        ++$error;
                        --$k;
                    };
                }
                    $error
                and unshift(@$act,undef);
            };

                @$act > 1
            and do {
                $nbrr += @$act - 2;
                ($$act[0] > 0 ? $nbsr : $nbrr) += 1;
                push(@{$$conflicts{FORCED}{DETAIL}{$stateno}{LIST}},
                    map { [ $term, $_ ] } splice(@$act,1));
            };
        }

            $nbsr
        and do {
            $$conflicts{FORCED}{TOTAL}[0]+=$nbsr;
            $$conflicts{FORCED}{DETAIL}{$stateno}{TOTAL}[0]+=$nbsr;
        };

            $nbrr
        and do {
            $$conflicts{FORCED}{TOTAL}[1]+=$nbrr;
            $$conflicts{FORCED}{DETAIL}{$stateno}{TOTAL}[1]+=$nbrr;
        };

    }

	$conflicts
}

###############################
# Make default reduce actions #
###############################
sub _SetDefaults {
	my($states)=@_;

    for my $state (@$states) {
        my($actions)=$$state{ACTIONS};
        my(%reduces,$default,$nodefault);

            exists($$actions{''})
        and do {
            $$actions{''}[0] = -$$actions{''}[0];
			++$nodefault;
        };

		#shift error token => no default
            exists($$actions{error})
        and $$actions{error}[0] > 0
        and ++$nodefault;

        for my $term (keys(%$actions)) {

			$$actions{$term}=$$actions{$term}[0];

                (   not defined($$actions{$term})
                 or $$actions{$term} > 0
                 or $nodefault)
            and next;

            push(@{$reduces{$$actions{$term}}},$term);
        }

			keys(%reduces) > 0
		or	next;

        $default=( map { $$_[0] }
                   sort { $$b[1] <=> $$a[1] or $$b[0] <=> $$a[0] }
                   map { [ $_, scalar(@{$reduces{$_}}) ] }
                   keys(%reduces))[0];

        delete(@$actions{ @{$reduces{$default}} });
        $$state{ACTIONS}{''}=$default;
    }
}

sub _LALR {
	my($grammar,$states) = @_;
	my($conflicts,$inconsistent);

    $inconsistent = _ComputeLA($grammar,$states);

    $conflicts = _SolveConflicts($grammar,$states,$inconsistent);
    _SetDefaults($states);

	$conflicts
}


1;

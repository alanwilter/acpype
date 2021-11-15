#
# Module Parse::Yapp::Output
#
# (c) Copyright 1998-2001 Francois Desarmenien, all rights reserved.
# (see the pod text in Parse::Yapp module for use and distribution rights)
#
package Parse::Yapp::Output;
@ISA=qw ( Parse::Yapp::Lalr );

require 5.004;

use Parse::Yapp::Lalr;
use Parse::Yapp::Driver;

use strict;

use Carp;

sub _CopyDriver {
	my($text)='#Included Parse/Yapp/Driver.pm file'.('-' x 40)."\n";
		open(DRV,$Parse::Yapp::Driver::FILENAME)
	or	die "BUG: could not open $Parse::Yapp::Driver::FILENAME";
	$text.="{\n".join('',<DRV>)."}\n";
	close(DRV);
	$text.='#End of include'.('-' x 50)."\n";
}

sub Output {
    my($self)=shift;

    $self->Options(@_);

    my($package)=$self->Option('classname');
    my($head,$states,$rules,$tail,$driver);
    my($version)=$Parse::Yapp::Driver::VERSION;
    my($datapos);
    my($text)=$self->Option('template') ||<<'EOT';
####################################################################
#
#    This file was generated using Parse::Yapp version <<$version>>.
#
#        Don't edit this file, use source file instead.
#
#             ANY CHANGE MADE HERE WILL BE LOST !
#
####################################################################
package <<$package>>;
use vars qw ( @ISA );
use strict;

@ISA= qw ( Parse::Yapp::Driver );
<<$driver>>

<<$head>>

sub new {
        my($class)=shift;
        ref($class)
    and $class=ref($class);

    my($self)=$class->SUPER::new( yyversion => '<<$version>>',
                                  yystates =>
<<$states>>,
                                  yyrules  =>
<<$rules>>,
                                  @_);
    bless($self,$class);
}

<<$tail>>
1;
EOT

	$driver='use Parse::Yapp::Driver;';

        defined($package)
    or $package='Parse::Yapp::Default';

	$head= $self->Head();
	$rules=$self->RulesTable();
	$states=$self->DfaTable();
	$tail= $self->Tail();

		$self->Option('standalone')
	and	$driver=_CopyDriver();

	$text=~s/<<(\$.+)>>/$1/gee;

	$text;
}

1;

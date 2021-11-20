package Chemistry::Obj;
$VERSION = 0.35;
# $Id: Obj.pm,v 1.28 2005/05/20 19:01:05 itubert Exp $
use 5.006;

use strict;
use Carp;

=head1 NAME

Chemistry::Obj - Abstract chemistry object

=head1 SYNOPSIS

    package MyObj;
    use base "Chemistry::Obj";
    Chemistry::Obj::accessor('color', 'flavor');

    package main;
    my $obj = MyObj->new(name => 'bob', color => 'red');
    $obj->attr(size => 42);
    $obj->color('blue');
    my $color = $obj->color;
    my $size = $obj->attr('size');

=head1 DESCRIPTION

This module implements some generic methods that are used by L<Chemistry::Mol>,
L<Chemistry::Atom>, L<Chemistry::Bond>, L<Chemistry::File>, etc.

=head2 Common Attributes

There are some common attributes that may be found in molecules, bonds, and
atoms, such as id, name, and type. They are all accessed through the methods of
the same name. For example, to get the id, call C<< $obj->id >>; to set the id,
call C<< $obj->id('new_id') >>.

=over 4

=item id

Objects should have a unique ID. The user has the responsibility for uniqueness
if he assigns ids; otherwise a unique ID is assigned sequentially.

=item name

An arbitrary name for an object. The name doesn't need to be unique.

=item type

The interpretation of this attribute is not specified here, but it's typically 
used for bond orders and atom types.

=item attr

A space where the user can store any kind of information about the object.  The
accessor method for attr expects the attribute name as the first parameter, and
(optionally) the new value as the second parameter. It can also take a hash or
hashref with several attributes. Examples:

    $color = $obj->attr('color');
    $obj->attr(color => 'red');
    $obj->attr(color => 'red', flavor => 'cherry');
    $obj->attr({color => 'red', flavor => 'cherry'});

=cut

sub attr {
    my $self = shift;
    my ($attr) = @_;
    if (ref $attr eq 'HASH') {
        $self->{attr} = { %$attr };
    } elsif (@_ == 1) {
        return $self->{attr}{$attr};
    } elsif (@_ == 0) {
        return {%{$self->{attr}}};
    } else {
        while (@_ > 1) {
            $attr = shift;
            $self->{attr}{$attr} = shift;
        }
    }
    $self;
}

=back

=head1 OTHER METHODS

=over

=item $obj->del_attr($attr_name)

Delete an attribute.

=cut

sub del_attr {
    my $self = shift;
    my $attr = shift;
    delete $self->{attr}{$attr};
}

# A generic class attribute set/get method generator
sub accessor {
    my $pkg = caller;
    no strict 'refs';
    for my $attribute (@_) {
        *{"${pkg}::$attribute"} =
          sub {
              my $self = shift;
              return $self->{$attribute} unless @_;
              $self->{$attribute} = shift;
              return $self;
          };
    }
}

sub print_attr {
    my $self = shift;
    my ($indent) = @_;
    my $ret = '';
    
    for my $attr (keys %{$self->{attr}}) {
        $ret .= "$attr: ".$self->attr($attr)."\n";
    }
    $ret and $ret =~ s/^/"    "x$indent/gem;
    $ret;
}

my $N = 0; # atom ID counter
sub nextID { "obj".++$N; }
sub reset_id { $N = 0; }

=item $class->new(name => value, name => value...) 

Generic object constructor. It will automatically call each "name" method with
the parameter "value". For example,

    $bob = Chemistry::Obj->new(name => 'bob', attr => {size => 42});

is equivalent to

    $bob = Chemistry::Obj->new;
    $bob->name('bob');
    $bob->attr({size => 42});

=cut

sub new {
    my $class = shift;
    my %args = @_;
    my $self = bless {
        id => $class->nextID,
        #$class->default_args, 
    }, ref $class || $class;
    $self->$_($args{$_}) for (keys %args);
    return $self;
}

#sub default_args { (id => shift->nextID) }

=back

=head1 OPERATOR OVERLOADING

Chemistry::Obj overloads a couple of operators for convenience.

=over

=cut

use overload 
    '""' => "stringify",
    'cmp' => "obj_cmp",
    '0+', => "as_number",
    fallback => 1,
    ;

=item ""

The stringification operator. Stringify an object as its id. For example, If an
object $obj has the id 'a1', print "$obj" will print 'a1' instead of something
like 'Chemistry::Obj=HASH(0x810bbdc)'. If you really want to get the latter,
you can call C<overload::StrVal($obj)>. See L<overload> for details.

=cut

sub stringify {
    my $self = shift;
    $self->id;
}

sub as_number {
    $_[0];
}

=item cmp

Compare objects by ID. This automatically overloads C<eq>, C<ne>, C<lt>, C<le>,
C<gt>, and C<ge> as well. For example, C<$obj1 eq $obj2> returns true if both
objects have the same id, even if they are different objects with different
memory addresses. In contrast, C<$obj1 == $obj2> will return true only if
C<$obj1> and C<$obj2> point to the same object, with the same memory address.

=cut

sub obj_cmp {
    my ($a, $b) = @_;
    no warnings;

    return $a->{id} cmp $b->{id};
}

=back

=cut

accessor(qw(name type));

sub id {
    my $self = shift;
    return $self->{id} unless @_;
    if ($self->{parent}) {
        my $new_id = shift;
        my $old_id = $self->{id};
        $self->{id} = $new_id;
        $self->{parent}->_change_id($old_id, $new_id);
    } else {
        $self->{id} = shift;
    }
}

# this is an experimental method and shouldn't be used!
sub use {
    my ($pack, $module, @args) = @_;
    $pack = ref $pack || $pack;
    my $args = @args ? "(@args)" : '';
    eval "package $pack; use $module $args";
}

1;

=head1 VERSION

0.35

=head1 SEE ALSO

L<Chemistry::Atom>, L<Chemistry::Bond>, L<Chemistry::Mol>

The PerlMol website L<http://www.perlmol.org/>

=head1 AUTHOR

Ivan Tubert-Brohman E<lt>itub@cpan.orgE<gt>

=head1 COPYRIGHT

Copyright (c) 2005 Ivan Tubert-Brohman. All rights reserved. This program is
free software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut


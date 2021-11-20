package Statistics::Descriptive::Smoother;

use strict;
use warnings;

use Carp;

our $VERSION = '3.0604';

sub instantiate {
    my ($class, $args) = @_;

    my $method     = delete $args->{method};
    my $coeff      = delete $args->{coeff} || 0;
    my $ra_samples = delete $args->{samples};
    my $ra_data    = delete $args->{data};
 
    if ($coeff < 0 || $coeff > 1) {
        carp("Invalid smoothing coefficient C $coeff\n");
        return;
    }
    if (@$ra_data < 2) {
        carp("Need at least 2 samples to smooth the data\n");
        return;
    }
    $method = ucfirst(lc($method));
    my $sub_class = __PACKAGE__."::$method";
    eval "require $sub_class";
    die "No such class $sub_class: $@" if $@;

    return $sub_class->_new({
        data       => $ra_data,
        samples    => $ra_samples,
        count      => scalar @$ra_data,
        coeff      => $coeff,
    });
}

sub get_smoothing_coeff { $_[0]->{coeff} }

sub set_smoothing_coeff {
    my ($self, $coeff) = @_;
    
    if ($coeff < 0 || $coeff > 1) {
        carp("Invalid smoothing coefficient C $coeff\n");
        return;
    }

    $self->{coeff} = $coeff;
    return 1;
}

1;

__END__

=head1 NAME

Statistics::Descriptive::Smoother - Base module for smoothing statistical data

=head1 SYNOPSIS

  use Statistics::Descriptive::Smoother;
  my $smoother = Statistics::Descriptive::Smoother->instantiate({
           method   => 'exponential',
           coeff    => 0.5,
           data     => [1, 2, 3, 4, 5],
           samples  => [110, 120, 130, 140, 150],
    });
  my @smoothed_data = $smoother->get_smoothed_data();

=head1 DESCRIPTION

This module provide methods to smooth the trend of a series of statistical data.

The methods provided are the C<Exponential> and the C<Weighted Exponential> (see respectively
C<Statistics::Descriptive::Smoother::Exponential> and C<Statistics::Descriptive::Smoother::Weightedexponential>
for more details).

This class is just a factory that will instantiate the object to perform the
chosen smoothing algorithm.

=head1 METHODS

=over 5

=item Statistics::Descriptive::Smoother->instantiate({});

Create a new Smoother object.

This method require several parameters:

=over 5

=item method

Method used for the smoothing. Allowed values are: C<exponential> and C<weightedexponential>

=item coeff

Smoothing coefficient. It needs to be in the [0;1] range, otherwise undef will be reutrned.
C<0> means that the series is not smoothed at all, while C<1> the series is universally equal to the initial unsmoothed value.

=item data

Array ref with the data of the series. At least 2 values are needed to smooth the series, undef is returned otherwise.

=item samples

Array ref with the samples each data value has been built with. This is an optional parameter since it is not used by all the
smoothing algorithm.

=back 

=item $smoother->get_smoothing_coeff();

Returns the smoothing coefficient.

=item $smoother->set_smoothing_coeff(0.5);

Set the smoothing coefficient value. It needs to be in the [0;1] range, otherwise undef will be reutrned.

=back 

=head1 AUTHOR

Fabio Ponciroli

=head1 COPYRIGHT

Copyright(c) 2012 by Fabio Ponciroli.

=head1 LICENSE

This file is licensed under the MIT/X11 License:
http://www.opensource.org/licenses/mit-license.php.

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

=cut

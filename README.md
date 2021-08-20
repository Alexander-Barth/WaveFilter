# WaveFilter

A filter is developed that removes inertia-gravity waves from the initial conditions of the linear shallow-water equations. </p>

These scripts implement the filtering procedure described in:

Filtering inertia-gravity waves from the initial conditions of the linear shallow-water equations

Alexander Barth, Jean-Marie Beckers, Aida Alvera-Azc&aacute;rate, Robert H. Weisberg, 2007, Ocean Modelling, Volume 19, Issues 3-4, 2007, Pages 204-218, doi: [10.1016/j.ocemod.2007.06.007](https://doi.org/10.1016/j.ocemod.2007.06.007)

Also available via [research gate](https://www.researchgate.net/publication/248496593_Filtering_inertia-gravity_waves_from_the_initial_conditions_of_the_linear_shallow_water_equations) and [orbi](http://hdl.handle.net/2268/4266).


## Download

The latest version of this program can be obtained at https://github.com/Alexander-Barth/WaveFilter

This program is released under the terms of the GNU General Public License version 2 or later.

## Requirements


The code works in [GNU Octave](http://www.octave.org) or for Matlab.
If you don't have Matlab installed, I recommend you to use GNU Octave.

The code has been tested with Octave 2.9.12 and Matlab 2013a on Linux.

## Installation


Extract the archive and add the script directory to your path by:


    addpath('/absolute/path/to/WaveFilter');

with obvious substitution.

## Testing


Run filter_test.

``` matlab
filter_test
```

You should see:

```
Testing filter_vort: OK
Testing filter_var:  OK
```


## Documentation


See the paper: "Filtering inertia-gravity waves from the initial conditions of the linear shallow-water equations", [10.1016/j.ocemod.2007.06.007](https://doi.org/10.1016/j.ocemod.2007.06.007)  and `help filter_vort` and `help filter_var` in Octave or Matlab.


## Changes log

* version 1.0: Initial release. Support for curvilinear grids.
* version 1.1: Several performance optimizations.



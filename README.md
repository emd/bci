Python tools for retrieving DIII-D bi-color interferometer (bci) signals.


Installation:
=============


... on GA's Iris cluster:
-------------------------
Package management is cleanly handled on Iris via
[modules](https://diii-d.gat.com/diii-d/Iris#Environment_modules).
The `bci` package has a corresponding modulefile
[here](https://github.com/emd/modulefiles).

To use the `bci` package, change to the directory
you'd like to download the source files to and
retrieve the source files from github by typing

    $ git clone https://github.com/emd/bci.git

The created `bci` directory defines the
package's top-level directory.
The modulefiles should be similarly cloned.

Now, at the top of the corresponding
[modulefile](https://github.com/emd/modulefiles/blob/master/bci),
there is a TCL variable named `bci_root`;
this must be altered to point at the
top-level directory of the cloned `bci` package.
That's it! You shouldn't need to change anything else in
the modulefile. The `bci` module can
then be loaded, unloaded, etc., as is discussed in the
above-linked Iris documentation.

The modulefile also defines a series of automated tests
for the `bci` package. Run these tests at the command line
by typing

    $ test_bci

If the tests return "OK", the installation should be working.
(Currently, no automated tests are implemented).


... elsewhere:
--------------
Define an environmental variable `$bci_path` specifying
the appropriate MDSplus server's
[tree-path definitions](http://www.mdsplus.org/index.php?title=Documentation:Tutorial:RemoteAccess&open=51668177299325667246079&page=Documentation%2FThe+MDSplus+tutorial%2FRemote+data+access+in+MDSplus)
by, for example, adding the following to your `.bashrc`

    $ export bci_path='atlas.gat.com::/data/usershots/~t\;atlas.gat.com::/data/shots/~t/~f~e/~d~c\;atlas.gat.com::/data/orphans/\;atlas.gat.com::/data/models/~t'

Now, change to the directory you'd like to download the source files to
and retrieve the source files from github by typing

    $ git clone https://github.com/emd/bci.git

Change into the `bci` top-level directory by typing

    $ cd bci

For accounts with root access, install by running

    $ python setup.py install

For accounts without root access (e.g. a standard account on GA's Venus
cluster), install locally by running

    $ python setup.py install --user

To test your installation, run

    $ nosetests tests/

If the tests return "OK", the installation should be working.


Use:
====
The *single-pass* phase measured by
the DIII-D bi-color interferometer (BCI)
can be readily retrieved via the
the `bci.signal.Signal` class. For example, use:

```python
import bci

shot = 169572
tlim = [0, 2]        # [tlim] = s

sig_V2 = bci.signal.Signal(shot, chord='V2', beam='CO2', tlim=tlim)

```

to retrieve the phase signal from the `'V2'` `'CO2'` beam.
Valid values of `chord` are `{'V1', 'V2', 'V3', 'R0'}` and
valid values of `beam` are `{'CO2', 'HeNe'}`.
If `beam` is `'CO2'`, the `vibration_subtracted` keyword
can also be set to `True` to remove vibrational contributions
to the CO2-measured phase; that is

```python
sig_V2_no_vib = bci.signal.Signal(
    shot, chord='V2', beam='CO2', tlim=tlim,
    vibration_subtracted=True)

```

Vibrational contributions to the CO2-measured phase
are typically small for frequencies above 10 kHz.

The autospectral density of V2-measured phase
can then be computed and easily visualized using the
[random_data package](https://github.com/emd/random_data).
Specifically,

```python
import random_data as rd

# Spectral-estimation parameters
Tens = 5e-3         # Ensemble time length, [Tens] = s
Nreal_per_ens = 10  # Number of realizations per ensemeble

# Compute autospectral density
asd_V2 = rd.spectra.AutoSpectralDensity(
    sig_V2.x, Fs=sig_V2.Fs, t0=sig_V2.t0,
    Tens=Tens, Nreal_per_ens=Nreal_per_ens)

asd_V2.plotSpectralDensity(
    flim=[50e3, 400e3],
    xlabel='$t \, [\mathrm{s}]$',
    ylabel='$f \, [\mathrm{Hz}]$')

```

![autospectral_density_V2](https://raw.githubusercontent.com/emd/bci/master/figs/autospectral_density_V2.png)

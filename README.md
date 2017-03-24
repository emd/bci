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

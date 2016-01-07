Simply unpack the tar file, enter the source directory and type:

    make

Move the `nw` executable to somewhere in your path.

Move the `data/mdm78.mat` Dayhoff matrix to a suitable location
(e.g. `/usr/local/lib/data`) and set the environment variable `DATADIR` to
point to that directory.

Under C-shell:

    setenv DATADIR /usr/local/lib

Under Bourne shell:

    DATADIR=/usr/local/lib
    export DATADIR

Type:

    nw -h

to get help on command line options.

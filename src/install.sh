#!/bin/bash

meson setup builddir --prefix=$FOAM_USER_LIBBIN --libdir=.
meson install -C builddir/

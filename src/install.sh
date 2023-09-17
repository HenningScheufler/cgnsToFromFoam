#!/bin/bash

meson setup builddir --prefix=$PWD/../deps --libdir=lib
meson install -C builddir/

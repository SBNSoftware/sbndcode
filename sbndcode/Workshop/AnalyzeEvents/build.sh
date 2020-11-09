#!/bin/bash

cd $MRB_BUILDDIR
mrbsetenv && mrb i -j4 && mrbslp
cd -

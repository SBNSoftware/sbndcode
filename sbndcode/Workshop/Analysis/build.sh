#!/bin/bash

cd $MRB_BUILDDIR
mrbsetenv && mrb i -j3 && mrbslp
cd -

#!/bin/bash

cd $MRB_BUILDDIR
mrbsetenv && mrb i -j64 && mrbslp
cd -

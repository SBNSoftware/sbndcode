#!/bin/bash

cd $MRB_BUILDDIR
mrbsetenv && mrb i -j2 && mrbslp
cd -

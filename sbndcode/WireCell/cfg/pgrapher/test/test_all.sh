#!/bin/bash

mydir=$(dirname $(readlink -f $BASH_SOURCE))
pgdir=$(dirname $mydir)
cfgdir=$(dirname $pgdir)

errors=""
for jsfile in $mydir/test*.jsonnet
do
    dat="$(jsonnet -J $cfgdir $jsfile 2>/dev/null)"
    if [ -z "$dat" -o "$?" != "0" ] ; then
        errors="$jsfile $errors"
    fi
done

if [ -n "$errors" ] ; then
    echo "Failed tests:"
    for one in $errors; do
        echo -e "\t$one"
    done
    exit 1
fi

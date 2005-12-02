#!/bin/bash
for file in *.f; do
    sed -e "/c___c/,/c___c/ d" $file > temp.f
    cat copyright.txt temp.f > $file
done
rm -f temp.f
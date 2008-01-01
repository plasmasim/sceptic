#!/bin/bash
# Attach the copyright and version information at the top of files.
# If $1= -i use individual arguments. If $1 = anything else. Use that.
if [ ".$1" == . ] ; then 
    echo Getting the overall sceptic version from CVS.
    VERSION=`cvs log sceptic.F | grep head: | sed -e "s/head://"`
fi
for file in *.f *.F ; do
    echo -n "Versioning file $file   "
    cp $file $file.bak
    sed -e "/c___c/,/c___c/ d" $file > temp.tmp
    cat >$file <<EOF
c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___
c
c This code is copyright(c) (2003-8) Ian H Hutchinson hutch@psfc.mit.edu
c with major contributions by Leonardo Patacchini. 
c
c  It may be used freely with the stipulation that any scientific or
c scholarly publication concerning work that uses the code must give an
c acknowledgement referring to the papers I.H.Hutchinson, Plasma Physics
c and Controlled Fusion, vol 44, p 1953 (2002), vol 45, p 1477 (2003).
c  The code may not be redistributed except in its original package.
c
c No warranty, explicit or implied, is given. If you choose to build
c or run the code, you do so at your own risk.
EOF
if [ ".$1" != . ] ; then 
    if [ "$1" == "-i" ] ; then
	VERSION=`cvs log $file | grep head: | sed -e "s/head://"`
    else
	VERSION="$1" 
    fi
fi
echo Version: $VERSION
echo "c     Version: $VERSION   `date`" >> $file
echo c >> $file
echo c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___ >> $file
    cat temp.tmp >> $file
    rm $file.bak
done
rm -f temp.tmp

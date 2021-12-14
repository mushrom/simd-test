#!/bin/sh

if [ -z "$OPT" ]; then
	export OPT="0";
fi

if [ -z "$BITS" ]; then
	export BITS="24"
fi

case `uname -m` in
	i586|i686)
		export CFLAGS="-O$OPT -msse2 -DUSESSE -DBITS=$BITS"
		;;

	x86_64)
		export CFLAGS="-O$OPT -mavx -DUSESSE -DUSEAVX -DBITS=$BITS"
		;;

	*)
		export CFLAGS="-O$OPT -mfpu=neon -DUSENEON -DBITS=$BITS"
esac

g++ -c quads quads.cpp $CFLAGS
g++ -o quads quads.o   $CFLAGS

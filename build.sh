#!/bin/sh

case `uname -m` in
	i586|i686)
		g++ -o quads quads.cpp -O0 -msse -DUSESSE
		;;

	x86_64)
		g++ -o quads quads.cpp -O0 -mavx -DUSESSE -DUSEAVX
		;;

	*)
		g++ -o quads quads.cpp -O0 -mfpu=neon -DUSENEON
esac

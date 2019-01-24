# test program for entire udouble package
utest: udouble.h utest.cc
	g++ -o utest -Wall -ansi -pedantic -static utest.cc -lm

# test package for only that part of the udouble package that
# is featured in the article in the C/C++ Users Journal 3/96
udmstest:
	g++ -o udmstest -Wall -ansi -pedantic udmstest.cc -lm

clean:
	-rm -f core *.o *.lout *.lob a.out

tar: udouble.tar

udouble.tar: udms.h udmstest.cc udouble.h utest.cc makefile
	tar cf udouble.tar udms.h udmstest.cc udouble.h utest.cc makefile

z: udouble.tar.Z

Z: udouble.tar.Z

udouble.tar.Z: udouble.tar
	compress udouble.tar

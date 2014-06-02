
#CC=clang
#CXX=clang++
CC=gcc
CXX=g++
CFLAGS=-g -W -Wall -O0 -DDEBUG 
#CFLAGS= -O3  -fomit-frame-pointer  -march=native

INCCDS=./libsdsl/include/sdsl
INCDIVSUF=./libsdsl/include

all: ./libsdsl/lib/libsdsl.a libsdsl/lib/libdivsufsort.a fmbuild fmcount fmlocate fmrecover fmextract 

./libsdsl/lib/libsdsl.a: 
	./sdsl/install.sh
	
./libsdsl/lib/libdivsufsort.a: 
	./sdsl/install.sh
	
libfmindex.a:  util.o  FM.o 
	ar rcs libfmindex.a util.o FM.o
	
fmbuild: fmbuild.cpp ./libsdsl/lib/libdivsufsort.a ./libsdsl/lib/libsdsl.a libfmindex.a
	$(CXX) -I $(INCCDS) -I $(INCDIVSUF) $(CFLAGS) -o fmbuild fmbuild.cpp libfmindex.a ./libsdsl/lib/libdivsufsort.a ./libsdsl/lib/libsdsl.a 

fmcount: fmcount.cpp ./libsdsl/lib/libdivsufsort.a ./libsdsl/lib/libsdsl.a libfmindex.a
	$(CXX) -I $(INCCDS) -I $(INCDIVSUF) $(CFLAGS) -o fmcount fmcount.cpp libfmindex.a ./libsdsl/lib/libdivsufsort.a ./libsdsl/lib/libsdsl.a 
	
fmlocate: fmlocate.cpp ./libsdsl/lib/libdivsufsort.a ./libsdsl/lib/libsdsl.a libfmindex.a
	$(CXX) -I $(INCCDS) -I $(INCDIVSUF) $(CFLAGS) -o fmlocate fmlocate.cpp libfmindex.a ./libsdsl/lib/libdivsufsort.a ./libsdsl/lib/libsdsl.a

fmextract: fmextract.cpp ./libsdsl/lib/libdivsufsort.a ./libsdsl/lib/libsdsl.a libfmindex.a
	$(CXX) -I $(INCCDS) -I $(INCDIVSUF) $(CFLAGS) -o fmextract fmextract.cpp libfmindex.a ./libsdsl/lib/libdivsufsort.a ./libsdsl/lib/libsdsl.a

fmrecover: fmrecover.cpp ./libsdsl/lib/libdivsufsort.a ./libsdsl/lib/libsdsl.a libfmindex.a
	$(CXX) -I $(INCCDS) -I $(INCDIVSUF) $(CFLAGS) -o fmrecover fmrecover.cpp libfmindex.a ./libsdsl/lib/libdivsufsort.a ./libsdsl/lib/libsdsl.a


	
# pattern rule for all objects files
%.o: %.c *.h
	$(CXX) -I $(INCCDS) -I $(INCDIVSUF) -c $(CFLAGS) $< -o $@

%.o: %.cpp *.h
	$(CXX) -I $(INCCDS) -I $(INCDIVSUF) -c $(CFLAGS) $< -o $@

clean:
	rm -f fmbuild fmcount fmlocate fmextract fmrecover libfmindex.a *~ *.o ;
	
cleanall:
	rm -f fmbuild fmcount fmlocate fmextract fmrecover libfmindex.a *~ *.o ; rm -rf libsdsl; cd sdsl/build; ./clean.sh ; cd ../..;

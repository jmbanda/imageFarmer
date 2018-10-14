iDist : main.o idistance.o helper.o
	g++ -Wall main.o helper.o idistance.o -o iDist
	\mv iDist ../run/iDist

main.o: main.cpp helper.h idistance.h 
	g++ -c main.cpp

helper.o: helper.cpp helper.h
	g++ -c helper.cpp

idistance.o: idistance.cpp idistance.h definitions.h
	g++ -c idistance.cpp

clean:
	\rm *.o *~ iDist

tar:
	tar cfv idistance.tar main.cpp helper.h helper.cpp idistance.h idistance.cpp definitions.h makefile readme.txt


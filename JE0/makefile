EXEC = thermoBlade

CLIB = Flib.a -L./lib -Flib -I./include						

$(EXEC): src/solvers/infoBlade/infoBlade.cpp
	g++ -c lib/blockMeshGenerator/blockMeshGenerator.cpp
	g++ -c lib/ODE_solver/solver1.cpp
	g++ -c lib/pathFindeer/pathFindeer.cpp
	gcc -c third_party/sqlite-amalgamation-3460100/shell.c
	gcc -c third_party/sqlite-amalgamation-3460100/sqlite3.c


	ar rcs Flib.a blockMeshGenerator.o solver1.o pathFindeer.o shell.o sqlite3.o

	g++ -o $@ $^ $(CLIB)		

	rm blockMeshGenerator.o solver1.o pathFindeer.o Flib.a shell.o sqlite3.o											

clean:
	rm -f $(EXEC)
.PHONY: clean	

clear:
	rm blockMeshGenerator.o solver1.o pathFindeer.o Flib.a shell.o sqlite3.o infoBlade.o										
.PHONY: clear
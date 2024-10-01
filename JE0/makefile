EXEC = thermoBlade

DEBUG = debug

CLIB = Flib.a -L./lib -Flib -I./include						

$(EXEC): src/solvers/simBlade/simBlade.cpp
	g++ -c lib/blockMeshGenerator/blockMeshGenerator.cpp
	g++ -c lib/ODE_solver/solver1.cpp
	g++ -c lib/pathFindeer/pathFindeer.cpp
	gcc -c third_party/sqlite-amalgamation-3460100/shell.c
	gcc -c third_party/sqlite-amalgamation-3460100/sqlite3.c
	g++ -c src/solvers/infoBlade/infoBlade.cpp
	g++ -c src/solvers/thermoBlade/thermoBlade.cpp
	g++ -c src/solvers/aeroBlade/aeroBlade.cpp

	ar rcs Flib.a blockMeshGenerator.o solver1.o pathFindeer.o shell.o sqlite3.o infoBlade.o thermoBlade.o aeroBlade.o

	g++ -O2 -march=native -o $@ $^ $(CLIB)		

	rm blockMeshGenerator.o solver1.o pathFindeer.o Flib.a shell.o sqlite3.o infoBlade.o thermoBlade.o aeroBlade.o			

$(DEBUG): src/solvers/simBlade/simBlade.cpp
	g++ -c lib/blockMeshGenerator/blockMeshGenerator.cpp
	g++ -c lib/ODE_solver/solver1.cpp
	g++ -c lib/pathFindeer/pathFindeer.cpp
	gcc -c third_party/sqlite-amalgamation-3460100/shell.c
	gcc -c third_party/sqlite-amalgamation-3460100/sqlite3.c
	g++ -c src/solvers/infoBlade/infoBlade.cpp
	g++ -c src/solvers/thermoBlade/thermoBlade.cpp
	g++ -c src/solvers/aeroBlade/aeroBlade.cpp

	ar rcs Flib.a blockMeshGenerator.o solver1.o pathFindeer.o shell.o sqlite3.o infoBlade.o thermoBlade.o aeroBlade.o

	g++ -fsanitize=address -g -o $@ $^ $(CLIB)		

	rm blockMeshGenerator.o solver1.o pathFindeer.o Flib.a shell.o sqlite3.o infoBlade.o thermoBlade.o aeroBlade.o
 
clean:
	rm -f $(EXEC)
.PHONY: clean	

clear:
	rm blockMeshGenerator.o solver1.o pathFindeer.o Flib.a shell.o sqlite3.o infoBlade.o thermoBlade.o aeroBlade.o					
.PHONY: clear
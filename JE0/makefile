EXEC = thermoBlade

CLIB = Flib.a -L./lib -Flib -I./include

$(EXEC): src/solvers/thermoBlade/thermoBlade.cpp
	g++ -c lib/blockMeshGenerator/blockMeshGenerator.cpp
	g++ -c lib/ODE_solver/solver1.cpp
	g++ -c lib/pathFindeer/pathFindeer.cpp

	ar rcs Flib.a blockMeshGenerator.o solver1.o pathFindeer.o

	g++ -o $@ $^ $(CLIB)		

	rm blockMeshGenerator.o solver1.o pathFindeer.o Flib.a													

clean:
	rm -f $(EXEC)
.PHONY: clean	
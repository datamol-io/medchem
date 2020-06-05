# Location for the installation

.PHONY: supporting_libraries
.PHONY: Molecule
.PHONY: test
.PHONY: all

all: supporting_libraries Molecule

supporting_libraries:
	cd lilly/supporting_libraries && make install

Molecule:
	cd lilly/Molecule && make install

install: supporting_libraries Molecule
	@echo "No separate installation required" >&2

clean:
	cd lilly/supporting_libraries && make clean
	cd lilly/Molecule && make clean

uninstall:
	cd lilly/supporting_libraries && make uninstall
	cd lilly/Molecule && make uninstall

all:
	pip3 install -e .
doxygen:
	doxygen doc/Doxyfile

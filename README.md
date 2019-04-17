
# VHEGEN

[VHEGEN: A Vibronic Hamiltonian Expansion GENerator for trigonal and tetragonal polyatomic systems](placeholder_link).

* VHEGEN is a Python package capable of generating the Hamiltonian operator expressions for all possible Jahn-Teller and pseudo Jahn-Teller problems of trigonal and tetragonal symmetry, for up to two vibrational modes, at an arbitrarily high order of expansion. VHEGEN may be used as an executable script or as an importable package.

## Compatibility 

Python 2.7+/3.5+ in Linux, Windows, and macOS are supported.

## Dependencies

VHEGEN requires base Python 2.7+/3.5+ and Python computer algebra system module Sympy. It is highly recommended the system also has a LaTeX/TeX distribution installed, such as TeX Live. A LaTeX/TeX installation will allow VHEGEN to execute `pdflatex` on the output TEX file and generate a LaTeX-typeset output PDF.

Required:

* [Python 2.7+/3.5+](https://www.python.org/)
* [sympy](https://www.sympy.org)
* [numpy](https://www.numpy.org/)
* [a TeX distribution](http://www.tug.org/interest.html#free)

If pip is installed, the Python dependencies for VHEGEN can be installed via `pip install -r requirements.txt` in the main package directory.

# User Guide

## Configuration

By default, the program settings are read from `config.cfg` in the base `vhegen` directory. This file must be edited manually if desired. Alternatively, the user may specify a custom configuration path to be read by typing `python vhegen.py --c=directory/filename`.

Configuration file arguments:

* `"input"`
	* `"static"`
		* The problem parameters and filename are handed to VHEGEN as additional arguments. When `input=static`, the program must be called with additional arguments:
			* `"--sym="`
				* specifies the problem's symmetry point group classification in static input mode.
			* `"--states"`
				* specifies the electronic states involved in the problem in static input mode.
			* `"--modes"`
				* specifies the vibrational modes involved in the problem in static input mode.
			* `"--o"`
				* specifies the desired order(s) of expansion in static input mode.
			* `"--f"`
				* specifies the name given to output files in static input mode. If not specified, the output filename will default to `"output"`.

		* Example: To specify 3rd-to-6th order expansions of the (E+A1)⊗(e+a2) problem in C3v symmetry in static input mode, type `"python vhegen.py --sym=c3v --states=E+A1 --modes=e+a2 --o=3,6 --f=example"`.

	* `"dynamic"`
		* The problem parameters and filename are entered through dynamic input prompts after calling the program by typing `"python vhegen.py"` (with `--c` argument if a custom config path is desired). Thus, any of the additional arguments mentioned above for `input=static` will be ignored when in dynamic input mode.
		* Type `"exit"` at any input prompt to exit dynamic input and cancel the program run.
		* Type `"help"` at an input prompt to receive information about valid entries.
* `"e_coords"`
	* `"cart"`
		* Describe cartesian
	* `"pol"`
		* Describe polar

## Input arguments

* Symmetry
	* The problem's symmetry must be specified by a conventional trigonal or tetragonal point group classification. Case insensitive.

* Electronic states (vibrational modes)
	* Electronic state(s) (vibrational mode(s)) are specified by standard Mulliken convention for irreducible representations under the symmetry point group of interest. Case insensitive. The number of electronic states (vibrational modes) must be either 1 or 2. If 2 states (modes) are specified, split the states by either `"+"` or `","`, including zero spaces.

* Orders of expansion
	* Order of Hamiltonian matrix element expansion may be any non-negative integer. If multiple orders of expansion are desired, one can specify a range of expansion orders m-n by typing `"m,n"`, where m must be less than n.

* Filename
	* The output filename may be any string. It is recommended to avoid special characters and spaces. This input may be left blank to default the output filename to `"output"`.

	* WARNING: Specifying a filename already used (including `"output`) in the output directory `vhegen/outputs` will automatically overwrite the previous outputs of that name. Ensure you are not overwriting output files you wish to keep before executing.

## Output files

* `.tex` output

* `.log` output

* `.pdf` output

## VHEGEN as a package

The VHEGEN as a package section.

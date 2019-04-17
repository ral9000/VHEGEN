# VHEGEN: A Vibronic Hamiltonian Expansion GENerator for trigonal and tetragonal polyatomic systems

[future link to paper](placeholder).

* VHEGEN (V-ibronic H-amiltonian E-xpansion GEN-erator) is a Python package capable of symbolically generating arbitrarily high order expansion formulas for vibronic Hamiltonians. The program covers all >6000 Jahn-Teller and pseudo-Jahn-Teller problems in trigonal and tetragonal symmetries with up to two vibrational modes.

## Compatibility 

Python 2.7+/3.5+ in Linux, Windows, and macOS are supported.

## Dependencies

* [Python 2.7+/3.5+](https://www.python.org/)
* [sympy](https://www.sympy.org)
* [numpy](https://www.numpy.org/)
* [a TeX distribution](http://www.tug.org/interest.html#free)

VHEGEN requires base Python 2.7+/3.5+ and external libraries Sympy and Numpy. It is highly recommended the system also has a TeX distribution installed, such as TeX Live. A TeX installation will allow VHEGEN to execute `pdflatex` on the output `.tex` file and generate a LaTeX-typeset output `.pdf` file. If `pip` is installed, the Python dependencies for VHEGEN can be installed via `pip install -r requirements.txt` in the main package directory.

# Using VHEGEN procedurally

VHEGEN is called procedurally by executing `python vhegen.py` in the main package directory. The program settings for procedural use are read from the configuration file `config.cfg` in the same directory, the contents of which are described below.

## Configuration

 In the configuration file, the user specifies the five parameters: `input`, `pdf_out`, `log_out`, `e_coords`, and `basis`. 
 
 | Parameter  | Values              |
|------------|---------------------|
| `input`    | `static`, `dynamic` |
| `pdf_out`  | `true`, `false`     |
| `log_out`  | `true`, `false`     |
| `e_coords` | `pol`, `cart`       |
| `basis`    | `complex`, `real`   |
 
 Parameter `input` is used to specify either static or dynamic input mode (discussed in the next section). When set to `true`(`false`), parameter `pdf_out` enables(disables) application of the TeX `pdflatex` command to produce a `.pdf` output file. The output file contains the input vibronic interaction, the matrix form of the vibronic Hamiltonian, and all explicit matrix element expansions typeset via LaTeX. Thus, `pdf_out` should only be set to `true` if a TeX distribution is installed on the system. Parameter `log_out` similarly enables or disables output of a text file, containing all relevant information used in the expansion process, including: the independent matrix elements and their symmetry eigenvalues, and the root expansion formulas along with their constraints. All expansions in Sympy readable syntax can be found in the `log` file. Parameter `e_coords` is used to specify the coordinate system used for expressing e-type vibrational modes. When `e_coords` is set to `pol`, the expansions will be kept in polar coordinates as they were originally constructed. When `e_coords` is set to `cart`, the expansions are converted to cartesian coordinates. Both polar and cartesian coordinate expansions may be included in the final output by setting `e_coords=both`. The `basis` parameter defines which basis to use for vibronic Hamiltonians involving E-type electronic states. When set to `complex`/`real`, the output vibronic Hamiltonian matrix elements are for the complex/real E component states. Both complex and real representations can be output by setting `basis=both`.

## Input

* Symmetry
	* The problem's symmetry must be specified by a conventional trigonal or tetragonal point group classification. Case insensitive.

## Output

* `.tex` output

* `.log` output

* `.pdf` output

# Using VHEGEN as a package

The VHEGEN as a package section.

# Contributing



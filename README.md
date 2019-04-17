# VHEGEN: A Vibronic Hamiltonian Expansion GENerator for trigonal and tetragonal polyatomic systems

[future link to paper](placeholder).

* VHEGEN (V-ibronic H-amiltonian E-xpansion GEN-erator) is a Python package capable of symbolically generating arbitrarily high order expansion formulas for vibronic Hamiltonians. The program covers all >6000 Jahn-Teller and pseudo-Jahn-Teller problems in trigonal and tetragonal symmetries with up to two vibrational modes.

## Scope 

VHEGEN can currently generate vibronic Hamiltonians for all (P)JT problems of trigonal and tetragonal symmetries including up to two vibrational modes. Therefore it covers all unimodal and bimodal problems in point groups:
* C3
* C3v
* C3h
* D3
* D3h
* D3d
* C4
* S4
* C4v
* C4h
* D4
* D4h
* D2d

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
 
 Parameter `input` is used to specify either static or dynamic input mode (discussed in the next section). When set to `true`(`false`), parameter `pdf_out` enables(disables) application of the TeX `pdflatex` command to produce a `.pdf` output file. The output file contains the input vibronic interaction, the matrix form of the vibronic Hamiltonian, and all explicit matrix element expansions typeset via LaTeX. Thus, `pdf_out` should only be set to `true` if a TeX distribution is installed on the system. Parameter `log_out` similarly enables or disables output of a text file, containing all relevant information used in the expansion process, including the independent matrix elements and their symmetry eigenvalues, and the root expansion formulas along with their constraints. All expansions in Sympy readable syntax can be found in the `log` file. Parameter `e_coords` is used to specify the coordinate system used for expressing e-type vibrational modes. When `e_coords` is set to `pol`, the expansions will be kept in polar coordinates as they were originally constructed. When `e_coords` is set to `cart`, the expansions are converted to cartesian coordinates. Both polar and cartesian coordinate expansions may be included in the final output by setting `e_coords=both`. The `basis` parameter defines which basis to use for vibronic Hamiltonians involving E-type electronic states. When set to `complex`/`real`, the output vibronic Hamiltonian matrix elements are for the complex/real E component states. Both complex and real representations can be output by setting `basis=both`.

## Input

The procedural input to VHEGEN consists of specifying: point group symmetry; irreps of electronic states; irreps of vibrational modes; order(s) of expansion; and an output filename. There exist two input approaches, namely the "dynamic input" and "static input" modes, selected using the `input` parameter in `config.cfg`.

Both modes of input follow the same general rules for input syntax. All input parameters are case-insensitive, except the output filename. Symmetries are specified by their usual point group classification. Electronic states and vibrational modes are specified by standard Mulliken symbols for their irreducible representations. Single primes in Mulliken symbols are denoted by an apostrophe `'`, and double primes are denoted by two apostrophes `''` or a quotation mark `"`. For pJT problems that involve two states or bimodal problems that involve two vibrational modes, either a plus sign `+` or a comma `,` can be used to separate the two specified states or modes. When specifying the orders of expansion, a range of orders can be stated by providing non-negative integers as lower and upper inclusive bounds separated by a comma, e.g., `0,6`. If expansion at a single order is desired, then only an integer is keyed in. Output filename should avoid any operating system dependent illegal characters. If the filename parameter is unspecified or left empty, the output filename will default to value `output`.

### Dynamic

When dynamic input mode is specified in the configuration file, the user enters problem parameters and output filename dynamically through terminal prompts after executing `python vhegen.py`. Any additional arguments made when calling the program will be ignored if set to dynamic input. The user will be re-prompted at each stage if an invalid parameter is entered. The user can receive lists of valid inputs by entering `list`, and can quit the program by entering `exit`.

### Static

When input mode is set to static, a vibronic problem is specified in-line with the execution of the program via additional arguments. The following additional arguments must be specified in the static input mode: `--sym` for point group symmetry, `--states` for electronic state(s), `--modes` for vibrational mode(s), and `--o` for order(s) of expansion. Specification of a filename is done by optional argument `--f`. The general syntax for specifying the additional arguments when executing the program is `--arg=value`. For example, execution of problem E x (e+a) in C3v symmetry at the third to sixth orders is accomplished in static input mode by `python vhegen.py --sym=C3V --states=E --modes=e+a1 --o=3,6`.

## Output

* `.tex` output

* `.log` output

* `.pdf` output

# Using VHEGEN as a package

The VHEGEN as a package section.

# Contributing



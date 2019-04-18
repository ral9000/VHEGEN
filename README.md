# VHEGEN: A vibronic Hamiltonian expansion generator for trigonal and tetragonal polyatomic systems

VHEGEN (V-ibronic H-amiltonian E-xpansion GEN-erator) is a Python package capable of symbolically generating arbitrarily high order expansion formulas for Jahn-Teller and pseudo-Jahn-Teller Hamiltonians in vibrational coordinates. 

## Scope 

VHEGEN covers all bimodal Jahn-Teller and pseudo-Jahn-Teller problems in trigonal and tetragonal symmetries, i.e) point groups: C3, C3v, C3h, D3, D3h, D3d, C4, S4, C4v, C4h, D4, D4h, and D2d. Unimodal problems are treated as special cases of their bimodal analogues and hence also covered.

## Compatibility 

VHEGEN has been tested with Python 2.7 and 3.7 on Linux, macOS, and Windows 10 operating systems.

## Dependencies

* [Python 2.7/3.7](https://www.python.org/)
* [SymPy](https://www.sympy.org)
* [NumPy](https://www.numpy.org/)
* [a TeX distribution](http://www.tug.org/interest.html#free)

VHEGEN requires base Python 2 or 3 and external libraries SymPy and NumPy. It is highly recommended the system also has a TeX distribution installed, such as TeX Live. A TeX installation will allow VHEGEN to execute `pdflatex` on the output `.tex` file and generate a LaTeX-typeset output `.pdf` file. If `pip` is installed, the Python dependencies for VHEGEN can be installed via `pip install -r requirements.txt` in the main package directory.

## Installation

Installation of the VHEGEN package is as simple as uncompressing the master `VHEGEN` folder to the user's preferred write-allowed directory. This primes the program for procedural use as an executable script, as well as for importing as a package local to the installation directory.

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

When called procedurally, outputs generated by VHEGEN are found in subdirectory `VHEGEN/outputs`. VHEGEN may produce three output files upon completion of a problem:

* `.tex` file: Always produced. This file contains all final matrix element expansions in the basis, coordinate system, and orders of expansion(s) specified, the matrix form of the vibronic Hamiltonian, and a count of free parameters needed to be fitted. All information is typeset to a compilable LaTeX document.

* `.log` file: Produced if `log_out=true` in `config.cfg`. This is a text file containing auxiliary information regarding the matrix element expansion process, including the independent matrix elements, their symmetry eigenvalues, their root formulas, and the appropriate contraints. It also contains all final matrix element expansions in Sympy syntax.

* `.pdf` file: Produced if `pdf_out=true` in `config.cfg`. This file is a read-friendly version of the `.tex` output file compiled via `pdflatex`.

# Using VHEGEN as a package

Herein we describe important methods and attributes of the VHEGEN class, made accessible by importing VHEGEN as a package via `import vhegen`. 

## Initialization

To initialize an instance of the VHEGEN class, a dictionary containing all vibronic problem parameters for the instance must be specified as its only argument. The input may be prepared by the `prepare_input` function from the `input` module, which takes the five parameters: point group symmetry, the electronic state(s), the vibrational mode(s), order(s) of expansion, and an optional filename -- all given as strings following the syntax rules for procedural input. E.g.) initializing an instance of VHEGEN for the (E+A)x(e+a) problem in C4 symmetry at 12th order is shown below.
~~~~
import vhegen as vhe

params = vhe.inp.prepare_input(sym='C4', states='E+A', modes='e+a', orders='12', filename='output')

vhegen_instance = vhe.VHEGEN(params)
~~~~

## Generating matrix element expansions

Below are the methods which must be sequentially called to generate the full matrix element expansions for a specified vibronic problem. For more details about the methods described, please see Section 5.3 in the associated paper.

| Method                              | Description                                                                                                                                                                                                                       |
|-------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `set_e_coords(e_coordinate_system)` | Defines the coordinate system for problems involving e-type vibrational modes. Allowed arguments are `cart`, `pol`, or `both`.                                                                                                    |
| `set_basis(E_component_basis)`      | Defines the E state component basis for problems involving E-type states. Allowed arguments are `complex`, `real`, and `both`.                                                                                                    |
| `get_eigenvals()`                   | Obtains unique matrix elements and their symmetry eigenvalues by initializing attribute `eigenvals`.                                                                                                                              |
| `get_matrix_form()`                 | Obtain matrix product form of the vibronic Hamiltonia, and stores in attribute `matrix`.                                                                                                                                          |
| `get_formulas()`                    | Performs lookup of root expansion formulas along with the required constraints, storing them in attributes `formulas` and `constraints` respectively.                                                                             |
| `get_expansions()`                  | Generate term-by-term expansions at all specified orders, storing them in a nested dictionary in attribute `expansions`. If `basis` is set to `both`, expansions in the real basis will be stored in attribute `real_expansions`. |


## The `auto` method

The methods discussed in the previous subsection can be immediately performed after initialization of a VHEGEN instance by the `auto` method. This method sequentially calls all mandaory processes discussed above to generate the `expansions` attribute including all matrix elements. It will also transform to the real E component basis if applicable. Below is an example script of how one would ggenerate the expanded matrix elements for a specific vibronic Hamiltonian at the desired orders.
~~~~
import vhegen as vhe 

params = vhe.inp.prepare_input(sym=`D4h', states=`Eg+A1u', modes=`eg+b2u', orders=`0,10')

vhegen_instance = vhe.VHEGEN(params)

vhegen_instance.auto()

print(vhegen_instance.expansions)
~~~~

## Expansion output

TBC

## Fitting of expansion coefficients

This subsection will discuss the parameter fitting process.

# Contributing

This section will discuss how others can contribute to the VHEGEN code.

# Authors

* Robert A. Lang (University of Toronto)
* Riley J. Hickman (University of Toronto)
* Tao Zeng (Carleton University)

# Amber force field test suite

This test suite attempts to cover as many residues as possible in the force
fields we support.  To run the tests, use `test_amber.py`:

```
usage: test_amber.py [-h] [--sander] [--openmm-amber] [--openmm-ffxml]
                     [--absolute-tolerance kcal/mol]
                     [--relative-tolerance ratio] [--debug-exceptions]
                     [--dump]
                     [--openmm-platform {automatic,Reference,CPU,OpenCL}]
                     [--ffxml-directory path]
                     [test_path ...]

Test script for Amber force field conversion to OpenMM FFXML format

positional arguments:
  test_path             Test directories containing test specifications
                        (test.json) to run

options:
  -h, --help            show this help message and exit
  --sander              Compute energies with AmberTools using sander
  --openmm-amber        Compute energies with OpenMM using OpenMM
                        AmberPrmtopFile
  --openmm-ffxml        Compute energies with OpenMM using OpenMM ForceField
  --absolute-tolerance kcal/mol
                        Absolute tolerance for matching energies (default
                        0.0002 kcal/mol)
  --relative-tolerance ratio
                        Relative tolerance for matching energies (default
                        0.0002)
  --debug-exceptions    Treat exceptions as fatal errors and display
                        tracebacks rather than marking tests as failures
  --dump                Dump input files for debugging
  --openmm-platform {automatic,Reference,CPU,OpenCL}
                        Platform for OpenMM to use, or "automatic" to let
                        OpenMM choose (default Reference)
  --ffxml-directory path
                        Path to OpenMM force field files (default
                        ../../openmmforcefields/ffxml/amber)
```

Each subdirectory in the `cases/` directory corresponds to a given test topology
that can be run with a number of force fields.  Passing `cases/*` as the
`test_path` argument to `test_amber.py` in your shell (*i.e.*, selecting all
subdirectories) runs the full test suite.

## Limitations

* ff03ua (the Amber united-atom force field) is not currently tested due to LEaP
  failing to parameterize many simple systems.  It is not clear that there is an
  easy way around this.
* Currently failing due to problems with the conversion:
  * ff03 endcaps require manual modification
  * Problems with inconsistent improper peripheral atom ordering
  * Other tests that need to be investigated
* Remaining to be implemented:
  * Tests for water and ions
  * Tests for GLYCAM
  * Tests with combinations of kinds of macromolecules

## Information for test case developers

Because not every residue is supported by every force field of a given kind in
the Amber force field distribution, setting up the test cases is not easily
fully automated.  However, some tools are provided to make this task easier if
it is necessary to expand the test suite to new force fields in the future.
`scan_leaprc.py` inspects the force fields in the AmberTools installation that
it looks for based on the setting of the `AMBERHOME` environment variable:

```
usage: scan_leaprc.py [-h] (--list-leaprc-files | --categorize-residues)

Searches for residues in Amber force fields in an AmberTools installation (set
$AMBERHOME or make it your current working directory)

options:
  -h, --help            show this help message and exit
  --list-leaprc-files   Lists all leaprc files to standard output
  --categorize-residues
                        Reads a list of leaprc files from standard input and
                        categorizes residues as to their membership in these
                        files to standard output
```

The `--categorize-residues` option will list residues along with what force
fields support them, to make it easier to design a single test case for several
residues that can be run against several force fields without extra effort.
You will have to manually match the appropriate set of leaprc files with the
corresponding FFXML files based on the definitions in the YAML files run by the
main conversion script in the containing folder `/amber/`.

To actually build the test cases in `cases/`, run `generate_cases.py`.  This is
not required to run the tests as the case files are distributed with this
repository, but it is only necessary if you want to add a new test case or
modify an existing one (*e.g.*, add a new force field to it).  See the source
code for `generate_cases.py` for examples.  In addition to specifying a set of
chain sequences, you may have to specify some special LEaP commands to build
branched chains, or work around idiosyncracies in residue definitions.  Consult
the AmberTools manual for more details on how to use LEaP.

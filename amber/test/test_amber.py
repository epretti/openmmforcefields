import argparse
import enum
import itertools
import numpy
import openmm.app
import os
import shutil
import subprocess
import sys
import tempfile
import tomllib

# The AMBER value can be found in unitio.c of the LEaP source code.
# The OpenMM value is from SimTKOpenMMRealType.h.
AMBER_VACUUM_PERMITTIVITY = openmm.unit.AVOGADRO_CONSTANT_NA * openmm.unit.elementary_charge ** 2 / (4 * numpy.pi * 18.2223 ** 2 * openmm.unit.kilocalorie_per_mole * openmm.unit.angstrom)
OPENMM_VACUUM_PERMITTIVITY = 8.8541878128e-12 * openmm.unit.farad / openmm.unit.meter
AMBER_ELECTROSTATIC_SCALE = AMBER_VACUUM_PERMITTIVITY / OPENMM_VACUUM_PERMITTIVITY

DEFAULT_ABSOLUTE_TOLERANCE = 1e-4 # kcal/mol
DEFAULT_RELATIVE_TOLERANCE = 1e-4
DEFAULT_PERTURB_DISTANCE = 0.1 # Å
DEFAULT_OPENMM_PLATFORM = "Reference"
DEFAULT_FFXML_DIRECTORY = "../../openmmforcefields/ffxml/amber"

PADDING_DISTANCE = 5.0 # Å
MINIMUM_CUTOFF = 10.0 # Å
HARMONIC_STRENGTH = 10.0
SYSTEM_OPTIONS = dict(rigidWater=False, removeCMMotion=False)

ENERGY_PRINT_UNIT = openmm.unit.kilocalorie_per_mole
DO_COLOR = sys.stdout.isatty()

class ForceGroup(enum.Enum):
    """
    Force groups for decomposing contributions to the total energy.
    """

    BONDS_UB = 1
    ANGLES = 2
    TORSIONS = 3
    CMAP = 4
    NONBONDED = 5

# For each force group, the Amber labels and the OpenMM force types.
FORCE_GROUP_DATA = {
    ForceGroup.BONDS_UB: (("BOND", "UB"), (openmm.HarmonicBondForce,)),
    ForceGroup.ANGLES: (("ANGLE",), (openmm.HarmonicAngleForce,)),
    ForceGroup.TORSIONS: (("DIHED", "IMP"), (openmm.PeriodicTorsionForce,)),
    ForceGroup.CMAP: (("CMAP",), (openmm.CMAPTorsionForce,)),
    ForceGroup.NONBONDED: (("VDWAALS", "EELEC", "1-4 NB", "1-4 EEL"), (openmm.NonbondedForce, openmm.CustomNonbondedForce, openmm.CustomBondForce)),
}

AMBER_ELECTROSTATIC_GROUPS = ("EELEC", "1-4 EEL")

def main():
    """
    Main entry point for test script.  Parses command-line arguments and runs.
    """

    parser = argparse.ArgumentParser(description="Test script for Amber force field conversion to OpenMM FFXML format")

    parser.add_argument(
        "test_path",
        nargs="*",
        help="Test specifications (in a TOML file) to run",
    )

    # Flags to compute energies using particular methods.
    parser.add_argument(
        "--sander",
        action="store_true",
        help="Compute energies with AmberTools using sander",
    )
    parser.add_argument(
        "--openmm-amber",
        action="store_true",
        help="Compute energies with OpenMM using OpenMM AmberPrmtopFile",
    )
    parser.add_argument(
        "--openmm-ffxml",
        action="store_true",
        help="Compute energies with OpenMM using OpenMM ForceField",
    )

    # Flags to control energy comparison.
    parser.add_argument(
        "--absolute-tolerance",
        type=float,
        default=DEFAULT_ABSOLUTE_TOLERANCE,
        metavar="kcal/mol",
        help=f"Absolute tolerance for matching energies (default {DEFAULT_ABSOLUTE_TOLERANCE} kcal/mol)",
    )
    parser.add_argument(
        "--relative-tolerance",
        type=float,
        default=DEFAULT_RELATIVE_TOLERANCE,
        metavar="ratio",
        help=f"Relative tolerance for matching energies (default {DEFAULT_RELATIVE_TOLERANCE})",
    )

    # Flags to control adding replicates involving random coordinate
    # perturbation.
    parser.add_argument(
        "--perturb-replicates",
        type=int,
        default=0,
        metavar="count",
        help="Repeat tests with positions perturbed from minimum energy this many times",
    )
    parser.add_argument(
        "--perturb-distance",
        type=float,
        default=DEFAULT_PERTURB_DISTANCE,
        metavar="Å",
        help=f"Perturb positions by up to this distance (default {DEFAULT_PERTURB_DISTANCE} Å)",
    )
    parser.add_argument(
        "--perturb-seed",
        type=int,
        metavar="seed",
        help="Random seed for initializing configurations and repeating tests with perturbed positions",
    )

    # Miscellaneous options.
    parser.add_argument(
        "--debug-exceptions",
        action="store_true",
        help="Treat exceptions as fatal errors rather than marking tests as failures",
    )
    parser.add_argument(
        "--dump",
        action="store_true",
        help="Dump input files for debugging",
    )
    parser.add_argument(
        "--openmm-platform",
        default=DEFAULT_OPENMM_PLATFORM,
        choices=["automatic"] + [openmm.Platform.getPlatform(index).getName() for index in range(openmm.Platform.getNumPlatforms())],
        help=f"Platform for OpenMM to use, or \"automatic\" to let OpenMM choose (default {DEFAULT_OPENMM_PLATFORM})",
    )

    # Specify directories where parameter files can be found.
    parser.add_argument(
        "--ffxml-directory",
        default=DEFAULT_FFXML_DIRECTORY,
        metavar="path",
        help=f"Path to OpenMM force field files (default {DEFAULT_FFXML_DIRECTORY})",
    )

    run(parser.parse_args())

def run(args):
    """
    Runs the specified tests with the specified options and exits.  The return
    code will be the number of failed tests.

    Parameters
    ----------
    args : argparse.Namespace
        Command-line arguments.
    """

    runner = TestRunner(
        args.sander,
        args.openmm_amber,
        args.openmm_ffxml,
        args.absolute_tolerance * openmm.unit.kilocalorie_per_mole,
        args.relative_tolerance,
        args.perturb_replicates,
        args.perturb_distance,
        numpy.random.default_rng(seed=args.perturb_seed),
        args.dump,
        args.openmm_platform,
        args.ffxml_directory,
    )

    test_results = []
    for test_path in args.test_path:
        print(f"Reading tests from test file {test_path!r}...")
        print()

        with open(test_path, "rb") as test_spec_file:
            test_specs = tomllib.load(test_spec_file)
            for test_spec in test_specs.get("test", []):
                test_name = test_spec.get("name", "<unnamed>")
                test_label = f"{test_path}:{test_name}"
                if args.debug_exceptions:
                    term_failure_count = runner.run_test(test_spec)
                    test_results.append((test_label, term_failure_count == 0, f"{term_failure_count} terms exceeded tolerance"))
                else:
                    try:
                        term_failure_count = runner.run_test(test_spec)
                    except Exception as error:
                        test_results.append((test_label, False, f"{type(error).__name__}: {error}"))
                    else:
                        test_results.append((test_label, term_failure_count == 0, f"{term_failure_count} terms exceeded tolerance"))

    if not test_results:
        print("No tests run (use --help for information on how to specify tests)")
        sys.exit()

    # Summarize test results.
    test_count = len(test_results)
    test_names, success_flags, _ = zip(*test_results)
    success_count = sum(success_flags)
    failure_count = test_count - success_count

    print(f"Ran {test_count} tests ({success_count} successes, {failure_count} failures)")
    test_name_width = max(map(len, test_names))
    for test_name, success_flag, test_message in test_results:
        formatted_line = f"{'Succeeded' if success_flag else 'Failed':9}  {test_message}"
        print(f"    {test_name:{test_name_width}}  {formatted_line if success_flag else color_message(formatted_line)}")

    sys.exit(failure_count)

class TestRunner:
    """
    Holds options used for running tests.

    Parameters
    ----------
    do_sander : bool
        Whether or not to compute energies with AmberTools.
    do_openmm_amber : bool
        Whether or not to compute energies with OpenMM using OpenMM
        AmberPrmtopFile.
    do_openmm_ffxml : bool
        Whether or not to compute energies with OpenMM using OpenMM ForceField.
    absolute_tolerance : openmm.unit.Quantity
        Absolute tolerance for matching energies.
    relative_tolerance : float
        Relative tolerance for matching energies.
    perturb_replicates : int
        Number of sets of perturbed coordinates to test.
    perturb_distance : float
        Maximum distance by which to perturb coordinates.
    perturb_rng: numpy.random.Generator
        The random number generator used to perturb coordinates.
    do_dump : bool
        Whether or not to dump input files for debugging.
    openmm_platform : str
        Name of an OpenMM Platform to use, or "automatic" to let OpenMM choose.
    ffxml_directory : str
        Directory containing OpenMM force field files.
    """

    def __init__(
        self,
        do_sander,
        do_openmm_amber,
        do_openmm_ffxml,
        absolute_tolerance,
        relative_tolerance,
        perturb_replicates,
        perturb_distance,
        perturb_rng,
        do_dump,
        openmm_platform,
        ffxml_directory,
    ):
        self.do_sander = do_sander
        self.do_openmm_amber = do_openmm_amber
        self.do_openmm_ffxml = do_openmm_ffxml
        self.absolute_tolerance = absolute_tolerance
        self.relative_tolerance = relative_tolerance
        self.perturb_replicates = perturb_replicates
        self.perturb_distance = perturb_distance
        self.perturb_rng = perturb_rng
        self.do_dump = do_dump
        self.openmm_platform = openmm_platform
        self.ffxml_directory = ffxml_directory

    def run_test(self, test_spec):
        """
        Runs the specified test with the stored options.

        Parameters
        ----------
        test_spec : dict
            A test case specification.

        Returns
        -------
        int
            The number of failed checks.
        """

        test_name = test_spec.get("name", "<unnamed>")
        chains = test_spec.get("chains", [])
        leaprcs = test_spec.get("leaprcs", [])
        ffxmls = test_spec.get("ffxmls", [])

        dump_prefix = "".join(character if character.isascii() and character.isalnum() else "_" for character in test_name)
        force_field = openmm.app.ForceField(*(os.path.join(self.ffxml_directory, ffxml) for ffxml in ffxmls))

        with tempfile.TemporaryDirectory() as temp_dir:
            # Generate initial coordinates with LEaP.
            with open(os.path.join(temp_dir, "test.leap"), "w") as leap_file:
                for leaprc in leaprcs:
                    print(f"source {leaprc}", file=leap_file)
                for chain_index, chain in enumerate(chains):
                    sequence = " ".join(chain)
                    print(f"test_chain_{chain_index} = sequence {{ {sequence} }}", file=leap_file)
                sequence = " ".join(f"test_chain_{chain_index}" for chain_index in range(len(chains)))
                print(f"test_chains = sequence {{ {sequence} }}", file=leap_file)
                print("saveAmberParm test_chains test.top test.crd", file=leap_file)
                print("savePdb test_chains test.pdb", file=leap_file)
                print("quit", file=leap_file)
            result = subprocess.run(["tleap", "-f", "test.leap"], cwd=temp_dir, capture_output=True)

            if self.do_dump:
                with open(f"{dump_prefix}_leap.log", "wb") as dump_file:
                    dump_file.write(result.stdout)
                    dump_file.write(result.stderr)
            
            if result.returncode:
                print(result.stdout.decode())
                print(result.stderr.decode())
                raise RuntimeError(f"tleap exited with code {result.returncode}")
        
            # Load PDB output.
            pdb = openmm.app.PDBFile(os.path.join(temp_dir, "test.pdb"))
            positions = numpy.array([position.value_in_unit(openmm.unit.angstrom) for position in pdb.positions])

            # We may have to patch up the PDB topology based on the prmtop if
            # non-standard residues were involved.
            prmtop = openmm.app.AmberPrmtopFile(os.path.join(temp_dir, "test.top"))
            pdb_atoms = list(pdb.topology.atoms())
            pdb_bonds = set()
            for bond in pdb.topology.bonds():
                index_1 = bond.atom1.index
                index_2 = bond.atom2.index
                pdb_bonds.add((min(index_1, index_2), max(index_1, index_2)))
            prmtop_bonds = set()
            for bond in prmtop.topology.bonds():
                index_1 = bond.atom1.index
                index_2 = bond.atom2.index
                prmtop_bonds.add((min(index_1, index_2), max(index_1, index_2)))
            if pdb_bonds - prmtop_bonds:
                raise ValueError("PDB topology has bonds not present in prmtop topology")
            for index_1, index_2 in sorted(prmtop_bonds - pdb_bonds):
                pdb.topology.addBond(pdb_atoms[index_1], pdb_atoms[index_2])

            # Space out chains from each other to prevent overlapping atoms.
            x_ranges = []
            for chain in pdb.topology.chains():
                x_positions = [positions[atom.index][0] for atom in chain.atoms()]
                x_ranges.append((min(x_positions), max(x_positions)))
            x_offsets = []
            for chain_index in range(len(x_ranges)):
                if chain_index:
                    x_offsets.append(PADDING_DISTANCE + x_offsets[chain_index - 1] + x_ranges[chain_index - 1][1] - x_ranges[chain_index][0])
                else:
                    x_offsets.append(0.0)
            for chain, x_offset in zip(pdb.topology.chains(), x_offsets):
                for atom in chain.atoms():
                    positions[atom.index, 0] += x_offset

            # Add some random perturbations to the positions.
            for index in range(positions.shape[0]):
                positions[index] += self.perturb_distance * random_in_unit_ball(self.perturb_rng)
            
            # Create an OpenMM system to generate an initial configuration by
            # minimizing the energy subject to a harmonic force that pulls
            # particles together.
            system = force_field.createSystem(pdb.topology, **SYSTEM_OPTIONS)
            force = openmm.CustomExternalForce(f"{HARMONIC_STRENGTH}*(x*x+y*y+z*z)")
            for particle_index in range(system.getNumParticles()):
                force.addParticle(particle_index)
            system.addForce(force)
            context = self._make_context_from_system(system)
            context.setPositions(positions * openmm.unit.angstrom)
            openmm.LocalEnergyMinimizer.minimize(context)
            coordinates_list = [context.getState(positions=True).getPositions(asNumpy=True)]

            # Generate sets of perturbed positions if desired.
            for _ in range(self.perturb_replicates):
                perturbed_positions = numpy.array(coordinates_list[0].value_in_unit(openmm.unit.angstrom))
                for particle_index in range(perturbed_positions.shape[0]):
                    perturbed_positions[particle_index] += self.perturb_distance * random_in_unit_ball(self.perturb_rng)
                context.setPositions(perturbed_positions * openmm.unit.angstrom)
                context.applyConstraints(context.getIntegrator().getConstraintTolerance())
                coordinates_list.append(context.getState(positions=True).getPositions(asNumpy=True))

            # Run tests and get results.
            print(f"Test {test_name!r}:")
            results_sets = {}
            if self.do_sander:
                print("    (Running sander)")
                results_sets["sander"] = self.get_sander_energies(temp_dir, coordinates_list, f"{dump_prefix}_sander" if self.do_dump else None)
            if self.do_openmm_amber:
                print("    (Running OpenMM reading Amber prmtop file)")
                results_sets["openmm_amber"] = self.get_openmm_amber_energies(temp_dir, coordinates_list, f"{dump_prefix}_openmm_amber" if self.do_dump else None)
            if self.do_openmm_ffxml:
                print("    (Running OpenMM creating system from FFXML)")
                results_sets["openmm_ffxml"] = self.get_openmm_ffxml_energies(temp_dir, coordinates_list, f"{dump_prefix}_openmm_ffxml" if self.do_dump else None, force_field, pdb.topology)
            print()

        if len(results_sets) < 2:
            raise ValueError("must specify at least 2 methods to compare")

        failure_count = 0
        for (name_1, results_1), (name_2, results_2) in itertools.combinations(results_sets.items(), 2):
            print(f"    Comparing {name_1} vs. {name_2}")

            # Print results for energies.
            print(f"        {'Energy error':20}  {'|ΔE| (kcal/mol)':20}  {'Relative':12}")
            for force_group in (None, *ForceGroup):
                difference_data = []
                for result_1, result_2 in zip(results_1, results_2):
                    energy_1 = result_1[force_group]
                    energy_2 = result_2[force_group]
                    difference_data.append((energy_1, energy_2, *self._compare_energies(energy_1, energy_2)))
                failure_count += self._format_energy_difference_data("TOTAL" if force_group is None else force_group.name, difference_data)
            print()

        return failure_count

    def get_sander_energies(self, temp_dir, coordinates_list, dump_prefix):
        """
        Computes energies with AmberTools using sander.

        Parameters
        ----------
        temp_dir : str
            Path to a temporary directory containing test files.
        coordinates_list : list
            A list of arrays of positions.  Energies will be evaluated for each
            set of coordinates.
        dump_prefix : str or None
            Prefix for names of input files to dump for debugging, or None to
            skip dumping.

        Returns
        -------
        list
            For each set of coordinates, a dictionary containing the energy for
            each ForceGroup (or None for the total potential energy).
        """

        # Make an input file for sander.  Choose a cutoff larger than the
        # diagonal of the axis-aligned bounding box of all coordinates (a cheap
        # way to ensure the cutoff will include all pairs).
        with open(os.path.join(temp_dir, "test.nml"), "w") as sander_file:
            cutoff = MINIMUM_CUTOFF
            for coordinates in coordinates_list:
                coordinates_angstrom = coordinates.value_in_unit(openmm.unit.angstrom)
                cutoff = max(cutoff, 0.01 + 1.01 * numpy.linalg.norm(numpy.amax(coordinates_angstrom, axis=0) - numpy.amin(coordinates_angstrom, axis=0)))
            print(f"\n&cntrl\nioutfm=0,ntb=0,cut={cutoff}\n/", file=sander_file)
        
        if dump_prefix is not None:
            shutil.copyfile(os.path.join(temp_dir, "test.nml"), f"{dump_prefix}.nml")
            shutil.copyfile(os.path.join(temp_dir, "test.top"), f"{dump_prefix}.top")
        
        results = []
        for coordinates in coordinates_list:
            # Write coordinates and call sander.
            with open(os.path.join(temp_dir, "test.crd"), "w") as inpcrd_file:
                print(f"\n{coordinates.shape[0]:6}", file=inpcrd_file)
                coordinates_flat = coordinates.value_in_unit(openmm.unit.angstrom).flatten()
                for line_index in range((coordinates_flat.size + 5) // 6):
                    print("".join(f"{value:12.7f}" for value in coordinates_flat[6 * line_index:6 * (line_index + 1)]), file=inpcrd_file)
            sander_result = subprocess.run(["sander", "-O", "-i", "test.nml", "-o", "test.log", "-p", "test.top", "-c", "test.crd", "-inf", "test.ene"], cwd=temp_dir, capture_output=True)
            
            if sander_result.returncode:
                try:
                    with open(os.path.join(temp_dir, "test.log"), "r") as log_file:
                        print(log_file.read())
                except FileNotFoundError:
                    pass
                raise RuntimeError(f"sander exited with code {sander_result.returncode}")

            # Read energy file.
            with open(os.path.join(temp_dir, "test.ene"), "r") as energy_file:
                # Adapted from OpenFF Interchange.
                # See https://github.com/openforcefield/openff-interchange/.
                energy_data = dict()
                energy_file.readline()
                energy_file.readline()
                for energy_line in energy_file:
                    for field_range in ((1, 24), (26, 50), (52, 79)):
                        energy_term = energy_line[field_range[0]:field_range[1]]
                        if "=" in energy_term:
                            energy_type, energy_value = energy_term.split("=")
                            energy_data[energy_type.strip()] = float(energy_value.strip())

            result = {}
            for force_group in ForceGroup:
                result[force_group] = sum(energy_data.get(label, 0.0) for label in FORCE_GROUP_DATA[force_group][0]) * openmm.unit.kilocalorie_per_mole
            result[None] = energy_data["EPtot"] * openmm.unit.kilocalorie_per_mole

            # Do correction for Amber's value of the electric constant.
            electrostatic_energy = sum(energy_data.get(label, 0.0) for label in AMBER_ELECTROSTATIC_GROUPS) * openmm.unit.kilocalorie_per_mole
            electrostatic_shift = electrostatic_energy * (AMBER_ELECTROSTATIC_SCALE - 1.0)
            result[ForceGroup.NONBONDED]+= electrostatic_shift
            result[None] += electrostatic_shift

            results.append(result)
        
        return results

    def get_openmm_amber_energies(self, temp_dir, coordinates_list, dump_prefix):
        """
        Computes energies with OpenMM using OpenMM AmberPrmtopFile.

        Parameters
        ----------
        temp_dir : str
            Path to a temporary directory containing test files.
        coordinates_list : list
            A list of arrays of positions.  Energies will be evaluated for each
            set of coordinates.
        dump_prefix : str or None
            Prefix for names of input files to dump for debugging, or None to
            skip dumping.

        Returns
        -------
        list
            For each set of coordinates, a dictionary containing the energy for
            each ForceGroup (or None for the total potential energy).
        """

        prmtop_path = os.path.join(temp_dir, "test.top")
        if dump_prefix is not None:
            shutil.copyfile(prmtop_path, f"{dump_prefix}.top")
        prmtop = openmm.app.AmberPrmtopFile(prmtop_path)
        return self._evaluate_openmm(prmtop.createSystem(**SYSTEM_OPTIONS), coordinates_list, dump_prefix, prmtop.topology)
    
    def get_openmm_ffxml_energies(self, temp_dir, coordinates_list, dump_prefix, force_field, topology):
        """
        Computes energies with OpenMM using OpenMM AmberPrmtopFile.

        Parameters
        ----------
        temp_dir : str
            Path to a temporary directory containing test files.
        coordinates_list : list
            A list of arrays of positions.  Energies will be evaluated for each
            set of coordinates.
        dump_prefix : str or None
            Prefix for names of input files to dump for debugging, or None to
            skip dumping.
        force_field : openmm.app.ForceField
            The force field to use to create a system.
        topology : openmm.app.Topology
            The topology to use to create a system.

        Returns
        -------
        list
            For each set of coordinates, a dictionary containing the energy for
            each ForceGroup (or None for the total potential energy).
        """

        return self._evaluate_openmm(force_field.createSystem(topology, **SYSTEM_OPTIONS), coordinates_list, dump_prefix, topology)

    def _evaluate_openmm(self, system, coordinates_list, dump_prefix, topology):
        """
        Evaluates energies for an OpenMM system.  This will modify the system.

        Parameters
        ----------
        system : openmm.System
            The system for which to create a context and perform the evaluation.
        coordinates_list : list
            A list of arrays of positions.  Energies will be evaluated for each
            set of coordinates.
        dump_prefix : str or None
            Prefix for names of input files to dump for debugging, or None to
            skip dumping.
        topology : openmm.app.Topology
            Topology for PDB dumping.

        Returns
        -------
        list
            For each set of coordinates, a dictionary containing the energy for
            each ForceGroup (or None for the total potential energy).
        """

        # Set force groups of all forces in the System.
        for force_index, force in enumerate(system.getForces()):
            for force_group in ForceGroup:
                if isinstance(force, FORCE_GROUP_DATA[force_group][1]):
                    force.setForceGroup(force_group.value)
                    break
            else:
                raise RuntimeError(f"unknown force {type(force).__name__}")

        if dump_prefix is not None:
            with open(f"{dump_prefix}.xml", "w") as dump_file:
                dump_file.write(openmm.XmlSerializer.serialize(system))
            with open(f"{dump_prefix}.pdb", "w") as pdb_file:
                openmm.app.PDBFile.writeFile(topology, coordinates_list[0], pdb_file)

        context = self._make_context_from_system(system)

        results = []
        for coordinates in coordinates_list:
            context.setPositions(coordinates)

            result = {}

            # Evaluate energies for the entire system.
            state = context.getState(getEnergy=True)
            result[None] = state.getPotentialEnergy()

            # Evaluate energies for each force group.
            for force_group in ForceGroup:
                state = context.getState(getEnergy=True, groups={force_group.value})
                result[force_group] = state.getPotentialEnergy()

            results.append(result)

        return results

    def _compare_energies(self, energy_1, energy_2):
        """
        Helper function for comparing energies.

        Parameters
        ----------
        energy_1 : openmm.unit.Quantity
            The first energy to compare.
        energy_2 : openmm.unit.Quantity
            The second energy to compare.

        Returns
        -------
        openmm.unit.Quantity, float
            The magnitude of the difference between energies, and the relative
            difference.  Both values will be non-negative.  The relative
            difference will use whatever energy value has the larger magnitude
            as a reference scale, and will be zero if both energies are zero.
        """

        absolute_difference = abs(energy_2 - energy_1)
        if energy_1 or energy_2:
            relative_difference = absolute_difference / max(abs(energy_1), abs(energy_2))
        else:
            relative_difference = 0.0
        return absolute_difference, relative_difference

    def _format_energy_difference_data(self, label, difference_data):
        """
        Prints a table of energy difference data.

        Parameters
        ----------
        label : str
            A label for the relevant force group.
        difference_data : list(tuple(openmm.unit.Quantity, openmm.unit.Quantity, openmm.unit.Quantity, float))
            For each set of coordinates tested, energies evaluated in two ways,
            followed by the output of _compare_energies().

        Returns
        -------
        bool
            Whether or not any test failures were detected.
        """

        energies_1, energies_2, absolute_differences, relative_differences = zip(*difference_data)
        failures = [
            absolute_difference > self.absolute_tolerance and relative_difference > self.relative_tolerance
            for absolute_difference, relative_difference in zip(absolute_differences, relative_differences)
        ]

        any_failures = any(failures)
        formatted_line = f"{max(absolute_differences).value_in_unit(ENERGY_PRINT_UNIT):20.12e}  {max(relative_differences):12.8f}"
        if any_failures:
            formatted_line = f"{color_message(formatted_line)}  E_1 (kcal/mol)        E_2 (kcal/mol)"
        print(f"        {label:20}  {formatted_line}")
        if any_failures:
            for energy_1, energy_2, failure in zip(energies_1, energies_2, failures):
                formatted_line = f"{energy_1.value_in_unit(ENERGY_PRINT_UNIT):20.12e}  {energy_2.value_in_unit(ENERGY_PRINT_UNIT):20.12e}"
                print(f"{' ' * 66}{color_message(formatted_line) if failure else formatted_line}")

        return any_failures
    
    def _make_context_from_system(self, system):
        """
        Creates an OpenMM Context from a System.

        Parameters
        ----------
        system : openmm.System
            The System from which to create a Context.
        
        Returns
        -------
        openmm.Context
            A Context with a dummy Integrator and the requested Platform.
        """

        # Dummy integrator should not be utilized except to get a default
        # constraint tolerance.
        integrator = openmm.VerletIntegrator(0.001)

        if self.openmm_platform == "automatic":
            return openmm.Context(system, integrator)
        else:
            return openmm.Context(system, integrator, openmm.Platform.getPlatformByName(self.openmm_platform))

def random_in_unit_ball(generator):
    """
    Generates a random point in the closed unit ball.

    Parameters
    ----------
    generator : numpy.random.Generator
        Generator used to produce random values.

    Returns
    -------
    numpy.ndarray
        3-vector of floating-point numbers uniformly randomly drawn from the
        closed unit ball.
    """

    # Generate by rejection sampling starting from an axis-aligned cube of side
    # length 2 centered at the origin.
    while True:
        vector = generator.uniform(-1, 1, 3)
        if vector @ vector <= 1:
            return vector

def color_message(message, ansi="1;38;5;216"):
    """
    Generates ANSI escape sequences to highlight a message with color or other
    effects.

    Parameters
    ----------
    message : str
        The message to color.
    ansi : str
        The ANSI formatting codes to apply.

    Returns
    -------
    str
        A colored message, if the standard output stream is a TTY; otherwise,
        the message unmodified.
    """

    return f"\x1b[{ansi}m{message}\x1b[0m" if DO_COLOR else message

if __name__ == "__main__":
    main()

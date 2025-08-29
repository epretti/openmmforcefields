#!/usr/bin/env python

"""
A helper script to generate test cases with LEaP and OpenMM by setting up
initial coordinates and minimizing using an existing OpenMM force field.  The
purpose of this step is just to provide reasonable sets of coordinates for
testing, so the force field used in test case generation need not correspond
with the final force fields used for running the test case.
"""

import json
import numpy
import openmm.app
import os
import subprocess
import tempfile

DEBUG = False
OFFSET_DISTANCE = 20.0 # Å
PADDING_DISTANCE = 5.0 # Å
FFXML_DIRECTORY = "../../openmmforcefields/ffxml/amber"
SYSTEM_OPTIONS = dict(rigidWater=False, removeCMMotion=False, nonbondedMethod=openmm.app.NoCutoff)
HARMONIC_STRENGTH = 10.0 # kJ/mol/nm^2
LANGEVIN_TEMPERATURE = 300.0 * openmm.unit.kelvin
LANGEVIN_FRICTION = 1 / (100 * 0.002 * openmm.unit.picosecond)
LANGEVIN_ERROR = 0.0001
LANGEVIN_STEPS = 1000
OPENMM_PLATFORM = "Reference"
REPLICATE_COUNT = 5
MINIMIZE_TOLERANCE = 10.0 * openmm.unit.kilojoule_per_mole / openmm.unit.nanometer
RNG = numpy.random.default_rng(0x7ed8fa99f7932a36)

class FF:
    DNA_BSC0 = (["oldff/leaprc.DNA.bsc0"], ["DNA.bsc0.xml"])
    DNA_BSC1 = (["leaprc.DNA.bsc1"], ["DNA.bsc1.xml"])
    DNA_OL15 = (["leaprc.DNA.OL15"], ["DNA.OL15.xml"])
    DNA_OL21 = (["leaprc.DNA.OL21"], ["DNA.OL21.xml"])
    LIPID_17 = (["oldff/leaprc.lipid17"], ["lipid17.xml"])
    LIPID_21 = (["leaprc.lipid21"], ["lipid21.xml"])
    MULTI_AM1 = (["leaprc.ffAM1"], ["ffAM1.xml"])
    MULTI_FF03 = (["oldff/leaprc.ff03"], ["ff03.xml"])
    MULTI_FF10 = (["oldff/leaprc.ff10"], ["ff10.xml"])
    MULTI_FF14IPQ = (["oldff/leaprc.ff14ipq"], ["ff14ipq.xml"])
    MULTI_FF14SB = (["oldff/leaprc.ff14SB"], ["ff14SB.xml"])
    MULTI_FF14SBREDQ = (["oldff/leaprc.ff14SB.redq"], ["ff14SB.redq.xml"])
    MULTI_FF94 = (["oldff/leaprc.ff94"], ["ff94.xml"])
    MULTI_FF96 = (["oldff/leaprc.ff96"], ["ff96.xml"])
    MULTI_FF98 = (["oldff/leaprc.ff98"], ["ff98.xml"])
    MULTI_FF99 = (["oldff/leaprc.ff99"], ["ff99.xml"])
    MULTI_FF99BSC0 = (["oldff/leaprc.ff99bsc0"], ["ff99bsc0.xml"])
    MULTI_FF99SB = (["oldff/leaprc.ff99SB"], ["ff99SB.xml"])
    MULTI_FF99SBILDN = (["oldff/leaprc.ff99SBildn"], ["ff99SBildn.xml"])
    MULTI_FF99SBNMR = (["oldff/leaprc.ff99SBnmr"], ["ff99SBnmr.xml"])
    MULTI_PM3 = (["leaprc.ffPM3"], ["ffPM3.xml"])
    PROTEIN_FB15 = (["leaprc.protein.fb15"], ["protein.fb15.xml"])
    PROTEIN_FF03R1 = (["leaprc.protein.ff03.r1"], ["protein.ff03.r1.xml"])
    PROTEIN_FF14SB = (["leaprc.protein.ff14SB"], ["protein.ff14SB.xml"])
    PROTEIN_FF14SBONLYSC = (["leaprc.protein.ff14SBonlysc"], ["protein.ff14SBonlysc.xml"])
    PROTEIN_FF15IPQ = (["leaprc.protein.ff15ipq"], ["protein.ff15ipq.xml"])
    PROTEIN_FF15IPQVAC = (["leaprc.protein.ff15ipq-vac"], ["protein.ff15ipq-vac.xml"])
    PROTEIN_FF19IPQ = (["leaprc.protein.ff19ipq"], ["protein.ff19ipq.xml"])
    PROTEIN_FF19SB = (["leaprc.protein.ff19SB"], ["protein.ff19SB.xml"])
    PROTEIN_PHOSAA10 = (["oldff/leaprc.ff99SB", "leaprc.phosaa10"], ["ff99SB.xml", "phosaa10.xml"])
    PROTEIN_PHOSAA14SB = (["oldff/leaprc.ff14SB", "leaprc.phosaa14SB"], ["ff14SB.xml", "phosaa14SB.xml"])
    PROTEIN_PHOSAA19SB = (["leaprc.protein.ff19SB", "leaprc.phosaa19SB"], ["protein.ff19SB.xml", "phosaa19SB.xml"])
    PROTEIN_PHOSFB18 = (["leaprc.protein.fb15", "leaprc.phosfb18"], ["protein.fb15.xml", "phosfb18.xml"])
    RNA_OL3 = (["leaprc.RNA.OL3"], ["RNA.OL3.xml"])
    RNA_ROC = (["leaprc.RNA.ROC"], ["RNA.ROC.xml"])
    RNA_YIL = (["leaprc.RNA.YIL"], ["RNA.YIL.xml"])
    # TODO: leaprc.GLYCAM_06j-1

def main():
    def get_test_name(name):
        print(f"Preparing test case {name!r}")
        name = os.path.join("cases", name)
        os.makedirs(name, exist_ok=True)
        name = os.path.join(name, "test")
        return name

    def exclude(ff_list, *ff_exclude):
        ff_list = ff_list.copy()
        for ff_item in ff_exclude:
            ff_list.remove(ff_item)
        return ff_list

    # Protein:

    protein_ff = [
        FF.PROTEIN_FF14SB,

        FF.MULTI_AM1,
        FF.MULTI_FF03,
        FF.MULTI_FF10,
        FF.MULTI_FF14IPQ,
        FF.MULTI_FF14SB,
        FF.MULTI_FF14SBREDQ,
        FF.MULTI_FF94,
        FF.MULTI_FF96,
        FF.MULTI_FF98,
        FF.MULTI_FF99,
        FF.MULTI_FF99BSC0,
        FF.MULTI_FF99SB,
        FF.MULTI_FF99SBILDN,
        FF.MULTI_FF99SBNMR,
        FF.MULTI_PM3,
        FF.PROTEIN_FB15,
        FF.PROTEIN_FF03R1,
        FF.PROTEIN_FF14SBONLYSC,
        FF.PROTEIN_FF15IPQ,
        FF.PROTEIN_FF15IPQVAC,
        FF.PROTEIN_FF19IPQ,
        FF.PROTEIN_FF19SB,
    ]

    # Standard amino acids and capping groups.
    for index, residues in enumerate([
        ("ALA", "ARG", "ASN"),
        ("ASP", "GLN", "GLU"),
        ("GLY", "HID", "HIE"),
        ("HIP", "ILE", "LEU"),
        ("LYS", "MET", "PHE"),
        ("PRO", "SER", "THR"),
        ("TRP", "TYR", "VAL"),
    ]):
        create_case(
            [(f"N{residue}", residue, residue, residue, f"C{residue}") for residue in residues],
            ["test_case = sequence { $chains }"],
            protein_ff, get_test_name(f"protein.{index + 1:02}"),
        )
    create_case(
        [("ACE", "CYS", "CYS", "CYS", "NME"), ("ACE", "CYS", "CYS", "CYS", "NHE")],
        # ff03's ACE and NME are broken (they don't terminate chains correctly).
        ["set chain_1 tail null", "test_case = sequence { $chains }"],
        protein_ff, get_test_name("protein.08"),
    )
    # Disulfide bridges.
    create_case(
        [("NCYX", "ALA", "CCYX"), ("NALA", "ALA", "CYX", "ALA", "CALA"), ("NALA", "ALA", "CYX", "ALA", "CALA")],
        ["test_case = sequence { $chains }", "bond test_case.1.SG test_case.6.SG", "bond test_case.3.SG test_case.11.SG"],
        exclude(protein_ff, FF.PROTEIN_FF03R1), get_test_name("protein.09"),
    )
    # Alternate protonation states.
    create_case(
        [("NALA", "ALA", "ASH", "GLH", "LYN", "ALA", "CALA"), ("NALA", "ALA", "ASH", "GLH", "LYN", "ALA", "CALA")],
        ["test_case = sequence { $chains }"],
        protein_ff, get_test_name("protein.10"),
    )
    create_case(
        [("NALA", "CYM", "CYM", "CYM", "CALA"), ("NALA", "CYM", "CYM", "CYM", "CALA")],
        ["test_case = sequence { $chains }"],
        exclude(protein_ff, FF.MULTI_FF14IPQ), get_test_name("protein.11"),
    )
    # Special residues only supported by certain IPQ force fields.
    create_case(
        [("ACE", "CYX", "NME", "MTB"), ("NNLE", "NLE", "NLE", "NLE", "CNLE"), ("NASH", "ALA", "CASH"), ("NCYM", "ALA", "CCYM"), ("NGLH", "ALA", "CGLH"), ("NALA", "ALA", "CLYN")],
        ["test_case = sequence { $chains }", "bond test_case.2.SG test_case.4.SG"],
        # Force OpenMM to write all of the CONECT records to non-standard
        # residues: LEaP gives them standard names, so OpenMM misses the bonds
        # to extra protons.  Renaming them is a workaround.
        [FF.PROTEIN_FF15IPQ, FF.PROTEIN_FF15IPQVAC, FF.PROTEIN_FF19IPQ], get_test_name("protein.12"), {"ASP": "XXX"},
    )
    # Hydroxyproline.
    create_case(
        [("NALA", "HYP", "HYP", "HYP", "CALA"), ("NALA", "HYP", "HYP", "HYP", "CALA")],
        ["test_case = sequence { $chains }"],
        [FF.PROTEIN_FF14SB, FF.MULTI_FF10, FF.MULTI_FF14SB, FF.MULTI_FF14SBREDQ, FF.PROTEIN_FF14SBONLYSC, FF.PROTEIN_FF19SB], get_test_name("protein.13"),
    )
    create_case(
        [("NALA", "ALA", "CHYP"), ("NALA", "ALA", "CHYP")],
        ["test_case = sequence { $chains }"],
        [FF.PROTEIN_FF14SB, FF.MULTI_FF14SB, FF.MULTI_FF14SBREDQ, FF.PROTEIN_FF14SBONLYSC, FF.PROTEIN_FF19SB], get_test_name("protein.14"),
    )
    # Phosphorylated residues.
    create_case(
        [("NALA", "ALA", "PTR", "S1P", "SEP", "T1P", "TPO", "Y1P", "ALA", "CALA"), ("NALA", "ALA", "PTR", "S1P", "SEP", "T1P", "TPO", "Y1P", "ALA", "CALA")],
        ["test_case = sequence { $chains }"],
        [FF.PROTEIN_PHOSAA14SB, FF.MULTI_FF10, FF.PROTEIN_PHOSAA10, FF.PROTEIN_PHOSAA19SB, FF.PROTEIN_PHOSFB18], get_test_name("protein.15"),
    )
    create_case(
        [("NALA", "H2D", "H2D", "H2D", "CALA"), ("NALA", "H2D", "H2D", "H2D", "CALA")],
        ["test_case = sequence { $chains }"],
        [FF.PROTEIN_PHOSAA14SB, FF.MULTI_FF10, FF.PROTEIN_PHOSAA19SB], get_test_name("protein.16"),
    )
    create_case(
        [("NALA", "ALA", "H1D", "H1E", "H2E", "ALA", "CALA"), ("NALA", "ALA", "H1D", "H1E", "H2E", "ALA", "CALA")],
        ["test_case = sequence { $chains }"],
        [FF.PROTEIN_PHOSAA14SB, FF.PROTEIN_PHOSAA19SB], get_test_name("protein.17"),
    )

    # DNA:

    dna_ff = [
        FF.DNA_OL15,

        FF.DNA_BSC0,
        FF.DNA_BSC1,
        FF.DNA_OL21,
        FF.MULTI_AM1,
        FF.MULTI_FF03,
        FF.MULTI_FF10,
        FF.MULTI_FF14IPQ,
        FF.MULTI_FF14SB,
        FF.MULTI_FF14SBREDQ,
        FF.MULTI_FF94,
        FF.MULTI_FF96,
        FF.MULTI_FF98,
        FF.MULTI_FF99,
        FF.MULTI_FF99BSC0,
        FF.MULTI_FF99SB,
        FF.MULTI_FF99SBILDN,
        FF.MULTI_FF99SBNMR,
        FF.MULTI_PM3,
        FF.RNA_ROC,
    ]

    create_case(
        [("DA5", "DA", "DA", "DC", "DG", "DT", "DA", "DA3"), ("DA5", "DA", "DA3"), ("DC5", "DC", "DC3"), ("DG5", "DG", "DG3"), ("DT5", "DT", "DT3"), ("DAN",), ("DCN",), ("DGN",), ("DTN",)],
        ["test_case = sequence { $chains }"],
        dna_ff, get_test_name("nucleic.01"),
    )

    # RNA:

    rna_new_ff = [
        FF.RNA_OL3,

        FF.MULTI_FF10,
        FF.MULTI_FF14IPQ,
        FF.MULTI_FF14SB,
        FF.MULTI_FF14SBREDQ,
        FF.RNA_ROC,
        FF.RNA_YIL,
    ]
    rna_old_ff = [
        FF.MULTI_FF99SB,

        FF.MULTI_AM1,
        FF.MULTI_FF03,
        FF.MULTI_FF94,
        FF.MULTI_FF96,
        FF.MULTI_FF98,
        FF.MULTI_FF99,
        FF.MULTI_FF99BSC0,
        FF.MULTI_FF99SBILDN,
        FF.MULTI_FF99SBNMR,
        FF.MULTI_PM3,
    ]

    create_case(
        [("A5", "A", "A", "C", "G", "U", "A", "A3"), ("A5", "A", "A3"), ("C5", "C", "C3"), ("G5", "G", "G3"), ("U5", "U", "U3")],
        ["test_case = sequence { $chains }"],
        rna_new_ff, get_test_name("nucleic.02"),
    )
    create_case(
        [("AN",), ("CN",), ("GN",), ("UN",), ("OHE", "A3"), ("OHE", "C3"), ("OHE", "G3"), ("OHE", "U3")],
        ["test_case = sequence { $chains }"],
        rna_new_ff, get_test_name("nucleic.03"),
    )

    # Lipids:

    # Cholesterol.
    create_case(
        [("CHL",), ("CHL",)],
        ["test_case = sequence { $chains }"],
        [FF.LIPID_17, FF.LIPID_21], get_test_name("lipid.01"),
    )
    # Standard lipids supported by Lipid17.
    create_case(
        [("LAL", "PC", "MY"), ("PA", "PE", "ST")],
        ["set chain_1 tail null", "test_case = sequence { $chains }"],
        [FF.LIPID_17, FF.LIPID_21], get_test_name("lipid.02"),
    )
    create_case(
        [("OL", "PH-", "AR"), ("DHA", "PS", "LAL")],
        ["set chain_1 tail null", "test_case = sequence { $chains }"],
        [FF.LIPID_17, FF.LIPID_21], get_test_name("lipid.03"),
    )
    create_case(
        [("MY", "PGR", "PA"), ("ST", "PGR", "OL")],
        ["set chain_1 tail null", "test_case = sequence { $chains }"],
        [FF.LIPID_17, FF.LIPID_21], get_test_name("lipid.04"),
    )
    # Sphingomyelin is only supported by Lipid21.  Note that only SA can be
    # attached as the 2nd substituent to SPM; there aren't other parameters.
    create_case(
        [("AR", "SPM", "SA"), ("DHA", "SPM", "SA")],
        ["set chain_1 tail null", "test_case = sequence { $chains }"],
        [FF.LIPID_21], get_test_name("lipid.05"),
    )
    # Enantiomers of phosphatidylglycerol.
    create_case(
        [("LAL", "PGR", "MY"), ("PA", "PGS", "ST")],
        ["set chain_1 tail null", "test_case = sequence { $chains }"],
        [FF.LIPID_21], get_test_name("lipid.06"),
    )

def create_case(chains, leap_commands, force_field_list, output_prefix, residue_name_replacements=None):
    """
    Creates a test case.  The LEaP commands should create an object test_case
    containing the topology for the test case from chains chain_1, chain_2, etc.
    """

    main_force_field = force_field_list[0]
    main_leaprc_list, main_ffxml_list = main_force_field

    with tempfile.TemporaryDirectory() as temp_dir:
        # Generate initial coordinates with LEaP.
        create_leap_input(temp_dir, main_leaprc_list, chains, leap_commands)
        try:
            run_leap(temp_dir)
        except Exception as exception:
            raise RuntimeError(main_force_field) from exception

        # Copy the LEaP input file for this first test case.
        with open(os.path.join(temp_dir, "test.leap"), "rb") as leap_in:
            with open(f"{output_prefix}.0.leap", "wb") as leap_out:
                leap_out.write(leap_in.read())

        # Load the prmtop generated by LEaP.
        try:
            prmtop = openmm.app.AmberPrmtopFile(os.path.join(temp_dir, "test.top"))
        except Exception as exception:
            raise RuntimeError(main_force_field) from exception

        # Load the PDB generated by LEaP.
        pdb = openmm.app.PDBFile(os.path.join(temp_dir, "test.pdb"))
        positions = numpy.array([position.value_in_unit(openmm.unit.angstrom) for position in pdb.positions])

    # We may have to patch up the PDB topology based on the prmtop if
    # non-standard residues were involved.  Find the bonds in the PDB and the
    # prmtop.
    pdb_atoms = list(pdb.topology.atoms())
    pdb_bonds = set()
    pdb_bond_indices = dict()
    for index, bond in enumerate(pdb.topology._bonds):
        index_1 = bond.atom1.index
        index_2 = bond.atom2.index
        pdb_bond_key = (min(index_1, index_2), max(index_1, index_2))
        pdb_bonds.add(pdb_bond_key)
        pdb_bond_indices[pdb_bond_key] = index
    prmtop_bonds = set()
    for bond in prmtop.topology.bonds():
        index_1 = bond.atom1.index
        index_2 = bond.atom2.index
        prmtop_bonds.add((min(index_1, index_2), max(index_1, index_2)))

    # Add any bonds in the prmtop that aren't in the PDB.
    for index_1, index_2 in sorted(prmtop_bonds - pdb_bonds):
        if DEBUG:
            print(f"Adding missing bond {pdb_atoms[index_1]} - {pdb_atoms[index_2]}")
        pdb.topology.addBond(pdb_atoms[index_1], pdb_atoms[index_2])

    # Complain about any bonds in the PDB that aren't in the prmtop.
    pdb_bonds_remove = []
    for index_1, index_2 in sorted(pdb_bonds - prmtop_bonds):
        if DEBUG:
            print(f"Found extra bond {pdb_atoms[index_1]} - {pdb_atoms[index_2]}")
        atom_1 = pdb_atoms[index_1]
        atom_2 = pdb_atoms[index_2]
        if atom_1.name == atom_2.name == "SG" and atom_1.residue.name == atom_2.residue.name == "CYS":
            # OpenMM might have incorrectly added a disulfide, so take it out.
            pdb_bonds_remove.append(pdb_bond_indices[index_1, index_2])
        else:
            raise ValueError(f"Invalid bond {pdb_atoms[index_1]} - {pdb_atoms[index_2]} in PDB file")
    for index in reversed(sorted(pdb_bonds_remove)):
        del pdb.topology._bonds[index]

    # Replace residue names (sometimes necessary to prevent OpenMM from getting
    # confused about standard bonds while writing residues).
    if residue_name_replacements is not None:
        for residue in pdb.topology.residues():
            replacement = residue_name_replacements.get(residue.name)
            if replacement is not None:
                residue.name = replacement

    # Space out chains from each other to prevent overlapping atoms.
    x_ranges = []
    for chain in pdb.topology.chains():
        x_positions = [positions[atom.index][0] for atom in chain.atoms()]
        x_ranges.append((min(x_positions), max(x_positions)))
    xyz_offsets = []
    for chain_index in range(len(x_ranges)):
        y_offset = RNG.uniform(-OFFSET_DISTANCE, OFFSET_DISTANCE)
        z_offset = RNG.uniform(-OFFSET_DISTANCE, OFFSET_DISTANCE)
        if chain_index:
            xyz_offsets.append((PADDING_DISTANCE + xyz_offsets[chain_index - 1][0] + x_ranges[chain_index - 1][1] - x_ranges[chain_index][0], y_offset, z_offset))
        else:
            xyz_offsets.append((0.0, y_offset, z_offset))
    for chain, xyz_offset in zip(pdb.topology.chains(), xyz_offsets):
        for atom in chain.atoms():
            positions[atom.index] += xyz_offset

    # Create an OpenMM system.
    try:
        openmm_force_field = openmm.app.ForceField(*(os.path.join(FFXML_DIRECTORY, ffxml) for ffxml in main_ffxml_list))
        system = openmm_force_field.createSystem(pdb.topology, **SYSTEM_OPTIONS)
    except Exception as exception:
        raise RuntimeError(main_force_field) from exception

    # Add a harmonic force to pull molecules together.
    harmonic_force = openmm.CustomExternalForce(f"{HARMONIC_STRENGTH}*(x^2+y^2+z^2)")
    for particle_index in range(system.getNumParticles()):
        harmonic_force.addParticle(particle_index)
    system.addForce(harmonic_force)

    # Delete bare charges in some Amber models that cause serious instabilities.
    # This is just for generating configurations automatically and will not
    # affect testing with the real Amber parameters later elsewhere.
    nonbonded_force, = (force for force in system.getForces() if isinstance(force, openmm.NonbondedForce))
    for index in range(system.getNumParticles()):
        charge, sigma, epsilon = nonbonded_force.getParticleParameters(index)
        if charge and not epsilon:
            nonbonded_force.setParticleParameters(index, 0.0, sigma, epsilon)
    for index in range(nonbonded_force.getNumExceptions()):
        index_1, index_2, charge, sigma, epsilon = nonbonded_force.getExceptionParameters(index)
        if charge and not epsilon:
            nonbonded_force.setExceptionParameters(index, index_1, index_2, 0.0, sigma, epsilon)

    # Generate an initial configuration by minimizing and integrating.
    integrator = openmm.VariableLangevinIntegrator(LANGEVIN_TEMPERATURE, LANGEVIN_FRICTION, LANGEVIN_ERROR)
    integrator.setRandomNumberSeed(RNG.integers(1 << 31))
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName(OPENMM_PLATFORM))
    context.setPositions(positions * openmm.unit.angstrom)
    print("    Generating initial coordinates...")
    openmm.LocalEnergyMinimizer.minimize(context, MINIMIZE_TOLERANCE)
    context.setVelocitiesToTemperature(LANGEVIN_TEMPERATURE, RNG.integers(1 << 31))
    integrator.step(LANGEVIN_STEPS)
    coordinates_list = [context.getState(positions=True).getPositions(asNumpy=True)]

    # Generate sets of perturbed positions if desired.
    while len(coordinates_list) < REPLICATE_COUNT:
        print("    Generating perturbed coordinates...")
        integrator.step(LANGEVIN_STEPS)
        coordinates_list.append(context.getState(positions=True).getPositions(asNumpy=True))

    # Write positions to a PDB.
    with open(f"{output_prefix}.pdb", "w") as pdb_out:
        openmm.app.PDBFile.writeHeader(pdb.topology, pdb_out)
        for index, coordinates in enumerate(coordinates_list):
            openmm.app.PDBFile.writeModel(pdb.topology, coordinates, pdb_out, index + 1)
        openmm.app.PDBFile.writeFooter(pdb.topology, pdb_out)

    # Reload the PDB for testing.
    test_pdb = openmm.app.PDBFile(f"{output_prefix}.pdb")

    # Process the other force fields.
    for index, force_field in enumerate(force_field_list):
        if not index:
            continue

        leaprc_list, ffxml_list = force_field

        with tempfile.TemporaryDirectory() as temp_dir:
            create_leap_input(temp_dir, leaprc_list, chains, leap_commands)
            try:
                run_leap(temp_dir)
            except Exception as exception:
                raise RuntimeError(force_field) from exception

            # Copy the LEaP input file for this test case.
            with open(os.path.join(temp_dir, "test.leap"), "rb") as leap_in:
                with open(f"{output_prefix}.{index}.leap", "wb") as leap_out:
                    leap_out.write(leap_in.read())

            # Make sure we can load the prmtop generated by LEaP.
            try:
                prmtop = openmm.app.AmberPrmtopFile(os.path.join(temp_dir, "test.top"))
            except Exception as exception:
                raise RuntimeError(force_field) from exception

            # Make sure we can load the FFXML files for this test case.
            try:
                openmm.app.ForceField(*(os.path.join(FFXML_DIRECTORY, ffxml) for ffxml in ffxml_list)).createSystem(test_pdb.topology, **SYSTEM_OPTIONS)
            except Exception as exception:
                raise RuntimeError(force_field) from exception

    # Write a manifest of the force fields for this test case.
    with open(f"{output_prefix}.json", "w") as json_file:
        json.dump(force_field_list, json_file, indent=4)

def create_leap_input(temp_dir, leaprc_list, chains, leap_commands):
    """
    Creates a LEaP input file `test.leap` in a temporary directory.
    """

    with open(os.path.join(temp_dir, "test.leap"), "w") as leap_file:
        for leaprc in leaprc_list:
            print(f"source {leaprc}", file=leap_file)
        for index, chain in enumerate(chains):
            chain_string = " ".join(chain)
            print(f"chain_{index + 1} = sequence {{ {chain_string} }}", file=leap_file)
        for leap_command in leap_commands:
            print(leap_command.replace("$chains", " ".join(f"chain_{index + 1}" for index in range(len(chains)))), file=leap_file)
        # Work around a LEaP bug that sometimes skips CMAPs (dumping and
        # reloading the structure we've built seems to normalize things).
        print("saveOff test_case test.off", file=leap_file)
        print("loadOff test.off", file=leap_file)
        print("saveAmberParm test_case test.top test.crd", file=leap_file)
        print("savePdb test_case test.pdb", file=leap_file)
        print("quit", file=leap_file)

def run_leap(temp_dir):
    """
    Calls LEaP on a LEaP input file created by `create_leap_input()`.
    """

    print("    Running LEaP...")
    result = subprocess.run(["tleap", "-f", "test.leap"], cwd=temp_dir, capture_output=True)
    stdout = result.stdout.decode()
    stderr = result.stderr.decode()
    failed = result.returncode or "Exiting LEaP: Errors = 0" not in stdout

    if failed or DEBUG:
        print("=" * 80, "LEaP output".center(80), "=" * 80, sep="\n")
        print(stdout)
        print("=" * 80)
    if failed:
        raise RuntimeError(f"LEaP exited with code {result.returncode}")

def random_in_unit_ball(rng):
    """
    Generates a random point in the closed unit ball.
    """

    # Generate by rejection sampling starting from an axis-aligned cube of side
    # length 2 centered at the origin.
    while True:
        vector = rng.uniform(-1, 1, 3)
        if vector @ vector <= 1:
            return vector

if __name__ == "__main__":
    main()

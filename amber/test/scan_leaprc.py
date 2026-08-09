#!/usr/bin/env python

"""
A helper script that finds all residues present in all Amber force fields in a
given AmberTools installation.  Useful for analyzing what residues need to be
tested with what sets of force fields, and what test cases to prepare.  Manual
interpretation of the output in order to prepare appropriate test cases is
necessary, so this script is provided as a standalone tool rather than
integrated into an automatic test case generator, but it may still be useful for
extending the test suite to new Amber force fields.
"""

import argparse
import functools
import glob
import json
import os
import sys
import warnings

def main():
    parser = argparse.ArgumentParser(description="Searches for residues in Amber force fields in an AmberTools installation (set $AMBERHOME or make it your current working directory)")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--list-leaprc-files", action="store_true", help="Lists all leaprc files to standard output")
    group.add_argument("--categorize-residues", action="store_true", help="Reads a list of leaprc files from standard input and categorizes residues as to their membership in these files to standard output")
    arguments = parser.parse_args()

    if arguments.list_leaprc_files:
        for leaprc_file in find_all_leaprc_files():
            print(leaprc_file)
    elif arguments.categorize_residues:
        print(json.dumps(categorize_residues(line for line in map(str.strip, sys.stdin) if line), indent=4))

def categorize_residues(leaprc_files):
    """
    Categorizes residues based on which LEaP input files they appear in.
    """

    # Construct a unique ordered list of LEaP input files provided.
    leaprc_files_ordered = {}
    for leaprc_file in leaprc_files:
        leaprc_files_ordered[leaprc_file] = len(leaprc_files_ordered)

    # Construct a unique ordered list of residues with the LEaP input files in
    # which they can be found.
    residues = {}
    for leaprc_file in leaprc_files_ordered:
        for residue in scan_leaprc_file(leaprc_file):
            residues.setdefault(residue, (len(residues), set()))[1].add(leaprc_file)

    # Make categories of residues based on the force fields supporting them.
    categories = {}
    for residue, (_, residue_leaprc_files) in residues.items():
        categories.setdefault(frozenset(residue_leaprc_files), set()).add(residue)

    # Sort and return results.
    return sorted((sorted(category, key=leaprc_files_ordered.get), sorted(category_residues, key=lambda residue: residues[residue][0])) for category, category_residues in categories.items())

@functools.cache
def get_leap_root():
    """
    Returns the root directory where LEaP stores data by reading the `AMBERHOME`
    environment variable.  If not set, uses the current working directory.
    """

    return os.path.realpath(os.path.join(os.environ.get("AMBERHOME", "."), "dat", "leap"))

def find_all_leaprc_files():
    """
    Returns all of the LEaP input files installed with AmberTools.
    """

    return sorted(glob.glob("**/leaprc.*", root_dir=os.path.join(get_leap_root(), "cmd"), recursive=True))

def find_leap_file(resource, start):
    """
    Finds an Amber library or Prep file by name `resource` starting in the
    directory `start` under the LEaP root directory, and descending into
    subdirectories if the file cannot be found.
    """

    start_path = os.path.join(get_leap_root(), start)
    try_path = os.path.join(start_path, resource)
    if os.path.exists(try_path):
        return os.path.realpath(try_path)

    for child_path in sorted(os.listdir(start_path)):
        if os.path.isdir(os.path.join(start_path, child_path)):
            try:
                return find_leap_file(resource, os.path.join(start, child_path))
            except FileNotFoundError:
                continue

    raise FileNotFoundError(resource)

def scan_leaprc_file(leaprc_path):
    """
    Returns all of the names of the entries from an Amber LEaP input file.
    """

    results = []
    results_set = set()

    with open(os.path.join(get_leap_root(), "cmd", leaprc_path), "r") as leaprc_file:
        for line in leaprc_file:
            fields = line.strip().split("#", maxsplit=1)[0].strip().split(maxsplit=1)
            if not fields:
                continue

            command = fields[0].lower()
            try:
                if command == "loadoff":
                    generator = scan_off_file(find_leap_file(fields[1], "lib"))
                elif command == "loadamberprep":
                    generator = scan_prep_file(find_leap_file(fields[1], "prep"))
                else:
                    continue
            except FileNotFoundError as error:
                warnings.warn(f"LEaP file {leaprc_path!r} references missing file {fields[1]}")
                continue

            for result in generator:
                if result in results_set:
                    warnings.warn(f"LEaP file {leaprc_path!r} includes multiple entries for {result!r}")
                else:
                    results.append(result)
                    results_set.add(result)

    return results

def scan_off_file(off_path):
    """
    Returns all of the names of the entries in an Amber library file.
    """

    with open(off_path, "r") as off_file:
        looking_for_index = True
        for line in off_file:
            if looking_for_index:
                if line.strip().lower().startswith("!!index"):
                    looking_for_index = False
            else:
                line = line.strip()
                if line.startswith("!"):
                    break
                if not (line.startswith("\"") and line.endswith("\"")):
                    raise ValueError(line)
                yield line[1:-1]

def scan_prep_file(prep_path):
    """
    Returns all of the names of the entries in an Amber Prep file.
    """

    with open(prep_path, "r") as prep_file:
        lines = iter(prep_file)
        next(lines)
        next(lines)

        while True:
            fields = next(lines).strip().split()
            if fields and fields[0].lower() == "stop":
                break
            next(lines)
            yield next(lines).strip().split()[0]
            while True:
                fields = next(lines).strip().split()
                if not fields:
                    continue
                if fields[0].lower() == "done":
                    break

if __name__ == "__main__":
    main()

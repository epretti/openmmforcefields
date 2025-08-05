COMBINED_FORCE_FIELDS = [
    (["oldff/leaprc.ff94"], ["ff94.xml"]),
    (["oldff/leaprc.ff96"], ["ff96.xml"]),
    (["oldff/leaprc.ff98"], ["ff98.xml"]),
    (["oldff/leaprc.ff99"], ["ff99.xml"]),
    (["oldff/leaprc.ff99SB"], ["ff99SB.xml"]),
    (["oldff/leaprc.ff99SBildn"], ["ff99SBildn.xml"]),
    (["oldff/leaprc.ff99SBnmr"], ["ff99SBnmr.xml"]),
    (["oldff/leaprc.ff99bsc0"], ["ff99bsc0.xml"]),
    (["oldff/leaprc.ff03"], ["ff03.xml"]),
    (["oldff/leaprc.ff10"], ["ff10.xml"]),
    (["oldff/leaprc.ff14SB.redq"], ["ff14SB.redq.xml"]),
    (["oldff/leaprc.ff14SB"], ["ff14SB.xml"]),
    (["oldff/leaprc.ff14ipq"], ["ff14ipq.xml"]),
    (["leaprc.ffAM1"], ["ffAM1.xml"]),
    (["leaprc.ffPM3"], ["ffPM3.xml"]),
]

DNA_FORCE_FIELDS = [
    (["oldff/leaprc.DNA.bsc0"], ["DNA.bsc0.xml"]),
    (["leaprc.DNA.bsc1"], ["DNA.bsc1.xml"]),
    (["leaprc.DNA.OL15"], ["DNA.OL15.xml"]),
    (["leaprc.DNA.OL21"], ["DNA.OL21.xml"]),
    (["leaprc.RNA.ROC"], ["RNA.ROC.xml"]),
]

DNA_NAMES = ["DA", "DC", "DG", "DT"]

RNA_FORCE_FIELDS = [
    (["leaprc.RNA.OL3"], ["RNA.OL3.xml"]),
    (["leaprc.RNA.YIL"], ["RNA.YIL.xml"]),
    (["leaprc.RNA.ROC"], ["RNA.ROC.xml"]),
]

RNA_NAMES = ["A", "C", "G", "U"]

RNA_NEW_NAMES_FORCE_FIELDS = [
    "oldff/leaprc.ff10",
    "oldff/leaprc.ff14SB.redq",
    "oldff/leaprc.ff14SB",
    "oldff/leaprc.ff14ipq",
    "leaprc.RNA.OL3",
    "leaprc.RNA.YIL",
    "leaprc.RNA.ROC",
]

PROTEIN_FORCE_FIELDS = [
    (["leaprc.protein.ff03.r1"], ["protein.ff03.r1.xml"]),
    (["leaprc.protein.ff03ua"], ["protein.ff03ua.xml"]),
    (["leaprc.protein.ff14SB"], ["protein.ff14SB.xml"]),
    (["leaprc.protein.ff14SBonlysc"], ["protein.ff14SBonlysc.xml"]),
    (["leaprc.protein.fb15"], ["protein.fb15.xml"]),
    (["leaprc.protein.ff15ipq"], ["protein.ff15ipq.xml"]),
    (["leaprc.protein.ff15ipq-vac"], ["protein.ff15ipq-vac.xml"]),
    (["leaprc.protein.ff19SB"], ["protein.ff19SB.xml"]),
    (["leaprc.protein.ff19ipq"], ["protein.ff19ipq.xml"]),
]

PROTEIN_NAMES = [
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HID",
    "HIE",
    "HIP",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
]

PHOSPHORYLATED_FORCE_FIELDS = [
    (["oldff/leaprc.ff99SB", "leaprc.phosaa10"], ["ff99SB.xml", "phosaa10.xml"], []),
    (["leaprc.protein.fb15", "leaprc.phosfb18"], ["protein.fb15.xml", "phosfb18.xml"], []),
    (["oldff/leaprc.ff14SB", "leaprc.phosaa14SB"], ["ff14SB.xml", "phosaa14SB.xml"], ["H1D", "H1E", "H2D", "H2E"]),
    (["leaprc.protein.ff19SB", "leaprc.phosaa19SB"], ["protein.ff19SB.xml", "phosaa19SB.xml"], ["H1D", "H1E", "H2D", "H2E"]),
    (["oldff/leaprc.ff10"], ["ff10.xml"], ["H2D"]),
]

PHOSPHORYLATED_NAMES = [
    "PTR",
    "S1P",
    "SEP",
    "T1P",
    "TPO",
    "Y1P",
]

def main():
    make_nucleic_test()
    make_protein_test()

    # Issues with lipids:
    # SA missing in place of PA (sphingosine has same topology as palmitoyl?  Or maybe not.  Or lipid17 doesn't have this.)
    # PGR and PGS are different, but we only end up including one of them
    # SPM 

    # TODO:
    #   proteins: linked cysteines
    #   glycans (note OME)
    #   HYP, NLN, OLP, OLS, and OLT with glycan linkage
    #   the zoo of water and ions
    #   lipids: modular
    #   lipids: composite (requires understanding CHARMMLIPID2AMBER)
    
    # TODO: phosphorylated
    # TODO: water and ions

def make_nucleic_test():
    with open("nucleic.toml", "w") as file:
        for leaprcs, ffxmls in COMBINED_FORCE_FIELDS + DNA_FORCE_FIELDS:
            write_test(file, f"{ffxmls[0]} DNA tetramers", [[f"{name}5", name, name, f"{name}3"] for name in DNA_NAMES], leaprcs, ffxmls)
            write_test(file, f"{ffxmls[0]} DNA monomers", [[f"{name}N"] for name in DNA_NAMES], leaprcs, ffxmls)
        for leaprcs, ffxmls in COMBINED_FORCE_FIELDS + RNA_FORCE_FIELDS:
            new_rna = leaprcs[0] in RNA_NEW_NAMES_FORCE_FIELDS
            rna_prefix = "" if new_rna else "R"
            write_test(file, f"{ffxmls[0]} RNA tetramers", [[f"{rna_prefix}{name}5", f"{rna_prefix}{name}", f"{rna_prefix}{name}", f"{rna_prefix}{name}3"] for name in RNA_NAMES], leaprcs, ffxmls)
            write_test(file, f"{ffxmls[0]} RNA nucleosides", [[f"{rna_prefix}{name}N"] for name in RNA_NAMES], leaprcs, ffxmls)
            if new_rna:
                # "Old" RNA force fields don't have the OHE to make nucleotides.
                write_test(file, f"{ffxmls[0]} RNA nucleotides", [["OHE", f"{rna_prefix}{name}3"] for name in RNA_NAMES], leaprcs, ffxmls)

def make_protein_test():
    with open("protein.toml", "w") as file:
        for leaprcs, ffxmls in COMBINED_FORCE_FIELDS + PROTEIN_FORCE_FIELDS:
            # Standard residues.
            protein_names = PROTEIN_NAMES
            if ffxmls[0] in ["protein.ff15ipq.xml", "protein.ff15ipq-vac.xml", "protein.ff19ipq.xml"]:
                # >=15 ipq force fields support norleucine.
                protein_names = protein_names + ["NLE"]
            for index in range(0, len(protein_names), 4):
                write_test(file, f"{ffxmls[0]} tetramers {index + 1}-{min(index + 4, len(protein_names))}", [[f"N{name}", name, name, f"C{name}"] for name in protein_names[index:index + 4]], leaprcs, ffxmls)
            
            # Capping groups.
            chains = [["ACE", "CALA"], ["NALA", "NME"]]
            if ffxmls[0] != "protein.ff03ua.xml":
                chains.append(["NALA", "NHE"])
            write_test(file, f"{ffxmls[0]} caps", chains, leaprcs, ffxmls)

            # Other protonation states.
            if ffxmls[0] != "protein.ff03ua.xml":
                chains = []
                for name in ["ASH", "CYM", "GLH", "LYN"]:
                    # CYM is missing from ff14ipq, terminal versions are only
                    # supported by >=15 ipq force fields, and NLYN isnt'
                    # supported by any Amber force field.
                    if name == "CYM" and ffxmls[0] == "ff14ipq.xml":
                        continue
                    use_c = ffxmls[0] in ["protein.ff15ipq.xml", "protein.ff15ipq-vac.xml", "protein.ff19ipq.xml"]
                    use_n = use_c and name != "LYN"
                    chain = []
                    chain.append(f"N{name}" if use_n else "NGLY")
                    chain.append(name)
                    chain.append(f"C{name}" if use_c else "CGLY")
                    chains.append(chain)
                write_test(file, f"{ffxmls[0]} protonation", chains, leaprcs, ffxmls)
        
        for leaprcs, ffxmls, extra_names in PHOSPHORYLATED_FORCE_FIELDS:
            write_test(file, f"{ffxmls[0]} phosphorylated", [["NGLY", name, "CGLY"] for name in PHOSPHORYLATED_NAMES + extra_names], leaprcs, ffxmls)

def write_test(file, name, chains, leaprcs, ffxmls):
    print("[[test]]", file=file)
    print(f"name = {name!r}", file=file)
    print(f"chains = {chains!r}", file=file)
    print(f"leaprcs = {leaprcs!r}", file=file)
    print(f"ffxmls = {ffxmls!r}", file=file)
    print(file=file)

if __name__ == "__main__":
    main()

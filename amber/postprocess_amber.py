import io
import os
import lxml.etree as etree

FFXML_DIR = "../openmmforcefields/ffxml/amber"

def main():
    """
    Apply miscellaneous fixes to generated FFXML files.  It is better to do this
    separately than to try to add yet another layer of complexity to the script
    `convert_amber.py`.  Eventually if this infrastructure can be replaced, and
    such fixes are still required, the approach used in the CHARMM conversion
    for patching can be employed.
    """

    # Fix ACE and NME in ff03 and ff03r1:
    for force_field in ("ff03", "protein.ff03.r1"):
        path = os.path.join(FFXML_DIR, f"{force_field}.xml")
        tree = etree.parse(path)
        find_and_remove_xml_elements(tree.getroot(), tree.findall("./Residues/Residue[@name='ACE']/ExternalBond[@atomName='HH31']"))
        find_and_remove_xml_elements(tree.getroot(), tree.findall("./Residues/Residue[@name='NME']/ExternalBond[@atomName='CH3']"))
        tree.write(path)

def find_and_remove_xml_elements(root, target_elements):
    root[:] = [child for child in root if child not in target_elements]
    for child in root:
        find_and_remove_xml_elements(child, target_elements)

if __name__ == "__main__":
    main()

# default bonds that can be cut in either directions (hence the set definition)
# key=resname value=list(set(atomname1, atomname2), ...)
cutable = {
    "ALA": [{"CA", "CB"}],
    "ARG": [{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD"}],
    "ASN": [{"CA", "CB"}],
    "ASP": [{"CA", "CB"}],
    "CYS": [{"CA", "CB"}],
    "CYM": [{"CA", "CB"}],
    "GLN": [{"CA", "CB"}, {"CB", "CG"}],
    "GLU": [{"CA", "CB"}, {"CB", "CG"}],
    "GLY": [],
    "HIS": [{"CA", "CB"}, {"CB", "CG"}],
    "HSD": [{"CA", "CB"}, {"CB", "CG"}],
    "HSE": [{"CA", "CB"}, {"CB", "CG"}],
    "HID": [{"CA", "CB"}, {"CB", "CG"}],
    "HIE": [{"CA", "CB"}, {"CB", "CG"}],
    "HSP": [{"CA", "CB"}, {"CB", "CG"}],
    "ILE": [{"CA", "CB"}, {"CB", "CG1"}, {"CB", "CG2"}, {"CG1", "CD"}],
    "LEU": [{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD1"}, {"CG", "CD2"}],
    "LYS": [{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD"}, {"CD", "CE"}],
    "MET": [{"CA", "CB"}, {"CB", "CG"}],
    "PHE": [{"CA", "CB"}, {"CB", "CG"}],
    "PRO": [],
    "SER": [{"CA", "CB"}],
    "THR": [{"CA", "CB"}, {"CB", "CG2"}],
    "TRP": [{"CA", "CB"}, {"CB", "CG"}],
    "TYR": [{"CA", "CB"}, {"CB", "CG"}],
    "VAL": [{"CA", "CB"}, {"CB", "CG1"}, {"CB", "CG2"}],
    "ADP": [{"C4'", "C5'"}],
    "ATP": [{"C4'", "C5'"}],
    "GDP": [{"C4'", "C5'"}],
    "GTP": [{"C4'", "C5'"}],
    "NAP": [{"C1", "C2"}, {"C11", "C12"}],
    "BGLC": [{"C5", "C6"}]  # oligosacharide residue
}
# default bonds that can be cut in only one direction, QM side goes first
# key=resname value=list(list(QMatom, MMatom), ...)
directional = {
    "ALA": [["CA", "C"], ["CA", "N"]],
    "ARG": [["CA", "C"], ["CA", "N"], ["CD", "NE"]],
    "ASN": [["CA", "C"], ["CA", "N"], ["CB", "CG"]],
    "ASP": [["CA", "C"], ["CA", "N"], ["CB", "CG"]],
    "CYS": [["CA", "C"], ["CA", "N"], ["CB", "CG"]],
    "CYM": [["CA", "C"], ["CA", "N"], ["CB", "CG"]],
    "GLN": [["CA", "C"], ["CA", "N"], ["CG", "CD"]],
    "GLU": [["CA", "C"], ["CA", "N"], ["CG", "CD"]],
    "GLY": [["CA", "C"], ["CA", "N"]],
    "HIS": [["CA", "C"], ["CA", "N"]],
    "HSD": [["CA", "C"], ["CA", "N"]],
    "HID": [["CA", "C"], ["CA", "N"]],
    "HSE": [["CA", "C"], ["CA", "N"]],
    "HIE": [["CA", "C"], ["CA", "N"]],
    "HSP": [["CA", "C"], ["CA", "N"]],
    "ILE": [["CA", "C"], ["CA", "N"]],
    "LEU": [["CA", "C"], ["CA", "N"]],
    "LYS": [["CA", "C"], ["CA", "N"], ["CE", "NZ"]],
    "MET": [["CA", "C"], ["CA", "N"], ["CG", "SD"], ["CE", "SD"]],
    "PHE": [["CA", "C"], ["CA", "N"]],
    "PRO": [["CA", "C"]],
    "SER": [["CA", "C"], ["CA", "N"], ["CB", "OG"]],
    "THR": [["CA", "C"], ["CA", "N"], ["CB", "OG1"]],
    "TRP": [["CA", "C"], ["CA", "N"]],
    "TYR": [["CA", "C"], ["CA", "N"]],
    "VAL": [["CA", "C"], ["CA", "N"]],
    "GTP": [["C1'", "N9"]],
    "NAP": []
}
# Cuts between residues are allowed by default, exceptions:
inter_residue_nocut = [{"N", "C"}]
# Formally charged atoms by default.
# Only used to suggest a charge for the QM region selected.
# list(list(resname, atomname), ...)
negative = [
    ["ASP", "OD2"], ["GLU", "OE2"], ["CTER", "OT2"], ["CYM", "SG"], ["OH", "O1"], ["CLA", "CLA"],
    ["ALA", "OT2"], ["ARG", "OT2"], ["ASN", "OT2"], ["ASP", "OT2"], ["CYS", "OT2"], ["GLN", "OT2"],
    ["GLU", "OT2"], ["GLY", "OT2"], ["HIS", "OT2"], ["HSD", "OT2"], ["HID", "OT2"], ["HSE", "OT2"],
    ["HIE", "OT2"], ["HSP", "OT2"], ["ILE", "OT2"], ["LEU", "OT2"], ["LYS", "OT2"], ["MET", "OT2"],
    ["PHE", "OT2"], ["PRO", "OT2"], ["SER", "OT2"], ["THR", "OT2"], ["TRP", "OT2"], ["TYR", "OT2"],
    ["VAL", "OT2"],
    ["GTP", "O1A"], ["GTP", "O1B"], ["GTP", "O1G"], ["GTP", "O2G"],
    ["GTP", "O1A"], ["GTP", "O1B"], ["GTP", "O1G"], ["GTP", "O2G"],
]
positive = [
    ["ARG", "NH1"], ["HSP", "NE2"], ["LYS", "NZ"], ["NTER", "N"], ["GLYP", "N"], ["LIT", "LIT"],
    ["CES", "CES"], ["POT", "POT"], ["RUB", "RUB"], ["SOD", "SOD"],
]
plus2 = [
    ["MG", "MG"], ["BAR", "BAR"], ["CAL", "CAL"], ["CD2", "CD"], ["ZN2", "ZN"]
]

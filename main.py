import argparse
import logging
from trimmer import Trimmer


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('structure', help="Structure file, e.g. PDB")
    parser.add_argument('topology', help="Corresponding topology file, e.g. PSF")
    parser.add_argument("selection", help="Central residue in segid:resid format e.g. PROA:43\n"
                                          "Optionally, specify an atom after a dash e.g. PROB:123-SG\n"
                                          "Multi-atom selection can be added with semicolon separator")
    parser.add_argument("-r", "--radius", default=5.0, type=float,
                        help="Radius for selection around said residue")
    parser.add_argument("-v", "--verbose", action="store_true", default=False,
                        help="Print verbose output.")
    parser.add_argument("--noh", action="store_true", default=False,
                        help="Disregard hydrogens for initial selection.")
    parser.add_argument("-o", "--output", default="selection",
                        help="Output suffix (selection)")
    parser.add_argument("-l", "--link", default="charmm",
                        help="Method of link atom placement (charmm|geometry)")
    parser.add_argument("-c", "--deleteCH", action="store_true", default=False,
                        help="Delete methane and ethane fragments from the selection.")
    parser.add_argument("-ms", "--multi-structure", nargs="+", default=[],
                        help="Additional structure to consider during selection")
    parser.add_argument("-mt", "--multi-topology", nargs="+", default=[],
                        help="Additional topologies for structures (eg. for different ligands)")
    args = parser.parse_args()
    if args.verbose:
        logging.basicConfig(encoding='utf-8', level=logging.DEBUG)
    else:
        logging.basicConfig(encoding='utf-8', level=logging.INFO)
    # end of parameters
    trimmer = Trimmer()
    trimmer.read_and_clean(args.topology, args.structure)
    for s in args.multi_structure:
        if len(args.multi_topology) == 0:
            trimmer.read_and_clean(args.topology, s, additional=True)
        elif len(args.multi_topology) == len(args.multi_structure):
            trimmer.read_and_clean(args.multi_topology[args.multi_structure.index(s)], s, additional=True)
        else:
            logging.error("Multi-structure selection requires multiple topologies, or omit for default.")
    trimmer.get_core(args.selection, args.radius, heavy=args.noh)
    if len(args.multi_structure) > 0:
        trimmer.merge_core(args.verbose)
    trimmer.build_selection(write=args.output, link=args.link, delete_CH=args.deleteCH)

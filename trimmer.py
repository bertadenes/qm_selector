from static import (cutable, directional, inter_residue_nocut,
                    positive, negative, plus2)
import logging
import numpy as np
import networkx as nx
import MDAnalysis as mda
from MDAnalysis.topology.guessers import guess_types


class Trimmer:
    """
    Class to make a chemically sensible trimmed structure centred around a protein residue for QM cluster or QM/MM
    calculations.

    Attributes
    structure (mda.Universe): MDanalysis universe of the whole system
    core (mda.AtomGroup): MDanalysis object of the initial selection
    selection (ndarray): indices of the selected atoms
    links (ndarray): indeces of cut bonds [0, :] contains inside atoms, [1, :] has the outside ones.
    """
    def __init__(self):
        """
        Instance variable placeholders.
        """
        self.structure = None
        self.multi_structure = []
        self.structure_files = []
        self.core_indices = None
        self.multi_core = []
        self.selection = None
        self.string_selection = None
        self.links = None
        self.radius = None
        self.centre = None
        self.ressel = None
        self.atom2fragment = None
        self.fragment2atom = None
        return

    def read_and_clean(self, top, structure, additional=False):
        """
        Tested with psf and prmtop topology files.
        It removes UB term bonds from TIP3 waters, to conform rdkit conversion.
        Similar actions may be needed in case of other topology specific problems.

        Args:
            top (str): path to topology file
            structure (str): path to structure file
            additional (bool): specify if this is the base structure or an extra

        Returns:
        """
        s = mda.Universe(top, structure)
        # elements are needed for rdkit
        guessed_elements = guess_types(s.atoms.names)
        s.add_TopologyAttr('elements', guessed_elements)
        # remove UB from TIP3 waters
        del_bonds = []
        waters = s.select_atoms("segid WAT? or segid SOLV")
        for r in waters.residues:
            del_bonds.append(
                [r.atoms.indices[np.where(r.atoms.names == "H1")], r.atoms.indices[np.where(r.atoms.names == "H2")]])
        del_bonds = np.array(del_bonds).reshape((len(del_bonds), 2))
        s.delete_bonds(del_bonds)
        if additional:
            self.multi_structure.append(s)
        else:
            self.structure = s
        self.structure_files.append(structure)
        return

    def get_core(self, selection, base_radius, write=False, heavy=False):
        """
        Defines the core selection around the user defined list.

        :param selection: definition by the segid plus residue number (e.g. PROA:43)
        :type selection: str
        :param base_radius: selection radius around the central residue in Angstroms
        :type base_radius: float
        :param write: writes a file of the selection if specified
        :type write: bool
        :param heavy: ignore H atoms in selection
        :type heavy: bool
        :return:
        :rtype:
        """
        self.radius = base_radius
        self.centre = selection
        atoms = selection.split(';')
        ressel = ""
        for atom in atoms:
            identifier = atom.split('-')
            segid, resnum = identifier[0].split(':')
            if ressel != "":
                ressel += " or "
            if len(identifier) == 2:
                ressel += f"(segid {segid} and resnum {resnum} and name {identifier[1]})"
            elif len(identifier) == 1:
                ressel += f"(segid {segid} and resnum {resnum})"
        self.ressel = ressel
        if heavy:
            defined = self.structure.select_atoms(ressel)
            heavy = self.structure.select_atoms('not element H')
            core = heavy.select_atoms(f"(around {base_radius} ({ressel})) or ({ressel})")
            self.core_indices = (core | defined).indices
            for s in self.multi_structure:
                sh = s.select_atoms('not element H')
                defined = s.select_atoms(ressel)
                core = sh.select_atoms(f"(around {base_radius} ({ressel})) or ({ressel})")
                self.multi_core.append((core | defined).indices)
        else:
            self.core_indices = self.structure.select_atoms(f"(around {base_radius} ({ressel})) or ({ressel})").indices
            for s in self.multi_structure:
                self.multi_core.append(s.select_atoms(f"(around {base_radius} ({ressel})) or ({ressel})").indices)
        if write:
            self.structure.atoms[self.core_indices].write(write)
        return

    def merge_core(self, debug=True):
        """
        Combine the selection from multiple structures.
        :param debug: verbose flag
        :type debug: bool
        :return:
        :rtype:
        """
        logging.info(f"{len(self.core_indices)} atoms in the base core selection")
        # self.core_indices = np.unique(np.concatenate([self.core_indices, np.concatenate(self.multi_core)])).tolist()
        core = []
        for i in self.core_indices:
            core.append(f"(atom {self.structure.atoms[i].segid} "
                        f"{self.structure.atoms[i].resid} {self.structure.atoms[i].name})")
        for i in range(len(self.multi_core)):
            for j in self.multi_core[i]:
                card = (f"(atom {self.multi_structure[i].atoms[j].segid} "
                        f"{self.multi_structure[i].atoms[j].resid} {self.multi_structure[i].atoms[j].name})")
                if card not in core:
                    core.append(card)
        for s in core:
            logging.debug(s)
        selection = " or ".join(core)
        self.string_selection = selection
        self.core_indices = self.structure.select_atoms(selection).indices
        logging.info(f"{len(self.core_indices)} atoms in the combined core selection")
        for i in range(len(self.multi_core)):
            self.multi_core[i] = self.multi_structure[i].select_atoms(selection).indices
        if debug:
            logging.debug("WARNING: this selection is not updated for multiple topologies.")
            j = 1
            for l in self.multi_core:
                logging.debug(f"merging in {len(l)} atoms from additional structure #{j}")
                j += 1
        return

    def extend_selection(self, delete_CH=True):
        """
        Performs an iterative extension of the core selection based on external rules. The process stops when it finds
        nothing new to be added. It also defined bonds broken.

        Args:
            delete_CH (bool): if True, removes small aliphatic bits (up to ethane).

        Returns:

        """
        extended = True
        to_be_linked = []
        selection = self.core_indices
        while extended:
            extended = False
            new_selection = selection
            for i in selection:
                for j in self.structure.atoms[i].bonded_atoms.indices:
                    if j not in selection:
                        resp = self.check_interface(self.structure.atoms[j], self.structure.atoms[i])
                        if resp["include"] is not None:
                            extended = True
                            new_selection = np.append(new_selection, resp["include"])
                        elif resp["link"] and (i, j) not in to_be_linked:
                            to_be_linked.append((i, j))
            selection = new_selection
        self.selection = np.unique(selection)
        to_be_linked = np.array(to_be_linked)
        if delete_CH:
            self.selection = self.clear_aliphatic(self.structure, self.selection)
        # cleaning up cut bonds
        logging.info("Links are to be made between:")
        delete_row = []
        for i in range(to_be_linked.shape[0]):
            if to_be_linked[i, 1] in self.selection:
                delete_row.append(i)
                logging.debug("Removing unnecessary link atom to " + self.structure.atoms[to_be_linked[i, 0]].__str__())
            elif to_be_linked[i, 0] not in self.selection:
                delete_row.append(i)
                logging.debug(
                    "QM atom removed (likely aliphatic): " + self.structure.atoms[to_be_linked[i, 0]].__str__())
            else:
                logging.info(f"QM atom {to_be_linked[i, 0]:d} and MM atom {to_be_linked[i, 1]:d}")
                logging.debug(self.structure.atoms[to_be_linked[i, 0]].__str__() + self.structure.atoms[
                    to_be_linked[i, 1]].__str__())
        self.links = np.delete(to_be_linked, delete_row, axis=0)
        logging.debug(self.selection)
        return

    def bridge_links(self):
        """
        Sometimes multiple links to the same MM atom are created, resulting in close link atoms.
        Desired behaviour is to include the MM atoms to the selection.
        :return: Whether changes were made to the core selection or not.
        :rtype: bool
        """
        MMatom, counts = np.unique(self.links[:, 1], return_counts=True)
        if np.max(counts) > 1:
            logging.debug("Continuing QM extension, as links to the same atom are identified.")
            self.core_indices = np.sort(np.concatenate([self.selection, MMatom[np.where(counts > 1)]]))
            return True
        return False

    def build_selection(self, write="selection", delete_CH=True, link="charmm"):
        """
        Top loop for QM extension. May include sub-loops depending on core selection updates.
        :param write: filename for the PDB file to write the selection. False is taken for none.
        :type write: String
        :return:
        :rtype:
        """
        iterate = True
        while iterate:
            self.extend_selection(delete_CH=delete_CH)
            iterate = self.bridge_links()
        if write:
            if link == "charmm":
                self.write_charmm_input(f"charmm_qmmm_{self.structure_files[0].split('.')[0]}_{write}.stream")
            selected = self.structure.atoms[self.selection]
            selected.write(f"{self.structure_files[0].split('.')[0]}_{write}.pdb")
            # mm = self.structure.atoms - selected
            # mm_charges = np.concatenate([mm.positions, mm.charges[:, None]], axis=1)
            # np.savetxt("mm_external.txt", mm_charges, delimiter=" ")
            cards = []
            for i in selected.indices:
                cards.append(f"(atom {self.structure.atoms[i].segid} "
                             f"{self.structure.atoms[i].resid} {self.structure.atoms[i].name})")
            if link == "geometry":
                self.add_links_to_pdb(f"{self.structure_files[0].split('.')[0]}_{write}.pdb")
            for i in range(len(self.multi_structure)):
                core = self.multi_structure[i].select_atoms(self.string_selection)
                final = self.multi_structure[i].select_atoms(' or '.join(cards))
                idx = (core | final).indices
                if delete_CH:
                    idx = self.clear_aliphatic(self.multi_structure[i], idx)
                selected = self.multi_structure[i].atoms[idx]
                # mm = self.multi_structure[i].atoms - selected
                # mm_charges = np.concatenate([mm.positions, mm.charges[:, None]], axis=1)
                selected.write(f"{self.structure_files[i+1].split('.')[0]}_{write}.pdb")
                if link == "geometry":
                    self.add_links_to_pdb(f"{self.structure_files[i+1].split('.')[0]}_{write}.pdb", i)
        else:
            if link == "charmm":
                self.write_charmm_input()
        return

    def add_links_to_pdb(self, filename, index=None):
        """
        Add link atoms to a PDB file written by MDA.
        :param filename: PDB filename to be modified
        :type filename: string
        :param index: structure index for processing multiple structures
        :type index: int
        :return:
        :rtype:
        """
        with open(filename, "r") as f:
            lines = f.readlines()
        last_atom = 0
        geom_written = False
        links_written = False
        with open(filename, "w") as f:
            for l in lines:
                if l.startswith("ATOM") or l.startswith("HETATM"):
                    last_atom = int(l[7:11])
                    geom_written = True
                elif geom_written and not links_written:
                    for link in self.links:
                        last_atom += 1
                        resname = self.structure.atoms[link[0]].resname
                        resid = self.structure.atoms[link[0]].resid
                        segid = self.structure.atoms[link[0]].segid
                        try:
                            cID = self.structure.atoms[link[0]].chainID
                        except mda.exceptions.NoDataError:
                            cID = "X"
                        if index is None:
                            qm2mm = self.structure.trajectory[0].positions[link[1]] - \
                                    self.structure.trajectory[0].positions[link[0]]
                            linkpos = self.structure.trajectory[0].positions[link[0]] + qm2mm * 0.7261
                            try:
                                cID = self.structure.atoms[link[0]].chainID
                            except mda.exceptions.NoDataError:
                                cID = "X"
                        else:
                            qmatom = self.multi_structure[index].select_atoms(f"atom {segid} {resid} "
                                                                              f"{self.structure.atoms[link[0]].name}")
                            mmatom = self.multi_structure[index].select_atoms(f"atom "
                                                                              f"{self.structure.atoms[link[1]].segid} "
                                                                              f"{self.structure.atoms[link[1]].resid} "
                                                                              f"{self.structure.atoms[link[1]].name}")
                            qm2mm = mmatom.positions[0] - qmatom.positions[0]
                            linkpos = self.multi_structure[index].trajectory[0].positions[qmatom.atoms.indices] + \
                                      qm2mm * 0.7261
                        linkpos = linkpos.reshape(3)
                        f.write(f"ATOM  {last_atom:5d}  QH  {resname:3s} {cID:1s}{resid:4d}    "
                                f"{linkpos[0]:8.3f}{linkpos[1]:8.3f}{linkpos[2]:8.3f}{1.0:6.2f}{0.0:6.2f}      "
                                f"{segid[:4]:4s} H  \n")
                    links_written = True
                f.write(l)
        return

    def write_charmm_input(self, filename="charmm_qmmm.inp", trim=True):
        """
        Writes a selection for QM region in charmm format, as well as linkatom definitions.
        Makes a suggestion to the formal charge of the QM region. Atom formal charge guessing is not yet implemented
        in MDanalysis, therefore it is based on atoms names in CHARMM36m. If you use another convention, change the
        corresponding lists.
        It also defines two selections for non-PBC QM/MM, if trim==True:
        fixsel: 25 A around the residue centred
        flexsel: 20 A around the residue centred

        Returns:

        """
        selected = self.structure.atoms[self.selection]
        with open(filename, "w") as f:
            f.write("* generated during automated trimming\n")
            f.write("* \n")
            f.write(f"! core selection is {self.radius} around {self.centre}\n")
            f.write(f"! defining {selected.n_atoms} QM atoms\n")
            if trim:
                centre = self.structure.select_atoms(self.ressel)
                # remove atoms outside 25 A
                f.write(f"define fixsel select .byres. ((segid {centre.segids[0]} .and. resi {centre.resids[0]})"
                        f".around. 25.0) show end\n")
                f.write(f"delete atom select .not. fixsel end\n")
                f.write("\nset CHG ?CGTOT\n")
            # add link atoms first, so selection is correct afterward
            for i in range(self.links.shape[0]):
                qmatom = self.structure.atoms[self.links[i, 0]]
                mmatom = self.structure.atoms[self.links[i, 1]]
                f.write(f"addlinkatom QQH{i + 1:d}"
                        f" {qmatom.segid:s} {qmatom.resid:d} {qmatom.name:s}"
                        f" {mmatom.segid:s} {mmatom.resid:d} {mmatom.name:s}\n")
            # QM selection
            # has to split it to avoid long line errors in charmm
            c = 1
            j = 0
            formal_charge = 0
            for i in range(selected.n_atoms):
                if j == 0:
                    f.write(f"define qm{c} sele -\n")
                a = selected.atoms[i]
                f.write(f"(atom {a.segid:6s}{a.resid:5d} {a.name:6s})")
                if i + 1 == selected.n_atoms or j == 49:
                    f.write(" - \n")
                else:
                    f.write(" .or. - \n")
                j += 1
                if [a.resname, a.name] in positive:
                    formal_charge += 1
                elif [a.resname, a.name] in negative:
                    formal_charge -= 1
                elif [a.resname, a.name] in plus2:
                    formal_charge += 2
                if j == 50:
                    f.write("end\n")
                    c += 1
                    j = 0
            if j != 0:
                f.write("end\n")
            f.write("define qm sele -\n")
            f.write(" .or. ".join([f"qm{i + 1}" for i in range(c)]))
            f.write(" show end\n")
            f.write(f"! the formal charge of the above selection is guessed: {formal_charge}\n")
            # link atom positioning
            for i in range(self.links.shape[0]):
                qmatom = self.structure.atoms[self.links[i, 0]]
                mmatom = self.structure.atoms[self.links[i, 1]]
                f.write("lonepair colinear scaled -.7261 -\n")
                f.write(f"sele type QQH{i + 1:d} show end"
                        f" sele atom {qmatom.segid:s} {qmatom.resid:d} {qmatom.name:s} show end"
                        f" sele atom {mmatom.segid:s} {mmatom.resid:d} {mmatom.name:s} end\n")
            # setting link atom mass
            f.write("\nscalar mass show select type qq* show end\n")
            f.write("scalar mass set  1.008000 select type qq* show end\n")
            f.write("scalar mass show select type qq* show end\n")
            if trim:
                # selecting atoms around 20 A, fix the rest
                f.write(f"define flexsel select .byres. ((segid {centre.segids[0]} .and. resi {centre.resids[0]})"
                        f".around. 20.0) show end\n")
                f.write("cons fix select none end\n")
                f.write("cons fix select ( .not. flexsel ) show end\n")
            # writing the QM selection out
            f.write("open write unit 10 card name ./qm_region.pdb\n")
            f.write("write coor pdb unit 10 card sele qm .or. type qq* end\n")
            f.write("\n")
        return

    @staticmethod
    def clear_aliphatic(structure, selection):
        """
        Takes the instance selection and removes separate aliphatic parts up to 2 heavy atoms which may remain from
        the geometric/topological selection. It relies on conversion to rdkit Mol format.

        Returns:

        """
        selected = structure.atoms[selection]
        to_remove = []
        # graph approach
        dm = np.linalg.norm(selected.positions[None, :] - selected.positions[:, None], axis=2)
        graph = nx.from_numpy_array(dm <= 2)
        for c in nx.connected_components(graph):
            hetatom = False
            c_count = 0
            for i in c:
                if selected.elements[i] not in ["C", "H"]:
                    hetatom = True
                elif selected.elements[i] == "C":
                    c_count += 1
            if not hetatom and c_count < 3:
                to_remove += list(c)
        return np.delete(selection, to_remove)

    @staticmethod
    def check_interface(outatom, inatom):
        """
        Processing function of an interface bond in protein context. Topology object are from MDanalysis.

        Args:
            outatom (mda.Atom): atom bonded to inatom, but outside the selection of interest.
            inatom (mda.Atom): atom in the selection with a dangling valance.

        Returns:
            dict: context, defining atoms to be included or link to be defined
            :param inatom:
            :type inatom:
            :return:
            :rtype:
        """
        resp = {"include": None, "link": False}
        if outatom.residue == inatom.residue:
            try:
                if {inatom.name, outatom.name} in cutable[inatom.resname]:
                    logging.debug("bond cut: " + outatom.__str__() + inatom.__str__())
                    resp["link"] = True
                elif [inatom.name, outatom.name] in directional[inatom.resname]:
                    logging.debug("bond cut: " + outatom.__str__() + inatom.__str__())
                    resp["link"] = True
                else:
                    logging.debug("including " + outatom.__str__())
                    resp["include"] = outatom.index
            except KeyError:
                logging.debug(f"residue {inatom.resname.__str__()} has no cutable bond registered")
                resp["include"] = outatom.index
        else:
            if {inatom.name, outatom.name} in inter_residue_nocut:
                logging.debug("extending to residue" + outatom.__str__())
                resp["include"] = outatom.index
        return resp

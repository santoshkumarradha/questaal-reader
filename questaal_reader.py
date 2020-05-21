from dotmap import DotMap
from pymatgen import Structure
import numpy as np
import os
ry2ev = 13.605662285137
from math import ceil
from itertools import chain, product
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.electronic_structure.core import Spin, Orbital
from pymatgen.core.lattice import Lattice
#from sumo.io.questaal import QuestaalSite, labels_from_syml
import numpy as np
from pymatgen.core.structure import Structure
import os
import warnings
import re
_bohr_to_angstrom = 0.5291772
_ry_to_ev = 13.605693009


class reader:
    def __init__(self, output_file="output"):
        self.fname = output_file
        self.data = get_data(self.fname)
        self.iterations = make_iterations(self.data)
        # self.natoms=get_nbas(self.data)
        # self.nbas=self.natoms

        #---finding and setting up structure
        self.species = get_species(self.data)
        self.atoms = self.species
        self.structure = make_structure(self.data)

        #---setting energy data
        self.niter = len(self.iterations)
        set_iteration_energy(self.data, self.iterations)
        get_charges(self.data, self.iterations)

        #---setting energy data
        self.energy = self.iterations[-1].ehf
        self.ehf = self.iterations[-1].ehf
        self.ehk = self.iterations[-1].ehk
        #---setting band data in eV
        set_band_data(self.data, self.iterations)
        self.gap = self.iterations[-1].gap
        self.valance_band_max = self.iterations[-1].valance_band_max
        self.conduction_band_min = self.iterations[-1].conduction_band_min

    def get_variables(self):
        lst = {
            "data": "raw data string",
            "Iterations":
            "Iteration object, further contains other data about iterations",
            "structure": "calculations structure (returns in pymatgen format",
            "energy": "total energy in (ev) of the final iteration (ehf)",
            "ehf": "final energy ehf",
            "ehk": "final energy ehk",
            "atoms": "species name in lmf codes",
            "gap": "final band gap from given k mesh (eV)",
            "valance_band_max": "valance band max energy in code(eV)",
            "conduction_band_min": "conduction band min in code (eV)"
        }
        for i in lst:
            print(i, "-", lst[i])


def get_data(fname="output"):
    try:
        with open(fname, 'r') as f:
            data = f.read()
        return data
    except IOError:
        print("file not found")


def get_lines(data, key, num_lines=0, return_index=False):
    '''
    get the num_lines along with line matching "key" in data text
    if return_index returns the index
    '''
    index = [i for i, s in enumerate(data.splitlines()) if key in s]
    values = [data.splitlines()[i:i + num_lines + 1] for i in index]
    index_vals = [list(range(i, i + num_lines + 1)) for i in index]
    values = [list(filter(lambda name: name.strip(), i)) for i in values]
    # if return_index:
    #   return values,index_vals
    # else:
    return values


def make_structure(data):
    screen_charge, charge_dict = get_final_charge(data)

    return Structure(lattice=get_lattice(data),
                     species=get_z(data, get_nbas(data)),
                     coords=get_atomicpos(data, frac=True),
                     coords_are_cartesian=False,
                     charge=screen_charge,
                     site_properties=charge_dict)


def get_nbas(data):
    ''' 
    Extract number of atoms
    '''
    try:
        nbas = int(
            get_lines(data, "nbas")[0][0].partition('nbas = ')[2].split()[0])
    except ValueError:
        print("unable to fing number of atoms")
    return nbas


def get_species(data):
    '''
    return species type for all atoms present
    '''
    natoms = get_nbas(data)
    return [
        ''.join(i[0].split()[-1].split(':')[::-1]).lower()
        for i in get_lines(data, "   species  ")[:natoms]
    ]


def make_iterations(data):
    '''
    make an iteration object to hold information
    '''
    iterations = []
    for i in range(get_iter(data)):
        dummy = DotMap()
        dummy.niter = i
        iterations.append(dummy)
    return iterations


def get_iter(data):
    '''
    returns total num of iterations
    '''
    return int(
        get_lines(
            data, "iteration",
            return_index=False)[-1][0].split("iteration ")[-1].split()[0])


def get_z(data, natoms):
    '''get the Z's of atoms'''
    z = []
    for j in range(natoms):
        z.append(
            float(
                get_lines(data,
                          " site  {}  z=".format(j + 1))[0][0].split()[3]))
    z = np.array(z).astype(np.int)
    return z


def get_atomicpos(data, frac=True):
    '''
    data: string with output
    frac: bool returns frac coordinates 
    '''
    natoms = get_nbas(data)
    if frac:
        start = 5
        frac_pos = [
            i.split()[start:start + 3] for i in get_lines(
                data, "pos (Cartesian coordinates)", num_lines=natoms)[0][1:]
        ]
        frac_pos = np.array(frac_pos).astype(np.float)
        return frac_pos
    else:
        start = 3
        cart_pos = [
            i.split()[start:start + 3] for i in get_lines(
                data, "pos (Cartesian coordinates)", num_lines=natoms)[0][1:]
        ]
        cart_pos = np.array(cart_pos).astype(np.float) / 1.8897259885789
        return cart_pos


def get_lattice(data):
    '''gets the lattice vector'''
    lat = get_lines(data, "Plat", num_lines=5)
    plat = np.array([i.split()[:3] for i in lat[0][1:4]]).astype(np.float)
    alat = float(lat[0][-1].split("alat = ")[1].split()[0])
    lattice = alat * plat / 1.8897259885789
    return lattice


def get_energy(data):
    '''gets the energy at each iteration'''
    ehf = [i[0].split()[2].split("=")[-1] for i in get_lines(data, " nit=")]
    ehk = [i[0].split()[3].split("=")[-1] for i in get_lines(data, " nit=")]
    ehf = np.array(ehf).astype(np.float) * ry2ev
    ehk = np.array(ehk).astype(np.float) * ry2ev
    return ehf, ehk


def set_iteration_energy(data, iterations):
    '''set energy to iteration object'''
    ehf, ehk = get_energy(data)
    for j in range(len(iterations)):
        iterations[j].ehf = ehf[j]
        iterations[j].ehk = ehf[j]
        iterations[j].energy = ehf[j]


def get_charges(data, iterations=None):
    '''
    get charge data for each iterations
    '''
    if iterations is None:
        iterations = make_iterations(data)
    for j in range(len(iterations)):
        iter_i = iterations[j].niter
        key = "charges:       old"
        natoms = get_nbas(data)
        iter_data_txt = get_lines(data, key, natoms + 1,
                                  return_index=False)[iter_i]
        natoms = get_nbas(data)
        species = get_species(data)
        iteration_charge_data = DotMap()
        charge = DotMap()
        niter = get_iter(data)
        for i in range(natoms + 1):
            line_data = iter_data_txt[1:][i].split()
            if i == 0:
                name = "smooth"
                tmp = 0
            else:
                name = species[i - 1]
                tmp = 1
            iteration_charge_data[name].old_charge = float(line_data[tmp + 1])
            iteration_charge_data[name].new_charge = float(line_data[tmp + 2])
            iteration_charge_data[name].screened_charge = float(line_data[tmp +
                                                                          3])
            iteration_charge_data[name].rms_charge = float(line_data[tmp + 4])
            iteration_charge_data[name].diff_charge = float(line_data[tmp + 5])
            iteration_charge_data[name].charge = float(line_data[tmp + 2])
        iterations[j].charge = iteration_charge_data

    if iterations is None:
        return iterations


def get_final_charge(data):
    natoms = get_nbas(data)
    key = "charges:       old"
    natoms = get_nbas(data)
    charge_data = get_lines(data, key, natoms + 1, return_index=False)[-1]
    screen_charge = float(charge_data[1:][0].split()[2])
    final_charges = [float(i.split()[3]) for i in charge_data[1:][1:]]
    charge_dict = {"charge": final_charges}
    return screen_charge, charge_dict


def get_band_data(data):
    '''get band gap data'''
    vals = [[i[0].split()[2], i[0].split()[5], i[0].split()[11]]
            for i in get_lines(data, "gap")]
    vals = np.array(vals).astype(float)
    valance_band_max = vals.T[0]
    conduction_band_min = vals.T[1]
    gap = vals.T[2]
    return valance_band_max, conduction_band_min, gap


def set_band_data(data, iterations):
    '''set to iterations data'''
    valance_band_max, conduction_band_min, gap = get_band_data(data)
    for i in range(len(iterations)):
        iterations[i].valance_band_max = valance_band_max[i] * ry2ev
        iterations[i].conduction_band_min = conduction_band_min[i] * ry2ev
        iterations[i].gap = gap[i]


#------------------- Helper for reading bnds file. make sure site file, bnds and syml file are present
# -------- and


# Class ported from Sumo
class QuestaalSite(object):
    def __init__(self,
                 nbas,
                 vn=3.,
                 io=15,
                 alat=1.,
                 xpos=True,
                 read='fast',
                 sites=None,
                 plat=[1, 0, 0, 0, 1, 0, 0, 0, 1]):
        sites = sites or []
        if nbas != len(sites):
            raise AssertionError()
        if len(plat) != 9:
            raise AssertionError()

        if read != 'fast':
            raise Exception("Algebraic expressions not supported, use 'fast'")
        if io != 15:
            warnings.warn(
                "Only site.ext format 15 supported at present \n if things dont work That might be the problem"
            )

        self.nbas, self.vn, self.io, self.alat = nbas, vn, io, alat
        self.xpos, self.read, self.sites, self.plat = xpos, read, sites, plat

        is_empty = re.compile('E\d*$')
        empty_sites = [
            site for site in sites
            if is_empty.match(site['species']) is not None
        ]
        self.nbas_empty = len(empty_sites)

    @property
    def structure(self):
        # Get lattice vectors in Angstrom
        lattice = Lattice(self.plat)
        lattice = Lattice(lattice.matrix * self.alat * _bohr_to_angstrom)

        # Get corresponding lists of species and positions by making a list of
        # pairs and unpacking with zip
        if self.xpos:
            species_coords = [(site['species'], site['pos'])
                              for site in self.sites]
            species, coords = zip(*species_coords)

            return Structure(lattice,
                             species,
                             coords,
                             coords_are_cartesian=False)
        else:
            species_coords = [
                (site['species'],
                 [x * self.alat * _bohr_to_angstrom for x in site['pos']])
                for site in self.sites
            ]
            species, coords = zip(*species_coords)

            return Structure(lattice,
                             species,
                             coords,
                             coords_are_cartesian=True)

    @classmethod
    def from_file(cls, filename):
        with open(filename, 'rt') as f:
            lines = f.readlines()

        header = lines[0]
        sites = [line for line in lines if line[0] not in '#%']

        # Some of the header info does not use '=' so handle separately
        header_items = header.strip().split()
        if header_items[0] != '%' or header_items[1] != 'site-data':
            raise AssertionError()

        xpos = True if 'xpos' in header_items else False
        read = 'fast' if 'fast' in header_items else False

        header_clean = ' '.join(x for x in header_items
                                if x not in ('%', 'site-data', 'xpos', 'fast'))

        tags = re.findall(r'(\w+)\s*=', header_clean)  # Find tags
        # Split on tags to find tag parameters
        tag_data = re.split(r'\s*\w+\s*=\s*', header_clean)[1:]
        tag_dict = dict(zip(tags, tag_data))

        vn = float(tag_dict['vn']) if 'vn' in tag_dict else 3.
        io = int(tag_dict['io']) if 'io' in tag_dict else 15.
        nbas = int(tag_dict['nbas']) if 'nbas' in tag_dict else 15.
        alat = float(tag_dict['alat']) if 'alat' in tag_dict else 1.
        plat = ([float(x) for x in tag_dict['plat'].split()]
                if 'plat' in tag_dict else [1, 0, 0, 0, 1, 0, 0, 0, 1])

        # Convert sites to structured format
        # (If needed, could support other 'io' options here)
        sites = [{
            'species': site.split()[0],
            'pos': [float(x) for x in site.split()[1:4]]
        } for site in sites]

        return cls(nbas,
                   vn=vn,
                   io=io,
                   alat=alat,
                   xpos=xpos,
                   read=read,
                   sites=sites,
                   plat=plat)


def labels_from_syml(syml_file):
    labels = {}

    with open(syml_file, 'r') as f:
        lines = f.readlines()
    for line in lines:
        npts, x1, y1, z1, x2, y2, z2, *label_text = line.split()
        if len(label_text) < 3:
            pass
        else:
            kpt1 = tuple(map(float, (x1, y1, z1)))
            kpt2 = tuple(map(float, (x2, y2, z2)))

            label_text = ' '.join(label_text)  # Undo previous split
            label1_label2 = label_text.split(' to ')
            if len(label1_label2) != 2:
                raise ValueError("Not clear how to interpret labels from "
                                 "this line: {}".format(line))
            label1, label2 = label1_label2
            labels.update({label1: kpt1, label2: kpt2})
    return labels


def get_bands(fname):
    '''
    returns a pymatgen BandStructureSymmLine object for easy plotting.
    '''
    filenames = [fname]
    bnds_file = filenames[0]
    ext = bnds_file.split('.')[-1]
    bnds_folder = os.path.join(bnds_file, os.path.pardir)

    site_file = os.path.abspath(
        os.path.join(bnds_folder, 'site.{}'.format(ext)))

    if os.path.isfile(site_file):
        site_data = QuestaalSite.from_file(site_file)
        bnds_lattice = site_data.structure.lattice
        alat = site_data.alat
    else:
        raise IOError('Site file {} not found: '
                      'needed to determine lattice'.format(site_file))

    syml_file = os.path.abspath(
        os.path.join(bnds_folder, 'syml.{}'.format(ext)))
    if os.path.isfile(syml_file):
        bnds_labels = labels_from_syml(syml_file)
    else:
        bnds_labels = {}

    with open(bnds_file, 'r') as f:
        kpoints = []

        # Read heading, get metadata and check no orbital projections used
        nbands, efermi, n_color_wts, *_ = f.readline().split()

        if int(n_color_wts) > 0:
            raise NotImplementedError(
                "Band data includes orbital data: "
                "this format is not currently supported.")

        nbands, efermi = int(nbands), float(efermi)
        eig_lines = ceil(nbands / 10)

        # Check if first two kpts are the same: if so, assume there are two
        # spin channels
        _ = f.readline()
        kpt1 = list(map(float, f.readline().split()))

        for line in range(eig_lines):  # Skip over the eigenvalues
            _ = f.readline()  # for now: re-read file later
        kpt2 = list(map(float, f.readline().split()))
        if len(kpt1) != 3 or len(kpt2) != 3:
            raise AssertionError()

        if kpt1 == kpt2:
            spin_pol = True
        else:
            spin_pol = False

    def _read_eigenvals(f, nlines):
        lines = [f.readline() for i in range(nlines)]
        # This statement is very "functional programming"; read it
        # backwards.  List of split lines is "flattened" by chain into
        # iterator of values; this is fed into map to make floats and
        # stored to a list
        return list(map(float, chain(*(line.split() for line in lines))))

    with open(bnds_file, 'r') as f:
        _ = f.readline()
        if spin_pol:  # Need to read two spin channels
            block_nkpts = int(f.readline().strip()) // 2
            eigenvals = {Spin.up: [], Spin.down: []}
        else:
            block_nkpts = int(f.readline().strip())
            eigenvals = {Spin.up: []}

        while block_nkpts > 0:  # File should be terminated with a 0
            for i in range(block_nkpts):
                kpoint = list(map(float, f.readline().split()))
                kpoints.append(np.array(kpoint) / (alat * _bohr_to_angstrom))

                eigenvals[Spin.up].append(_read_eigenvals(f, eig_lines))

                if spin_pol:
                    spin_down_kpoint = list(map(float, f.readline().split()))
                    if spin_down_kpoint != kpoint:
                        raise AssertionError(
                            "File interpreted as spin-polarised, but this"
                            " kpoint only has one entry: {}".format(kpoint))
                    eigenvals[Spin.down].append(_read_eigenvals(f, eig_lines))

            block_nkpts = int(f.readline().strip())
            if spin_pol:
                block_nkpts = block_nkpts // 2
    eigenvals = {
        key: np.array(data).T * _ry_to_ev
        for key, data in eigenvals.items()
    }
    efermi *= _ry_to_ev
    if os.path.isfile(site_file):
        site_data = QuestaalSite.from_file(site_file)
        bnds_lattice = site_data.structure.lattice
        alat = site_data.alat
    coords_are_cartesian = True
    labels = bnds_labels
    if coords_are_cartesian:
        for label, coords in labels.items():
            labels[label] = np.array(coords) / (alat * _bohr_to_angstrom)
    else:
        for label, coords in labels.items():
            labels[label] = np.dot(
                coords,
                bnds_lattice.reciprocal_lattice_crystallographic.matrix)
    return BandStructureSymmLine(
        kpoints,
        eigenvals,
        bnds_lattice.reciprocal_lattice_crystallographic,
        efermi,
        labels,
        coords_are_cartesian=True)

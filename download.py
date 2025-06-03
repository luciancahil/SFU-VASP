from mp_api.client import MPRester
from pymatgen.io.vasp import Poscar

api_key = "6RbbnmuRQ49U7m01pOyZ29LPU1EIIJHv"  # Your Materials Project API key
material_id = "mp-23445"

with MPRester(api_key) as mpr:
    structure = mpr.get_structure_by_material_id(material_id)
    poscar = Poscar(structure)
    poscar.write_file("POSCAR")
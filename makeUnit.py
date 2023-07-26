from pymatgen.core.structure import Structure
from pymatgen.core.periodic_table import Element
from pathlib import Path

import shutil
from string import Template
a_str = "Cs"
b_strs = ["Pb", "Sn", "Ge"]
x_strs= ["I", "Br", "Cl"]
pot_dir = "/depot/amannodi/data/Pseudopotentials/POTCARs/"

for b_str in b_strs:
    for x_str in x_strs:
        b = Element(b_str)
        x = Element(x_str)
        pervoskite_base = Structure.from_file("refDIR/POSCAR")
        pervoskite_base.replace_species({Element("Pb"): b, Element("I"): x})
        new_structure = pervoskite_base
        bulk_dir = "runDIR/"+a_str + b_str + x_str + "3/bulk"
        Path(bulk_dir).mkdir(parents=True, exist_ok=True)
        new_structure.to(fmt='poscar', filename=(bulk_dir + "/POSCAR"))

        shutil.copy("refDIR/inFilesBulk/INCAR", bulk_dir + "/INCAR")
        shutil.copy("refDIR/inFilesBulk/KPOINTS", bulk_dir + "/KPOINTS")
        shutil.copy("refDIR/inFilesBulk/job.sh", bulk_dir + "/job.sh")

        with open(bulk_dir + "/job.sh", 'r') as file:
               lines = file.readlines()
        lines[5]= "#SBATCH --job-name "+ a_str +b_str +x_str +"3_bulk"
        with open(bulk_dir + "/job.sh", 'w') as file:
            file.writelines(lines)

        pot_paths = [pot_dir + a_str, pot_dir + b_str , pot_dir + x_str]
        with open(bulk_dir+'/POTCAR','wb') as potcar:
               for pot in pot_paths:
                      with open(pot,'rb') as source:
                             shutil.copyfileobj(source, potcar)

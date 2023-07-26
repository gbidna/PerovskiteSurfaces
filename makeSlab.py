from pymatgen.core.surface import *
from pymatgen.core.structure import Structure
import pathlib
import shutil

pot_dir = "/depot/amannodi/data/Pseudopotentials/POTCARs/"

a_spec = "Cs"
b_spec = ["Pb", "Sn", "Ge"]
x_spec = ["Cl", "I", "Br"]

for b_atom in b_spec:
    for x_atom in x_spec:
        mainPath = "runDIR/" + a_spec + b_atom  +  x_atom  + "3/"
        unit  = Structure.from_file(mainPath + "bulk/CONTCAR")

        unit_dict = unit.as_dict(0)
        a = unit_dict['lattice']['matrix'][0][0]
        unit_cells = 6

        slabs = generate_all_slabs(
            structure=unit,
            max_index=1,
            min_slab_size=(unit_cells+1),
            min_vacuum_size=20/a,
            center_slab =True,
            symmetrize=True,
            in_unit_planes=True,
            lll_reduce=True)

        n = len(slabs)
        prev_miller =(0,0,0)

        for i in range(0, n):
            slab = slabs[i]
            comp = slab.composition.as_dict()
            comp_keys = comp.keys()
            comp_name = ''.join(comp_keys)+ "3"
            print(comp_name)
            m_index = slab.miller_index
            print(m_index)
            if m_index == prev_miller:
                runName =''.join([str(ele) for ele in m_index]) + "Term2"
            else:
                runName = ''.join([str(ele) for ele in m_index]) + "Term1"
            dir_path = mainPath + runName
            pathlib.Path(dir_path).mkdir(parents=True, exist_ok=True)
            prev_miller = m_index
            slab.to(fmt='poscar', filename=(dir_path + "/POSCAR"))
            shutil.copy("refDIR/inFilesSlab/INCAR", dir_path + "/INCAR")
            shutil.copy("refDIR/inFilesSlab/KPOINTS", dir_path + "/KPOINTS")
            shutil.copy("refDIR/inFilesSlab/job.sh", dir_path + "/job.sh")
            with open(dir_path + "/job.sh", 'r') as file:
                lines = file.readlines()
                lines[5]= "#SBATCH --job-name "+ comp_name + runName
            with open(dir_path + "/job.sh", 'w') as file:
                file.writelines(lines)
            pot_paths = [pot_dir + x for x in comp_keys]

            with open(dir_path+'/POTCAR','wb') as potcar:
                for pot in pot_paths:
                    with open(pot,'rb') as source:
                        shutil.copyfileobj(source, potcar)

from pymatgen.core.surface import *
from pymatgen.core.periodic_table import Element
from pymatgen.core.sites import Site
from pymatgen.core.structure import Structure, Molecule
from pymatgen.analysis.adsorption import *
from  pathlib import Path
import shutil
pot_dir = "/depot/amannodi/data/Pseudopotentials/POTCARs/"
a_atom = "Cs"
b_spec = ["Pb", "Sn", "Ge"]
x_spec = ["Cl", "I", "Br"]
terms = ["111Term1", "111Term2"]
refDIR ="/depot/amannodi/data/MCHP_Database/Reference_compounds/"
i_comp = 0
text_size = 36
for b_atom in b_spec:
    for x_atom in x_spec:
        for term in terms:
            miller = tuple(map(int, term[:3]))
            comp_name =a_atom+b_atom+x_atom+"3"
#            blk =Structure.from_file("runDIR/"+comp_name+"/bulk/CONTCAR"
            clean_path = "runDIR/"+comp_name+"/"+term+"/"
            struct=Structure.from_file(clean_path+"CONTCAR")
            struct.make_supercell([2,2,1])
            slab =Slab(lattice=struct.lattice,
                species=struct.species_and_occu,
                coords=struct.frac_coords,
                miller_index=miller,
                oriented_unit_cell=struct,
                shift=0,
                scale_factor=np.eye(3, dtype=int),
                site_properties=struct.site_properties)
            search = AdsorbateSiteFinder(slab, list(miller))
            x=search.find_adsorption_sites()
            site_type=["hollow", "ontop", "bridge"]
            atoms =["H","O"]
            for site in site_type:
                for atom in atoms:
                    mol =Molecule([atom], [[0,0,0]])
                    slab = search.add_adsorbate(molecule=mol, ads_coord=x[site][0])
                    new_path ="runDIR/"+atom+"_adsorb/"+comp_name+"+"+term+"/"+site+"/"
                    Path(new_path).mkdir(parents=True, exist_ok=True)
                    slab.to(fmt='poscar', filename=new_path +"POSCAR")
                    comp_keys = slab.composition.as_dict().keys()
                    shutil.copy("refDIR/inFilesSlab/INCAR", new_path + "/INCAR")
                    shutil.copy("refDIR/inFilesSlab/KPOINTS", new_path + "/KPOINTS")
                    shutil.copy("refDIR/inFilesSlab/job.sh", new_path + "/job.sh")
                    with open(new_path + "/job.sh", 'r') as file:
                        lines = file.readlines()
                        lines[5]= "#SBATCH --job-name "+ comp_name + "-"+term+site
                    with open(new_path + "/job.sh", 'w') as file:
                        file.writelines(lines)
                    pot_paths = [pot_dir + x for x in comp_keys]
                    with open(new_path+'/POTCAR','wb') as potcar:
                        for pot in pot_paths:
                           with open(pot,'rb') as source:
                               shutil.copyfileobj(source, potcar)

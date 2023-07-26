import numpy as np
from pymatgen.io.vasp.outputs import Outcar
from pymatgen.core.structure import Structure
import matplotlib.pyplot as plt
import matplotlib as mpl
a_spec = "Cs"
b_spec = ["Ge", "Sn", "Pb"]
x_spec = ["I", "Br", "Cl"]
terms = ["100Term1", "100Term2", "110Term1", "110Term2", "111Term1", "111Term2"]
refDIR ="/depot/amannodi/data/MCHP_Database/Reference_compounds/"
i_comp = 0
decomps =[]
for x_atom in x_spec:
    x_dcomp=[]
    for b_atom in b_spec:
        i_comp=i_comp+1
        mainPath = "runDIR/" + a_spec + b_atom  +  x_atom  + "3/"
        bulk_energy = Outcar(mainPath+"bulk/OUTCAR").final_fr_energy

        AX_path = refDIR+"AX_halides/"+a_spec+x_atom+"/PBE/"
        AX_energy = Outcar(AX_path +"OUTCAR.gz").final_fr_energy
        n_ax=Structure.from_file(AX_path +"POSCAR").composition.as_dict()[x_atom]
        AX_energy = AX_energy /n_ax

        BX2_path =  refDIR+"BX2_halides/" + b_atom + x_atom + "2/PBE/"
        BX2_energy = Outcar(BX2_path + "OUTCAR.gz").final_fr_energy
        n_bx2 = Structure.from_file(BX2_path +"POSCAR").composition.as_dict()[b_atom]
        BX2_energy = BX2_energy/n_bx2

        A_path =refDIR + "Elemental_std/" +a_spec+"/PBE/"
        A_energy = Outcar(A_path+"OUTCAR.gz").final_fr_energy
        n_ABULK = Structure.from_file(A_path +"POSCAR").composition.as_dict()[a_spec]
        A_energy = A_energy/n_ABULK
     
        B_path =refDIR + "Elemental_std/" +b_atom+"/PBE/"
        B_energy = Outcar(B_path+"OUTCAR.gz").final_fr_energy
        n_BBULK = Structure.from_file(B_path +"POSCAR").composition.as_dict()[b_atom]
        B_energy = B_energy/n_BBULK

        X_path =refDIR + "Elemental_std/" +x_atom+"/PBE/"
        X_energy = Outcar(X_path+"OUTCAR.gz").final_fr_energy
        n_XBULK = Structure.from_file(X_path +"POSCAR").composition.as_dict()[x_atom]
        X_energy = X_energy/n_XBULK
        print("Cs"+b_atom+x_atom+"3")
        delhax = AX_energy - A_energy - X_energy 
        delhbx2 = BX2_energy - B_energy - 2*X_energy
        delhabx3 = bulk_energy - A_energy - B_energy -3*X_energy
        decomp = bulk_energy - AX_energy - BX2_energy
        x_dcomp.append(decomp)
        print(decomp)
    decomps.append(x_dcomp)

fig, ax = plt.subplots()
im = ax.imshow(decomps, cmap=mpl.colormaps['plasma'])

ax.set_xticks(np.arange(len(b_spec)), labels=b_spec, size=20)
ax.set_yticks(np.arange(len(x_spec)), labels=x_spec, size=20)
bar = plt.colorbar(im)

bar.set_label("$\Delta H_{rxn}$(eV)", loc="center", size=30)



plt.xlabel("B Site Species", fontsize=30)
plt.ylabel("X Site Species", fontsize=30)
fig.tight_layout()
#ax.set_title("Decomposition Energy", fontsize=35)
plt.show()

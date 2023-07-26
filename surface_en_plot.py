import numpy as np 
from pymatgen.io.vasp.outputs import Outcar
from pymatgen.core.structure import Structure
import matplotlib.pyplot as plt



a_spec = "Cs"
b_spec = ["Pb", "Sn", "Ge"]
x_spec = ["Cl", "I", "Br"]
terms = ["100Term1", "100Term2", "110Term1", "110Term2", "111Term1", "111Term2"]
refDIR ="/depot/amannodi/data/MCHP_Database/Reference_compounds/"
i_comp = 0
text_size = 25
for b_atom in b_spec:
    for x_atom in x_spec:
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
        
        delhax = AX_energy - A_energy - X_energy 
        delhbx2 = BX2_energy - B_energy - 2*X_energy
        delhabx3 = bulk_energy - A_energy - B_energy -3*X_energy
        if i_comp == 1:
           axis = plt.subplot(1,9,1)
           axis.set_ylabel("Surface Energy(eV)", fontdict={'size':text_size*0.85})
           axis.yaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))

        else:
           plt.subplot(1,9,i_comp, sharey=axis)
           plt.tick_params('y', labelleft=False)
        print("Cs" +b_atom+x_atom+"3")

        for term in terms:
            mainPath = "runDIR/" + a_spec + b_atom  +  x_atom  + "3/"
            struct = Structure.from_file(mainPath + term +"/CONTCAR")
            slab_En = Outcar(mainPath+term+"/OUTCAR").final_fr_energy
            comp = struct.composition.as_dict()
            notFound  = True
            n_bulk=0
            while notFound:
                temp_a = comp[a_spec]-1
                temp_b = comp[b_atom]-1
                temp_x = comp[x_atom]-3
                if temp_a<0 or temp_b<0 or temp_x<0:  
                    notFound = False
                else:
                    comp[a_spec] = temp_a
                    comp[b_atom] = temp_b
                    comp[x_atom] = temp_x
                    n_bulk = n_bulk + 1

            a = np.array(struct.as_dict()['lattice']['matrix'][0])
            b = np.array(struct.as_dict()['lattice']['matrix'][1])
            area = np.linalg.norm(np.cross(a, b))
           # X rich plot(d_mu_X=0)
            #begin = delhabx3-delhax
            #end =  delhbx2
            #d_mu_B  = np.linspace(begin, end, num=2)
          
            #surface_energy_x_rich = [0.5/area*(slab_En - n_bulk*bulk_energy - comp[a_spec]*(delhabx3 - dB +A_energy) -  comp[b_atom]*(dB + B_energy)) for dB in d_mu_B]
            #surface_energy =surface_energy_x_rich

            #B rich plot(d_mu_B=0)
            begin = (delhabx3-delhax)/2
            end = delhbx2/2
            d_mu_X  = np.linspace(begin, end, num=2)
            surface_energy_B_rich = [0.5/area*(slab_En - n_bulk*bulk_energy - comp[a_spec]*(delhabx3 - 2*dX +A_energy) -  comp[x_atom]*(dX + X_energy)) for dX in d_mu_X]
            surface_energy= surface_energy_B_rich

            term_title = r'(%s) $ A_{%i} B_{%i} X_{%i}$' % (term[:3], int(comp[a_spec]), int(comp[b_atom]), int(comp[x_atom]))
            print(term_title)
            plt.plot(d_mu_X, surface_energy, label = term_title)
            plt.xticks(np.around([begin,end],decimals=2))
        plt.title(r'$%s %s %s_{3}$' % (a_spec, b_atom, x_atom), fontsize=text_size*0.65)
        plt.suptitle("B Rich Growth Conditions", fontsize=text_size)
        plt.xlabel(r'$\Delta\mu_{'+x_atom+'}$(eV)', fontsize=text_size*0.55)
plt.rcParams.update({'font.size': text_size*0.65})
plt.legend(bbox_to_anchor=(1.05, 0.5), loc="center left", title="Surface Termination")
plt.show()
#get buldk energy
#get excess atoms in slab


#get binary compound formation energies
#get elemental energy formation energy 

#calculate surface energies for all surfaces in allowed domain
#identify lowest surface energy at each chemical potentila point 
#color that region to signify the surface most stable there

#goal ==calculate surface phase diagram 

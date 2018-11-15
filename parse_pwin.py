import numpy as np
import f90nml

def parse_pwin_nmlist(filein):

    nml = f90nml.read(filein)

    return nml

class cards_pwin():

    def __init__(self):

        self.atomic_species={}
        self.atomic_positions={}
        self.cell_parameters=[]
        self.k_points=[]
        self.unit_pos = "alat"
        self.unit_cellp = "bohr"
        self.ktype = "tpiba"
        self.nk = 0 

    def read(self, pwin, ntyp, nat):

        with open(pwin, "r") as f:
            lines = f.readlines()

        for i in range(0, len(lines)):

            if (lines[i].strip()=="ATOMIC_SPECIES"): 
                
                k=0                
                for j in range(i+1, len(lines)):
                    if not lines[j].strip():
                        continue
                    else:
                        k+=1
                        ll = lines[j].split() 
                        self.atomic_species[ll[0]] = {"mass":ll[1], "pseudopot":ll[2]}
                        if (k==ntyp):
                            break

            if ("ATOMIC_POSITIONS" in lines[i]):

                ap_optionlist = ["alat", "bohr", "angstrom", "crystal_sg", "crystal"]
                apu = filter(lambda x: x in lines[i], ap_optionlist)
                if not apu:
                    self.unit_pos = "alat"
                else:
                    self.unit_pos = apu[0]

                k=0
                for j in range(i+1, len(lines)):
                    if not lines[j].strip():
                        continue
                    else:
                        k+=1
                        ll = lines[j].split()
                        self.atomic_positions[k] = {"element":ll[0], "position":np.asarray([ll[1], ll[2], ll[3]],dtype="float")}
                        if (k==nat):
                            break

            if ("CELL_PARAMETERS" in lines[i]):

                cp_optionlist = ["alat", "bohr", "angstrom"]
                cpu = filter(lambda x: x in lines[i], cp_optionlist) 
                if not cpu:
                    self.unit_cellp = "bohr"
                else:
                    self.unit_cellp = cpu[0]

                k=0
                for j in range(i+1, len(lines)):
                    if not lines[j].strip():
                        continue
                    else:
                        k+=1
                        self.cell_parameters.append(np.asarray(lines[j].split()))
                        if (k==3):
                            break

            if ("K_POINTS" in lines[i]):

                k_optionlist = ["automatic", "gamma", "tpiba_b" , "tpiba_c", "tpiba", "crystal_b", "crystal_c", "crystal"]
                kt = filter(lambda x: x in lines[i], k_optionlist)
                if not kt:
                    self.ktype = "tpiba"
                else:
                    self.ktype = kt[0]

                if (self.ktype == "automatic"):
                    for j in range(i+1, len(lines)):
                        if not lines[j].strip():
                            continue
                        else: 
                            self.k_points = np.asarray(lines[j].split(),dtype="int") 
                            break
                elif (self.ktype != "gamma") and (self.ktype != "automatic"):
                    for j in range(i+1, len(lines)):
                        if not lines[j].strip():
                            continue
                        else:
                            self.nk = int(lines[j].split()[0])
                            l_nk = j
                            break
                    
                    self.k_points = [] 
                    k=0
                    for j in range(l_nk+1, len(lines)):
                        if not lines[j].strip():
                            continue
                        else:
                            k+=1
                            self.k_points.append(np.asarray(lines[j].split(),dtype="float"))
                            if (k==self.nk):
                                break

#            if ("OCCUPATIONS" in lines[i]):


#            if ("CONSTRAINTS" in lines[i]):


#            if ("ATOMIC_FORCES" in lines[i]):



        return 

    def write(self,fout):
        """
        Writes the cards in the output file fout
        """
         
        return 

if __name__ == '__main__':
    cd = cards_pwin()
    cd.read('scf.in',1,2)
    print cd.atomic_species
    print cd.atomic_positions
    print cd.cell_parameters
    print cd.nk
    print cd.k_points

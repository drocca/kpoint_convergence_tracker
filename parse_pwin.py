import numpy as np
import f90nml
import sys

def parse_pwin_nmlist(filein):
    """
    This function uses the f90nml library to parse the namelists in the input of pw.x
    For example, nml["system"]["ecutwfc"] will contain the cutoff
    """
    nml = f90nml.read(filein)

    return nml

class cards_pwin():
    """ Class defined to store the information contained in the cards of the pw.x input """

    def __init__(self):

        # Example: .atomic_species["Si"]["mass"] will give the atomic mass of Si in the file
        # .atomic_species["Si"]["pseudopot"] the pseudopotential
        self.atomic_species={} 

        # Example: .atomic_positions[2]["element"] will give the type of the second element in the list of atoms
        # .atomic_positions[2]["position"] its position (numpy array)
        self.atomic_positions={}

        # If cell_parameters is present it will contain a list with 3 numpy arrays (the 3 cell parameters)
        self.cell_parameters=[]

        # If K_POINTS automatic it stores a 6 component 1D array with grid and shift 
        self.k_points_auto=np.asarray([1, 1, 1, 0, 0, 0], dtype="int")

        # If K_POINTS is not automatic or gamma
        # this list will contain nk 4 component 1D arrays with k-points and weights 
        self.k_points_list=[]

        self.unit_pos = "alat"
        self.unit_cellp = "bohr"
        self.ktype = "tpiba"
        self.nk = 0 

    def read(self, pwin, ntyp, nat):
        """ 
        Reads from pwin and set up the attributes of the class cards_pwin from the input 
        pwin: pw.x input file
        ntyp: Number of types of atoms
        nat: Number of atoms    
        """

        # Open pw.x input file
        with open(pwin, "r") as f:
            lines_tot = f.readlines()

        # Eliminates empty lines
        lines= []
        for i in range(0, len(lines_tot)):
            if not lines_tot[i].strip():
                continue
            else:
                lines.append(lines_tot[i].strip())
            

        for i in range(0, len(lines)):

            if ("ATOMIC_SPECIES" in lines[i]): 
                
                # Creates a dictionary with atomic species
                # Example: {"Si":{"mass": "28.0855", "pseudopot": "Si.UPF"}}
                for j in range(i+1, i+1+ntyp):
                    ll = lines[j].split() 
                    self.atomic_species[ll[0]] = {"mass":ll[1], "pseudopot":ll[2]}

            if ("ATOMIC_POSITIONS" in lines[i]):

                # Determines the units; if none is provided it uses alat 
                ap_optionlist = ["alat", "bohr", "angstrom", "crystal_sg", "crystal"]
                apu = filter(lambda x: x in lines[i], ap_optionlist)
                if not apu:
                    self.unit_pos = "alat"
                else:
                    self.unit_pos = apu[0]

                # Creates a dictionary with atomic positions and type
                # For example for the second atom in the list we could have {2: "element":"Si", "position":[0.25 0.25 0.25]}
                for counter, j in enumerate(range(i+1, i+1+nat)):
                    ll = lines[j].split()
                    self.atomic_positions[counter+1] = {"element":ll[0], "position":np.asarray([ll[1], ll[2], ll[3]],dtype="float")}

            if ("CELL_PARAMETERS" in lines[i]):

                # Determines the units
                cp_optionlist = ["alat", "bohr", "angstrom"]
                cpu = filter(lambda x: x in lines[i], cp_optionlist) 
                if not cpu:
                    self.unit_cellp = "bohr"
                else:
                    self.unit_cellp = cpu[0]

                # Set up a list with 3 arrays corresponding to the cell parameters
                for j in range(i+1, i+4):
                    self.cell_parameters.append(np.asarray(lines[j].split(),dtype="float"))

            if ("K_POINTS" in lines[i]):

                # Determines if k-points is automatic, gamma, or explicit list of k-points 
                # If list of k-points, it determines the units
                k_optionlist = ["automatic", "gamma", "tpiba_b" , "tpiba_c", "tpiba", "crystal_b", "crystal_c", "crystal"]
                kt = filter(lambda x: x in lines[i], k_optionlist)
                if not kt:
                    self.ktype = "tpiba"
                else:
                    self.ktype = kt[0]

                # If k-points is automatic, the grid and shift are stored
                if (self.ktype == "automatic"):

                    self.k_points_auto = np.asarray(lines[i+1].split(),dtype="int") 

                # If not automatic or gamma the full list of k-points is stored 
                elif (self.ktype != "gamma") and (self.ktype != "automatic"):
                    self.nk = int(lines[i+1].split()[0])
                    
                    k=0
                    for j in range(i+2, i+2+self.nk):
                        self.k_points_list.append(np.asarray(lines[j].split(),dtype="float"))

            # These cards are not (yet) implemented
            if ("OCCUPATIONS" in lines[i]):
                sys.exit("Parsing of OCCUPATIONS not implemented")

            if ("CONSTRAINTS" in lines[i]):
                sys.exit("Parsing of CONSTRAINTS not implemented")

            if ("ATOMIC_FORCES" in lines[i]):
                sys.exit("Parsing of ATOMIC_FORCES not implemented")


        return 

    def write(self,fout):
        """
        Writes the cards in a new input file called fout
        """

        fout.write("ATOMIC_SPECIES"+"\n")
        for key in self.atomic_species:
            fout.write(key+" "+self.atomic_species[key]["mass"]+" "+self.atomic_species[key]["pseudopot"]+"\n")

        fout.write("ATOMIC_POSITIONS"+" "+"{"+self.unit_pos+"}"+"\n")
        for i in range(1, len(self.atomic_positions)+1): 
            pos = self.atomic_positions[i]["position"]
            fout.write(self.atomic_positions[i]["element"]+" "+('%12.8f %12.8f %12.8f'%tuple(pos))+"\n")

        if len(self.cell_parameters)>0: 
            fout.write("CELL_PARAMETERS"+" "+"{"+self.unit_cellp+"}"+"\n")
            fout.write(('%12.8f %12.8f %12.8f'%tuple(self.cell_parameters[0])).strip()+"\n")
            fout.write(('%12.8f %12.8f %12.8f'%tuple(self.cell_parameters[1])).strip()+"\n")
            fout.write(('%12.8f %12.8f %12.8f'%tuple(self.cell_parameters[2])).strip()+"\n")

        fout.write("K_POINTS"+" "+"{"+self.ktype+"}"+"\n")
        
        if (self.ktype == "automatic"):
            grd = [str(x) for x in self.k_points_auto]
            fout.write(grd[0]+" "+grd[1]+" "+grd[2]+" "+grd[3]+" "+grd[4]+" "+grd[5])
        elif (self.ktype != "gamma") and (self.ktype != "automatic"):
            fout.write(str(self.nk)+"\n")
            for i in range(0, self.nk):  
                fout.write(('%12.8f %12.8f %12.8f %12.8f'%tuple(self.k_points_list[i])).strip()+"\n")
       
        return

    def set_kgrid_auto(self,grid_new):
        """ 
        Set up the k-points to automatic 
        grid_new: numpy array with the size of the grid and the shift
        Example: [2 2 2 0 0 0] 
        """
        self.ktype = "automatic"
        self.k_points_auto = grid_new         
        return 

class pwout_scf():
    """
    Class defined to store some information from the output of pw.x
    (a small subset of the information contained in this file) 
    """

    def __init__(self):

        # total energy
        self.tot_energy = 0.0

        # number of electrons; it's useful because the tot energy is an extensive quantity
        self.n_elec = 0

        # Forces
        self.forces = []

        # Stress tensor
        self.stress = []

        # Reciprocal lattice axes (useful to study k-point convergence) 
        self.rec_axes = []
        
    def read(self, pwout, nat):
        """
        Reads from pwout and setup the attributes of the class pwout_scf 
        pwout: output of a scf calculation of pw.x
        nat: number of atoms
        """

        # Open pw.x output file (scf calculation)
        with open(pwout, "r") as f:
            lines_tot = f.readlines()

        # Eliminates empty lines
        lines= []
        for i in range(0, len(lines_tot)):
            if not lines_tot[i].strip():
                continue
            else:
                lines.append(lines_tot[i].strip())


        for i in range(0, len(lines)):

            if (lines[i][0] == "!"):
                ll = lines[i].split() 
                self.tot_energy = np.float(ll[4])

            if ("number of electrons" in lines[i]):
                ll = lines[i].split()
                self.n_elec = np.float(ll[4]) 

            if ("Forces" in lines[i]):
                for counter, j in enumerate(range(i+1, i+1+nat)):
                    ll = lines[j].split()
                    self.forces.extend([ll[6], ll[7], ll[8]])

            if ("total   stress" in lines[i]):
                for j in range(i+1, i+4):
                    ll = lines[j].split()
                    self.stress.extend([ll[0], ll[1], ll[2]])

            if ("reciprocal axes" in lines[i]):
                for j in range(i+1, i+4):
                    ll = lines[j].split()
                    self.rec_axes.append(np.asarray([ll[3], ll[4], ll[5]],dtype='float'))


if __name__ == '__main__':

    # This is just for a quick testing 

    # pw.x input file
    name_pwin = 'scf.in'

    # parsing namelists
    pw_nml = parse_pwin_nmlist(name_pwin)

    # parsing cards
    n_typ = pw_nml["system"]["ntyp"]
    n_at = pw_nml["system"]["nat"]

    pw_cards = cards_pwin()
    pw_cards.read(name_pwin,n_typ,n_at) 
 
    # writing back to a file
    with open("scf_parsed.in", 'w') as fout:
        pw_nml.write(fout)
        pw_cards.write(fout) 

    # pw.x output file of a scf calculation
    name_scfout = 'scf.out'
 
    # parsing output
    scfout = pwout_scf() 
    scfout.read(name_scfout,n_at) 

    print scfout.tot_energy
    print scfout.n_elec
    print scfout.forces
    print scfout.stress 
    print scfout.rec_axes

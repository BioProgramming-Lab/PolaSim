import numpy as np


kon_A = 1.1
kon_B = 1.1

koff_A = 0.1
koff_B = 0.1

koff_A2B = 1.0
koff_B2A = 1.0

function_radius_A2B = 1.0
function_radius_B2A = 1.0

diffusion_coefficient_A = 0.01
diffusion_coefficient_B = 0.01

production_rate_A = 1.0
production_rate_B = 1.0

degradation_rate_A = 1.0
degradation_rate_B = 1.0

ellipsoid_axis_a = 3.0
ellipsoid_axis_b = 1.0
ellipsoid_axis_c = 1.0

class molecule:
    def __init__(self, type:str, phi:float, theta:float, OnMembrane:bool):
        self.type = type
        self.phi = phi
        self.theta = theta
        self.OnMembrane = OnMembrane
    
    def On_membrane_diffuse(self, dt:float):
        if self.OnMembrane:
            # Calculate the diffusion distance using a Gaussian distribution
            diffusion_distance = np.random.normal(loc=0.0, scale=np.sqrt(2 * diffusion_coefficient_A * dt))

            # Convert spherical coordinates to Cartesian coordinates
            x = ellipsoid_axis_a * np.sin(self.theta) * np.cos(self.phi)
            y = ellipsoid_axis_b * np.sin(self.theta) * np.sin(self.phi)
            z = ellipsoid_axis_c * np.cos(self.theta)

            # Apply the diffusion distance to the Cartesian coordinates
            x += diffusion_distance * np.random.normal()
            y += diffusion_distance * np.random.normal()
            z += diffusion_distance * np.random.normal()

            # Convert back to spherical coordinates
            r = np.sqrt((x / ellipsoid_axis_a)**2 + (y / ellipsoid_axis_b)**2 + (z / ellipsoid_axis_c)**2)
            self.theta = np.arccos(z / (r * ellipsoid_axis_c))
            self.phi = np.arctan2(y / ellipsoid_axis_b, x / ellipsoid_axis_a)
        else:
            pass

class cell:
    def __init__(self, size_a:float, size_b:float, size_c:float, num_init_A_molecule:int, num_init_B_molecule:int):
        self.size_a = size_a
        self.size_b = size_b
        self.size_c = size_c
        self.num_cyto_A_molecule = num_init_A_molecule
        self.num_cyto_B_molecule = num_init_B_molecule
        self.mem_molecules = []
    
    def cyto2mem(self, number:int, type:str):
        for i in range(number):
            phi = np.random.uniform(0, 2*np.pi)
            theta = np.random.uniform(0, np.pi)
            self.mem_molecules.append(molecule(type, phi, theta, True))
        if type == 'A':
            self.num_cyto_A_molecule -= number
        else:
            self.num_cyto_B_molecule -= number

    def mem2cyto(self, index:int):
        if self.mem_molecules[index].type == 'A':
            self.num_cyto_A_molecule += 1
        else:
            self.num_cyto_B_molecule += 1
        del self.mem_molecules[index]

    def koff_cal(self, index:int):
        basal_koff = 0.0
        if self.mem_molecules[index].type == 'A':
            basal_koff = koff_A
        else:
            basal_koff = koff_B
        for i in range(len(self.mem_molecules)):
            if i != index and self.mem_molecules[i].type != self.mem_molecules[index].type:
                # Convert spherical coordinates to Cartesian coordinates for both molecules
                x1 = ellipsoid_axis_a * np.sin(self.mem_molecules[index].theta) * np.cos(self.mem_molecules[index].phi)
                y1 = ellipsoid_axis_b * np.sin(self.mem_molecules[index].theta) * np.sin(self.mem_molecules[index].phi)
                z1 = ellipsoid_axis_c * np.cos(self.mem_molecules[index].theta)

                x2 = ellipsoid_axis_a * np.sin(self.mem_molecules[i].theta) * np.cos(self.mem_molecules[i].phi)
                y2 = ellipsoid_axis_b * np.sin(self.mem_molecules[i].theta) * np.sin(self.mem_molecules[i].phi)
                z2 = ellipsoid_axis_c * np.cos(self.mem_molecules[i].theta)

                # Calculate the Euclidean distance between the two points in Cartesian coordinates
                distance = np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
                if self.mem_molecules[i].type == 'A':
                    if distance < function_radius_A2B:
                        basal_koff += koff_A2B
                else:
                    if distance < function_radius_B2A:
                        basal_koff += koff_B2A
        return basal_koff
    
    def update_status(self, dt:float):
        for i in range(len(self.mem_molecules)-1, -1, -1):
            if np.random.uniform() < self.koff_cal(i) * dt:
                self.mem2cyto(i)
            else:
                self.mem_molecules[i].On_membrane_diffuse(dt)
        number_A2mem = np.random.poisson(kon_A * self.num_cyto_A_molecule * dt)
        number_B2mem = np.random.poisson(kon_B * self.num_cyto_B_molecule * dt)
        if number_A2mem > self.num_cyto_A_molecule:
            number_A2mem = self.num_cyto_A_molecule
        if number_B2mem > self.num_cyto_B_molecule:
            number_B2mem = self.num_cyto_B_molecule
        self.cyto2mem(number_A2mem, 'A')
        self.cyto2mem(number_B2mem, 'B')
        # print("mem ", len(self.mem_molecules))
        # print ("cyto ", self.num_cyto_A_molecule, self.num_cyto_B_molecule)
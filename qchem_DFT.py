import numpy as np

'''
A simple script for the extraction of unperturbed properties
from a DFT calculation made with QChem software.
Here the input requests in QChem .inp file:
$rem
METHOD           CAM-B3LYP
CIS_N_ROOTS      5
CIS_SINGLETS     true
CIS_TRIPLETS     false
SYMMETRY         FALSE
SYM_IGNORE       TRUE
RPA              FALSE
BASIS            6-311++G(2d,2p)
STS_MOM          true
ESP_CHARGES      2
STATE_ANALYSIS true
MOLDEN_FORMAT true
NTO_PAIRS 3
ESP_GRID -3
$end
'''


'''
Insert below:
- filename  = .log filename
- n_at      = nubmer of atoms
- n_es      = number of excited states
'''

filename = 'gua_n_CAM_esp_norelax.log'

n_at = 16
n_es = 5
n_tr = int( (n_es)*(n_es - 1) / 2 )


'''
Key lines to extract properties
'''
gs_energy_line = 'Total energy in the final basis set'
es_energy_line = 'Total energy for state'
energies = []

gs_perm_dipole_line = 'Electron Dipole Moments of Ground State'
es_perm_dipole_line = 'Electron Dipole Moments of Singlet Excited State'
gs_es_dipole_line =  'Transition Moments Between Ground and Singlet Excited States'
es_es_dipole_line = 'Transition Moments Between Singlet Excited States'

list_mom = []
dipole_moms = np.zeros((n_es + 1, n_es + 1,3))

gs_esp_charges = 'Merz-Kollman ESP Net Atomic Charges'
es_esp_charges = 'ESP charges for excited states'

gs_esp = []
esp = np.zeros((n_at, n_es + 1))


'''
file opening and extraction
'''

with open(filename, 'r') as log:

    for line in log:
        
        # **** Energy section ****

        if gs_energy_line in line:

            # print(line)
            
            energies.append(float(line.split()[8]))
            
            # print(energies)

        if es_energy_line in line:

            # print(line)

            energies.append(float(line.split()[5]))
            
            # print(energies)


        # **** Dipole moments section ****

        if gs_perm_dipole_line in line:

            # print(line)
            
            for i in range(3):
                next(log)    
            
            list_mom = next(log).split()

            dipole_moms[0, 0, :] = [float(i) for i in list_mom[1:4]]

        if es_perm_dipole_line in line:

            # print(line)

            for i in range(3):
                next(log)
            
            for i in range(n_es):
                
                list_mom = next(log).split()

                dipole_moms[i+1, i+1, :] = [float(i) for i in list_mom[1:4]]
        
        if gs_es_dipole_line in line:

            # print(line)

            for i in range(3):
                next(log)
            
            for i in range(n_es):
                
                list_mom = next(log).split()

                dipole_moms[0, i+1, :] = \
                    dipole_moms[i+1, 0, :] = [float(i) for i in list_mom[2:5]]
            
        if es_es_dipole_line in line:

            # print(line)

            for i in range(3):
                next(log)
            
            for i in range(n_tr):
                
                list_mom = next(log).split()

                h, l = int(list_mom[0]), int(list_mom[1])

                dipole_moms[h, l, :] = \
                    dipole_moms[l, h, :] = [float(j) for j in list_mom[2:5]]

        
        # **** ESP charges section ****

        if gs_esp_charges in line:
            
            for i in range(3):
                next(log)

            for i in range(n_at):
                
                gs_esp.append(next(log).split()[2])

            # print(gs_esp)

            esp[:, 0] = [float(i) for i in gs_esp]

        if es_esp_charges in line:

            next(log)

            for i in range(n_at):
                
                n_at_esp = next(log).split()
                
                # print(n_at_esp)

                esp[i, 1:] = [float(j) for j in n_at_esp[1:]]
                
            esp = esp.T

            # print(esp)

log.close()

print('------------------------------------------------------------------- \n')
print('Energy of states: \n', energies)
print('------------------------------------------------------------------- \n')
print('Electric dipole moment matrix: \n', dipole_moms)
print('------------------------------------------------------------------- \n')
print('ESP charges of states: \n', esp)
print('------------------------------------------------------------------- \n')

'''
Saving unperturbed properties to files 
'''

np.savetxt('energies.pmm', energies, fmt='%.6f')

with open('dipole_matrix.pmm', 'w') as dip_out_file:

    for i in range(n_es + 1):
        for j in range(n_es + 1):
            
            dip_out_file.write('{} {} {:.6f} {:.6f} {:.6f} \n'.format(i, j, \
                dipole_moms[i, j, 0], dipole_moms[i, j, 1], dipole_moms[i, j, 2]))

dip_out_file.close()

np.savetxt('esp_matrix.pmm', esp, fmt='%.6f')
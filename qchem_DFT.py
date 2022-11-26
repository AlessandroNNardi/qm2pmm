import numpy as np


'''
Insert filename below
'''
filename = 'DFT/gua_n_CAM_esp_norelax.log'

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


with open(filename, 'r') as log:

    for line in log:

        # Energy section

        if gs_energy_line in line:

            # print(line)
            
            energies.append(float(line.split()[8]))
            print(energies)

        if es_energy_line in line:

            # print(line)

            energies.append(float(line.split()[5]))
            
            # print(energies)


        # Dipole moments section

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

log.close()

print(energies)
print(dipole_moms)
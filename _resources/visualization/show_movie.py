#!/usr/bin/python3 
"""
    Given a file containing a stream of occupation numbers,
    generate an animated GIF movie (and compress the resulting file
    by using a reduced colour palette). 
"""
import numpy as np
import matplotlib.pyplot as plt 
import re 
import imageio 
import subprocess

DIR="./"

infile_dn = DIR+'Fock_samples_ncpu00000_dn.dat' 
infile_up = DIR+'Fock_samples_ncpu00000_up.dat' 


def inputline2mat(line):
    BSS_sign = float(line.split()[0])           # type I sign (sign problem in the Hamiltonian)
    Re_phase = float(line.split()[1])           # type II sign (sampling of pseudofermion density matrix)
    Im_phase = float(line.split()[2])           # Re_phase = real_part(e^{i phase}), Im_phase = imag_part(e^{i phase})
    reweight_factor = float(line.split()[3])    
    occ_vector = np.array([int(b) for b in line.split()[4:]])

    n = len(occ_vector)
    l = int(np.sqrt(float(n)))
    M = occ_vector.reshape([l,l])

    return BSS_sign, Re_phase, Im_phase, reweight_factor, M

# create animation from image list 
image_list = []

counter = 0
skip = 0
max_img = 120

print("Combining snapshots for spin up and down.")
print("No. of snapshots:", max_img)
print("Please wait...")

with open(infile_up, 'r') as fh_up:
    with open(infile_dn, 'r') as fh_dn:
        while (True):
            if (counter < max_img):
                counter += 1 
            else:
                break
            print("image counter=", counter)

            line_up = fh_up.readline(); line_dn = fh_dn.readline()
            if ((line_up == None) or (line_dn == None)):
                break
            # skip comments
            if ((not re.match(r'^#', line_up)) and (not re.match(r'^#', line_dn))):
                BSS_sign_up, Re_phase_up, Im_phase_up, weight_up, M_up = inputline2mat(line_up)
                BSS_sign_dn, Re_phase_dn, Im_phase_dn, weight_dn, M_dn = inputline2mat(line_dn)

                mask_double = (M_up == 1) & (M_dn == 1)
                mask_hole = (M_up == 0) & (M_dn == 0)        
                mask_single_up = (M_up == 1) & (M_dn == 0)
                mask_single_dn = (M_up == 0) & (M_dn == 1)

                # Comment: It would be better to represent the state space on a site 
                #    by single characters. Unfortunately, the broadcasting operation 
                #    which replaces each entry of the occupation vector by the three RBG channels
                #    works only for numpy arrays, not for dictionaries. 
                #
                # SINGLE_UP = 'u'
                # SINGLE_DN = 'd'
                # DOUBLE = 'D'
                # HOLE = 'h'

                SINGLE_UP = 0
                SINGLE_DN = 1
                DOUBLE = 2
                HOLE = 3

                M = np.array(M_up, dtype=np.int8)
                M[mask_single_up] = SINGLE_UP
                M[mask_single_dn] = SINGLE_DN
                M[mask_double] = DOUBLE
                M[mask_hole] = HOLE

                # See comment above. 
                # palette = {'u': np.array([0.0, 0.1, 0.9], dtype=np.float32), # single up spin
                #            'd': np.array([0.0, 0.4, 0.9], dtype=np.float32), # single down spin
                #            'D': np.array([1.0, 1.0, 1.0], dtype=np.float32), # doubly occupied site
                #            'h': np.array([0.0, 0.0, 0.0], dtype=np.float32)  # empty site 
                #             }

                palette = np.array([ np.array([0.0, 0.1, 0.9], dtype=np.float32), # single up spin
                                     np.array([0.0, 0.4, 0.8], dtype=np.float32), # single down spin
                                     np.array([1.0, 1.0, 1.0], dtype=np.float32), # doubly occupied site
                                     np.array([0.0, 0.0, 0.0], dtype=np.float32) ])  # empty site 
                            

                RGB_img = palette[M]
                
                if (counter > skip):
                    fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(21, 7))
                    ax1.set_title("SPIN UP \n sampling sign = %d, reweighting = %f" % (Re_phase_up, weight_up))
                    ax2.set_title("SPIN DOWN \n sampling sign = %d, reweighting = %f" % (Re_phase_dn, weight_dn))
                    ax3.set_title("Combined (UP, DOWN)\n sign = %d, reweighting = %f\nBSS_sign = %d"\
                         % (Re_phase_up*Re_phase_dn, weight_up*weight_dn, BSS_sign_up*BSS_sign_dn))
                    ax1.imshow(M_up, origin="lower")
                    ax2.imshow(M_dn, origin="lower")
                    ax3.imshow(RGB_img, origin="lower")
                    plt.show()



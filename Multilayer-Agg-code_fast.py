import os
import numpy as np
import random as ran
import time
import cProfile
import multiprocessing as mp
from concurrent.futures import ThreadPoolExecutor
from numba import njit


def get_user_input():
    while True:  # choose monomer size
        try:
            r_mm = float(input('Enter the size of monomers: '))
        except ValueError:
            # Not a valid number
            print('You must enter a valid number')
            r_mm = float(input('Enter the size of monomers: '))
        else:
            # No error; stop the loop
            break

    while True:  # choose layers
        try:
            l = int(input('number of layers: '))
        except ValueError:
            # Not a valid number
            print('You must enter a valid number')
            l = int(input('number of layers: '))
        else:
            # No error; stop the loop
            break

    layers = np.empty(l, dtype=object)

    for i in range(l):
        n = int(input('how many particles in layer' + ' ' + str(i) + ':?'))
        layers[i] = [0] * n

    return r_mm, l, layers


def initialize_output_file(n_m, l, r_c, Re, Im):
    Wavelength = 2450 / 1000
    Re = 3.4080000000000004
    Im = 0.02024772462765312
    
    save_path = '/Users/raomorusupalli/Documents/UNI/Honours/project/GMM/'
    name_of_file = 'aggregate_NIR' + str(n_m) + '_' + str(l)
    completeName = os.path.join(save_path, name_of_file + ".k")

    with open(completeName, 'w') as f:
        f.write(str(Wavelength) + '\n')
        f.write(str(n_m + 1) + '\n')
        f.write(" ".join([str(0.), str(0.), str(0.), str(r_c), str(Re), str(Im), '\n']))

    return completeName

@njit
def get_random_position(i, r_c, r_m):
    pos = np.random.rand(3)
    pos[0] = 2.0 * pos[0] - 1.0
    r = np.sqrt(1 - pos[0] ** 2)

    pos[1] = 2.0 * np.pi * pos[1]
    pos[2] = r * np.cos(pos[1])
    pos[1] = r * np.sin(pos[1])

    fac = (r_c + (i + 1) * r_m)
    pos *= fac

    return pos

@njit
def check_overlap(posvec, actp, r_m, idx):
    for i in range(idx):
        dist = np.linalg.norm(actp[i] - posvec)
        if dist < 2 * r_m:
            return True
    return False

def create_particles_in_layer(args):
    layer_index, layer, r_c, r_m = args
    particles = np.empty((len(layer), 3), dtype=np.float64)
    idx = 0

    while idx < len(layer):
        posvec = get_random_position(layer_index, r_c, r_m)
        overlap = check_overlap(posvec, particles, r_m, idx)

        if not overlap:
            particles[idx] = posvec
            idx += 1

    return particles

if __name__ == '__main__':

    def create_dust_particles(r_mm, l, layers, Re, Im):
        r_c = 1
        r_m = r_mm / 1e3
        n_m = int(sum(len(x) for x in layers))

        completeName = initialize_output_file(n_m, l, r_c, Re, Im)

        with ThreadPoolExecutor() as executor:
            args = [(i, layers[i], r_c, r_m) for i in range(l)]
            layer_particles = list(executor.map(create_particles_in_layer, args))
        actp = [particle for particles in layer_particles for particle in particles]


        with open(completeName, 'a') as f:
            for val in actp:
                f.write('{} {} {}'.format(val[0], val[1], val[2]) + ' ' + str(r_m) + ' ' + str(Re) + ' ' + str(Im) + '\n')
            f.close()

        return actp

def main():
    r_mm, l, layers = get_user_input()
    Re = 3.4080000000000004
    Im = 0.02024772462765312
    start_time = time.time()
    create_dust_particles(r_mm, l, layers, Re, Im)
    print("--- %s seconds ---" % (time.time() - start_time))
    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    cProfile.run('main()')

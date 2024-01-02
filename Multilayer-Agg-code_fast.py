import os
import numpy as np
import random as ran
import time
import sys
import cProfile
import multiprocessing as mp
from concurrent.futures import ThreadPoolExecutor
from numba import njit, prange

def get_user_input():
    n_m = int(sys.argv[1])
    r_m = int(sys.argv[2])
    total_monomers = int(sys.argv[3])

    monomers_per_layer = total_monomers // r_m
    remaining_monomers = total_monomers % r_m

    layers = np.empty(r_m, dtype=object)
    for i in range(r_m):
        num_monomers = monomers_per_layer + (1 if i < remaining_monomers else 0)
        layers[i] = [0] * num_monomers
#remaining monomers are added to the first layer something to potentially think about ?
    return n_m, r_m, layers

#def get_wavelength():
#    while True:  # ask for wavelength
#        try:
#            w = input('Enter the wavelength in nm (default 870): ')
#            if w == '':
#                w = 870  # default value if user just presses Enter
#            else:
#                w = float(w)
#        except ValueError:
#            print('You must enter a valid number or press Enter for the default value.')
#        else:
#            break
#    return w / 1000  # return wavelength in um

#def get_user_input():
#    while True:  # choose monomer size
#        try:
#            r_mm = float(input('Enter the size of monomers: '))
#        except ValueError:
#            # Not a valid number
#            print('You must enter a valid number')
#            r_mm = float(input('Enter the size of monomers: '))
#        else:
#            # No error; stop the loop
#            break
#
#    while True:  # choose layers
#        try:
#            l = int(input('number of layers: '))
#        except ValueError:
#            # Not a valid number
#            print('You must enter a valid number')
#            l = int(input('number of layers: '))
#        else:
#            # No error; stop the loop
#            break
#this is for the grid
#    n_m = int(sys.argv[1])
#    r_m = int(sys.argv[2])
#    total_monomers = int(sys.argv[3])
#    total_monomers = int(os.environ['num_monomers']) 
#    monomers_per_layer = total_monomers // l
#    remaining_monomers = total_monomers % l
#
#    layers = np.empty(l, dtype=object)
#    for i in range(l):
#        num_monomers = monomers_per_layer + (1 if i < remaining_monomers else 0)
#        layers[i] = [0] * num_monomers
# uncomment for running individual and not on an automated grid 
#    layers = np.empty(l, dtype=object)

#    for i in range(l):
#        n = int(input('how many particles in layer' + ' ' + str(i) + ':?'))
#        layers[i] = [0] * n
#
#    return r_mm, l, layers


def initialize_output_file(total_particles, l, r_c, Re, Im):
    Wavelength = 870/ 1000
#    Re = 3.4080000000000004
#    Im = 0.02024772462765312
# to ignore the central sphere we make it a vaccum 
    Re = 1
    Im = 0   
    save_path = os.path.expanduser('~/runs/Agg_models/input')
    name_of_file = 'aggregate_NIR' + str(total_particles) + '_' + str(l)
    completeName = os.path.join(save_path, name_of_file + ".k")

    with open(completeName, 'w') as f:
        f.write(str(Wavelength) + '\n')
        f.write(str(total_particles + 1) + '\n')
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

def calculate_spacing_factor(r_mm, l, layers):
    # You can use any formula that takes into account the input parameters
    # This is just an example; you may need to adjust the formula
    # based on your specific requirements and observations
    average_monomers = sum(len(layer) for layer in layers) / l
    return 1.0 + 0.01 * average_monomers

#@njit
#def check_overlap(posvec, actp, r_m, idx, inflate_distance):
#    for i in range(idx):
#        dist = np.linalg.norm(actp[i] - posvec)
#        if dist <  (2 * r_m * inflate_distance):
#            return True
#    return False
#taken out inflate distance 
@njit(parallel=True)
def check_overlap(posvec, actp, r_m, idx, spacing_factor):
    inflated_distance_squared = (2 * r_m * spacing_factor ) ** 2
    for i in prange(idx):
        dist_squared = np.sum((actp[i] - posvec) ** 2)
        if dist_squared < inflated_distance_squared :
            return True
    return False

def create_particles_in_layer(args, prev_layer=None, next_layer=None):
    layer_index, layer, r_c, r_m, spacing_factor = args
    #spacing_factor = 1.09
    effective_radius = r_m * spacing_factor
    layer_surface_area = surface_area_of_layer(layer_index, r_c, effective_radius)
    max_num_monomers = int(layer_surface_area // (4 * np.pi * effective_radius ** 2))
    num_mon = min(len(layer), max_num_monomers)
    
    particles = np.empty((num_mon, 3), dtype=np.float64)
    idx = 0

    while idx < num_mon:
        posvec = get_random_position(layer_index, r_c, r_m)
        inflate_distance = spacing_factor
        overlap = check_overlap(posvec, particles, r_m, idx, inflate_distance)

        if not overlap and prev_layer is not None:
            overlap = check_overlap(posvec, prev_layer, r_m, len(prev_layer), inflate_distance)

        if not overlap and next_layer is not None:
            overlap = check_overlap(posvec, next_layer, r_m, len(next_layer), inflate_distance)
      
        if not overlap:
            particles[idx] = posvec
            idx += 1
            print(f"Layer {layer_index}: Added monomer {idx}/{num_mon}")
        else:
            while overlap:
                posvec = get_random_position(layer_index, r_c, r_m)
#                inflate_distance += 0.02
                overlap = check_overlap(posvec, particles, r_m, idx, inflate_distance)

                if not overlap and prev_layer is not None:
                    overlap = check_overlap(posvec, prev_layer, r_m, len(prev_layer), inflate_distance)

                if not overlap and next_layer is not None:
                    overlap = check_overlap(posvec, next_layer, r_m, len(next_layer), inflate_distance)

    return particles, num_mon

def surface_area_of_layer(layer_index, r_c, r_m):
    layer_radius = r_c + (layer_index + 1) * r_m
    return 4 * np.pi * layer_radius ** 2


def total_monomer_surface_area(n_monomers, r_m):
    return n_monomers * 4 * np.pi * r_m ** 2


def extrapolate_refractive_index(wavelength):
    # Read in the data from the file
    data = pd.read_csv('draine_optics.dat', header=None, delimiter=r"\s+", engine='python')
    data.columns = ['Wav', 'Real', 'Img']

    f = interpolate.interp1d(np.log(data.Wav), np.log(data.Img), fill_value='extrapolate')
    g = interpolate.interp1d(np.log(data.Wav), np.log(data.Real), fill_value='extrapolate')

    im = np.exp(f(np.log(wavelength)))
    rel = np.exp(g(np.log(wavelength)))
    
    return rel, im

if __name__ == '__main__':


    def main():
        r_mm, l, layers = get_user_input()
#        Wavelength = get_wavelength()
#        Re, Im = extrapolate_refractive_index(Wavelength)
        Re = 3.4080000000000004
        Im = 0.02024772462765312
        r_c = 1
        r_m = r_mm / 1e3
        n_m = int(sum(len(x) for x in layers)) 
        spacing_factor = calculate_spacing_factor(r_mm, l, layers)
        spacing_factor = 1.09 
        start_time = time.time()
        with ThreadPoolExecutor() as executor:
            args = [(i, layers[i], r_c, r_m, spacing_factor) for i in range(l)]
#            layer_particles_and_counts = list(executor.map(create_particles_in_layer, args))
            layer_particles_and_counts = []
            for i, layer_args in enumerate(args):
                prev_layer = layer_particles_and_counts[i-1][0] if i > 0 else None
                next_layer = None if i == l - 1 else layers[i + 1]
                particles, count = create_particles_in_layer(layer_args, prev_layer, next_layer)
                layer_particles_and_counts.append((particles, count))
        layer_particles = [x[0] for x in layer_particles_and_counts]
        actual_counts = [x[1] for x in layer_particles_and_counts]

        for layer_index, count in enumerate(actual_counts):
            if count < len(layers[layer_index]):
               print(f"Physical limit reached. Actual number of monomers used in layer {layer_index}: {count}")
        actp = [particle for particles in layer_particles for particle in particles]
        total_particles = sum(actual_counts)
        print(f"TOTAL_PARTICLES:{total_particles}")
        completeName = initialize_output_file(total_particles, l, r_c, Re, Im)
        with open(completeName, 'a') as f:
            for val in actp:
                f.write('{} {} {}'.format(val[0], val[1], val[2]) + ' ' + str(r_m) + ' ' + str(Re) + ' ' + str(Im) + '\n')
            f.close()


        print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()

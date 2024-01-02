import os 
import sys

def initialize_output_file(monomer_size, Re=3.4080000000000004, Im=0.02024772462765312):
    Wavelength = 870 / 1000
# uncomment to creategrains of simillar volume to aggreagte    
#    r_c = monomer_size * (num_layers * spacing_factor)

    save_path = os.path.expanduser('~/runs/Agg_models/input/single_run') 
    name_of_file = 'monomer_NIR' + '_' + str(monomer_size)
    completeName = os.path.join(save_path, name_of_file + ".k")

    with open(completeName, 'w') as f:
        f.write(str(Wavelength) + '\n')
        f.write('1\n')  # Only one monomer
        f.write(" ".join([str(0.0), str(0.0), str(0.0), str(monomer_size) , str(Re), str(Im), '\n']))

def main():
    monomer_size = float(sys.argv[1])  # From command line argument
#    num_layers = int(sys.argv[2])  # From command line argument    
#    spacing_factor = 1.09
    initialize_output_file(monomer_size)

if __name__ == "__main__":
    main()


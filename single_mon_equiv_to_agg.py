import os
import sys

def initialize_output_file(monomer_size, num_layers, spacing_factor, Re=3.4080000000000004, Im=0.02024772462765312):
    Wavelength = 870 / 1000
    r_c = 1 + monomer_size/1000 * (num_layers * spacing_factor)

    save_path = os.path.expanduser('~/runs/Agg_models/input/single_run')

    # Create directory if it doesn't exist
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    # Modify this line to match the filename pattern in the Bash script
    name_of_file = 'monomer_NIR' + '_' + str(monomer_size) + '_' + str(num_layers)
    
    completeName = os.path.join(save_path, name_of_file + ".k")

    try:
        with open(completeName, 'w') as f:
            f.write(str(Wavelength) + '\n')
            f.write('1\n')  # Only one monomer
            f.write(" ".join([str(0.0), str(0.0), str(0.0), str(r_c), str(Re), str(Im), '\n']))
    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)

def main():
    if len(sys.argv) < 3:
        print("Usage: python script.py <monomer_size> <num_layers>")
        sys.exit(1)
      
    monomer_size = int(sys.argv[1])  # From command line argument
    num_layers = int(sys.argv[2])  # From command line argument
    
    spacing_factor = 1.09
    initialize_output_file(monomer_size, num_layers, spacing_factor)

if __name__ == "__main__":
    main()


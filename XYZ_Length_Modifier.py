"""
This code extracts the lines in the first N timesteps from an XYZ file (where N is specified
by user) and writes them to an output file.

The format of the XYZ file can be seen in my repository. 

Author: Henry Yip/DrOreilleNehCinq
"""

# Need to be an integer
Num_required = int(input("Number of timesteps you want to extract: "))

# Need to end in .xyz
Input_File_Name = str(input("Number of input file: "))

# Need to end in .xyz
Output_File_Name = str(input("Number of output file: "))

def extract_first_N_terms(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        term_count = 0
        while term_count < Num_required:
            try:
                # Reading the number of planets
                num_planets = infile.readline().strip()
                if not num_planets:
                    break  
                outfile.write(num_planets + '\n')

                # Reading the "Point" line (which correspond to the total number of timesteps in input file) 
                point_line = infile.readline().strip()
                outfile.write(point_line + '\n')

                # Read the Planets Coordinates (Position and Velocity)
                for _ in range(int(num_planets)):
                    atom_line = infile.readline().strip()
                    outfile.write(atom_line + '\n')

                term_count += 1
            except Exception as e:
                print(f"Error while processing the file: {e}")
                break

extract_first_N_terms(Input_File_Name, Output_File_Name)

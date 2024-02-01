#!/bin/bash

# Define the path to the executable
executable="exec/a1"

# Define the input files
input_files=("model/bunny.obj" "model/sphere.obj" "model/teapot.obj")

# Define the resolutions
resolutions=(32 64)

# Loop over each input file
for input_file in "${input_files[@]}"; do
    # Extract the base name of the file without the extension
    base_name=$(basename "$input_file" .obj)

    # Loop over each resolution
    for resolution in "${resolutions[@]}"; do
        # Construct the output file name
        output_file="out/${base_name}_${resolution}.obj"

        # Run the executable with the current input file, output file, and resolution, q?
        $executable "$input_file" "$output_file" "$resolution" "0"
    done
done

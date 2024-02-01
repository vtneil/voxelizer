#!/bin/bash

# Define the path to the executable
executables=("exec/voxelizer_bvh" "exec/voxelizer_cpu")

# Define the input files
input_file="model/bunny.obj"

# Define the resolutions
resolutions=(8 16 32 64)

# Loop over each resolution
for resolution in "${resolutions[@]}"; do
    # Loop over each input file
    for executable in "${executables[@]}"; do
        # Extract the base name of the file without the extension
        base_name=$(basename "$input_file" .obj)

        # Construct the output file name
        output_file="bench_out/${executable##*/}_${base_name}_${resolution}.obj"

        # Run the executable with the current input file, output file, and resolution, q?
        $executable "$input_file" "$output_file" "$resolution" "0"
    done
done

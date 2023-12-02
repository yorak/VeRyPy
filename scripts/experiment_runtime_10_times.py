import subprocess
import re
import numpy as np
from collections import defaultdict
from sys import argv

# Command to run the algorithm
command = argv[1]
outfile = argv[2]
print("VeRyPy command to run is:", command)
print("Store log to file:", outfile)


# Dictionaries to store runtimes and solution lengths for each algorithm
runtimes = defaultdict(list)
solution_lengths = defaultdict(list)

# Repeat the process 10 times
for i in range(10):
    # Execute the command and capture the output
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, text=True)
    output_lines = result.stdout.split('\n')

    # Save the output to a file
    with open(outfile, "a") as file:
        file.write(f"##########\n")
        file.write(f"Run {i+1}:\n")
        file.write(f"##########\n")
        file.write(result.stdout) 
        file.write("\n\n")

    for line in output_lines:
        print(line)
        # Extract the algorithm name
        algorithm_match = re.search(r'Solving khadija with (.+)', line)
        if algorithm_match:
            current_algorithm = algorithm_match.group(1)

        # Extract the runtime
        runtime_match = re.search(r'solution in ([0-9.]+) s', line)
        if runtime_match:
            runtimes[current_algorithm].append(float(runtime_match.group(1)))

        # Extract the solution length
        length_match = re.search(r'SOLUTION LENGTH: ([0-9.]+)', line)
        if length_match:
            solution_lengths[current_algorithm].append(float(length_match.group(1)))

# Calculate mean and standard deviation for each algorithm
for algorithm in runtimes:
    mean_runtime = np.mean(runtimes[algorithm])
    stddev_runtime = np.std(runtimes[algorithm])
    mean_length = np.mean(solution_lengths[algorithm])
    stddev_length = np.std(solution_lengths[algorithm])

    print(f"Algorithm: {algorithm}")
    print(f"Mean Runtime: {mean_runtime} s, StdDev Runtime: {stddev_runtime} s")
    print(f"Mean Solution Length: {mean_length}, StdDev Solution Length: {stddev_length}")
    print("\n")

# Replace the placeholder command with the actual command to execute your algorithm.


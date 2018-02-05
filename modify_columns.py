import numpy as np

# Reading from file

"""Change the filename here"""
filename = "associationruletestdata.txt"

# Temporary file for modifying the data.
temp_file = "tempfile.txt"

file = np.genfromtxt(filename, dtype=str, delimiter="\t")

file = np.matrix(file)
(rows, cols) = file.shape

# Modifying the data
index = 0;
with open(filename, 'r')  as main_file:
    with open(temp_file, "w") as write_file:
        for line in main_file:
            line = line.split()
            index = index + 1;
            # print len(line)
            for i in range(cols):
                if (i + 1 != cols):
                    line[i] = "G" + str(i + 1) + "_" + line[i]
            # print line
            length = len(line)
            if (length == cols + 1):
                disease = line[cols - 1] + " " + line[cols]
            else:
                disease = line[cols - 1]
            # print disease
            row = ""
            for i in range(cols - 1):
                row = row + line[i] + "\t"
            row = row + disease + "\n"
            write_file.write(row)

file1 = np.genfromtxt(temp_file, dtype=str, delimiter='\t')
file1 = np.matrix(file1)
(rows1, cols1) = file1.shape

assert (rows,cols) == (rows1,cols1)

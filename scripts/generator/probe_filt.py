# MIT License
#
# Copyright (c) 2025 CTLab-ITMO
#
# Authors: Daniil Smutin, Aleksandr Serdiukov, Vitalii Dravgelis, Artem Ivanov,
# Aleksei Zabashta, Sergey Muravyov, and the CTLab-ITMO university team.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


# imports
import sys
import numpy as np

# functions
def hash_check(str_x, hash_x, remove_arr, remove_hash):
        if not np.isin(hash_x, remove_hash):
            if not np.isin(str_x, remove_arr):
                return True
        return False
    
true_file, remove_file, output, max_out = sys.argv[1:]
max_out = int(max_out)
print(" -- load files")

# load remove list, hash and reorder
remove_arr = open(remove_file, "r").read().splitlines()
remove_hash = [hash(_) for _ in remove_arr]
print(" -- hashing done")

# init
results = []
analyze = 0
include = False
counter = 0

# read files
for line_str in open(true_file):
    if analyze == 0:
        line = line_str.split()
        analyze = int(line[6])

        str_x = line[0]
        hash_x = hash(str_x)

        # add multilocus hit check
        if hash_check(str_x, hash_x, remove_arr, remove_hash):
            #print("include ", str_x)
            include = True
            counter += 1
        else:
            #print("remove ", str_x)
            include = False

        if counter >= max_out:
                break
    else:
        analyze -= 1
    
    if include:
        results.append(line_str)

        
out = open(output, "w").write("".join(results))

print(" - done")
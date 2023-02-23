#!/usr/bin/env python3
#---------------------------------------------------------------------------
#    To sort the bed file second position in ascending order 
#---- * Usage : python3 sort_python.py -i {input_file} -o {output_file} -----*
#--------------------------------------------------------------------------
import argparse

""" Handle arguments """
parser = argparse.ArgumentParser(
    description="Sort the bed file"
)

parser.add_argument("-i", "--input_file",type=argparse.FileType('r'),
    help="input bed file for each sample",
)

parser.add_argument("-o", "--output_file", type=argparse.FileType('w'), 
      help="Output file"
      )

args = parser.parse_args()
bed_file = args.input_file
out_file = args.output_file


#create the function to sort the integer value for 2nd position
def sort_bed_file(file,out):
    lines = file.readlines()
    
    for i in range(len(lines) - 1):
        for j in range(i+1, len(lines)):
            if int(lines[i].split()[1]) > int(lines[j].split()[1]):
                    lines[i], lines[j] = lines[j], lines[i]
    
    out.writelines(lines)

sort_bed_file (bed_file,out_file)

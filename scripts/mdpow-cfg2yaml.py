#! /bin/python

import sys
import yaml
import os

# check extension

if len(sys.argv) < 3:
    print("Need input file and output file. Exiting...")
    exit(1)

if not sys.argv[1].split(".")[1] == "cfg":
    print("Please include *file*.cfg to convert to *file*.yml. Exiting...")
    exit(1)

print("Converting %s.cfg ..." % sys.argv[1].split(".")[0])

full_file = []

with open(sys.argv[1],'r') as f:
    for line in f:
        newline = line.strip()
        if not newline.startswith("#") and not newline == "":
            full_file.append(newline)

heads = []

for x in range(len(full_file)):
    if full_file[x].startswith("["):
        heads.append(x)

sections = []

for x in range(len(heads)-1):
    sections.append(full_file[heads[x]:heads[x+1]])

sections.append(full_file[heads[-1]:])

yaml_formatting = {}

for x in sections:
    options = x[1:]
    split = [[str(y.split("=")[0].strip()),str(y.split("=")[1].strip()).strip()] for y in options]
    yaml_formatting[x[0][1:-1]]=split

for x in yaml_formatting:
    d = dict((key,value) for (key,value) in yaml_formatting[x])
    yaml_formatting[x] = d

with open("result.yml",'w') as f:
    f.write(yaml.dump(yaml_formatting,default_flow_style=False))

with open("result.yml",'r') as infile, open(sys.argv[-1],'w') as outfile:
    data = infile.read()
    data = data.replace("'","")
    outfile.write(data)

os.remove("result.yml")

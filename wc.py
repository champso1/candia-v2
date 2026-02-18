import os


srcfilelist = os.listdir("src/")
srcfilelist.remove("fortran")
includefilelist = os.listdir("include/Candia-v2/")

filelistdict = {
    "src/":srcfilelist,
    "include/Candia-v2/":includefilelist
}
totalcount = 0

for prefix,filelist in filelistdict.items():
    for _f in filelist:
        filename = prefix+_f
        with open(filename, 'r') as f:
            count = len(f.readlines())
            totalcount += count
            print(f"File {filename} has {count} lines")

print(f"Total lines: {totalcount}")

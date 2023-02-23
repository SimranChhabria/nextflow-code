import glob, os

files = glob.glob('*.fastq.gz')

for file in files:
    file = file.strip()
    filearr = file.split('_')
    if(filearr[1].startswith('S')):
        del filearr[1]
        print(filearr[1])
    seqID = filearr[0]
    later = filearr[1:]
    #print (seqID)
    later = ".".join(later)
    #print(later)

    cmd = 'mv '+file+" "+seqID+".merged."+later
    
    #cmd = 'mv '+file+" "+".".join(filearr)
    os.system(cmd)

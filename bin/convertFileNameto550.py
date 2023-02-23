import glob, os

files = glob.glob('*.fastq.gz')

for file in files:
    file = file.strip()
    filearr = file.split('_')
    
    #assuming format is fixed
    seqID = filearr[0:3]
    remainder = filearr[3:]

    seqID_stitched = "-".join(seqID)
    remaining_stitched = "_".join(remainder)

    cmd = 'mv '+file+" "+seqID_stitched+"_"+remaining_stitched
    os.system(cmd)

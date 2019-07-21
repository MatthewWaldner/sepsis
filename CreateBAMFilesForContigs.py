import os
import sys

refNameBase = sys.argv[1]
assemblyGraph = sys.argv[2]
PairedReadsR1 = sys.argv[3]
PairedReadsR2 = sys.argv[4]
outDirectory = sys.argv[5]

print(refNameBase)
print(assemblyGraph)
print(PairedReadsR1)
print(PairedReadsR2)
print(outDirectory)

minimap2_Location = ''
samtools_Location = ''

os.system("python make_fasta_from_fastg.py -g "+assemblyGraph+" -o "+outDirectory+"/"+refNameBase+".nodes.fasta")

if (minimap2_Location != '' and samtools_Location != ''):
    os.system(minimap2_Location+'/minimap2 -ax sr '+outDirectory+'/'+refNameBase+'.nodes.fasta '+PairedReadsR1+' '+PairedReadsR2+' | '+samtools_Location+'/samtools view -buS - > '+outDirectory+'/'+refNameBase+'.reads_minimap2_pe.bam')
    os.system(samtools_Location+'/samtools sort '+outDirectory+'/'+refNameBase+'.reads_minimap2_pe.bam > '+outDirectory+'/'+refNameBase+'.reads_minimap2_pe.sort.bam')
    os.system(samtools_Location+'/samtools index '+outDirectory+'/'+refNameBase+'.reads_minimap2_pe.sort.bam')


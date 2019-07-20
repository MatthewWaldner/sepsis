import sys
import os
import re

inputFile = sys.argv[1]
inputFileBasename = sys.argv[2]
outFileName = sys.argv[3]

#requires GNU-sed, not mac sed
os.system("gsed '1~4 s/@/@"+inputFileBasename+"_/g' "+inputFile+" > "+outFileName)

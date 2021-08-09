#!/usr/bin/env python
"""
	The code to identify the graphlet signatures of all nodes in a network
	using the ORCA method.
	
	The script does the formatting of the given LEDA file to the format that 
	ORCA code can run. Then it reformats the output of the ORCA into ndump2 format.
	
	Run as:
		./count.py <LEDA_formatted_network>.gw
		
		Outputs: <LEDA_formatted_network>.ndump2 in the same folder with the input file
	
	Implemented by:
		Omer Nebil Yaveroglu
		05.11.2013 - 16:15
"""

import os
import sys
import pdb

"""
	Helper functions are defined here
"""
# Read the network in LEDA format
# Read LEDA formatted file into a network
def readLeda(networkFile):
	nodeList = []
	edgeList = []
	
	fRead = open(networkFile, 'r')
	
	mode = 0
	
	for line in fRead:
		if mode == 0 and line.startswith('|'):
			mode = 1
		
		if mode == 1 and line.startswith('|'):
			nodeList.append(line.strip().strip('|').strip('{').strip('}'))
		elif mode == 1 and not line.startswith('|'):
			mode = 2
		elif mode == 2 and line.strip().endswith('}|'):
			splitted = line.strip().split(' ')
			edgeList.append( [ int(splitted[0]) - 1, int(splitted[1]) - 1] )
	
	fRead.close()
	
	return (nodeList, edgeList)

# Write the network in a format that is ready to get executed by ORCA
def writeORCA(edgeList, nodeCount, outputFile):

	fWrite = open(outputFile, 'w')
	
	fWrite.write(str(nodeCount) + ' ' + str(len(edgeList)) + '\n')
	
	for edge in edgeList:
		fWrite.write(str(edge[0]) + ' ' + str(edge[1]) + '\n')

	fWrite.close()

# Read the temporary ndump2 file and create the original one
def formatNdump2(tempNdump2File, originalNdump2, nodeList):
	# Read the temporary result file
	fRead = open(tempNdump2File, 'r')
	fWrite = open(originalNdump2, 'w')

	lineID = 0
	for line in fRead:
		fWrite.write(str(nodeList[lineID]) + ' ' + line)
		
		lineID += 1

	fRead.close()	
	fWrite.close()
	

"""
	Main code starts here
"""
if __name__ == "__main__":
	netFileName = sys.argv[1]
	
	if not netFileName.endswith('.gw'):
		print 'ERROR: The network file should be in LEDA format!'
		exit(0)
	
	# Read the LEDA formatted network	
	(nodeList, edgeList) = readLeda(netFileName)
	
	# Write in ready to ORCA counting format
	outputFileName = netFileName.rsplit('.', 1)[0] + '_orca.txt'
	writeORCA(edgeList, len(nodeList), outputFileName)
	
	
	# Run the ORCA graphlet counting code with the resulting file
	tempNdump2File = netFileName.rsplit('.', 1)[0] + '_temp.ndump2'
	cmd = outputFileName + ' ' + tempNdump2File
	cmd = cmd.replace("l ","l\ ")
	cmd = cmd.replace("a ","a\ ")
	cmd = './orca 5 ' + cmd
	os.system(cmd)
	
	# Format the temp file to original format
	originalNdump2 = netFileName.rsplit('.', 1)[0] + '.ndump2'
	formatNdump2(tempNdump2File, originalNdump2, nodeList)
	
	# Counting finished remove the temp network file
	os.remove(outputFileName)
	os.remove(tempNdump2File)

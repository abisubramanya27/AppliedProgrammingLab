#Importing Libraries
import sys

def NetListReader(filePath):
	#Will contain all the sentences read from the file
	fullList = None

	try:
		with open(filePath,'r') as f:
			fullList = f.readlines()
			fullList = [s.strip() for s in fullList]
			#Checking if '.circuit' and '.end' is there in the Netlist
			if('.circuit' not in fullList or '.end' not in fullList) :
				print('Invalid Circuit Definition - No .circuit and/or .end found')
				sys.exit()
	#File Handling Errors
	except (IOError,OSError):
		print(f'Unable to Find/Open {filePath}')
		sys.exit()

	#Resultant Net List we need
	NetList = []
	l = -1
	for i in range(len(fullList)):
		sentence = fullList[i]
		if sentence == '.circuit':
			#Checking that previous '.circuit' ended with a '.end' before starting another
			if l != -1:
				print('Invalid Circuit Definition - .circuit already there')
				sys.exit()
			l = i
		elif sentence == '.end':
			#Checking there was a '.circuit' corresponding to the present '.end'
			if l == -1:
				print('Invalid Circuit Definition - No matching .circuit')
				sys.exit()
			l = -1
		else:
			if l != -1:
				NetList.append(sentence)

	#Checking for no Orphaned '.circuit'
	if l != -1:
		print('Invalid Circuit Definition - No matching .end')
		sys.exit()

	#Removing the comments as we are not going to process them
	NetList = [s.split('#')[0].strip() for s in NetList]
	#Removing extra white spaces
	NetList = [' '.join(s.split()) for s in NetList]

	#First Assignment - Printing in reverse
	for sentence in NetList[::-1] :
		words = [word for word in sentence.split()]
		for word in words[::-1] :
			print(word,end = ' ')
		print()


#Main Function
if __name__ == '__main__':
	if(len(sys.argv) != 2) :
		print('\nUsage: %s <inputfile>' % sys.argv[0])
		sys.exit()

	filePath = sys.argv[1]
	NetListReader(filePath)

#Importing Libraries
import sys
import numpy as np
import cmath as cm

#Class Definitions of Circuit Elements

#n1 and n2 are the nodes at two ends of the element
#Value is the value (in appropriate units) of the element
class Resistor():
	def __init__(self,**kwargs):
		self.name = kwargs.get('name',None)
		self.value = kwargs.get('value',0)
		self.n1 = kwargs.get('n1',0)
		self.n2 = kwargs.get('n2',0)

class Inductor():
	def __init__(self,**kwargs):
		self.name = kwargs.get('name',None)
		self.value = kwargs.get('value',0)
		self.n1 = kwargs.get('n1',0)
		self.n2 = kwargs.get('n2',0)

class Capacitor():
	def __init__(self,**kwargs):
		self.name = kwargs.get('name',None)
		self.value = kwargs.get('value',0)
		self.n1 = kwargs.get('n1',0)
		self.n2 = kwargs.get('n2',0)

class Ind_Vol_Source():		#Independent Voltage Source
	def __init__(self,**kwargs):
		self.name = kwargs.get('name',None)
		self.value = kwargs.get('value',0)
		self.n1 = kwargs.get('n1',0)
		self.n2 = kwargs.get('n2',-1)
		self.type = kwargs.get('type','dc')
		self.phase = kwargs.get('phase',float(0))

class Ind_Cur_Source():		#Independent Current Source
	def __init__(self,**kwargs):
		self.name = kwargs.get('name',None)
		self.value = kwargs.get('value',0)
		self.n1 = kwargs.get('n1',0)
		self.n2 = kwargs.get('n2',0)
		self.type = kwargs.get('type','dc')
		self.phase = kwargs.get('phase',float(0))

class VCVS():		#Voltage Controlled Voltage Source
	def __init__(self,**kwargs):
		self.name = kwargs.get('name',None)
		self.value = kwargs.get('value',0)
		self.n1 = kwargs.get('n1',0)
		self.n2 = kwargs.get('n2',-1)
		self.n3 = kwargs.get('n3',0)
		self.n4 = kwargs.get('n4',0)

class VCCS():		#Voltage Controlled Current Source
	def __init__(self,**kwargs):
		self.name = kwargs.get('name',None)
		self.value = kwargs.get('value',0)
		self.n1 = kwargs.get('n1',0)
		self.n2 = kwargs.get('n2',0)
		self.n3 = kwargs.get('n3',0)
		self.n4 = kwargs.get('n4',0)

class CCVS():		#Current Controlled Voltage Source
	def __init__(self,**kwargs):
		self.name = kwargs.get('name',None)
		self.value = kwargs.get('value',0)
		self.n1 = kwargs.get('n1','GND')
		self.n2 = kwargs.get('n2','something')
		self.n3 = kwargs.get('n3','GND')
		self.n4 = kwargs.get('n4','GND')

class CCCS():		#Current Controlled Current Source
	def __init__(self,**kwargs):
		self.name = kwargs.get('name',None)
		self.value = kwargs.get('value',0)
		self.n1 = kwargs.get('n1','GND')
		self.n2 = kwargs.get('n2','GND')
		self.n3 = kwargs.get('n3','GND')
		self.n4 = kwargs.get('n4','GND')

#Dictionary that stores nodes and branches number
code = 1  #The code increments as we add elements to the dictionary
hash = {'GND':0}

#The AC frequency if ac source is used
ACfreq =  None

#Incidence List
adjList = []

#Mx = b
M = None
b = None

no_voltage_sources = 0

#Circuit Elements present in the given circuit
Elements = []

def NetListReader(filePath):

	#########Start of Assignment 1#########
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

	global ACfreq
	count = 0
	l = -1
	for i in range(len(fullList)):
		sentence = fullList[i]
		if sentence == '.circuit':
			#Checking that previous '.circuit' ended with a '.end' before starting another
			if count > 0 or l != -1:
				print('Invalid Circuit Definition - .circuit already there')
				sys.exit()
			l = i
		elif sentence == '.end':
			#Checking if there was a '.circuit' corresponding to the present '.end'
			if l == -1:
				print('Invalid Circuit Definition - No matching .circuit')
				sys.exit()
			l = -1
			count += 1
		elif sentence.split()[0] == '.ac':
			#Checking if .circuit and .end block is over
			if count <= 0:
				print('Invalid Circuit Definition - No matching .circuit')
				sys.exit()
			#The last value in the sentence is the frequency
			#Storing 2*pi*f = omega directly
			ACfreq = cm.pi*2*float(sentence.split()[-1])
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

    #########End of Assignment 1#########

	global hash
	global code
	global Elements,no_voltage_sources

	branches = []

	for i in range(len(NetList)):
		sentence = NetList[i]
		words = sentence.split()

		#Checking proper format of circuit element definitions
		if len(words) < 4:
			print('Invalid circuit element definition!!')
			sys.exit()

		name = words[0]
		#Creating a new node incase the node doesn't previously exist
		try:
			n1 = hash[words[1]]
		except KeyError:
			hash[words[1]] = code
			n1 = code
			code += 1
		try:
			n2 = hash[words[2]]
		except KeyError:
			hash[words[2]] = code
			n2 = code
			code += 1

		#After checking circuit element definition formats, adding all the elements to "Element" list
		if name[0].upper() == 'R':
			if len(words) != 4:
				print('Check the format of resistor element definition!!')
				sys.exit()
			try:
				value = float(words[-1])
			except (TypeError,ValueError):
				print('Invalid Resistor Value!!')
				sys.exit()
			Elements.append(Resistor(value = value,name = name,n1 = n1,n2 = n2))

		elif name[0].upper() == 'L':
			if len(words) != 4:
				print('Check the format of inductor element definition!!')
				sys.exit()
			if ACfreq == None:
				print('Cannot handle Inductor in DC circuit!!')
				sys.exit()
			try:
				value = float(words[-1])
			except (TypeError,ValueError):
				print('Invalid Inductor Value!!')
				sys.exit()
			Elements.append(Inductor(value = value,name = name,n1 = n1,n2 = n2))

		elif name[0].upper() == 'C':
			if len(words) != 4:
				print('Check the format of capacitor element definition!!')
				sys.exit()
			if ACfreq == None:
				print('Cannot handle Capacitor in DC circuit!!')
				sys.exit()
			try:
				value = float(words[-1])
			except (TypeError,ValueError):
				print('Invalid Capacitor Value!!')
				sys.exit()
			Elements.append(Capacitor(value = value,name = name,n1 = n1,n2 = n2))

		elif name[0].upper() == 'V':
			if (words[3].lower() not in ('dc','ac')) or (words[3].lower() == 'dc' and len(words) != 5) or (words[3].lower() == 'ac' and len(words) != 6):
				print('Check the format of voltage source definition!!')
				sys.exit()
			if words[3].lower() == 'ac' and ACfreq == None:
				print('Cannot handle AC source in DC circuit!!')
				sys.exit()
			elif words[3].lower() == 'ac':
				value = words[4]
				phase = words[-1]
			if words[3].lower() == 'dc' and ACfreq != None:
				print('Cannot handle DC source in AC circuit!!')
				sys.exit()
			elif words[3].lower() == 'dc':
				value = words[4]
				phase = '0'
			#Assuming input phase to be in degrees
			try:
				value = float(value)
				phase = float(phase)
			except (TypeError,ValueError):
				print('Invalid Voltage Value/Phase!!')
				sys.exit()
			#Cannot find individual currents incase many voltage sources are in parallel, so effectively one voltage source is enough
			if hash.get(str(n1)+'-'+str(n2),None) == None and hash.get(str(n2)+'-'+str(n1),None) == None:
				hash[str(n1)+'-'+str(n2)] = code
				branches.append(str(n1)+'-'+str(n2))
				no_voltage_sources += 1
				#If AC circuit, peak-to-peak voltage is given, max voltage is half of it
				if words[3].lower() == 'ac':
					value *= .5
				Elements.append(Ind_Vol_Source(value = value,phase = phase,type = words[3].lower(),name = name,n1 = n1,n2 = n2))

		elif name[0].upper() == 'I':
			if (words[3].lower() not in ('dc','ac')) or (words[3].lower() == 'dc' and len(words) != 5) or (words[3].lower() == 'ac' and len(words) != 6):
				print('Check the format of current source definition!!')
				sys.exit()
			if words[3].lower() == 'ac' and ACfreq == None:
				print('Cannot handle AC source in DC circuit!!')
				sys.exit()
			elif words[3].lower() == 'ac':
				value = words[4]
				phase = words[-1]
			if words[3].lower() == 'dc' and ACfreq != None:
				print('Cannot handle DC source in AC circuit!!')
				sys.exit()
			elif words[3].lower() == 'dc':
				value = words[4]
				phase = '0'
			#Assuming input phase to be in degrees
			try:
				value = float(value)
				phase = float(phase)
			except (TypeError,ValueError):
				print('Invalid Current Value/Phase!!')
				sys.exit()
			#If AC circuit, peak-to-peak current is given, max current is half of it
			if words[3].lower() == 'ac':
				value *= .5
			Elements.append(Ind_Cur_Source(value = value,phase = phase,type = words[3].lower(),name = name,n1 = n1,n2 = n2))
	    
	    #Since we dont deal with dependent sources for now
		else:
			print('Invalid circuit element!!')
			sys.exit()

	#Now we add the number code for branches atlast, because we want the nodes first then the branches to ease forming equations
	for branch in branches:
		hash[branch] = code
		code += 1

	
def solveNetList():
	global M,b,adjList,Elements
	#len(hash) - say,N is the total number of variables we need to figure out. So M is NxN and b is Nx1 matrix
	#Mx = b is the equation, we need x
	M = np.zeros(shape = (len(hash),len(hash)),dtype = complex)
	b = np.zeros(shape = len(hash),dtype = complex)
	
	#hash contains the code (number starting from 0) for all variables that we need to figure out - node voltages and branch currents across voltage sources
	#adjList contains len(hash) lists of elements which are incident to the node or are part of the branch
	adjList = list([] for _ in range(len(hash)))
	for element in Elements:
		if element.name[0] == 'V':
			#We added branch current across voltage source from n1 to n2 node only, and not the opposite way 
			adjList[hash[str(element.n1)+'-'+str(element.n2)]].append(element)
		adjList[element.n1].append(element)
		adjList[element.n2].append(element)

	for i in range(len(adjList)):
		for ele in adjList[i]:
			#First equation for V0 node will be V0 = 0
			if i == 0:
				M[i][0] += 1
				b[i] += 0+0j
			#Last (no_voltage_sources) equations will be for voltage across voltage sources
			elif i + no_voltage_sources >= len(adjList):
				n1,n2 = ele.n1,ele.n2
				M[i][n1] += 1
				M[i][n2] -= 1
				value = ele.value
				if ele.type == 'ac':
					value = cm.rect(ele.value,(cm.pi/180)*ele.phase)
				b[i] += value
			#After V0,next (len(adjList)-no_voltage_sources-1) questions will be nodal equations
			#For the passive elements 1/impedance must be added at origin node and subtracted at ending node
			elif ele.name[0].upper() in ('L','R','C') :
				value = (1/ele.value)
				if ele.name[0].upper() == 'L':
					value = (1/(1j*ACfreq*ele.value))
				elif ele.name[0].upper() == 'C':
					value = (1j*ACfreq*ele.value)
				if i == ele.n1:
					M[i][ele.n1] += value
					M[i][ele.n2] -= value
				else:
					M[i][ele.n2] += value
					M[i][ele.n1] -= value
			#For current source the constant current value must be added/subtracted to b matrix based on its polarity
			#Taking n1 as start and n2 as end node for current source
			elif ele.name[0].upper() == 'I':
				value = ele.value
				if ele.type == 'ac':
					value = cm.rect(ele.value,(cm.pi/180)*ele.phase)
				if i == ele.n1:
					b[i] -= value
				else:
					b[i] += value
			#For voltage source the current (stored as variable) through it must be added/subtracted based on its polarity
			#Taking n1 as positive and n2 as negative node for current source
			else:
				if i == ele.n1:
					M[i][hash[str(ele.n1)+'-'+str(ele.n2)]] += 1
				else:
					M[i][hash[str(ele.n1)+'-'+str(ele.n2)]] -= 1

	if(np.linalg.det(M) == 0):
		print('Circuit/Equations not quite correct. Check again!!')
		sys.exit()

	#We now have the M and b matrices. Lets solve for x
	x = np.linalg.solve(M,b)
	
	#######Printing the values of the unknows

	tmp = 0
	#Generating list of keys to get node name corresponding to node's number code
	keys_val = list(hash.items())
	#Sorting by value so that the keys are in order of increasing number
	keys_val.sort(key = lambda x: x[-1])
	keys = [k for (k,_) in keys_val]
	for k in keys:
		if ACfreq == None:
			val = x[hash[k]].real
			if k.find('-') != -1:
				n1,n2 = k.split('-')
				print(f"I in branch {keys[int(n1)]}-{keys[int(n2)]} : {val:.5E}")
			else:
				print(f"V at node {k} : {val:.5E}")
		else:
			val = cm.polar(x[hash[k]])
			#Outputting phase in degrees
			if k.find('-') != -1:
				n1,n2 = k.split('-')
				print(f"I in branch {keys[int(n1)]}-{keys[int(n2)]} : {val[0]:.5E} L({val[1]*(180/cm.pi):.5f})")
			else:
				print(f"V at node {k} : {val[0]:.5E} L({val[1]*(180/cm.pi):.5f})")
		tmp += 1


#Main Function
if __name__ == '__main__':
	if(len(sys.argv) != 2) :
		print('\nUsage: %s <inputfile>' % sys.argv[0])
		sys.exit()

	filePath = sys.argv[1]
	#######Extracts the information about circuit and elements from the Netlist file
	NetListReader(filePath)

	#######Formas the M and b matrices to solve Mx = b equation and print the x matrix
	solveNetList()

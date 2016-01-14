import xml.etree.ElementTree as ET

def setParameter( root, name, value ):
	for child in root.iter('Parameter'):
		if(child.attrib['name']==name):
			child.attrib['value']=str(value)

def set( fname, name, value ):
	tree = ET.parse(fname)
	#
	root = tree.getroot()
	#
	setParameter( root, name, value )
	#for child in root.iter('Parameter'):
		#if(child.attrib['name']==name):
			#child.attrib['value']=str(value)
	#	
	tree.write(fname)

def setIDS( fname ):
	tree = ET.parse(fname)
	#
	root = tree.getroot()
	#
	i = 0
	for child in root.iter('Parameter'):
		child.attrib['id']=str(i)
		i = i+1
	for child in root.iter('ParameterList'):
		child.attrib['id']=str(i)
		i = i+1
	#	
	tree.write(fname)

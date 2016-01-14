import os
from numpy import linspace
import numpy as np
from math import pi
from platform_paths import *
import xml.etree.ElementTree as ET

tree = ET.parse('plSMGIn.xml')

root = tree.getroot()

exe = 'peri_navierXML'


os.chdir( exe_path )
os.system( 'make '+exe+' -j4' )




case_path = ['','','','']

lambdas = [ 0.11, 0.21, 0.31 ]

case_path[0] = '/case_study2'
if not os.path.exists( data_path+case_path[0] ):
	os.mkdir( data_path+case_path[0] )
for lam in lambdas:
	case_path[1] = '/ST_'+str(lam)
	if not os.path.exists( data_path+case_path[0]+case_path[1] ):
		os.mkdir( data_path+case_path[0]+case_path[1] )
	os.chdir( data_path+case_path[0]+case_path[1] )
	os.system(' rm ./* -r -v  ')
	for child in root.iter('Parameter'):
		if(child.attrib['name']=='alpha2'):
			child.attrib['value']=str(2.*3.1416*lam*200.)
	tree.write('parameterIn.xml')
	os.system( exe_pre(4,' -W 21:00 ') + exe_path+exe + ' > output ' )


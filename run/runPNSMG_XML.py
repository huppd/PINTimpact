import os
from platform_paths import *
import xml.etree.ElementTree as ET

tree = ET.parse('plSMGIn.xml')

root = tree.getroot()

exe = 'peri_navierXML'


os.chdir(EXE_PATH)
os.system('make '+exe+' -j4')


CASE_PATH = ['']*4

lambdas = [0.11, 0.21, 0.31]

CASE_PATH[0] = '/case_study2'
if not os.path.exists(DATA_PATH+CASE_PATH[0]):
    os.mkdir(DATA_PATH+CASE_PATH[0])
for lam in lambdas:
    CASE_PATH[1] = '/ST_'+str(lam)
    if not os.path.exists(DATA_PATH+CASE_PATH[0]+CASE_PATH[1]):
        os.mkdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
    os.chdir(DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
    os.system(' rm ./* -r -v  ')
    for child in root.iter('Parameter'):
        if(child.attrib['name']=='alpha2'):
            child.attrib['value']=str(2.*3.1416*lam*200.)
    tree.write('parameterIn.xml')
    os.system(exe_pre(4,' -W 21:00 ') + EXE_PATH+exe + ' > output ')

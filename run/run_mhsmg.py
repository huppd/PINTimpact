""" runner for peri_naverXML deprecated"""
import os
import xml.etree.ElementTree as ET

import platform_paths as pp

TREE = ET.parse('plSMGIn.xml')

ROOT = TREE.getroot()

EXE = 'peri_navierXML'


os.chdir(pp.EXE_PATH)
os.system('make '+EXE+' -j4')


CASE_PATH = ['']*4

LAMBDAS = [0.11, 0.21, 0.31]

CASE_PATH[0] = '/case_study2'
if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]):
    os.mkdir(pp.DATA_PATH+CASE_PATH[0])
for lam in LAMBDAS:
    CASE_PATH[1] = '/ST_'+str(lam)
    if not os.path.exists(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1]):
        os.mkdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
    os.chdir(pp.DATA_PATH+CASE_PATH[0]+CASE_PATH[1])
    os.system(' rm ./* -r -v  ')
    for child in ROOT.iter('Parameter'):
        if child.attrib['name'] == 'alpha2':
            child.attrib['value'] = str(2.*3.1416*lam*200.)
    TREE.write('parameterIn.xml')
    os.system(pp.exe_pre(4, ' -W 21:00 ') + pp.EXE_PATH+EXE + ' > output ')

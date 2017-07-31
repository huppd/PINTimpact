""" manipulates XML files for mains"""
import xml.etree.ElementTree as ET


def set_parameter(root, name, value):
    """ sets parameter to value """
    for child in root.iter('Parameter'):
        if child.attrib['name'] == name:
            child.attrib['value'] = str(value)


def set_insublist(root, sublist_name, parameter_name, pvalue):
    """ set paramter in sublist """
    for child in root.iter('ParameterList'):
        if child.attrib['name'] == sublist_name:
            for cchild in child.iter('Parameter'):
                if cchild.attrib['name'] == parameter_name:
                    cchild.attrib['value'] = str(pvalue)


def set_ids(fname):
    """ sets id correct """
    tree = ET.parse(fname)
    #
    root = tree.getroot()
    #
    i = 0
    for child in root.iter('Parameter'):
        child.attrib['id'] = str(i)
        i = i+1
    for child in root.iter('ParameterList'):
        child.attrib['id'] = str(i)
        i = i+1
    #
    tree.write(fname)


def set_output(fname):
    """ removes all 'Output Stream' and 'MyPID' """
    tree = ET.parse('parameterOut.xml')
    #
    root = tree.getroot()
    # remove Output Stream
    parents = root.findall('.//Parameter[@name="Output Stream"]...')
    for parent in parents:
        for child in parent.iter():
            if child.attrib['name'] == 'Output Stream':
                parent.remove(child)
    # remove My PID
    parents = root.findall('.//Parameter[@name="MyPID"]...')
    for parent in parents:
        for child in parent.iter():
            if child.attrib['name'] == 'MyPID':
                parent.remove(child)
    #
    tree.write(fname)
    set_ids(fname)

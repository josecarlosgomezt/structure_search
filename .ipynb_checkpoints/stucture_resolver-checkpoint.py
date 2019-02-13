import pubchempy as pbc
import pandas as pd
import urllib3



# Initialize urllib3
http = urllib3.PoolManager()
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)


# Get Inchi from cactus
def inchi_from_cactus(identifier):
    url = (f'https://cactus.nci.nih.gov/chemical/structure/{identifier}/inchi')
    
    try:
        response = http.request('GET', url)
    except:
        return False
    if "Bad" in str(response.data):
        return False
    if "found" in str(rlesponse.data):
        return False
    inchi = str(response.data.decode("UTF-8"))
    if not 'InChI' in inchi:
        return False
    return inchi


# Get Inchi from drugbankID
def inchi_from_drugbank(identifier):

    url = (f'https://www.drugbank.ca/structures/small_molecule_drugs/'
           f'{identifier}.inchi') 
    try:
        response = http.request('GET', url)
    except:
        return False

    if "Bad" in str(response.data):
        return False
    if "found" in str(response.data):
        return False
    inchi = str(response.data.decode("UTF-8"))
    if not 'InChI' in inchi:
        return False
    return inchi

# Get Inchi from pubchempy
def inchi_from_pubchem(identifier):
    try:
        comp =  pbc.get_compounds(identifier, namespace='name')
        inchi = comp[0].inchi
    except Exception as e:
        return False
    if not 'InChI' in inchi:
        return False
    return inchi
    


# Main function. Return a DataFrame with found structures
def add_inchis(frame):
    names = False
    CAS = False
    drugbank = False
    inchis = []
    if 'name' in frame.columns:
        names = frame['name']
    if 'CAS' in frame.columns:
        CAS = frame['CAS']
    if 'DrugbankID' in frame.columns:
        drugbank = frame['DrugbankID']
    
    for i in range(len(frame)):
        print('Molecule: ' + str(i))
        inchi = False
        if len(CAS) > 0:
            'Trying CAS'
            inchi = inchi_from_cactus(CAS[i])
        if not inchi  and len(names) > 0:
            print('Trying name')
            inchi = inchi_from_cactus(names[i])
            if not inchi:
                inchi = inchi_from_pubchem(names[i])
        if not inchi and len(drugbank) > 0:
            print('Trying drugbank')
            inchi = inchi_from_drugbank(drugbank[i])
        print(inchi)
        print(f'Molecule {str(i)} finished.\n')
        inchis.append(inchi)
        
    frame['InChI'] = inchis
    return frame

        
import pubchempy as pbc
import pandas as pd
import urllib3


http = urllib3.PoolManager()
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)



def inchi_from_cactus(identifier):

    '''
    
    Info
    ----

    Function to get InChIs from CASRN or name or smile using cactus web service
    
    Parameters
    ----------
    
    identifier: str  
        String that can be a compound name or a CAS Registry number or smile
       
    Returns
    -------
    
    str
        InChI
    
    Example
    -------        

    inchi = InChi_from_cactus('c1ccccc1')
         
    '''
    
    if str(identifier) == 'nan':
        identifier = 'error'
    else:
        identifier = identifier
    
    url = (f'https://cactus.nci.nih.gov/chemical/structure/{identifier}/inchi')

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



def inchi_from_drugbank(identifier):

    '''
    
    Info
    ----

    Function to get InChIs from Drug Bank ID using drugbank webservice
    
    Parameters
    ----------
    
    identifier: str  
        Drug Bank ID
       
    Returns
    -------
    
    str
        InChI
    
    Example
    -------        

    inchi = inchi_from_drugbank('DB11558')
         
    '''

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


def inchi_from_pubchem(identifier):
    
    '''
    
    Info
    ----

    Function to get InChIs from CASRN or name or smile using pubchem web service
    To use this function, one needs to install this : pip install pubchempy
    
    Parameters
    ----------
    
    identifier: str  
        String that can be a compound name or a CAS Registry number or smile
       
    Returns
    -------
    
    str
        InChI
    
    Example
    -------        

    inchi = inchi_from_pubchem('toluene')
         
    '''   

    try:
        comp =  pbc.get_compounds(identifier, namespace='name')
        inchi = comp[0].inchi
    except Exception as e:
        return False
    if not 'InChI' in inchi:
        return False
    return inchi
    


def add_inchis(frame, name = None, CASRN = None, DBID= None):

    '''
    
    Info
    ----

    Main Function to add InChI column to pandas dataframe from name, CASRN,
    smiles or DrugBank ID provided. This function calls other functions:
    - inchi_from_pubchem()
    - inchi_from_drugbank()
    - inchi_from_cactus()
    
    Parameters
    ----------
    
    frame: pd.DataFrame  
        Dataframe containing whole information(compound name or a CAS Registry number or smile)
    name: str
        Compound name or smile column
    CASRN: str
        CAS registry number
    DBID: str
        Drug Bank ID
           
    Returns
    -------
    
    frame: pd.DataFrame
        Dataframe with an inchi column. It will appear False if there is no inchi.
    
    Example
    -------        
    
    frame = add_inchis(frame)
    checking not found structures: frame[frame['InChI']==False]
    
    '''
       
    if name != None:
        names = frame[name]
    else:
        names = []
    
    if CASRN != None:
        CAS = frame[CASRN]
    else:
        CAS = []
    
    if DBID != None:
        drugbank = frame[DBID]
    else:
        drugbank = []
        
    inchis = []
 
    for i in range(len(frame)):
        inchi = False
        if len(CAS) > 0:
            inchi = inchi_from_cactus(CAS[i])
        if not inchi and len(names) > 0:
            inchi = inchi_from_cactus(names[i])
            if not inchi:
                inchi = inchi_from_pubchem(names[i])
        if not inchi and len(drugbank) > 0:
            inchi = inchi_from_drugbank(drugbank[i])
        inchis.append(inchi)
    frame['InChI'] = inchis
    return frame
    


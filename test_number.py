
import pytest


#from pathlib import Path
from structure_resolver import *
import pubchempy as pbc
import pandas as pd
import urllib3

#csv = Path("test_data/test_data.csv")


def structure_retrieval():
    """test compound retrieval"""

    frame = pd.read_csv("./test_data/test_data2.csv", sep="\t")
    frame = add_inchis(frame, name='name', CASRN='CAS',
                        DBID='DrugbankID')
    print(frame.shape)
    return(len(frame))


def test_number():
    assert structure_retrieval() > 10

"""
"""


from IIRrational.v2 import (
    data2filter,
    ResultsAid,
)

#from IIRrational.RDF import RationalDiscFilter
#from IIRrational.RDF2MRF import RDF2MRF
#from IIRrational.MRF import MultiReprFilterZ
from IIRrational.testing import IIRrational_data

#monkeypatch in the msurrogate annotations
ResultsAid._msurrogate_MT             = False
data2filter._msurrogate_MT        = True

__all__ = [
    'data2filter',
    'IIRrational_data',
]

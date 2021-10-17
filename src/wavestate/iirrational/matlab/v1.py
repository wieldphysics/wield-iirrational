"""
"""


from IIRrational.v1 import (
    data2filter,
    FitAid,
)

from IIRrational.testing import IIRrational_data

#monkeypatch in the msurrogate annotations
FitAid._msurrogate_MT             = False
data2filter._msurrogate_MT        = True

__all__ = [
    'data2filter',
    'IIRrational_data',
]

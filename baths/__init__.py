import numpy as np
### Importing baths
from Feom.baths.debye_bath import Debye_bath
from Feom.baths.debye_aaa import Debye_aaa              # NOTE: This is not implemented yet
from Feom.baths.debye_cothpoles import Debye_cothpoles  # NOTE: This is not implemented yet
### Mapping of potential names to potential classes
bath_mapping = {
'debye': Debye_bath,
'debyeAAA': Debye_aaa,  
'debyeCothpoles': Debye_cothpoles,
}

### Function to select potential class from a given string
def getbath(bath_name: str):
    if bath_name not in bath_mapping:
        raise ValueError(f'Invalid bath name: {bath_name}')
    return bath_mapping[bath_name]  # Return bath class

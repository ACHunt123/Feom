import numpy as np
### Importing baths
from Heom.baths.debye_bath import Debye_bath
### Mapping of potential names to potential classes
bath_mapping = {
'debye': Debye_bath,
}

### Function to select potential class from a given string
def getbath(bath_name: str):
    if bath_name not in bath_mapping:
        raise ValueError(f'Invalid bath name: {bath_name}')
    return bath_mapping[bath_name]  # Return bath class

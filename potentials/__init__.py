import numpy as np
### Importing potential classes
from Heom.potentials.harm_oscillator import Harmonic_oscillator
from Heom.potentials.spin_boson import Spin_boson

### Mapping of potential names to potential classes
potential_mapping = {
'harmonic': Harmonic_oscillator,
'spinboson': Spin_boson
}

### Function to select potential class from a given string
def getpotential(potential_name: str):
    if potential_name not in potential_mapping:
        raise ValueError(f'Invalid potential name: {potential_name}')
    return potential_mapping[potential_name]  # Return potential class

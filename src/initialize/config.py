import json
import numpy as np
from dataclasses import dataclass
from typing import Optional

# --- 1. The Helper to handle Numpy/Complex JSON ---
class ScientificEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, complex):
            return {"re": obj.real, "im": obj.imag}
        return super().default(obj)

def scientific_decoder(dct):
    if "re" in dct and "im" in dct:
        return complex(dct["re"], dct["im"])
    return dct

# --- 2. The Configuration Object ---
@dataclass
class SimConfig:
    system: dict
    bath: dict
    params: dict
    terminator: Optional[dict] = None

    def save(self, filename):
        """Saves this config to a clean JSON file."""
        with open(filename, 'w') as f:
            json.dump(self.__dict__, f, cls=ScientificEncoder, indent=4)

    @classmethod
    def load(cls, filename):
        """Loads from JSON and reconstructs the object."""
        with open(filename, 'r') as f:
            data = json.load(f, object_hook=scientific_decoder)
        
        # Optional: Convert lists back to arrays here if strict types are needed
        # e.g., data['system']['hamiltonian'] = np.array(data['system']['hamiltonian'])
        
        return cls(**data)
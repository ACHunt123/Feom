"""
config_definitions.py

Defines the required arguments and default structures for the HEOM solver.
"""

# -----------------------------------------------------------------------------
# SYSTEM REQUIREMENTS
# -----------------------------------------------------------------------------
REQUIRED_SYS_MANUAL = {
    'H_mat',   # The system Hamiltonian matrix (numpy array)
    's_mat',   # The system coupling operator matrix (numpy array)
}

# -----------------------------------------------------------------------------
# BATH REQUIREMENTS
# -----------------------------------------------------------------------------
REQUIRED_BATH_MANUAL = {
    'C_ks',      # the coefficients behind the exponential terms
    'gam_ks',   # the exponential decay terms
    'zeta',     # the coeficient behind the delta function
}

# -----------------------------------------------------------------------------
# MISC PARAMS REQUIREMENTS + DEFAULTS
# -----------------------------------------------------------------------------
REQUIRED_PARAMS_MANUAL = {
    'dt',       # Time step
    'tmax',     # Max simulation time (same units as dt)
    'L',        # Hierarchy depth
}
DEFAULT_PARAMS_MANUAL = {
            'noSIA': 0,         # Whether to use SIA or not
            'print_ADOs': 0,    # Whether to output all of the ADOs (very memory intensive)
        }

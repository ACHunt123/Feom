"""
config_requirements.py

Defines the required arguments and default structures for the HEOM solver.
"""

# -----------------------------------------------------------------------------
# SYSTEM REQUIREMENTS
# -----------------------------------------------------------------------------
REQUIRED_SYS_MANUAL = {
    'H_mat',   # The system Hamiltonian matrix (numpy array)
    'V_ks_plus',   # The system coupling operator matrices for the sigma=+ (list of numpy arrays)
}

# -----------------------------------------------------------------------------
# BATH REQUIREMENTS
# -----------------------------------------------------------------------------
REQUIRED_BATH_MANUAL = {
    'C_ks_plus',      # the coefficients behind the exponential terms
    'C_ks_mnus',      # the coefficients behind the exponential terms
    'gam_ks_plus',   # the exponential decay terms
    'gam_ks_mnus',   # the exponential decay terms
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
            'noSIA': False,         # Whether to use SIA or not
            'print_ADOs': False,    # Whether to output all of the ADOs (very memory intensive)
            'IT_RWAterms':False      # Whether to do the Ishizaki-Tanimura on the RWA inputs
        }

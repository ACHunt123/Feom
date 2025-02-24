#!/usr/bin/env python
# File: type_dicts.py
"""Dictionaries of types."""
import heom.states as states
import heom.potentials as pots
import heom.simulations as sims

potentials = {
        "PotEff": pots.PotEff,
        0: pots.PotHO,
        1: pots.PotTani,
        3: pots.PotPoly,
        "hydoh": pots.PotHydOH,
        "harmonic": pots.PotHarmonic,
        "anharmonic": pots.PotAnharmonic,
        "quartic": pots.PotQuartic,
        "harmoquartic": pots.PotHarmoquartic,
        "champagne": pots.PotChampagne,
        "morse": pots.PotMorse,
        "dw1": pots.PotDW1,
        "dw2": pots.PotDW2,
        "twoLsys":pots.PotTwoLsys,
        }

states = {
        "wavepacket": states.StateWP,
        "tanimura": states.StateTani,
        "thermal": states.StateThermEq,
        "loaded": states.StateLoaded,
        "thermal_aren": states.StateThermEq_aren,
        }

simulations = {
        "highT_qm": sims.SimQHEOM,
        "highT_cl": sims.SimCHEOM,
        "lowT_qm": sims.SimQHEOM_Mats,
        "lowT_cl": sims.SimCHEOM_Mats,
        "lowT_qm_ee": sims.SimQHEOM_Mats_EE,
        "heomc": sims.SimQHEOM_Mats_EE_heomc,
        #"pyheom": sims.SimQHEOM_pyheom,
        #5: sims.SimQHEOM_Mats_EE_scaled,
        #6: sims.SimCHEOM_1side,
        }

bath_types = {
        "none",
        "debye",
        "debye_correction",
        "debye_correction2",
        "debye_correction3_cl",
        "white",
        }

heom_switches = {
        "none",
        "potential",
        "tcf_qq",
        "kubo_qq",
        "tcf_q2q2",
        "kubo_q2q2",
        }

a_ren_types = {
        "none",
        "gamma*eta",
        "eta",
        }

truncation_types = {
        "zero",
        "anchor",
        }

scaling_types = {
        "none",
        "shi",
        }

pruning_types = {
        "none",
        "pruning",
        "radical",
        }

sop_decomposition_types = {
        "matsubara",
        "pade_N-1/N",
        "pade_N/N",
        "pade_N+1/N",
        "nbead",
        }

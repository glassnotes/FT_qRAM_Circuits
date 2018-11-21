# -*- coding: utf-8 -*-                                                            
#                                                                                  
# double_qram_circuits.py: qRAM circuits for Cartesian-product address spaces. 
#                                                                                  
# © 2018 Olivia Di Matteo (odimatte@uwaterloo.ca)                                  
#                                                                                  
# This file is part of the project FT_qRAM_Circuits.                                      
# Licensed under MIT License. 

from circuit import *
from surface_code import *

# Double QRAM circuits
# Note: In the worst case, the layers of T-depth from the two individual QRAMs don't line up,
#       in which case the T-depth would be the sum of the T-depths of both parts. In the best case,
#       they all line up, in which case the T-depth would be the max of both parts. In what's below,
#       we assume the worst. The normal depth is the max in either case. 
#
# We could also use the T-depth 1 implementation of the Toffoli here, but this would only save 2 layers
# of T depth overall, which isn't much especially for large memories

class Double_QRAM(qRAMCircuit):
    def __init__(self, n, q, k, b1, b2):
        super().__init__(n, q, k, b1, b2)

        self.params["name"] = "Double_QRAM"

        # We don't need to consider k > q or k < q cases here; you can find another combination of  
        # q, b1, b2, etc. because of all the symmetry (e.g. interchange upper and lower blocks)
        # For the normal case, we can't really do anything (that I immediately see) to uncompute without
        # running all the MPMCTs backwards. 
        self.params["n_qubits"] = 2*n + 1
        self.params["depth"] = 2 * ( max(pow(2,b1)*depth(k), pow(2,b2)*depth(n-k)) ) + 10
        self.params["t_count"] = 2 * ( pow(2,b1)*t_c(k) + pow(2,b2)*t_c(n-k) ) + 7
        self.params["t_depth"] = 2 * ( pow(2,b1)*t_d(k) + pow(2,b2)*t_d(n-k) ) + 3
        self.params["h_count"] = 2 * ( pow(2,b1)*h_c(k) + pow(2,b2)*h_c(n-k) ) + 2
        self.params["cnot_count"] = 2 * ( pow(2,b1)*cnot_c(k) + pow(2,b2)*cnot_c(n-k) ) + 7 
        self.params["cliffords"] = self.params["h_count"] + self.params["cnot_count"]


class Double_QRAM_Tier1Parallel(qRAMCircuit):
    def __init__(self, n, q, k, b1, b2):
        super().__init__(n, q, k, b1, b2)

        self.params["name"] = "Double_QRAM_Tier1Parallel"
        # For any parallel tiers, we can copy the outputs down to the even-superposition state
        # to eliminate the need to uncompute, at the cost of a couple CNOTs. But then, for the
        # non-parallel tier, we still need to run that layer backwards, so it would be kinda
        # pointless to do this.

        # Parallelize top tier + 1 for output + remaining n - k + all ancillaes + final output
        self.params["n_qubits"] = (k+1)*pow(2,b1) + 1 + (n-k) + pow(2,b1)*(k-1) + (n-k-1) + 1 + 1
        self.params["depth"] =  max(2*b1 + depth(k) + 1, pow(2,b2)*depth(n-k)) + max(2*b1 + 1, pow(2, b2)*depth(n-k))+ 10
        self.params["t_count"] = pow(2,b1)*t_c(k) + 2*pow(2,b2)*t_c(n-k) + 7
        self.params["t_depth"] = t_d(k) + 2*pow(2,b2)*t_d(n-k) + 3
        self.params["h_count"] = pow(2,b1)*h_c(k) + 2*pow(2,b2)*h_c(n-k) + 2
        # 2*(fanout top layer + parity + 1) + MPMCTs top layer + 2*MPMCTs bottom layer + 7 for Toffoli
        self.params["cnot_count"] = 2*(k*(pow(2,b1)-1) + pow(2,b1)-1 + 1) + pow(2,b1)*cnot_c(k) + 2*pow(2, b2)*cnot_c(n-k) + 7
        self.params["cliffords"] = self.params["h_count"] + self.params["cnot_count"]


class Double_QRAM_Tier2Parallel(qRAMCircuit):
    def __init__(self, n, q, k, b1, b2):
        super().__init__(n, q, k, b1, b2)

        self.params["name"] = "Double_QRAM_Tier2Parallel"

        self.params["n_qubits"] = k + 1 + (n-k+1)*pow(2,b2) + 1 + (k-1) + pow(2,b2)*(n-k) + 1
        self.params["depth"] =  max(pow(2,b1)*depth(k), 2*b2 + depth(n-k) + 1) + max(pow(2,b1)*depth(k), 2*b2 + 1) + 10
        self.params["t_count"] = 2*pow(2,b1)*t_c(k) + pow(2,b2)*t_c(n-k) + 7
        self.params["t_depth"] = 2*pow(2,b1)*t_d(k) + t_d(n-k) + 3
        self.params["h_count"] = 2*pow(2,b1)*h_c(k) + pow(2,b2)*h_c(n-k) + 2
        self.params["cnot_count"] = 2 * ( pow(2,b1)*cnot_c(k) + (n-k+1)*(pow(2,b2)-1) + 1) + pow(2,b2)*cnot_c(n-k) + 7
        self.params["cliffords"] = self.params["h_count"] + self.params["cnot_count"]


class Double_QRAM_Parallel(qRAMCircuit):
    def __init__(self, n, q, k, b1, b2):
        super().__init__(n, q, k, b1, b2)

        self.params["name"] = "Double_QRAM_Parallel"

        # Use the even-parity superposition trick for both tiers; need to compute all MPMCTs only once
        self.params["n_qubits"] = (k+1)*pow(2, b1) + 1 + (n-k+1)*pow(2, b2) + 1 + pow(2,b1)*(k-1) + pow(2,b2)*(n-k) + 1
        self.params["depth"] = max(b1 + depth(k) + b1 + 1, b2 + depth(n-k) + b2 + 1) + 10 + max(b1 + 1, b2 + 1)
        self.params["t_count"] = pow(2,b1)*t_c(k) + pow(2,b2)*t_c(n-k) + 7
        self.params["t_depth"] = t_d(k) + t_d(n-k) + 3
        self.params["h_count"] = pow(2,b1)*h_c(k) + pow(2,b2)*h_c(n-k) + 2
        self.params["cnot_count"] = pow(2,b1)*cnot_c(k) + pow(2,b2)*cnot_c(n-k) + 2*((k+1)*(pow(2,b1)-1) + (n-k+1)*(pow(2,b2)-1) + 2) + 7
        self.params["cliffords"] = self.params["h_count"] + self.params["cnot_count"]
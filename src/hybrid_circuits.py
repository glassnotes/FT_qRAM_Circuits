# -*- coding: utf-8 -*-                                                            
#                                                                                  
# hybrid_circuits.py: Hybrid qRAM circuits and parallelized versions. 
#                                                                                  
# © 2018 Olivia Di Matteo (odimatte@uwaterloo.ca)                                  
#                                                                                  
# This file is part of the project FT_qRAM_Circuits.                                      
# Licensed under MIT License. 

from circuit import *

class Hybrid(qRAMCircuit):
    """ Hybrid circuits: control first on the output of first k bits, then use
        outputs to control on valid address on the last n - k bits. Perform both
        "tiers" of the circuit in series with no parallelization. """
    def __init__(self, n, q, k):
        super().__init__(n, q, k)

        self.params["name"] = "Hybrid"

        if k < q: # Worst case: 2^k k-controlled, 2^q (n-k+1)-controlled
            self.params["n_qubits"] = n + pow(2,k) + 1 + max(k-1, n-k) + 1
            self.params["depth"] = 2 * ( pow(2,k)*depth(k) + pow(2,q)*depth(n-k+1) ) + 1
            self.params["t_count"] = 2 * ( pow(2,k)*t_c(k) + pow(2,q)*t_c(n-k+1) )
            self.params["t_depth"] = 2 * ( pow(2,k)*t_d(k) + pow(2,q)*t_d(n-k+1) )
            self.params["h_count"] = 2 * ( pow(2,k)*h_c(k) + pow(2,q)*h_c(n-k+1) )
            self.params["cnot_count"] = 2 * (pow(2,k)*cnot_c(k) + pow(2,q)*cnot_c(n-k+1)) + 1
        else: # Worst case: 2^q k-controlled, 2^q (n-k+1)-controlled (no common substrings on the bits)
            self.params["n_qubits"] = n + pow(2,q) + 1 + max(k-1, n-k) + 1
            self.params["depth"] = 2 * ( pow(2,q)*depth(k) + pow(2,q)*depth(n-k+1) ) + 1
            self.params["t_count"] = 2 * ( pow(2,q)*t_c(k) + pow(2,q)*t_c(n-k+1) )
            self.params["t_depth"] = 2 * ( pow(2,q)*t_d(k) + pow(2,q)*t_d(n-k+1) )
            self.params["h_count"] = 2 * ( pow(2,q)*h_c(k) + pow(2,q)*h_c(n-k+1) )
            self.params["cnot_count"] = 2 * (pow(2,q)*cnot_c(k) + pow(2,q)*cnot_c(n-k+1)) + 1

        self.params["cliffords"] = self.params["h_count"] + self.params["cnot_count"]


class Hybrid_Tier1Parallel(qRAMCircuit):
    """ Hybrid circuit, but with first tier done in parallel, and second tier in series. """
    def __init__(self, n, q, k):
        super().__init__(n, q, k)

        self.params["name"] = "Hybrid_Tier1Parallel"

        if k < q:
            self.params["n_qubits"] = (k+1)*pow(2,k) + n - k + 1 + max(n-k, pow(2,k)*(k-1)) + 1
            self.params["depth"] = 2 * ( k + depth(k) + pow(2,q)*depth(n-k+1) ) + 1
            self.params["t_count"] = 2 * ( pow(2,k)*t_c(k) + pow(2,q)*t_c(n-k+1) )
            self.params["t_depth"] = 2 * ( t_d(k) + pow(2,q)*t_d(n-k+1) )
            self.params["h_count"] = 2 * ( pow(2,k)*h_c(k) + pow(2,q)*h_c(n-k+1) )
            self.params["cnot_count"] = 2 * ( k*(pow(2,k)-1) + pow(2,k)*cnot_c(k) + pow(2,q)*cnot_c(n-k+1) ) + 1
        else: # Worst case: 2^q k-controlled, 2^q (n-k+1)-controlled (no common substrings on the bits)
            self.params["n_qubits"] = (k+1)*pow(2,q) + n - k + 1 + max(n-k, pow(2,q)*(k-1)) + 1
            self.params["depth"] = 2 * ( q + depth(k) + pow(2,q)*depth(n-k+1) ) + 1
            self.params["t_count"] = 2 * ( pow(2,q)*t_c(k) + pow(2,q)*t_c(n-k+1) )
            self.params["t_depth"] = 2 * ( t_d(k) + pow(2,q)*t_d(n-k+1) )
            self.params["h_count"] = 2 * ( pow(2,q)*h_c(k) + pow(2,q)*h_c(n-k+1) )
            self.params["cnot_count"] = 2 * ( k*(pow(2,q)-1) + pow(2,q)*cnot_c(k) + pow(2,q)*cnot_c(n-k+1) ) + 1

        self.params["cliffords"] = self.params["h_count"] + self.params["cnot_count"]


class Hybrid_Tier2Parallel(qRAMCircuit):
    """ Hybrid circuit, but with second tier done in parallel, and first tier in series. """
    def __init__(self, n, q, k):
        super().__init__(n, q, k)

        self.params["name"] = "Hybrid_Tier2Parallel"

        if k < q:
            # In the worst case, one of the 2^k outputs must do 2^(q-1) + 1 of the n-k+1-controlled gates,
            # and all the remaining ones must do 1; assume then that one of the registers must perform
            # enough CNOTs to copy down to 2^(q-1) + 1 registers.
            self.params["n_qubits"] = k + (n-k+2)*pow(2,q) + max(k-1, pow(2,q)*(n-k)) + 1
            self.params["depth"] = 2 * (pow(2,k)*depth(k) + q + depth(n-k+1) + q) + 1
            self.params["t_count"] = 2 * ( pow(2,k)*t_c(k) + pow(2,q)*t_c(n-k+1) )
            self.params["t_depth"] = 2 * ( pow(2,k)*t_d(k) + t_d(n-k+1) )
            self.params["h_count"] = 2 * ( pow(2,k)*h_c(k) + pow(2,q)*h_c(n-k+1) )
            # CNOTs for first 2^q k + CNOTs to copy outputs + CNOTs to copy (n-k) in parallel + CNOTs for 2^q (n-k+1) controlled + CNOTs to merge outputs
            self.params["cnot_count"] = 2 * ( pow(2,k)*cnot_c(k) + pow(2,q-1) + (n-k)*(pow(2,q)-1) + pow(2,q)*cnot_c(n-k+1) + pow(2,q) - 1) + 1
        else: # Worst case: 2^q k-controlled, 2^q (n-k+1)-controlled (no common substrings on the bits)
            self.params["n_qubits"] = k + (n-k+2)*pow(2,q) + max(k-1, pow(2,q)*(n-k)) + 1
            self.params["depth"] = 2 * (pow(2,q)*depth(k) + q + depth(n-k+1) + q) + 1
            self.params["t_count"] = 2 * ( pow(2,q)*t_c(k) + pow(2,q)*t_c(n-k+1) )
            self.params["t_depth"] = 2 * ( pow(2,q)*t_d(k) + t_d(n-k+1) )
            self.params["h_count"] = 2 * ( pow(2,q)*h_c(k) + pow(2,q)*h_c(n-k+1) )
            self.params["cnot_count"] = 2 * (  pow(2,q)*cnot_c(k) + pow(2,q-1) + (n-k)*(pow(2,q)-1) + pow(2,q)*cnot_c(n-k+1) + pow(2,q) - 1) + 1

        self.params["cliffords"] = self.params["h_count"] + self.params["cnot_count"]


class Hybrid_Parallel(qRAMCircuit):
    """ Hybrid circuit with both tiers done in parallel (one after the other, as the
        results of the second tier depend on those of the first. """
    def __init__(self, n, q, k):
        super().__init__(n, q, k)

        self.params["name"] = "Hybrid_Parallel"

        if k < q: # Similar worst case as above; one of the 2^k top-tier outputs must be copied down to 2^(q-1) + 1 fresh outputs
            self.params["n_qubits"] = k*pow(2,k) + (n-k+2)*pow(2,q) + max(pow(2,k)*(k-1), pow(2,q)*(n-k)) + 1
            self.params["depth"] = 2 * (q + depth(k) + q + depth(n-k+1) + q) + 1
            self.params["t_count"] = 2 * ( pow(2,k)*t_c(k) + pow(2,q)*t_c(n-k+1) )
            self.params["t_depth"] = 2 * ( t_d(k) + t_d(n-k+1) )
            self.params["h_count"] = 2 * ( pow(2,k)*h_c(k) + pow(2,q)*h_c(n-k+1) )
            self.params["cnot_count"] = 2 * ( k*(pow(2,k)-1) + pow(2,k)*cnot_c(k) + pow(2,q-1) + (n-k)*(pow(2,q)-1) + pow(2,q)*cnot_c(n-k+1) ) + 1
        else: 
            self.params["n_qubits"] = k*pow(2,q) + (n-k+2)*pow(2,q) + max(pow(2,q)*(k-1), pow(2,q)*(n-k)) + 1
            self.params["depth"] = 2 * (q + depth(k) + q + depth(n-k+1) + q) + 1
            self.params["t_count"] = 2 * ( pow(2,q)*t_c(k) + pow(2,q)*t_c(n-k+1) )
            self.params["t_depth"] = 2 * ( t_d(k) + t_d(n-k+1) )
            self.params["h_count"] = 2 * ( pow(2,q)*h_c(k) + pow(2,q)*h_c(n-k+1) )
            self.params["cnot_count"] = 2 * ( k*(pow(2,q)-1) + pow(2,q)*cnot_c(k) + pow(2,q-1) + (n-k)*(pow(2,q)-1) + pow(2,q)*cnot_c(n-k+1) ) + 1

        self.params["cliffords"] = self.params["h_count"] + self.params["cnot_count"]
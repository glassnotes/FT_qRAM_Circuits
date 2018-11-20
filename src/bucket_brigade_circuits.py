# -*- coding: utf-8 -*-                                                            
#                                                                                  
# bucket_brigade_circuits.py: Bucket brigade circuits and parallelized versions. 
#                                                                                  
# © 2018 Olivia Di Matteo (odimatte@uwaterloo.ca)                                  
#                                                                                  
# This file is part of the project FT_qRAM_Circuits.                                      
# Licensed under MIT License. 

from circuit import *

class BucketBrigade(qRAMCircuit):
    """ Vanilla bucket brigade circuit, as presented in Fig. 9 of Srinivasan et al NJP paper.

        Note to self: updated resource counts p. 153 of nb.
    """
    def __init__(self, n):
        super().__init__(n)

        self.params["name"] = "BucketBrigade"

        self.params["n_qubits"] = n + pow(2,n+1) + 5
        self.params["depth"] = 45*pow(2,n) + 2*n - 60
        self.params["t_count"] = 21*pow(2,n) - 28
        self.params["t_depth"] = 3*pow(2,n) - 4
        self.params["h_count"] = 4*pow(2,n) - 6
        self.params["cnot_count"] = 50*pow(2,n) - 66
        self.params["cliffords"] = self.params["h_count"] + self.params["cnot_count"]


class BucketBrigadeParallel(qRAMCircuit):
    """ Parallelize all Toffolis when possible; more logical qubits, but huge savings in depth. 
       
        In this circuit, we use a special technique with the Toffoli target register, by preparing it
        in the superposition state of all even parities. This eliminates the need to uncompute
        the parallel cascade of Toffolis after doing the parity calculation and copying the value down
        to the output qubit, saving us a layer of T-depth and 2^n Toffolis.

        Creating the even parity superposition is done by applying a Hadamard to (n-1) of the qubits,
        and then computing their parity onto the last one. This can be done in logarithmic depth.
        Since this incurs no extra T-count, and as it only needs to be done once when doing 
        multiple sequential queries, we ignore the resources needed for this here as their effect
        is negligible.

        Note to self: updated resource counts p. 155 of nb.
    """
    def __init__(self, n):
        super().__init__(n)
        
        self.params["name"] = "BucketBrigadeParallel"
        
        self.params["n_qubits"] = 8*pow(2,n) 
        self.params["depth"] = pow(n, 2) + 35*n - 15 
        self.params["t_count"] = 21*pow(2,n) - 28 
        self.params["t_depth"] = 2*n - 1
        self.params["h_count"] = 6*pow(2,n) - 8
        self.params["cnot_count"] = 60*pow(2,n) - 2*n - 74
        self.params["cliffords"] = self.params["h_count"] + self.params["cnot_count"]
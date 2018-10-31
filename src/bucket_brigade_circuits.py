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
    """ Vanilla bucket brigade circuit, as presented in Fig. 9 of Srinivasan et al NJP paper."""
    def __init__(self, n):
        super().__init__(n)

        self.params["name"] = "BucketBrigade"

        self.params["n_qubits"] = n + pow(2,n+1) + 6
        self.params["depth"] = 56*pow(2,n) + 2*n - 53
        self.params["t_count"] = 28*(pow(2,n)-1)
        self.params["t_depth"] = 4*(pow(2,n)-1)
        self.params["h_count"] = 4*(pow(2,n)-1)
        self.params["cnot_count"] = 62*pow(2,n) - 59
        self.params["cliffords"] = self.params["h_count"] + self.params["cnot_count"]


class BucketBrigadeParallel(qRAMCircuit):
    """ Parallelize all Toffolis when possible; more logical qubits, but huge savings in depth. """
    def __init__(self, n):
        super().__init__(n)
        
        self.params["name"] = "BucketBrigadeParallel"
        
        self.params["n_qubits"] = 8*pow(2,n) 
        self.params["depth"] = pow(n,2) + 37*n + 3 
        self.params["t_count"] = 28*(pow(2,n)-1) 
        self.params["t_depth"] = 2*n 
        self.params["h_count"] = 8*(pow(2,n)-1) 
        self.params["cnot_count"] = 78*pow(2,n) - 2*n - 75
        self.params["cliffords"] = self.params["h_count"] + self.params["cnot_count"]


class BucketBrigadeParallelEvenParities(qRAMCircuit):
    """ In this circuit, we do a special trick with the Toffoli target register, by preparing it
        in the superposition state of all even parities. This eliminates the need to uncompute
        the parallel cascade of Toffolis after doing the parity calculation and copying the value down
        to the output qubit, saving us a layer of T-depth and 2^n Toffolis.

        However, creating the superposition of all even parities incurs an exponential blowup
        in circuit depth (not T-depth). If we are doing exponentially many queries, or working
        with a small n, this is not a big deal. It may become significant for large n, though, and 
        maybe not be worthwhile when performing only one or two queries.
    """
    def __init__(self, n):
        super().__init__(n)
        
        self.params["name"] = "BucketBrigadeParallel"
        
        self.params["n_qubits"] = 8*pow(2,n) 
        self.params["depth"] = pow(2, n) + 0.5*(pow(n,2) + 13*n) + 3 
        self.params["t_count"] = 21*pow(2,n) - 28 
        self.params["t_depth"] = 2*n - 1
        self.params["h_count"] = 7*pow(2,n) - 9 
        self.params["cnot_count"] = 61*pow(2,n) - 2*n - 75
        self.params["cliffords"] = self.params["h_count"] + self.params["cnot_count"]
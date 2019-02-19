# -*- coding: utf-8 -*-                                                            
#                                                                                  
# lks_circuits.py: Improved qRAM circuits presented in https://arxiv.org/pdf/1812.00954.pdf
#                  (Low, Kliuchnikov, Schaeffer 2018). Improves the bounds on our basic 
#                  highly-sequential circuit, but only gives big-O notation instead of
#                  actual quantities. It also does not provide Clifford counts, but we typically
#                  assume these will not affect the overall cost much anyways.
#                                                                                  
# © 2019 Olivia Di Matteo (odimatteo@triumf.ca)                                  
#                                                                                  
# This file is part of the project FT_qRAM_Circuits.                                      
# Licensed under MIT License. 

from numpy import ceil, log

from circuit import *

class SelSwap(qRAMCircuit):
    """ Line 3 of Table II. Notation is a little bit different, but essentially, 
        we will feed it n such that N = 2^n is the number of slots in the memory, 
        lambda (l) provides a space-depth tradeoff (chosen from [1, N]), and b is the size
        of the number returned (for consistency with our circuits, we choose b = 1 since
        we are always returning only a single bit, but I guess this could be higher too.)
    """
    def __init__(self, n, l, b):
        super().__init__(n)

        self.params["name"] = "SelSwap"

        self.params["n_qubits"] = b * l + 2 * n 
        self.params["t_count"] = 4 * ceil((2 ** n) / l) + 8 * b * l
        self.params["t_depth"] = (2 ** n) / l + log(l)
        self.params["cliffords"] = 0
       

class SelSwapDirty(qRAMCircuit):
    """ Like the other SelSwap circuit but this one allows b * l of the input
        qubits to be 'dirty', and the circuit returns them to their initial state
        after the query. This is nice for alternating algorithm/qRAM query circuits
        because we don't have to add as many ancillae, we can just use the idle 
        algorithm qubits. It seems like here there is a small tradeoff in that we are
        adding extra qubits and T-count, but the fact that they can be dirty makes up
        for it.
    """
    def __init__(self, n, l, b):
        super().__init__(n)

        self.params["name"] = "SelSwapDirty"

        self.params["n_qubits"] = b * (l + 1) + 2 * n 
        self.params["t_count"] = 8 * ceil((2 ** n) / l) + 32 * b * l
        self.params["t_depth"] = (2 ** n) / l + log(l)
        self.params["cliffords"] = 0


class SelectV(qRAMCircuit):
    """ Based on the circuits in Appendix G4 of https://arxiv.org/pdf/1711.10980.pdf.
        This is a more specialized circuit that is meant to traverse all possible 
        configurations of the address variables and apply a different unitary to each.
        In our case, this would just be a NOT, or an identity, I suppose, depending on
        the values stored in memory. 

        Note: currently unclear what the T-depth of this guy is. Based on the image at
        the bottom of p.47, it would seem like the T-depth is the same as the T-count,
        because they are all applied in sequence. Need to check this further.
    """
    def __init__(self, n):
        super().__init__(n)
    
        self.params["name"] = "SelectV"

        self.params["n_qubits"] = n + (n - 1)  # n - 1 ancillae
        self.params["depth"] = -1 
        self.params["t_count"] = 7.5*pow(2,n) + 6*n - 28 
        self.params["t_depth"] = 7.5*pow(2,n) + 6*n - 28 
        self.params["h_count"] = 2*pow(2,n) + 2*n - 8 
        self.params["cnot_count"] = 7.5*pow(2,n) + 6*n - 26
        self.params["cliffords"] = self.params["h_count"] + self.params["cnot_count"]
# -*- coding: utf-8 -*-                                                            
#                                                                                  
# circuit.py: Parent classes for circuits.
#                                                                                  
# © 2018 Olivia Di Matteo (odimatte@uwaterloo.ca)                                  
#                                                                                  
# This file is part of the project FT_qRAM_Circuits.                                      
# Licensed under MIT License. 

# Resources for an MPMCT with n controls
def depth(n):
    return 28*n - 60

def t_c(n):
    return 12*n - 20

def t_d(n):
    return 4*(n-2)

def h_c(n):
    return 4*n - 6

def cnot_c(n):
    return 24*n - 40

class Circuit():
    """ Plain old circuit. We can initialize with gate counts only.
        Parameter names should be sufficiently descriptive. 
        Depth refers to the *full* depth of the circuit, both Clifford and T.
    """
    def __init__(self, n_qubits, t_count, t_depth, h_count, cnot_count, depth):
        # Store all variables in a dictionary, this makes it easier to dump to a CSV after. 
        self.params = {}

        # Logical qubit and gate counts
        self.params["n_qubits"] = n_qubits
        self.params["t_count"] = t_count
        self.params["t_depth"] = t_depth
        self.params["h_count"] = h_count
        self.params["cnot_count"] = cnot_count
        self.params["cliffords"] = h_count + cnot_count # h_count + cnot_count
        self.params["depth"] = depth

class qRAMCircuit(Circuit):
    def __init__(self, n, q=0, k=0, b1=0, b2=0):
        """ Parent circuit class for a qRAM circuit.

            :param int n: The number of address bits. Used for all circuits.
            :param int q: Indicates that there are 2^q 1s stored in the memory. 
                Used for all circuits except bucket brigade.
            :param int k: Hybrid splitting parameter. Does a set of multi-controlled
                gates conditioned on the first k, then next the last n - k. Valid only
                for Hybrid and Double qRAM circuits.
            :param int b1: For Double qRAM circuits only. Indicates the number of addresses
                (2^b1) in the first half of a Cartesian product address space.
            :param int b2: Same as b1, but for the second half of the address space.
        """
        # Initialize gate counts as empty for now, we will fill them in the child classes.
        super().__init__(0, 0, 0, 0, 0, 0)

        # Store the type of circuit; will be set by child classes
        self.params["name"] = "" 

        # Circuit specification parameters
        self.params["n"] = n
        self.params["q"] = q
        self.params["k"] = k
        self.params["b1"] = b1
        self.params["b2"] = b2
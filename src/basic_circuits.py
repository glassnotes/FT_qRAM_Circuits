from circuit import *

class LargeWidthSmallDepth(qRAMCircuit):
    """ Parallelize all the MPMCTs at the cost of a lot of extra logical qubits,
        but significantly decreased circuit/T depth """
    def __init__(self, n, q):
        super().__init__(n, q)

        self.params["name"] = "LargeWidthSmallDepth"

        self.params["n_qubits"] = n*pow(2,q+1) + 1
        self.params["depth"] = q + depth(n) + q + 2 + q 
        self.params["t_count"] = pow(2,q)*t_c(n)
        self.params["t_depth"] = t_d(n)
        self.params["h_count"] = pow(2,q)*h_c(n)
        self.params["cnot_count"] = 2*(n+1)*(pow(2,q)-1) + pow(2,q)*cnot_c(n) + 2 
        self.params["cliffords"] = self.params["h_count"] + self.params["cnot_count"]


class SmallWidthLargeDepth(qRAMCircuit):
    """ Run all the MPMCTs in sequence; few logical qubits, but massive depth. """
    def __init__(self, n, q):
        super().__init__(n, q)
        
        self.params["name"] = "SmallWidthLargeDepth"

        self.params["n_qubits"] = 2*n 
        self.params["depth"] = pow(2,q)*depth(n) 
        self.params["t_count"] = pow(2,q)*t_c(n)
        self.params["t_depth"] = pow(2,q)*t_d(n)
        self.params["h_count"] = pow(2,q)*h_c(n)
        self.params["cnot_count"] = pow(2,q)*cnot_c(n) 
        self.params["cliffords"] = self.params["h_count"] + self.params["cnot_count"]
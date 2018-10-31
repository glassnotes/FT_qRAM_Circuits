# -*- coding: utf-8 -*-                                                            
#                                                                                  
# surface_code.py: Surface code class and resource estimation. 
#                                                                                  
# © 2018 Olivia Di Matteo (odimatte@uwaterloo.ca)                                  
#                                                                                  
# This file is part of the project FT_qRAM_Circuits.                                      
# Licensed under MIT License.                                                      

import numpy as np

class SurfaceCode():
    def __init__(self, p_in=1e-4, p_g=1e-5, t_sc=200e-9):
        """ Initial a surface code with given physical parameters.

            :param float p_in: Input failure probability of magic states
            :param float p_g: Default failure probability of a gate
            :param float t_sc: Time for a surface code cycle.

        """
        self.params = { 
                    "p_in" : p_in, 
                    "p_g" : p_g,
                    "t_sc" : t_sc
                    }

    def defect_compute_resources(self, circuit, eps=1.0, num_factories="t_width"):
        """ Performs resource estimates using defect-based surface codes in the style of
                - Our old SHA analysis (https://link.springer.com/chapter/10.1007/978-3-319-69453-5_18)
                - Austin's PRA from 2012 (https://link.aps.org/doi/10.1103/PhysRevA.86.032324)

            :param float eps: Specifies ratio of distillation / logical error. Default 1 meaning they are the same.
            :param int num_factories: Number of factories to perform the estimate with.
                Default value of "t_width" will perform an estimate where the number of factories
                is computed according to a 'T-width', T_w = self.T_c / T_d, which is an *assumed* average
                of how many T-gates we need to perform in any given layer of T-depth.
        """
        if num_factories != "t_width":
            return NotImplemented

        resources = { # The set of values we will compute
                "distill_reqd" : True,
                "distill_distances" : [],
                "distill_factories" : 0,
                "distill_sim_states" : 0,
                "distill_cycles_per_state" : 0,
                "distill_total_cycles" : 0,
                "distill_logical_q_per_factory" : 0,
                "distill_logical_q" : 0,
                "distill_phys_q" : 0,
                "distill_time" : 0,
                "clifford_distance" : 0,
                "clifford_phys_q" : 0,
                "clifford_time" : 0,
                "clifford_total_avg_cycles" : 0,
                "total_time" : 0,
                "total_cost" : 0,
                "safety_factor" : 0
              }

        resources["eps"] = eps

        # First check if distillation is required and then find the distances
        p_out = 1. / circuit.params["t_count"] # Maximum failure prob of the state distillation
        d = [] # Distances of the code

        p = p_out

        if p >= self.params["p_in"]:
            resources["distill_reqd"] = False
        else:
            while p <= self.params["p_in"]:
                d_i = 1
                p_i = p

                while True:
                    if 192*d_i*pow(100*self.params["p_g"], (d_i+1.)/2) < (resources["eps"]*p_i) / (1 + resources["eps"]):
                        break
                    else:
                        d_i += 1

                d.append(d_i)

                p = pow(p_i / (35 * (1 + resources["eps"])), 1./3)

        if resources["distill_reqd"]:
            resources["distill_distances"] = d

            # Number of distillation cycles is the sum of distances * 10
            resources["distill_cycles_per_state"] = sum(d) * 10

            # Total number of cycles is the number of cycles for a layer of Ts times number of layers
            resources["distill_total_cycles"] = resources["distill_cycles_per_state"] * circuit.params["t_depth"]

            # Number of input states consumed - 16 for first layer, 16 x 15 for second
            # layer (since the 16th qubit can be reused from the previous layer, etc.
            resources["distill_logical_q_per_factory"] = 16 * pow(15, len(d) - 1)

            # The `bottom' layer of the distillation circuit has the largest footprint
            # i.e. number of physical qubits
            bottom_footprint = resources["distill_logical_q_per_factory"] * np.ceil(2.5 * 1.25 * pow(d[-1], 2))

            # Physical qubits for the middle footprint, if applicable
            # The number of magic states that we can produce simultaneous is the ratio of
            # the layers; if there is only 1 layer, we can only produce 1 state at a time.
            if len(d) > 1:
                middle_footprint = (resources["distill_logical_q_per_factory"]/ 15.) * np.ceil(2.5 * 1.25 * pow(d[-2], 2))
                resources["distill_sim_states"] = np.floor(bottom_footprint / middle_footprint)
            else:
                resources["distill_sim_states"] = 1

            # Number of physical qubits is the number of qubits in the bottom-most layer
            resources["distill_phys_q"] = bottom_footprint

            # Compute the T-width, i.e. the number of T gates per layer of T-depth
            # If we need to do more per layer than we can produce in the distillery,
            # we need to add additional distilleries
            T_w = np.ceil(circuit.params["t_count"] * 1.0 / circuit.params["t_depth"])
            if T_w > resources["distill_sim_states"]:
                resources["distill_factories"] = int(np.ceil(T_w / resources["distill_sim_states"]))
            else:
                resources["distill_factories"] = 1
            resources["distill_phys_q"] *= resources["distill_factories"]


            resources["distill_logical_q"] = resources["distill_factories"] * resources["distill_logical_q_per_factory"]

            # Compute the total amount of time it will take to distill all needed magic states
            resources["distill_time"] = self.params["t_sc"] * resources["distill_cycles_per_state"] * circuit.params["t_count"] / (resources["distill_sim_states"] * resources["distill_factories"])

        # -----------------------------
        # For the Clifford surface code
        # -----------------------------
        p_cliff = 1. / circuit.params["cliffords"]
        d_cliff = 1
        while pow(self.params["p_in"] / 0.0125, (d_cliff + 1.)/2) > p_cliff:
            d_cliff += 1
        resources["clifford_distance"] = d_cliff
        resources["clifford_phys_q"] = circuit.params["n_qubits"] * np.ceil(2.5 * 1.25 * (d_cliff ** 2))

        # Compute the number of cycles required to implement the Clifford portion
        # CNOTs take 2 cycles, Hadamards d cycles

        # First compute the weighted average number of cycles based on the number of
        # CNOTs and hadamards per layer
        f_cnot = circuit.params["cnot_count"] * 1./ circuit.params["cliffords"] # Fraction of CNOTs
        f_h = circuit.params["h_count"] * 1./ circuit.params["cliffords"] # Fraction of Hadamards
        cliffords_per_depth = (1. * circuit.params["cliffords"]) / (circuit.params["n_qubits"] * circuit.params["depth"]) # Cliffords per qubit per layer

        # Average cycles per layer of depth
        avg_cycles = cliffords_per_depth * (f_cnot * 2 + f_h * d_cliff)
        resources["clifford_total_avg_cycles"] = avg_cycles * circuit.params["depth"]
        resources["clifford_time"] = self.params["t_sc"] * resources["clifford_total_avg_cycles"]

        # Now compute the total costs
        total_logical_qubits = circuit.params["n_qubits"] + resources["distill_logical_q"]
        total_cycles = resources["distill_total_cycles"] +resources["clifford_total_avg_cycles"]

        resources["total_cost"] = np.log2(total_logical_qubits * total_cycles)
        resources["total_time"] = max(resources["clifford_time"], resources["distill_time"])
        resources["total_phys_q"] = resources["distill_phys_q"] + resources["clifford_phys_q"]

        # Return EVERYTHING
        return {"mode" : "defect", **resources, **circuit.params, **self.params}


    def lattice_compute_resources(self, circuit, num_factories=1, safety_factor=99, d1=15, d2=23):
        return NotImplemented
        

    def out(self):
        """ Somewhat pretty output of relevant surface code parameters.
        """
        print("Distillation parameters:")
        print("--------------------------------")
        if self.params["distill_reqd"]:
            print("eps: " + str(self.params["eps"]))
            print("p_in: {:.4g}".format(self.params["p_in"]))
            print("p_g: {:.4g}".format(self.params["p_g"]))
            print("Distance(s): " + str(self.params["distill_distances"]))
            print("Logical qubits per factory: " + str(self.params["distill_logical_q_per_factory"]))
            print("Simultaneous magic states per factory : " + str(self.params["distill_sim_states"]))
            print("Number of factories: " + str(self.params["distill_factories"]))
            print("Total logical qubits: {:.4g}".format(self.params["distill_logical_q"]))
            print("Cycles per set of magic states: " + str(self.params["distill_cycles_per_state"]))
            print("Total cycles: {:.4g}".format(self.params["distill_total_cycles"]))
            print("Physical qubits required: {:.4g}".format(self.params["distill_phys_q"]))
            print("Time for distillation: {:.4g} s".format(self.params["distill_time"]))
            print("Time for distillation: {:.4g} y".format(self.params["distill_time"] / (60 * 60 * 24 * 356)))
        else:
            print("No distillation required.")

        print("")
        print("Clifford surface code parameters")
        print("--------------------------------")
        print("Distance: " + str(self.params["clifford_distance"]))
        print("Physical qubits: {:.4g}".format(self.params["clifford_phys_q"]))
        print("Total avg cycles: {:.4g}".format(self.params["clifford_total_avg_cycles"]))
        print("Clifford time: {:.4g} s".format(self.params["clifford_time"]))
        print("Clifford time: {:.4g} y".format(self.params["clifford_time"] / (60 * 60 * 24 * 356)))
        print()
        print("Totals:")
        print("--------------------------------")
        print("Total time: {:.4g} s".format(self.params["total_time"]))
        print("Total time: {:.4g} y".format(self.params["total_time"] / (60 * 60 * 24 * 356)))
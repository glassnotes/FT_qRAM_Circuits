# -*- coding: utf-8 -*-                                                            
#                                                                                  
# surface_code.py: Surface code class and resource estimation. 
#                                                                                  
# © 2019 Olivia Di Matteo (odimatte@uwaterloo.ca)                                  
#                                                                                  
# This file is part of the project FT_qRAM_Circuits.                                      
# Licensed under MIT License.
#
# Routines for lattice surgery were provided by Austin Fowler, and slightly
# reworked to fit within the SurfaceCode class.                                                     

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
            bottom_footprint = resources["distill_logical_q_per_factory"] * np.ceil(2.5 * 1.25 * pow(2*d[-1], 2))

            # Physical qubits for the middle footprint, if applicable
            # The number of magic states that we can produce simultaneous is the ratio of
            # the layers; if there is only 1 layer, we can only produce 1 state at a time.
            if len(d) > 1:
                middle_footprint = (resources["distill_logical_q_per_factory"]/ 15.) * np.ceil(2.5 * 1.25 * pow(2*d[-2], 2))
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
        if circuit.params["cliffords"] > 0:
            p_cliff = 1. / circuit.params["cliffords"]
            d_cliff = 1
            while pow(self.params["p_in"] / 0.0125, (d_cliff + 1.)/2) > p_cliff:
                d_cliff += 1
            resources["clifford_distance"] = d_cliff
            resources["clifford_phys_q"] = circuit.params["n_qubits"] * np.ceil(2.5 * 1.25 * pow(2*d_cliff, 2))

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
        total_cycles = resources["distill_total_cycles"] + resources["clifford_total_avg_cycles"]

        if total_cycles > 0: # No distillation, and no Clifford counts provided
            resources["total_cost"] = np.log2(total_logical_qubits * total_cycles)
            resources["total_time"] = max(resources["clifford_time"], resources["distill_time"])
            resources["total_phys_q"] = resources["distill_phys_q"] + resources["clifford_phys_q"]
        else:
            resources["total_cost"] = 0
            resources["total_time"] = 0
            resources["total_phys_q"] = 0

        # Return EVERYTHING
        return {"mode" : "defect", **resources, **circuit.params, **self.params}


    def lattice_compute_resources(self, circuit):
        """ 
            The contents of this method were generously provided in some scripts by Austin Fowler.
            I have adjusted them slightly to fit within the structure of the rest of my code (noted). 
            Thanks, Austin!!

            This routine currently allows for only a single distillation factory; it should therefore
            *not* be compared directly to the defect-based case, unless defect-based estimation is done
            with num_factories = 1, but this is currently not implemented.

            (Note: currently hasn't been cross-checked against other estimates)
        """
        resources = {}

        def p_logical(p, d):
            return 0.1 * ((100 * p) ** ((d + 1) // 2))

        def distance(p, p_L):
            d = 3
            while p_logical(p, d) > p_L:
                d += 2
            return d

        def distillation_15_1_p_out(p):
            return 35 * (p ** 3)

        def distillation_8_2_p_out(p):
            return 28 * (p ** 2)

        # (ODM) This internal method is the one that actually does the resource counts.
        def overhead(T_count, Toffoli_count, qubit_count, error_budget, p, surface_code_cycle_time, routing_overhead,
                     L0_distillation_code_distance, L1_distillation_code_distance, L2_distillation_code_distance):    
            # The same value as the characteristic gate error rate, by using the techniques of Ying Li (https://arxiv.org/abs/1410.7808).
            L0_state_injection_error = p
            
            # The chance of a logical error occurring within a lattice surgery unit cell at L0 distillation code distance.
            L0_topological_error_per_unit_cell = p_logical(p, L0_distillation_code_distance)
            
            # It takes approximately 100 L0 unit cells to get the injected state where it needs to be and perform the T gate.
            L0_topological_error_from_T0_gate = 100 * L0_topological_error_per_unit_cell
            
            # Chance of failure of a T gate performed with an injected T state.
            L0_total_T_error = L0_state_injection_error + L0_topological_error_from_T0_gate

            # The chance of a logical error occurring within a lattice surgery unit cell at L1 distillation code distance.
            L1_topological_error_per_unit_cell = p_logical(p, L1_distillation_code_distance)
            
            # The L1 T factory uses approximately 1000 L1 unit cells.
            L1_topological_error_from_factory = 1000 * L1_topological_error_per_unit_cell
            
            # It takes approx. 100 L1 unit cells to get the L1 state produced by the factory to where it needs to be and perform the T gate.
            L1_topological_error_from_T_gate = 100 * L1_topological_error_per_unit_cell
            
            # The chance that undetected errors in the L0 T states input into the L1 factory result in an undetected error in an L1 T state being output.
            L1_distillation_error = distillation_15_1_p_out(L0_total_T_error)
            
            # Chance of failure of a T gate performed with a T state produced from the L1 factory.
            L1_total_T_error = L1_topological_error_from_factory + L1_topological_error_from_T_gate + L1_distillation_error

            # If error already low enough, do no more distillation, do no catalyzation.
            if L1_total_T_error < error_budget / (T_count + 4 * Toffoli_count):
                # Total footprint of the L2 factory (either CCZ or catalyzed T)
                physical_qubits_for_factories = 4 * 8 * 2 * (L1_distillation_code_distance ** 2)

                # Total number of rounds of error correction needed to produce all the magic states the algorithm needs with a single factory.
                cycles_spent_producing_magic_states = 5.75 * (T_count + 4 * Toffoli_count) * L1_distillation_code_distance
            
                # Chance any factory anywhere in the algorithm fails.
                total_error_from_distillation = L1_total_T_error * (T_count + 4 * Toffoli_count)

            else:
                # The chance of a logical error occurring within a lattice surgery unit cell at L2 distillation code distance.
                L2_topological_error_per_unit_cell = p_logical(p, L2_distillation_code_distance)
            
                # The L2 CCZ factory and catalyzed T factory both use approximately 1000 L2 unit cells.
                L2_topological_error_from_factory = 1000 * L2_topological_error_per_unit_cell
            
                # The chance that undetected errors in the L1 T states input into the L2 factory result in an undetected error in an L2 T state being output.
                L2_distillation_error = distillation_8_2_p_out(L1_total_T_error)
            
                # Chance of failure of a CCZ or pair of Ts performed with a CCZ state or pair of T states produced by the L2 factory.
                L2_total_CCZ_or_2T_error = L2_topological_error_from_factory + L2_distillation_error
            
                if L2_total_CCZ_or_2T_error < error_budget / (((T_count + 1) // 2) + Toffoli_count):
                    # Each T needs 1 T state, a single CCZ state can be catalyzed into 2 T states.
                    catalyzations_required = (T_count + 1) // 2
            
                    # Number of CCZ states that must be produced by the CCZ factory (including CCZ states consumed to run catalysis).
                    CCZ_states_required = catalyzations_required + Toffoli_count
            
                    # Total footprint of the L2 factory (either CCZ or catalyzed T)
                    physical_qubits_for_factories = 6 * (4 * 8 * 2 * (L1_distillation_code_distance ** 2)) + 4 * 8 * 2 * (L2_distillation_code_distance ** 2)
            
                    # Determine whether the L1 or L2 code distance limits factory execution time.
                    characteristic_temporal_distance = max(2 * L1_distillation_code_distance + 1, L2_distillation_code_distance)
            
                    # Total number of rounds of error correction needed to produce all the magic states the algorithm needs with a single factory.
                    cycles_spent_producing_magic_states = (5 * CCZ_states_required + 1 * catalyzations_required) * characteristic_temporal_distance

                    # Chance any factory anywhere in the algorithm fails.
                    total_error_from_distillation = L2_total_CCZ_or_2T_error * CCZ_states_required
                
                else:
                    # The chance that undetected errors in the L1 T states input into the L2 factory result in an undetected error in an L2 T state being output.
                    L2_distillation_error = distillation_15_1_p_out(L1_total_T_error)
                
                    # Chance of failure of a T performed with T states produced by the L2 factory.
                    L2_total_T_error = L2_topological_error_from_factory + L2_distillation_error
            
                    # Total footprint of the L2 factory (either CCZ or catalyzed T)
                    physical_qubits_for_factories = 8 * (4 * 8 * 2 * (L1_distillation_code_distance ** 2)) + 4 * 8 * 2 * (L2_distillation_code_distance ** 2)

                    # Determine whether the L1 or L2 code distance limits factory execution time.
                    characteristic_temporal_distance = max(2 * L1_distillation_code_distance + 1, L2_distillation_code_distance)

                    # Total number of rounds of error correction needed to produce all the magic states the algorithm needs with a single factory.
                    cycles_spent_producing_magic_states = 5.75 * (T_count + 4 * Toffoli_count) * characteristic_temporal_distance
            
                    # Chance any factory anywhere in the algorithm fails.
                    total_error_from_distillation = L2_total_T_error * (T_count + 4 * Toffoli_count)

            # Maximum tolerable error from data qubits to still be under the error budget.
            remaining_error_budget_for_data = max(0, error_budget - total_error_from_distillation)

            # Algorithm execution rounds times the number of data qubits and their communication channels.
            total_data_unit_cells = (1 + routing_overhead) * qubit_count * cycles_spent_producing_magic_states
            
            # Safe target error rate per data round.
            target_error_per_data_round = remaining_error_budget_for_data / total_data_unit_cells
            
            # Code distance required to achieve the safe target error rate per data round.
            data_code_distance = distance(p, target_error_per_data_round)
            
            # The chance of a logical error occurring within a lattice surgery unit cell at the data code distance.
            topological_error_per_data_unit_cell = p_logical(p, data_code_distance)
            
            # Total number of data qubits, including communication channels.
            physical_qubits_for_data = (1 + routing_overhead) * qubit_count * 2 * (data_code_distance ** 2)
            
            # Chance of any data logical qubit suffering a logical error anywhere in the algorithm.
            total_data_error = total_data_unit_cells * topological_error_per_data_unit_cell

            # Chance of any logical error anywhere in the algorithm.
            total_error = total_error_from_distillation + total_data_error
            
            # Test whether algorithm can be implemented with the chosen two levels of distillation, current value of p, and target error_budget. If no, try
            # increasing the error budget or lowering the physical gate error rate p.
            under_budget = (total_error < error_budget)

            # Total physical qubits required to run algorithm.
            total_physical_qubits = physical_qubits_for_factories + physical_qubits_for_data 
            
            # Total time in seconds to run algorithm, assuming no routing or Clifford bottlenecks.
            execution_time = cycles_spent_producing_magic_states * surface_code_cycle_time

            if under_budget:
                return_values = [L1_distillation_code_distance, L2_distillation_code_distance, total_error, physical_qubits_for_factories,
                                 physical_qubits_for_data, total_physical_qubits, execution_time]
            else:
                return_values = [0, 0, total_error, 0, 0, 0, 0]

            return return_values


        # (ODM) The code below is what was originally a method called 'algorithm_metrics'; it called 
        # the overhead function. I've adjusted it to work with the parameters from the input circuit.

        # Acceptable chance of a top-level distillation failure or a data storage failure occurring at any point during the algorithm.
        error_budget = 0.1
        
        # Gate error rate that sets how quickly logical errors are suppressed. 0.1% means a factor of 10 suppression with each increase of d by 2.
        # (ODM) : using my own parameters here and below
        p = self.params["p_g"]
        
        # Time in seconds required to execute a single round of the surface code circuit detecting errors.
        surface_code_cycle_time = self.params["t_sc"]
        
        # Should be 50% (0.5) unless you know what you are doing. Additional space needed for moving magic states and data qubits around in order to perform operations.
        routing_overhead = 0.5

        # Should be 7 unless you know what you are doing. Code distance used for state injection.
        L0_distillation_code_distance = 7
        
        # Min should be 15 unless you know what you are doing. Code distance used for level 1 factories. Max can be arbitraryily large > 15, but shouldn't need to be more than 100.
        min_L1_distillation_code_distance = 15
        max_L1_distillation_code_distance = 51
        
        # Min should be 15 unless you know what you are doing. Code distance used for level 2 factories. Max can be arbitraryily large > 15, but shouldn't need to be more than 100.
        min_L2_distillation_code_distance = 15
        max_L2_distillation_code_distance = 51

        results = [0, 0, 0, 0, 0, 0, 0]
        
        for L1_distillation_code_distance in range(min_L1_distillation_code_distance, max_L1_distillation_code_distance + 1, 2):
            for L2_distillation_code_distance in range(min_L2_distillation_code_distance, max_L2_distillation_code_distance + 1, 2):
                # (ODM) Important note: I have expressed everything explicitly in terms of Ts and not Toffolis. Therefore I have set
                # the Toffoli count to 0; I have, however, kept it in the method in case this gets changed later.
                new_results = overhead(circuit.params["t_count"], 0, circuit.params["n_qubits"], error_budget, p, surface_code_cycle_time, routing_overhead,
                                       L0_distillation_code_distance, L1_distillation_code_distance,  L2_distillation_code_distance)
                if results[0] == 0:
                    results = new_results
                elif new_results[0] != 0 and results[-2] > new_results[-2]:
                    results = new_results

        # (ODM) make it a dictionary
        resources["distill_distances"] = [results[1], results[0]]
        resources["total_error"] = results[2]
        resources["distill_phys_q"] = results[3]
        resources["clifford_phys_q"] = results[4]
        resources["total_phys_q"] = results[5]
        resources["total_time"] = results[-1]

        return {"mode" : "defect", **resources, **circuit.params, **self.params}


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

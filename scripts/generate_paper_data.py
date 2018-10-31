# -*- coding: utf-8 -*-                                                            
#                                                                                  
# generate_paper_data.py: Generate data for all circuit families used in arXiv:[TBD]. 
#                                                                                  
# © 2018 Olivia Di Matteo (odimatte@uwaterloo.ca)                                  
#                                                                                  
# This file is part of the project FT_qRAM_Circuits.                                      
# Licensed under MIT License. 

import numpy as np
import csv
import sys

sys.path.append("../src")

from surface_code import SurfaceCode
from bucket_brigade_circuits import *
from basic_circuits import *
from hybrid_circuits import *
from double_qram_circuits import *

# Initialize the surface code
sc = SurfaceCode()

# Grab a sample header
keys = sc.defect_compute_resources(BucketBrigade(30)).keys()

# Make a giant CSV file of everything
# Run for a range of n, range of q
with open("double_hybrid_bucket_basic_circuits.csv", 'w') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=keys)
    writer.writeheader()

    for n in range(15, 37):
        # Bucket brigade circuits; depend only on n
        writer.writerow(sc.defect_compute_resources(BucketBrigade(n)))
        writer.writerow(sc.defect_compute_resources(BucketBrigadeParallel(n)))

        # Basic circuits depend on n and memory fullness q
        for q in range(n - 10, n):
            writer.writerow(sc.defect_compute_resources(LargeWidthSmallDepth(n, q)))
            writer.writerow(sc.defect_compute_resources(SmallWidthLargeDepth(n, q)))

            # MPMCT stuff only works for 4 controls or higher
            for k in range(4, n - 3): 
                if k >= q:
                    continue

                writer.writerow(sc.defect_compute_resources(Hybrid(n, q, k)))
                writer.writerow(sc.defect_compute_resources(Hybrid_Tier1Parallel(n, q, k)))
                writer.writerow(sc.defect_compute_resources(Hybrid_Tier2Parallel(n, q, k)))
                writer.writerow(sc.defect_compute_resources(Hybrid_Parallel(n, q, k)))
                
                # Choose b1, b2 for hybrid parallel
                for b1 in range(0, int(q/2) + 5): 
                    b2 = q - b1

                    writer.writerow(sc.defect_compute_resources(Double_QRAM(n, q, k, b1, b2)))
                    writer.writerow(sc.defect_compute_resources(Double_QRAM_Tier1Parallel(n, q, k, b1, b2)))
                    writer.writerow(sc.defect_compute_resources(Double_QRAM_Tier2Parallel(n, q, k, b1, b2)))
                    writer.writerow(sc.defect_compute_resources(Double_QRAM_Parallel(n, q, k, b1, b2)))
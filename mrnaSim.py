import os
import sys
print(sys.version)
import random
import numpy as np
import math
from collections import defaultdict
from collections import namedtuple
from scipy.stats import poisson
from scipy.stats import nbinom
from scipy.stats import norm
from scipy.stats import logistic
import matplotlib.pyplot as plt
from mrna_expression import mrnaDict
from mrna_expression  import MasterMRNA
from CellDef import Cell
from MRNA import MRNA

#
# This is a simulation of mRNA synthesis and degradation. 
# mRNA will start to be synthesized at a given time, and will be degraded at a constant rate.
# The simulation will run for a given number of time steps, and the mRNA count will be recorded at each time step.
# the mRNA will start to be recorded 10 days before the final day of the simulation and will be sampled at 15% to match 10x Genomics RNA-seq data.
# 
# The simulation will start with a given number of HSC cells.
# HSC cells will divide and differentiate into other cell types. The total number of HSC cells will be fairly constant
# HSC cells will be differentiated into other cell types in a stochastic manner.
# 
class Simulation:
    def __init__(self, num_cells=1000, num_days=30, degradation_rate=0.1, synthesis_rate=0.5):
        self.num_cells = num_cells
        self.num_days = num_days
        self.degradation_rate = degradation_rate
        self.synthesis_rate = synthesis_rate
        self.mRNA_counts = defaultdict(list)
        self.time_steps = np.arange(0, num_days + 1)

    def run(self):
        for day in self.time_steps:
            mRNA_count = self.simulate_mRNA(day)
            self.mRNA_counts[day].append(mRNA_count)

    def simulate_mRNA(self, day):
        # Simulate mRNA synthesis and degradation
        synthesized_mRNA = poisson.rvs(self.synthesis_rate * self.num_cells)
        degraded_mRNA = nbinom.rvs(1, 1 - np.exp(-self.degradation_rate)) * self.num_cells
        return max(0, synthesized_mRNA - degraded_mRNA) 



                    

cells = []
for i in range(100):
    cell = Cell(i, "stem")
    cells.append(cell)             
#
# simulate in 6 hour intervals
#
days = 100
for simCount in range(0,days * 24,6):
    day = int(simCount / 24)  # Convert hour to day (6 hours per step)
    hour = simCount % 24  # Get the hour of the day
    newCells = []
    for cell in cells:
        d = cell.sim(time=simCount, simMRNA=(day >= 85))
        # daughter if any
        if (d is not None):
            if d.cellType != "mature":
                newCells.append(d)
            else:
                del(d)
        # starting cell
        if cell.cellType != 'mature':
            newCells.append(cell)
        else:
            del(cell)
    cells = newCells
    
    if hour == 0:
        print(f"Day {day} hour {hour}: Total cells = {len(cells)}")
        d = defaultdict(int)
        #
        # summarize cell types
        #
        for cell in cells:
            d[cell.cellType] += 1
        sd = sorted(d.items(), key=lambda x: x[1], reverse=True)
        for cellType, count in sd:
            if count > 0:
                print(f"Cell Type: {cellType}, Count: {count}") 
#
# quick and dirty way to print the mRNA counts for each cell
#
# for cell in cells:
#     if cell.cellType == "my1":
#         d = defaultdict(int)
#         for m in cell.mrna.mrna:
#             d[m[1]] += 1
#         sd = sorted(d.items(), key=lambda x: x[1], reverse=True)
#         print(f"Cell ID: {cell.cellID}, cloneID = {cell.cloneID} Cell Type: {cell.cellType}, mRNA Count: {len(cell.mrna.mrna)}")
#         print(f"cell cycle stage: {cell.cellCycleStage}, division count: {cell.division_count}")
#         print("Top mRNA counts:")
#         for gene, count in sd:
#             if count > 0:
#                 print(f"  Gene: {gene}, Count: {count}")
#
# out gene expression matrix
#   
gene_expression_matrix = np.matrix(np.zeros((len(cells), len(MasterMRNA))), dtype=int)
cl = list(cells)
for i in range(len(cells)):
    cell = cells[i]
    gene_expression_matrix[i, :] = 0  # Initialize row for each cell
    ml = list(MasterMRNA.keys())
    for j in range(len(ml)):
        g = ml[j]
        # get mRNA counts for the gene in the cell    
        data = cell.mrna.get_mrna_counts(g)
        # if the gene is not in the cell, set data to None
        if data is not None:
            gene_expression_matrix[i, j] = data
    #
    # apply PCR bias
    #
    cx = np.linspace(-7, 5, 1000)
    cdf = logistic.cdf(cx, loc=-2.5, scale=2)
    for j in range(len(ml)):
        ri = random.randint(0,1000-1)
        gene_expression_matrix[i, j] = gene_expression_matrix[i, j] * cdf[ri] * 2
#   
# write out the gene expression matrix to a CSV file  
#
output_file = '/Users/markenstrom/Documents/BM5BasedModel/python/5State/gene_expression_matrix.csv'
with open(output_file, 'w') as f:
    # Write header
    f.write("CellID," + ",".join(MasterMRNA.keys()) + "\n")
    # Write data
    for i in range(len(cells)):
        cellTAG = f"{cl[i].cellType}_{cl[i].cellID}_{cl[i].cloneID}_{cl[i].cellCycleStage}_{cl[i].division_count}"
        f.write(f"{cellTAG}," + ",".join(map(str, gene_expression_matrix[i, :].A1)) + "\n")      

#
# out gene expression matrix of simulated RNA-Seq data
# 
def applyRandomSampleing(data):
    count = 0
    for _ in range(data):
        r = random.random()
        if r < 0.15:  # Simulate 15% sampling
            count += 1
    return count    
#
# rnaSEQ sampling bias
#
x = np.arange(1, 100)
p = (2.0 - np.log10(x))
samples = np.random.choice(p, size=1000)

gene_expression_matrix = np.matrix(np.zeros((len(cells), len(MasterMRNA))), dtype=int)
cl = list(cells)
for i in range(len(cells)):
    cell = cells[i]
    gene_expression_matrix[i, :] = 0  # Initialize row for each cell
    ml = list(MasterMRNA.keys())
    for j in range(len(ml)):
        g = ml[j]
        # get mRNA counts for the gene in the cell    
        data = cell.mrna.get_mrna_counts(g)
        if data is None:   
            data = 0
        # apply 15% random sampling to the data
        data = applyRandomSampleing(data)  
        gene_expression_matrix[i, j] = data
    # apply PCR bias
    #bias_factors = np.random.lognormal(mean=0.0, sigma=0.5, size=len(ml))
    bias_factors = np.random.choice(p, size=len(ml))
    gene_expression_matrix[i, :] = np.multiply(gene_expression_matrix[i, :], bias_factors)
#   
# write out the gene expression matrix to a CSV file  
#
output_file = '/Users/markenstrom/Documents/BM5BasedModel/python/5State/gene_rnaseq_matrix.csv'
with open(output_file, 'w') as f:
    # Write header
    f.write("CellID," + ",".join(MasterMRNA.keys()) + "\n")
    # Write data
    for i in range(len(cells)):
        cellTAG = f"{cl[i].cellType}_{cl[i].cellID}_{cl[i].cloneID}_{cl[i].cellCycleStage}_{cl[i].division_count}"
        f.write(f"{cellTAG}," + ",".join(map(str, gene_expression_matrix[i, :].A1)) + "\n") 




#
# clone distribution
# 
d = defaultdict(int)
for cell in cells:
    d[cell.cloneID] += 1
sd = sorted(d.items(), key=lambda x: x[1], reverse=True)
print("Clone distribution:")
for cloneID, count in sd:
    if count > 0:
        print(f"Clone ID: {cloneID}, Count: {count}")   

# for a in cells:
#     print(f"Cell ID: {a.cellID}, Clone ID: {a.cloneID}, Cell Type: {a.cellType}, Cell Cycle Stage: {a.cellCycleStage}, Division Count: {a.division_count}")
#     for i in a.history:
#         print(i)
    
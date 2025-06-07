from mrna_expression import mrnaDict
from mrna_expression import MasterMRNA
from mrna_expression import cycleDict
from mrna_expression import degrageDict
from scipy.stats import poisson
import numpy as np
#
# MRNA represents all the mRNA in a cell.
#
class MRNA:
    def __init__(self):
        self.mrna = []
        self.synthesis_rate = 0.5  # Rate of synthesis
        self.degradation_rate = 0.01 # from 0.05  # Rate of degradation
    #
    # synthesizeCycle is called once per 6 hour to synthesize mRNA for cell cycle genes
    #
    def synthesizeCycle(self, cellType, cellCycleStage, time):
        #
        # cell cycle runs hourly
        #
        cell_cycle_genes = {}
        if (cellCycleStage >= 1) and (cellCycleStage <= 4):
            # S phase genes
            cell_cycle_genes = cycleDict['S']
        elif (cellCycleStage >= 5) and (cellCycleStage <= 8):
            # G2M phase genes
            cell_cycle_genes == cycleDict['G2M']
        #
        # cell cycle genes
        #
        if len(cell_cycle_genes) > 0:
            for gene,mean_expression in cell_cycle_genes.items():
                mean_expression = float(mean_expression)/5.0
                # Simulate mRNA synthesis for each gene based on its mean expression
                synthesized_mRNA = poisson.rvs(mean_expression)
                #print(f"Synthesized {synthesized_mRNA} mRNA from {mean_expression:2.2f} for gene {gene} in cell type {cellType} at time {time}")
                # Store the synthesized mRNA with the time of synthesis
                for _ in range(synthesized_mRNA):
                    # Append a tuple of (time, gene, count) to the mRNA list
                    self.mrna.append((time, gene))
    #
    # main synthesis function runs onece per day
    #
    def synthesize(self, cellType, cellCycleStage,time):
        # Simulate mRNA synthesis at a given time
        #synthesized_mRNA = poisson.rvs(self.synthesis_rate)
        
        if cellType not in mrnaDict:
            print(f"Cell type {cellType} not found in mRNA dictionary.")
            return

        ml = mrnaDict[cellType]
        for gene, mean_expression in ml.items():
            mean_expression = float(mean_expression)/2.0
            # Simulate mRNA synthesis for each gene based on its mean expression
            synthesized_mRNA = poisson.rvs(mean_expression)
            #print(f"Synthesized {synthesized_mRNA} mRNA from {mean_expression:2.2f} for gene {gene} in cell type {cellType} at time {time}")
            # Store the synthesized mRNA with the time of synthesis
            for _ in range(synthesized_mRNA):
                # Append a tuple of (time, gene, count) to the mRNA list
                self.mrna.append((time, gene))


    def degrade(self,current_time):
        # Simulate mRNA degradation at a given time
        nextmRNA = []
        for time_rna in self.mrna:
            dgRate = degrageDict.get(time_rna[1], self.degradation_rate)
            dt = current_time - time_rna[0]
            p_degrade = 1 - np.exp(-dgRate * dt)
            #print(f"Time since synthesis: {dt}, Probability of degradation: {p_degrade}")
            degrade = np.random.binomial(1, p_degrade)
            if degrade > 0:
                pass
                #print(f"Degraded mRNA at time {time_rna[0]}: {time_rna}")
            else:
                nextmRNA.append(time_rna)
        self.mrna = nextmRNA
        #print(f"Remaining mRNA after degradation at time {current_time}: {len(self.mrna)}")
    def get_mrna_counts(self,gene):
        # Get the count of mRNA for a specific gene
        count = 0
        for time, g in self.mrna:
            if g == gene:
                count += 1  
        return count#
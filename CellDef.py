#
#
#
#
from MRNA import MRNA
import random
#
#
#
class Cell:
    CellID = 0
    def __init__(self, cloneID, cellType):
        self.cellID = Cell.CellID
        Cell.CellID += 1
        self.cloneID = cloneID
        self.cellType = cellType
        self.mrna = MRNA()
        self.division_count = 0
        self.cellCycleStage = 0  # Initial cell cycle stage
        self.doHistory = False  # Flag to indicate if history is being tracked
        self.history = []
        # differentiation targets and ratea
        self.stem_targets = ["my1", "er1"]
        self.my1_targets = ["my2", "proB"]
        self.my2_targets = ['mature'] #["cd14_mono"]
        self.er1_targets = ['mature'] #["er2", "er2","bas"]
        self.er2_targets = ["mature"]
        self.preB_targets = ["mature"]
        self.proB_targets = ['mature'] #["preB"]
        self.bas_targets = ["mature"]
        self.cd14_mono_targets = ["mature"]   
    def clone(self):
        # Create a clone of the cell with the same ID and type
        d =  Cell(self.cloneID, self.cellType)
        d.division_count = self.division_count
        d.cellCycleStage = 0
        # Copy mRNA from the original cell to the clone
        d.mrna.mrna = self.mrna.mrna.copy()
        return d
    def sim(self,time,simMRNA=False):
        day = int(time / 24)  # Convert hour to day (6 hours per step)
        hour = time % 24  # Get the hour of the day
        # simulate mRNA synthesis and degradation if simMRNA is True
        if simMRNA:
            self.mrna.degrade(current_time=time)
            self.mrna.synthesizeCycle(cellType=self.cellType, cellCycleStage=self.cellCycleStage, time=time)
            if hour == 0:
                self.mrna.synthesize(cellType=self.cellType, cellCycleStage=self.cellCycleStage, time=time)
        #
        # Track history if doHistory is True
        #
        if self.doHistory:
            target = 'MCM5'
            targetCount = 0
            for (time, gene) in self.mrna.mrna:
                if gene == target:
                    targetCount += 1
            self.history.append((time, targetCount))
        #
        # Simulate cell division and differentiation based on cell type 
        # 
        ret = None           
        if self.cellType == "stem":
            ret = self.sim_param(cycleRate=0.05, mainDiffRate=0.5, targets=self.stem_targets)
        elif self.cellType == "my1":   
            ret =  self.sim_param(cycleRate=0.25, mainDiffRate=0.55, targets=self.my1_targets)
        elif self.cellType == "my2":
            ret =  self.sim_param(cycleRate=0.25, mainDiffRate=0.66, targets=self.my2_targets)
        elif self.cellType == "er1":
            ret =  self.sim_param(cycleRate=0.25, mainDiffRate=0.55, targets=self.er1_targets)
        elif self.cellType == "er2":
            ret =  self.sim_param(cycleRate=0.25, mainDiffRate=0.66, targets=self.er2_targets)
        elif self.cellType == "preB":   
            ret =  self.sim_param(cycleRate=0.35, mainDiffRate=0.90, targets=self.preB_targets)  # much higher to keep # low
        elif self.cellType == "proB":
            ret =  self.sim_param(cycleRate=0.25, mainDiffRate=0.66, targets=self.proB_targets)
        elif self.cellType == "bas":
            ret =  self.sim_param(cycleRate=0.25, mainDiffRate=0.66, targets=self.bas_targets)
        elif self.cellType == "cd14_mono":
            ret =  self.sim_param(cycleRate=0.25, mainDiffRate=0.90, targets=self.cd14_mono_targets)
        else:
            print(f"Unknown cell type: {self.cellType}")
        
        return ret

        
    def stochastic_choice(self, targets):
        """
        Randomly choose a target from the list of targets.
        """
        return random.choice(targets)

        
    def sim_param(self,cycleRate,mainDiffRate,targets):
        """
        Simulate cell division and differentiation based on parameters.
        """
        d1 = None
        # Simulate cell division and differentiation
        #
        # cellCycleStage 0 is G1 phase
        # cellCycleStage 1 to 2 is S phase
        # cellCycleStage 3 to 4 is G2 phase
        if self.cellCycleStage == 0:
            # 25% chance to start cycle
            if random.random() < cycleRate:
                self.cellCycleStage = 1  # Move to S phase
        else:
            self.cellCycleStage += 1
        
        if self.cellCycleStage > 4: 
            # cycle complete, divide
            self.division_count += 1
            self.cellCycleStage = 0
            d1 = self.clone()
            # chance to differentiate
            if random.random() < mainDiffRate:
                ct = self.stochastic_choice(targets)
                self.cellType = ct
            else:
                # remain same cell type
                pass
            # daughter, chance to remain same cell type
            if random.random() < mainDiffRate:
                ct = self.stochastic_choice(targets)
                d1.cellType = ct
            else:
                # remain same cell type
                pass

        return d1





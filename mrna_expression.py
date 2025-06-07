#
# cell state mRNA expression table
#
MasterMRNA = {}
mrnaDict = {}
cycleDict = {}
degrageDict = {}

# s_genes = []
# g2m_genes = []
# fh = [x.strip() for x in open('s_genes.csv','r')]
# for g in fh[1:]:
#     s_genes.append(g)
#     MasterMRNA[g] = 1
# fh = [x.strip() for x in open('g2m_genes.csv','r')]
# for g in fh[1:]:
#     g2m_genes.append(g)
#     MasterMRNA[g] = 1   


# def cycle_gene_expression(cycle_stage):
#     """
#     Returns a list of genes associated with the cell cycle.
#     """
#     if cycle_stage == 0:  # G1 phase  
#         return ""
#     elif cycle_stage > 1 and cycle_stage < 3 :  # S phase 
#         return s_genes
#     else:
#         return g2m_genes


def buildDict(cellType):
    """
    Build a dictionary of mRNA expression levels for a given cell type.
    """
    fn = '/Users/markenstrom/Documents/BM5BasedModel/R_Model_estimate/bm6_' + cellType + '_expression.csv'
    fh = [x.strip().split(',') for x in open(fn,'r')]
    gd = {}
    for gene,count in fh:
        count = int(count)
        gd[gene] = count
        MasterMRNA[gene] = 1
    mrnaDict[cellType] = gd

buildDict('stem')
buildDict('my1')
buildDict('my2')
buildDict('er1')
buildDict('er2')
buildDict('preB')
buildDict('proB')
buildDict('bas')
buildDict('cd14_mono')



def buildCycleDict(ss):
    """
    Build a dictionary of mRNA expression levels for a given cell type.
    """
    fn = '/Users/markenstrom/Documents/BM5BasedModel/R_Model_estimate/bm6_' + ss + '_expression.csv'
    fh = [x.strip().split(',') for x in open(fn,'r')]
    gd = {}
    for gene,count in fh:
        count = int(count)
        gd[gene] = count
        MasterMRNA[gene] = 1
    cycleDict[ss] = gd


buildCycleDict('S')
buildCycleDict('G2M')

l = [x.strip().split(',') for x in open('/Users/markenstrom/Documents/BM5BasedModel/R_Model_estimate/rna_degrade.csv','r')]

for line in l:
    degradation_rate = float(line[1])
    gene = line[0]
    degrageDict[gene] = degradation_rate

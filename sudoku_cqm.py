# Some tests with the cqm solver
import dimod
import neal
from dwave.system import DWaveSampler, EmbeddingComposite, LeapHybridBQMSampler

# The objective is that the sum of all cells equals (1+2+...+8+9)*9 = 405
goal = 405

# Constraints are that we don't need the same number on the same line
# and not the same in one of the 9 squares

indexes = range(0, 9)
choices = [i+1 for i in indexes]
vars = [[[dimod.Binary(f'v_{i}_{j}_{k}') for k in choices] for j in indexes] for i in indexes]

cqm = dimod.ConstrainedQuadraticModel()
#cqm.add_constraint(sum([vars[i][j][k] * (k+1) for i in indexes for j in indexes for k in indexes]) == goal)
# constraints in single cell
for i in indexes:
    for j in indexes:
        cqm.add_discrete(sum([vars[i][j][k] for k in indexes]) == 1, label=f'cell_{i}_{j}')
# constraints on columns
for i in indexes:
    for k in indexes:
        cqm.add_constraint(sum([vars[i][j][k] for j in indexes]) == 1, label=f'column_{i}_k{k}')
# constraints on lines
for j in indexes:
    for k in indexes:
        cqm.add_constraint(sum([vars[i][j][k] for i in indexes]) == 1, label=f'line_{j}_k{k}')

# constraints on sub-squares
for si in range(0,3):
    for sj in range(0,3):
        for k in indexes:
            cqm.add_constraint(sum([vars[i][j][k] for i in range(si*3, si*3+3) for j in range(sj*3,sj*3+3)]) == 1, label=f'sub-square_{si}_{sj}_k{k}')

# Initial positions
sudoku = [
    [0, 0, 0,  0, 0, 0,  1, 9, 0],
    [2, 3, 0,  0, 0, 0,  6, 0, 0],
    [0, 0, 0,  2, 4, 0,  0, 0, 0],

    [0, 0, 0,  0, 0, 0,  9, 6, 0],
    [0, 0, 0,  1, 6, 0,  0, 7, 0],
    [0, 4, 8,  0, 7, 0,  0, 0, 0],

    [0, 0, 1,  0, 0, 3,  4, 0, 5],
    [0, 0, 9,  0, 0, 8,  0, 0, 0],
    [0, 0, 6,  0, 0, 5,  8, 0, 0],
]
for i in indexes:
    for j in indexes:
        print(sudoku[i][j], end=' ')
    print("")
for i in indexes:
    for j in indexes:
        if sudoku[i][j] != 0:
            cqm.add_constraint(vars[i][j][sudoku[i][j]-1] + sum(vars[i][j][k]*999 for k in indexes if k != sudoku[i][j]-1) == 1, label=f'initial_{i}_{j}')

bqm, inverter = dimod.cqm_to_bqm(cqm)

#sampler = neal.SimulatedAnnealingSampler()
#sampleset = sampler.sample(bqm)
sampler = LeapHybridBQMSampler()
sampleset = sampler.sample(bqm)

invert = inverter(sampleset.first.sample)

print(sampleset)
sample = sampleset.first
for i in indexes:
    for j in indexes:
        print([k for k in choices if invert[f'v_{i}_{j}_{k}'] == 1], end='')
    print("")
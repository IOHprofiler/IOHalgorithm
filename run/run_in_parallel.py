import os, itertools
import numpy as np
from mpi4py import MPI  

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

## settings for experiments
# The settings of
algorithm_name = [
  "ea", "ea2", "ea23", "llea",
  "rls", "rs", "ghc", "fga",
  "2ratega", "umda", "sa", "sars", "fea"
]
suite_name = [
  "wmodelonemax",
  "wmodelleadingones"
]
dir_folder = "./"
eval_budget = 10000
ind_runs = 10
rand_seed = 666
## The setting of w-model is defined in main.cpp
## Please make sure that problem_id is in the corresponding domain!
problem_str = ["1-10", "11-20", "21-30","30-40", "40-49"]
instance_str = "1-5"
dimension = "20,100"

if rank == 0:
    alg_p = list(itertools.product(algorithm_name,problem_str, suite_name))
    N = len(alg_p)

    r = N % size
    step = int((N - r) / size)
    data = [alg_p[(i * step) : ((i + 1) * step)] for i in range(size)]
    for i in range(r) :
      data[i].append(alg_p[step * size + i]) 
else:
    data = None

data = comm.scatter(data, root=0)

for alg_p in data:
    algorithm = alg_p[0]
    problem = alg_p[1]
    suite = alg_p[2]
    command = './main {0} {1} {2} {3} {4} {5} {6} {7} {8}'.format(algorithm,suite, problem, instance_str, dimension, dir_folder, ind_runs, eval_budget, rand_seed)
    os.system(command)
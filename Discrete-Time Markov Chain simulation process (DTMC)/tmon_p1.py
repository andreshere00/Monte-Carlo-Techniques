#!/usr/bin/env python
# coding: utf-8

# # **Montecarlo Techniques - Practice 1**
# 
# Andrés Herencia López-Menchero
# 
# MUTECI - 2023/2024
# 
# ## **Old car simulation process**

import numpy as np
import random

### Implement the simulation model of the old car problem in R or python.

def dtmc(N = 10**3, transitionMatrix = [[1, 0, 0, 0], [0.05, 0.50, 0.20, 0.25], [0.10, 0.75, 0, 0.15], [0, 0, 0, 1]], initialState = 1, finalState = 0):

    """
    Simulate a Discrete-Time Markov Chain (DTMC).

    Parameters:
    - N (int): The number of steps to simulate the Markov Chain. Default is 1000.
    - transitionMatrix (list of lists): The transition matrix representing the DTMC. Default is a 4x4 matrix.
    - initialState (int): The initial state from which the simulation starts. Default is state 1.
    - finalState (int): The absorption state, if any, at which the simulation should stop. Default is state 0.

    Returns:
    - t: A list which contains the number of steps from the initialState to the finalState in each iteration.
    """
    
    # list where we will save our states
    states = [initialState]
    np.random.seed(123)
    t = np.ones(N)  # simulation time

    # check if we have an absorption state
    absorption_state = []
    for tm in range(0,len(transitionMatrix)):
        for m in transitionMatrix[tm]:
            if (m == 1)&(transitionMatrix[tm].index(m) == tm):
                absorption_state.append(tm)

    # simulation
    for n in range(N):
        i = initialState
        while (i not in absorption_state):
            states.append(i)
            j = -1; # columns of our transition matrix (this is the 'next state' column)
            pj = 0; # transition probability from state 'i' to state 'j'.
            threshold = np.random.uniform(0, 1)
            while (pj < threshold):
                j += 1
                pj = pj + transitionMatrix[i][j]
            i = j
            t[n] += 1
    return(t)

help(dtmc)

# ### Approximate the absorption time from `Operativo` to `Siniestrado` through simulation (years).

round(np.mean(dtmc()),4)


# ### Approximate through simulation the probability of absorption from `Averiado` to `Vendido` (%).

N = 10**3
round((N/sum(dtmc(initialState = 2, finalState = 3)))*100,4)
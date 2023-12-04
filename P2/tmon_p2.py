import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def mm1_model(N=100, L=1/7, mu=1/5, seed = 12345):
    
    """
    Simulates an M/M/1 queue system and calculates some performance metrics.

    Args:
    - N (int): Number of time units to simulate. Default = 100.
    - L (float): Average arrival rate (customers per time unit). Default = 1/7.
    - mu (float): Average service rate (customers per time unit). Default = 1/5.
    - seed (int): Seed for saving the random state. Default = 12345.

    Returns:
    tuple: A tuple containing:
        - queue_model (pandas.DataFrame): A DataFrame that records the system's state at each time unit.
        - Ls (float): Average number of customers in the system.
        - Lq (float): Average number of customers in the queue.
        - Ws (float): Average time customers spend in the system.
        - Wq (float): Average time customers spend in the queue.

    Raises:
    Exception: An exception is raised if the arrival rate (L) is greater than the service rate (mu), indicating that the system
    has not reached the stationary state (process explosive).
    """
    
    if L/mu > 1:
        raise Exception("The system has not reached stationary state. You must redefine the parameters of your model.")
        
    np.random.seed(seed)
    arrival_time = np.random.exponential(scale=1/L)
    service_time = np.random.exponential(scale=1/mu)
    stay = 0; t = 0; queue = 0; service = 0; arrivals = 0;

    queue_model = pd.DataFrame({
        't': [t],
        'queue': [queue],
        'service': [service],
        'arrivals': [arrivals],
        'stay': [stay],
        'arrival time': [arrival_time],
        'service time': [service_time]
    })

    while min(arrival_time, service_time) <= N:
        if arrival_time <= service_time:
            t = arrival_time
            if service > 0: # client stays in the queue (C1Q1)
                queue += 1
            else: # a client can be served in this moment (C1Q0)
                service = 1
            arrivals += 1
            arrival_time = t + np.random.exponential(scale=1/L)
        else:                   
            t = service_time
            if queue > 0: # client is dispatched from the queue to the system (C0Q1)
                queue -= 1
                service_time = t + np.random.exponential(scale=1/mu)
            else: # no one in queue and previous service has finished, system at rest (C0Q0)         
                service = 0
                service_time = arrival_time + np.random.exponential(scale=1/mu)
        stay = min(arrival_time, service_time) - t

        new_register = pd.DataFrame({
            't': [t],
            'queue': [queue],
            'service': [service],
            'arrivals': [arrivals],
            'stay': [stay],
            'arrival time': [arrival_time],
            'service time': [service_time]
        })
        queue_model = pd.concat([queue_model, new_register], ignore_index=True)
    
    queue_time = []; wait_times = []
    for q in range(1,len(queue_model['queue'])):
        if queue_model['queue'][q] > queue_model['queue'][q-1]:
            queue_time.append(queue_model['t'][q])                                                                                               
        if queue_model['queue'][q] < queue_model['queue'][q-1]:
            wait_times.append(queue_model['t'][q] - queue_time.pop(0))
            
    # Average clients number on the system (Ls)
    Ls = sum(queue_model['service'])/t
    # Average clients number on the queue (Lq)
    Lq = sum(queue_model['queue'])/t
    # Average clients time on the system (Ws)
    Ws = np.mean(np.diff(queue_model['service time'].unique()))
    # Average clients time on the queue (Wq)
    Wq = np.mean(wait_times)
    
    return queue_model, Ls, Lq, Ws, Wq
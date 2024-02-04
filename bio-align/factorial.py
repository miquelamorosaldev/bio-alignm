import numpy as np

def factorial(n):
    array = np.arange(n+1, dtype=np.dtype("int64"))
    array[0] = 1

    for i in range(1, n+1):
        array[i] = array[i-1]*i

    return array[n]

assert factorial(6) == 720, f"n = {factorial(6)}"

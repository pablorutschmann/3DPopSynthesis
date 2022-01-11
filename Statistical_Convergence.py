import numpy as np

lst = [100,200,400,800,1000]
print(lst[:-3])


N_samples = 100000
ran = np.random.uniform(10,100,N_samples)

integral = (0.5 * 100**2 - 0.5 * 10**2) / (100-10)


print(np.mean(ran))
print(np.std(ran))
print(integral)
import Library as l
import matplotlib.pyplot as plt
import numpy as np

simulation_space = l.Simulation_Space()
data = [[] for _ in range(1,8)]
plt.figure()
for i in range(1,8):
    print(i)
    inf = l.Infection(50*i,0.1,.2)
    simulation_space.objects = []
    simulation_space.populate([l.Susceptible(tuple(np.random.uniform(-5.0,5.0,size=2))) for _ in range(100)])
    simulation_space.populate(([l.Infected((0,0),infection=inf)]))
    for frame in range(1000):
        l.basic_simulate(inf,simulation_space)
        data[i-1].append(l.get_data(simulation_space,data_type=['N_infected'])[0])
    plt.plot(range(1000),data[i-1],label='%s'%(50*(i-1)))

plt.legend()
plt.show()


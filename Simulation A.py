"""
Disease Spread Simulation
Written by Nathan Diggins
"""

import numpy as np
import matplotlib.pyplot as plt
import Library as s

"""
Base Functions
"""


def simulate(infection, simulation_space):
    """
    Simulates a single step in the simulation process.
    :param transmission_rate: This is the rate of transmission within the radius
    :param simulation_space: This is the simulation space within which to simulate.
    :param transmission_radius: The possible radius of transmission.
    :return: Return None
    """
    # To start, we check if objects are close enough to transmit.
    for transmitter in simulation_space.objects:
        # We now iterate over all of the objects in the simulation space
        for receiver in simulation_space.objects[simulation_space.objects.index(transmitter):]:
            # This matches each transmitter with each receiver exactly once.
            if s.distance(transmitter, receiver) < infection.transmittance_radius and transmitter in simulation_space and receiver in simulation_space:
                """
                Now we know if the two objects are within the necessary distance.
                """
                if type(transmitter) == s.Infected and type(
                        receiver) == s.Susceptible and np.random.uniform() <= infection.transmittance(s.distance(transmitter,receiver)):
                    print("Receiver %s was exposed"%(receiver))
                    receiver.expose(infection,s.distance(transmitter,receiver))
                elif type(transmitter) == s.Susceptible and type(
                        receiver) == s.Infected and np.random.uniform() <= infection.transmittance(s.distance(transmitter,receiver)):
                    print("Transmitter %s was exposed"%(transmitter))
                    transmitter.expose(infection,s.distance(transmitter,receiver))
                else:
                    pass
            else:
                pass
    """
    Now we have spread the disease, lets look at curing people.
    """

    infecteds = [i for i in simulation_space.objects if type(i) == s.Infected]
    for object in infecteds:
        if np.random.uniform() <= (infection.death_rate/infection.time):
            object.kill()
        elif object.infected_time >= infection.time:
            # Clearly the object needs to be cured.
            object.cure()
        else:
            object.infected_time += 1

    """
    Now lets move the objects
    """

    for object in simulation_space.objects:
        s.move(object)

    return None



sim_space = s.Simulation_Space()
sick = s.Infection(500,0.5,0.2)
sim_space.populate([s.Susceptible(tuple(np.random.uniform(0.0-(sim_space.size[0]/2.0),0.0+(sim_space.size[0]/2.0),2)),
                                  direction=np.random.uniform(-np.pi,np.pi)) for i in range(50)])
sim_space.populate(s.Infected((1,1)))

immune = 0
infs= 1
sus = 50
data = []
while infs != 0:
    print(len(sim_space.objects),immune,infs,sus)
    simulate(sick,sim_space)
    data.append(s.get_data(sim_space))
    immune = data[-1][2]
    infs = data[-1][1]
    sus = data[-1][3]
plt.figure()
plt.plot(range(len(data)),[i[0] for i in data],label='Total')
plt.plot(range(len(data)),[i[1] for i in data],label="Infected")
plt.plot(range(len(data)),[i[2] for i in data],label='Immune')
plt.plot(range(len(data)),[i[3] for i in data],label='Susceptible')
plt.legend()
plt.show()
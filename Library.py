"""

Disease spread modeling library
Written by: Nathan Diggins

"""

import numpy as np


class Simulation_Space:
    """
    This is the space in which the simulation will occur.
    """

    def __init__(self, resolution=1000, dimensions=(10, 10)):
        """
            Intiates the Simulation Space
            :param resolution: This is the number of points to display in the simulation space.
            :param dimensions: Size of the simulation space.
        """
        self.size = dimensions
        self.resolution = resolution
        self.space = np.meshgrid(np.linspace(-dimensions[0] / 2, dimensions[0] / 2, resolution),
                                 np.linspace(-dimensions[1] / 2, dimensions[1] / 2, resolution))
        self.objects = []

    def __contains__(self, other):
        """
        This tells if the object is within the simulation space.
        :param other: (either tuple or an object type)
        :return: True or False
        """
        if isinstance(other, Object):
            # The object is an object and we need to check the position of the object.
            if other in self.objects and abs(other.position[0]) <= self.size[0] / 2.0 and abs(other.position[1]) <= self.size[
                1] / 2.0:
                # The object is inside of the Simulation Space
                return True
            else:
                return False
        if isinstance(other, tuple):
            # This object is a tuple, we just have to check the position
            if abs(other[0]) <= self.size[0] / 2.0 and abs(other[1]) <= self.size[1] / 2.0:
                return True
            else:
                return False

    def populate(self, objects):
        """
        Populates the simulation space
        :param objects: objects to be added to the simulation
        :return: null
        """
        if type(objects) in [Object, Infected, Immune, Susceptible]:
            # First we add the object to the simulation space
            self.objects.append(objects)
            objects.sim_space = self  # Now we add a simulation space to the object.
        elif type(objects) == list:
            self.objects += objects
            for i in objects:
                i.sim_space = self
        else:
            raise TypeError

    def remove(self, object):
        """
                removes an object from the population list
                :param object: Type object, the object to be removed.
                :return: The object if capable, null if failed.
                """
        self.objects.remove(object)

    def __len__(self):
        return len(self.objects)

    def resize(self, dimensions):
        return self.__init__(self.resolution, dimensions)

    def plot_space(self,subplot,style='r-',linewidth=3.0,alpha=1.0):
        """
        This will graph the borders of the board.
        :param subplot: This is a matplotlib.pyplot.subplot object
        :return: Return none
        """
        subplot.plot([-1*(self.size[0]/2.0),self.size[0]/2.0,self.size[0]/2.0,-1*(self.size[0]/2.0),-1*(self.size[0]/2.0)],
                     [self.size[1]/2.0,self.size[1]/2.0,-1*self.size[1]/2.0,-1*self.size[1]/2.0,self.size[1]/2.0],
                     style,
                     linewidth=linewidth,alpha=alpha)
        return None

    def plot_objects(self, subplot, color_dictionary=None,object_size=2):
        if color_dictionary is None:
            color_dictionary = {Infected: 'r.', Immune: 'b.', Susceptible: 'g.'}
        for object in self.objects:
            subplot.plot([object.position[0]],[object.position[1]],color_dictionary[type(object)],markersize=object_size)
        return None



class Object:
    def __init__(self, position, sim_space=None, direction=0.0):
        """
        Defines the base class for objects within the simulation. These objects can be whatever type is intended by the
        user by definition through a higher class with inheritance to Object.
        :param position: Position in space (x,y). This should only ever be a 2D vector as the simulation space will
                            provide z axis positioning if ever necessary.
        :param sim_space: If the user wishes to avoid using SimulationSpace.populate(), they can auto populate through
                            this parameter.
        """
        self.direction = direction
        self.position = position
        self.x, self.y = position
        # Adding the object to a simulation space if specified
        if sim_space != None:
            # The user wishes to add a simulation space
            try:
                sim_space.populate(self)
                self.sim_space = sim_space
            except:
                raise UserWarning(
                    "Simulation Space %s doesn't exist." % (sim_space))  # todo: Possible to return string of var name?
                return None
        else:
            self.sim_space = None


class Infected(Object):
    def __init__(self, position, sim_space=None, direction=0.0, infection=None):
        """
        Initiates the infected class of objects. These are objects that can currently transmit the disease.
        :param position: Location of the object
        :param sim_space: The simulation space in which to place the object.
        :param direction: facing location of the objects.
        :param infection: The infection to give to the object.
        """
        Object.__init__(self, position, sim_space=sim_space, direction=direction)
        self.infection = infection
        self.infected_time = 0
    def __repr__(self):
        return "Infected"
    def cure(self):
        """
        Cures the object, transferring it to the state of immune.
        :return: Return Immune type object
        """
        cured = Immune(self.position, sim_space=self.sim_space, direction=self.direction)
        print("Curing %s" % (self))
        self.sim_space.remove(self)
        del self

    def kill(self):
        print("Killing %s" % (self))
        self.sim_space.remove(self)
        del self


class Immune(Object):
    def __init__(self, position, sim_space=None, direction=0.0):
        """
        Initiates the infected class of objects. These are objects that can currently transmit the disease.
        :param position: Location of the object
        :param sim_space: The simulation space in which to place the object.
        :param direction: facing location of the objects.
        """
        Object.__init__(self, position, sim_space=sim_space, direction=direction)
    def __repr__(self):
        return "Immune"

class Susceptible(Object):
    def __init__(self, position, direction=0.0, sim_space=None):
        Object.__init__(self, position, sim_space=sim_space, direction=direction)
    def __repr__(self):
        return 'Susceptible'
    def expose(self, infection, distance):
        """
        Infects the susceptible object.
        :param infection: The infection to give them
        :return: None
        """
        print("Exposing %s" % (self))
        if np.random.uniform() <= infection.transmittance(distance):
            print("Infecting %s" % (self))
            infected = Infected(self.position, sim_space=self.sim_space, direction=self.direction, infection=infection)
            self.sim_space.remove(self)
            del self
        else:
            pass


class Infection:
    """
        Defines an infection type
    """

    def __init__(self, infection_time, transmittance_radius, death_rate, transmittance_function=lambda x: np.e ** (x)):
        """
        Initiates the Infection class
        :param transmittance_function: This is a lambda function of radius which expresses the likelihood of infection at distance R
        :param transmittance_radius: the radius underwhich to calculated the trans. function
        :param death_rate: the death rate of the disease.
        """
        self.transmittance_radius = transmittance_radius
        self.transmittance = transmittance_function
        self.death_rate = death_rate
        self.time = infection_time


def distance(object1, object2):
    """
    Calculates the distance between two objects.
    :param object1:
    :param object2:
    :return: <Float> Distance
    """
    return np.sqrt((object1.x - object2.x) ** 2 + (object1.y - object2.y) ** 2)


def move(obj, speed_parameter=0.1, rotation_scale_parameter=0.20):
    """
    This function allows for the motion of objects.
    :param rotation_scale_parameter: This is the standard deviation in the gaussian used to produce the rotation effect.
    :param obj: This is the object to move
    :param speed_parameter: This is the step speed of the object.
    :return: Returns the new position of the object
    """
    # TODO: Allow for different choices based on type, I.E. infected's don't move as much and susceptible avoids contact.
    out_of_bounds = None
    while out_of_bounds is not False:
        """
        This attempts to constrain the motion of the particle to only the simulation space.
        """
        obj.direction = np.random.normal(obj.direction, rotation_scale_parameter)
        temp_position = (obj.position[0] + (speed_parameter * np.cos(obj.direction)),
                         obj.position[1] + (speed_parameter * np.sin(obj.direction)))
        obj.x, obj.y = obj.position
        if temp_position in obj.sim_space:
            out_of_bounds = False
            obj.position = temp_position
        else:
            out_of_bounds = True
    return obj.position


def get_data(simulation_space, data_type='all'):
    """
    Gets data from the simulation space.
    :param simulation_space: <type sim_space> simulation space.
    :param data_type: What data to gather.
    :return: returns a list of the data points
    """
    data = []
    if data_type == 'all':
        data_type = ['N_total', 'N_infected', 'N_cured', 'N_susceptible']
    for d_type in data_type:
        if d_type == 'N_infected':
            data.append(len([obj for obj in simulation_space.objects if type(obj) == Infected]))
        if d_type == 'N_cured':
            data.append(len([obj for obj in simulation_space.objects if type(obj) == Immune]))
        if d_type == 'N_susceptible':
            data.append(len([obj for obj in simulation_space.objects if type(obj) == Susceptible]))
        if d_type == 'N_total':
            data.append(len([obj for obj in simulation_space.objects if type(obj) == Immune or type(obj) == Infected]))
    return data

def basic_simulate(infection, simulation_space):
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
            if distance(transmitter, receiver) < infection.transmittance_radius and transmitter in simulation_space and receiver in simulation_space:
                """
                Now we know if the two objects are within the necessary distance.
                """
                if type(transmitter) == Infected and type(
                        receiver) == Susceptible and np.random.uniform() <= infection.transmittance(distance(transmitter,receiver)):
                    print("Receiver %s was exposed"%(receiver))
                    receiver.expose(infection,distance(transmitter,receiver))
                elif type(transmitter) == Susceptible and type(
                        receiver) == Infected and np.random.uniform() <= infection.transmittance(distance(transmitter,receiver)):
                    print("Transmitter %s was exposed"%(transmitter))
                    transmitter.expose(infection,distance(transmitter,receiver))
                else:
                    pass
            else:
                pass
    """
    Now we have spread the disease, lets look at curing people.
    """

    infecteds = [i for i in simulation_space.objects if type(i) == Infected]
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
        move(object)

    return None
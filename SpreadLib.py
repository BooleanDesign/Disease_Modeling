"""

Disease spread modeling library
Written by: Nathan Diggins

"""

import numpy as np
import matplotlib.pyplot as plt  # todo: potentially unnecessary


class Simulation_Space:
    """
    Generation and geographic system for carrying out various simulations
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

    def __contains__(self,
                     object):  # TODO: Revise this to give a geometric idea of object being in or out, not the list concept.
        """
        Tests for membership in the simulation space.
        :param object: Any form of object type.
        :return: Boolean
        """
        if isinstance(object,Object):
            """
            An object is being checked
            """
            if object in self.objects and False not in [abs(object.position[i]) < self.size[i] / 2.0 for i in
                                                        range(len(self.size))]:
                return True
            else:
                return False
        elif type(object) == tuple:
            """
            A point was inserted as the location point
            """
            if False not in [abs(object[i]) < self.size[i] / 2.0 for i in
                             range(len(self.size))]:
                return True
            else:
                return False

    def populate(self, objects):
        """
        Populates the simulation space
        :param objects: objects to be added to the simulation
        :return: null
        """
        if type(objects) in [Object,Infected,Immune,Susceptible]:
            self.objects.append(objects)
            objects.sim_space = self
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
        try:
            self.objects.remove(object)
            return object
        except:
            return None

    def __len__(self):
        return len(self.objects)

    def resize(self, dimensions):
        return self.__init__(self.resolution, dimensions)


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

    def __del__(self):
        """
        Deletes the object from memory and from the simulation space.
        :return: None
        """
        #TODO: Figure out why there are two calls to the method
        if self.sim_space is None:
            return None
        else:
            if self in self.sim_space:
                self.sim_space.remove(self)
                return None


class Infected(Object):
    """
    These are objects which are actually infected with the virus and contagious
    """

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

    def cure(self):
        """
        Cures the object, transferring it to the state of immune.
        :return: Return Immune type object
        """
        cured = Immune(self.position, sim_space=self.sim_space, direction=self.direction)
        if self.sim_space is not None:
            # We are going to delete from the simspace
            self.__del__()
            return None
        else:
            self.__del__()
            return cured


class Immune(Object):
    def __init__(self, position, sim_space=None, direction=0.0):
        """
        Initiates the infected class of objects. These are objects that can currently transmit the disease.
        :param position: Location of the object
        :param sim_space: The simulation space in which to place the object.
        :param direction: facing location of the objects.
        """
        Object.__init__(self, position, sim_space=sim_space, direction=direction)


class Susceptible(Object):
    """
    This is the uninfected class
    """

    def __init__(self, position, direction=0.0, sim_space=None):
        Object.__init__(self, position, sim_space=sim_space, direction=direction)

    def infect(self, infection):
        """
        Infects the susceptible object.
        :param infection: The infection to give them
        :return: None
        """
        infected = Infected(self.position, sim_space=self.sim_space, direction=self.direction, infection=infection)
        if self.sim_space is not None:
            Object.__del__(self)
            return None
        else:
            self.__del__()
            return infected

class Infection:
    """
    Defines an infection type
    """
    def __init__(self,transmission_rate,transmission_radius,sickness_time):
        self.transmittance = transmission_rate
        self.radius = transmission_radius
        self.time = sickness_time

"""
Simulation Functions
"""

def distance(object1,object2):
    """
    Calculates the distance between two objects.
    :param object1:
    :param object2:
    :return: <Float> Distance
    """
    return np.sqrt((object1.x-object2.x)**2+(object1.y-object2.y)**2)


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

def get_data(simulation_space,data_type='all'):
    """
    Gets data from the simulation space.
    :param simulation_space: <type sim_space> simulation space.
    :param data_type: What data to gather.
    :return: returns a list of the data points
    """
    data = []
    if data_type == 'all':
        data_type = ['N_total','N_infected','N_cured','N_susceptible']
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

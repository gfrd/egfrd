#!/usr/env python

import numpy

from _gfrd import (
    SphericalShell,
    SphericalShellContainer,
    CylindricalShell,
    CylindricalShellContainer
    )
from utils import *
from single import NonInteractionSingle
from multi import Multi
from gfrdbase import (
    get_closest_surface
    )

class ShellContainer(object):

    def __init__(self, world):

        # the containers hold the spherical and cylindrical shells respectively 
        # the containers is a cache for the shell objects in the domain and is mostly
        # read only. The only method that writes the containers is self.move_shell
        self.containers = [SphericalShellContainer(world.world_size,
                                                   world.matrix_size),
                           CylindricalShellContainer(world.world_size,
                                                     world.matrix_size)]

        self.world = world
        self.user_max_shell_size = numpy.inf    # Note: shell_size is actually the RADIUS of the shell

    def get_matrix_cell_size(self):
        return self.containers[0].cell_size     # cell_size is the width of the (cubic) cell

    # Here shell_size is actually the shell radius (and not shell diameter)
    def set_user_max_shell_size(self, size):
        self.user_max_shell_size = size

    def get_user_max_shell_size(self):
        return self.user_max_shell_size

    def get_max_shell_size(self):
        return min(self.get_matrix_cell_size() * .5 / SAFETY,
                   self.user_max_shell_size)

    def get_container(self, shell):
    # Private method
    # Returns the container that holds 'shell'
        if type(shell) is SphericalShell:
            return self.containers[0]
        elif type(shell) is CylindricalShell:
            return self.containers[1]

    def move_shell(self, shell_id_shell_pair):
    # This is the only method that writes to the shell containers
        shell = shell_id_shell_pair[1]
        container = self.get_container(shell)
        container.update(shell_id_shell_pair)

    def remove_shell(self, shell_id_shell_pair):
        shell_id = shell_id_shell_pair[0]
        shell = shell_id_shell_pair[1]
        container = self.get_container(shell)
        del container[shell_id]

    def distance(self, pos1, pos2):
        # Export the distance calculation
        return self.world.distance(pos1, pos2)

    def get_intruders(self, position, radius, ignore):
        # gets the intruders in a spherical volume of radius 'radius'?
        # TODO make this surface specific -> in 2D only need to check in cylinder,
        # not in sphere
        intruders = []                  # intruders are domains within radius
        closest_domain = None           # closest domain, excluding intruders.
        closest_distance = numpy.inf    # distance to the shell of the closest domain.

        seen = set(ignore)
        for container in self.containers:
            neighbors = container.get_neighbors(position)
            for n in neighbors:
                domain_id = n[0][1].did
                distance = n[1]
                if distance > radius:
                    if distance < closest_distance:
                        # This is domain (the first one for this 
                        # container) that has a shell that is more than 
                        # radius away from pos.  If it is closer than 
                        # the closest such one we found so far: store 
                        # it. Always break out of the inner for loop and 
                        # check the other containers.
                        closest_domain = domain_id
                        closest_distance = distance
                        break
                    else:
                        break
                elif domain_id not in seen:
                    # Since a domain can have more than 1 shell (multis 
                    # for example), and for each shell there is an entry 
                    # in the shell container, we make sure each domain 
                    # occurs only once in the returned list here.
                    seen.add(domain_id)
                    intruders.append(domain_id)

        return intruders, closest_domain, closest_distance

    def get_closest_obj(self, pos, domains, ignore=[], ignores=[]):
        # ignore: domain ids.
        # TODO It's bad to pass the domains to this function, but is currently
        # only way to be able to return a domain object (instead of domain_id)
        closest_domain = None
        closest_distance = numpy.inf

        for container in self.containers:
            result = container.get_neighbors(pos)

            for shell_id_shell_pair, distance in result:
                domain_id = shell_id_shell_pair[1].did

                if ((domain_id not in ignore) and (distance < closest_distance)):
                    domain = domains[domain_id]
                    closest_domain, closest_distance = domain, distance
                    # Found yet a closer domain. Break out of inner for 
                    # loop and check other containers.
#                    break

        surface, distance = get_closest_surface(self.world, pos, ignores)

        if distance < closest_distance:
            return surface, distance
        else:
            return closest_domain, closest_distance

    def get_neighbor_domains(self, position, domains, ignore=[]):
    # returns the neighboring domains and the distance to their nearest shell (Multi's can have multiple)

        dists = {}

        for container in self.containers:
            neighbors = container.get_neighbors(position)

            for shell_id_shell_pair, distance in neighbors:
                domain_id = shell_id_shell_pair[1].did

                if (domain_id not in ignore) and \
                   ((domain_id not in dists) or (distance < dists[domain_id])):

                    dists[domain_id] = distance

        neighbor_domains = [(domains[domain_id], distance) for domain_id, distance in dists.iteritems()]
        return neighbor_domains


    def filter_distance(self, neighbor_domains, radius):
        # returns all the domains that are less than 'radius' distance away
        intruders = [(domain, distance) for (domain, distance) in neighbor_domains if distance < radius]
        return intruders        

    def filter_partners(self, neighbor_domains):

        partners = []
        for domain, distance in neighbor_domains:
            if (isinstance (domain, NonInteractionSingle) and domain.is_reset()) or \
               isinstance (domain, Multi):
                partners.append((domain, distance))

        return partners

    def sort_domain_distance(self, neighbor_domains):
        return sorted(neighbor_domains, key=lambda domain_dist: domain_dist[1])

    def get_neighbor_surfaces(self, position, ignores=[]):
        # analogous to get_neighbor domains, this returns a list of (surface, distance) tuples

        # TODO return a list of neighboring surfaces
        surface, distance = get_closest_surface(self.world, position, ignores)
        if surface:
            return [(surface, distance), ]
        else:
            return []

    def get_neighbors_within_radius_no_sort(self, pos, radius, ignore=[]):
        # Get neighbor domains within given radius.
        # Note that the function returns a generator
        #
        # ignore: domain ids.
        #
        # Only returns neighbor domain_ids, not the distances towards their 
        # shells. Can for example be used to try to clear all objects 
        # from a certain volume.

        try:
            for container in self.containers:
                result = container.get_neighbors_within_radius(pos, radius)
                # result = [((shell_id_shell_pair), distance), ]
                # Since a domain can have more than 1 shell (multis for 
                # example), and for each shell there is an entry in the 
                # shell container, we make sure each domain occurs only once 
                # in the returned list here.
                for did in uniq(s[0][1].did for s in result):
                    if did not in ignore:
                        yield did
        except AttributeError:
            pass

    def get_dids_shells(self):
        did_map = {}
        shell_map = {}
        for container in self.containers:
            if self.world.world_size != container.world_size:
                raise RuntimeError('self.world.world_size != '
                                   'container.world_size')
            for shell_id, shell in container:
                did_map.setdefault(shell.did, []).append(shell_id)
                shell_map[shell_id] = shell

        return did_map, shell_map

    def get_total_num_shells(self):
        return sum(len(container) for container in self.containers)


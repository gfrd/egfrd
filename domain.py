#!/usr/env python

import numpy
import myrandom

__all__ = [
    'Domain',
    'ProtectiveDomain',
    ]


class Domain(object):
# The domain is the main unit of eGFRD. Single, Pairs and Multis are all domains
# with an own domain_id

    SINGLE_SHELL_FACTOR = 2.0

    def __init__(self, domain_id):
	self.domain_id = domain_id	# identifier for this domain object
	self.event_id = None		# identifier for the event coupled to this domain

	self.event_type = None

	self.last_time = 0.0
	self.dt = 0.0

    def initialize(self, time):		# this method needs to be overloaded in all subclasses
        '''
        Sets the local time variables for the Domain (usually to the current simulation time).
        Do not forget to reschedule this Domain after calling this method.

        For the NonInteractingSingle, the radius (shell size) should be shrunken to the actual
        radius of the particle.

        initialize allows you to reuse the same Domain (of the same class with the same particles)
	at a different time.
        '''
	pass

    def calc_ktot(self, reactionrules):
	# calculates the total rate for a list of reaction rules
        # The probability for the reaction to happen is proportional to 
        # the sum of the rates of all the possible reaction types.
	k_tot = 0
        for rr in reactionrules:
            k_tot += rr.k
	return k_tot

    def draw_reaction_rule(self, reactionrules):
	# draws a reaction rules out of a list of reaction rules based on their
	# relative rates
        k_array = numpy.add.accumulate([rr.k for rr in reactionrules])
        k_max = k_array[-1]

        rnd = myrandom.uniform()
        i = numpy.searchsorted(k_array, rnd * k_max)

        return reactionrules[i]

    def create_new_shell(self):		# needs to be overloaded in subclasses
	pass				# creates an appropriate shell object


class ProtectiveDomain(Domain):
# Protective Domains are the proper eGFRD domains. The shell associated with a Protective
# Domains really exclude other particles from the simulation. These Domains are Singles and Pairs
# which have only one shell. Multi's on the other hand, have multiple shells and can still leave
# their shell within a timestep, breaking the principle of a Protective Domain

    def __init__(self, domain_id):
	Domain.__init__(self, domain_id)

	self.shell_id = None
	self.shell = None
	self.num_shells = 1		# all protective domains have only one shell

    def get_shell_id_shell_pair(self):
	return (self.shell_id, self.shell)

    def set_shell_id_shell_pair(self, shell_id_shell_pair):
	self.shell_id = shell_id_shell_pair[0]
	self.shell = shell_id_shell_pair[1]

    shell_id_shell_pair = property(get_shell_id_shell_pair, set_shell_id_shell_pair)

    def get_shell_list(self):
	return [(self.shell_id, self.shell), ]

    shell_list = property(get_shell_list)

    def get_shell_radius(self):		# returns the radius of the shell
	return self.shell.shape.radius

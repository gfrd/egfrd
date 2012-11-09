#!/usr/env python

import sys
import math
import numpy

from array import *

import ConfigParser as CP

from _gfrd import (
    Event,
    EventScheduler,
    Particle,
    SphericalShell,
    CylindricalShell,
    DomainIDGenerator,
    ShellIDGenerator,
    DomainID,
    ParticleContainer,
    CuboidalRegion,
    SphericalSurface,
    CylindricalSurface,
    DiskSurface,
    PlanarSurface,
    Surface,
    _random_vector,
    Cylinder,
    Disk,
    Sphere,
    Plane,
    NetworkRulesWrapper,
    )

from gfrdbase import *
from gfrdbase import DomainEvent
import model
import _gfrd

import logging


__all__ = [ 'save_state', 'load_state' ]


log = logging.getLogger('ecell')


def save_state(simulator, filename):
    

    #### BURST ALL DOMAINS ####
    # This is an important precondition for this
    # save routine to work properly.
    # In particular running this without bursting
    # all the domains will mess up the EventScheduler
    # because it will reconstruct it with different
    # event IDs
    
    if __debug__:
        log.info('Attempting to save the state of the simulator: running burst_all_domains()')

    simulator.burst_all_domains()

    # Re-seed the random number generator
    # (with a randomly drawn seed)
    # This seed has to be saved and used to
    # (re)seed the newly setup simulator when
    # loading to ensure that we will reproduce
    # the same trajectory.
    new_seed = simulator.rng.uniform_int(0, 1000000)
    simulator.reset_seed(new_seed)

    # Reset counters
    N_structure_types = 0
    N_species         = 0
    N_rules           = 0

    #### DEFINE THE CONFIG PARSER OBJECT ####
    cp = CP.ConfigParser()
    cp.optionxform = str # for upper case option names

    # Now add the relevant info in separate sections
    #### MODEL, WORLD and SEED ####
    cp.add_section('MODEL')
    cp.set('MODEL', 'world_size', simulator.world.world_size)

    cp.add_section('WORLD')
    cp.set('WORLD', 'world_size', simulator.world.world_size)
    cp.set('WORLD', 'matrix_size', simulator.world.matrix_size)

    cp.add_section('SEED')
    cp.set('SEED', 'seed', new_seed)

    #### STRUCTURE TYPES ####
    if __debug__:
        log.info('Saving structure types...')

    structure_type_names = [] # will store the names in a lookup list
    for structure_type in simulator.get_structure_types():

        id_int = id_to_int(structure_type.id)

        if structure_type['name'] != '' :
            name = structure_type['name']

        elif id_int == 1 : # the std. structure type
            name = 'default'

        sectionname = 'STRUCTURETYPE_' + str(id_int)
        cp.add_section(sectionname)        

        cp.set(sectionname, 'id', id_int)
        cp.set(sectionname, 'name', name)

        # Remember name for later
        structure_type_names.append( (id_int, name) )

        # Update counter
        N_structure_types += 1

    #### SPECIES ####
    if __debug__:
        log.info('Saving particle species...')

    for species in simulator.get_species():

        id_int = id_to_int(species.id)
        sectionname = 'SPECIES_' + str(id_int)
        cp.add_section(sectionname)

        cp.set(sectionname, 'id', id_int)
        #cp.set(sectionname, 'name', species.name]) # TODO SpeciesInfo does not have this property for some reason...
        cp.set(sectionname, 'radius', species.radius)
        cp.set(sectionname, 'D', species.D)
        cp.set(sectionname, 'v', species.v)
        cp.set(sectionname, 'structure_type_id', id_to_int(species.structure_type_id))

        # Update counter
        N_species += 1

    #### RULES ####
    if __debug__:
        log.info('Saving rules network...')

    extracted_rules = []  # to avoid double-extraction
    for species0 in simulator.get_species():

        # Reactions with one reactant (unimolecular or unbinding)
        rules = simulator.network_rules.query_reaction_rule(species0)
        
        for rr in rules:

              if not rr.id in extracted_rules:

                  rate = rr.k

                  r0 = rr.reactants[0]
                  reactant_ids = ( id_to_int(r0), )

                  if len(rr.products) == 1:

                      rtype = 'unimolecular'
                      p0 = rr.products[0]
                      product_ids  = ( id_to_int(p0), )

                  elif len(rr.products) == 2:

                      rtype = 'unbinding'
                      p0 = rr.products[0]
                      p1 = rr.products[1]
                      product_ids  = ( id_to_int(p0), id_to_int(p1) )

                  else:
                      raise LoadSaveError('Unimolecular or unbinding reaction rule has more than two products. Something is utterly wrong.')

                  sectionname = 'REACTIONRULE_' + str(rr.id)
                  cp.add_section(sectionname)

                  cp.set(sectionname, 'type', rtype)
                  cp.set(sectionname, 'reactant_ids', reactant_ids)
                  cp.set(sectionname, 'product_ids', product_ids)
                  cp.set(sectionname, 'rate', rate)

                  extracted_rules.append(rr.id)

        # Reactions with two reactants (particle species)
        for species1 in list(simulator.get_species()) + list(simulator.get_structure_types()):
              
            rules = simulator.network_rules.query_reaction_rule(species0, species1)

            for rr in rules:

                if not rr.id in extracted_rules and not len(rr.products)==0:

                    rate = rr.k

                    # Check whether we have a surface binding reaction.
                    # Since only species1 loops over StructureTypes we only have to check
                    # whether species1 is a SstructureType or not.
                    if species1.id in [st.id for st in list(simulator.get_structure_types())]:
                        # We must construct an id list from the iterator because otherwise
                        # the comparison does not work
                        rtype = 'binding_surface'

                    else:
                        rtype = 'binding_particle'

                    # Make sure that the first reactant stores species1, which is the surface
                    # in case we have a surface reaction.
                    if rr.reactants[0] == species1.id:
                        r0 = rr.reactants[0]
                        r1 = rr.reactants[1]
                    else:
                        # reverse order
                        r0 = rr.reactants[1]
                        r1 = rr.reactants[0]

                    reactant_ids = ( id_to_int(r0), id_to_int(r1) )

                    # Check whether there is a legal number of products
                    if len(rr.products) == 0:

                        product_ids  = None
                        rtype = 'annihilation'

                    elif len(rr.products) == 1:

                        p0 = rr.products[0]
                        product_ids  = ( id_to_int(p0), )

                    else:
                        raise LoadSaveError('Binding reaction rule has %s product(s). Something is utterly wrong.' % str(len(rr.products)) )           

                    # Construct the section info
                    sectionname = 'REACTIONRULE_' + str(rr.id)
                    cp.add_section(sectionname)

                    cp.set(sectionname, 'type', rtype)
                    cp.set(sectionname, 'reactant_ids', reactant_ids)
                    cp.set(sectionname, 'product_ids', product_ids)
                    cp.set(sectionname, 'rate', rate)

                    extracted_rules.append(rr.id)

    #### STRUCTURES ####
    if __debug__:
        log.info('Saving structures...')

    for structure in simulator.get_structures():

        id_int  = id_to_int(structure.id)
        sid_int = id_to_int(structure.sid)
        parent_id_int = id_to_int(structure.structure_id)
            # the structure_id property of structures is the parent structure id

        # Find the name of the structure type
        for i, name in structure_type_names: ## FIXME replace by dictionary?
            if i == sid_int:
                st_name = name

        # Find the keyword for the geometric type of this structure
        for object_type, key in structure_keywords:
            if isinstance(structure, object_type):
                geometry_key = key

        sectionname = 'STRUCTURE_' + geometry_key + '_' + str(id_int)
        cp.add_section(sectionname)

        cp.set(sectionname, 'id', id_int)
        cp.set(sectionname, 'parent_id', parent_id_int)
        cp.set(sectionname, 'structure_type_id', sid_int)
        cp.set(sectionname, 'structure_type_name', st_name)
        cp.set(sectionname, 'name', structure.name)

        cp.set(sectionname, 'position', list(structure.shape.position))            

        if isinstance(structure, CuboidalRegion):
            
            cp.set(sectionname, 'unit_x', list(structure.shape.unit_x))
            cp.set(sectionname, 'unit_y', list(structure.shape.unit_y))
            cp.set(sectionname, 'unit_z', list(structure.shape.unit_z))
            cp.set(sectionname, 'half_extent', list(structure.shape.half_extent))

        elif isinstance(structure, SphericalSurface):

            cp.set(sectionname, 'radius', structure.shape.radius)

        elif isinstance(structure, CylindricalSurface):

            cp.set(sectionname, 'unit_z', list(structure.shape.unit_z))
            cp.set(sectionname, 'radius', structure.shape.radius)
            cp.set(sectionname, 'half_length', structure.shape.half_length)

        elif isinstance(structure, DiskSurface):

            cp.set(sectionname, 'unit_z', list(structure.shape.unit_z))
            cp.set(sectionname, 'radius', structure.shape.radius)

        elif isinstance(structure, PlanarSurface):

            cp.set(sectionname, 'unit_x', list(structure.shape.unit_x))
            cp.set(sectionname, 'unit_y', list(structure.shape.unit_y))
            cp.set(sectionname, 'unit_z', list(structure.shape.unit_z))
            cp.set(sectionname, 'half_extent', list(structure.shape.half_extent))


    #### STRUCTURE CONNECTIVITY ####
    if __debug__:
        log.info('Saving structure connectivity information...')

    for structure in simulator.get_structures():

        if isinstance(structure, PlanarSurface):

            id_int  = id_to_int(structure.id)

            sectionname = 'STRUCTURECONNECTION_' + str(id_int)
            cp.add_section(sectionname)

            # Retrieve the neighbor structure IDs for this structure
            for n in range(0, 4):

                nid = simulator.world.get_neighbor_id(structure, n)
                nid_int  = id_to_int(nid)

                cp.set(sectionname, 'neighbor_id_'+str(n), nid_int)

    #### PARTICLES ####
    if __debug__:
        log.info('Saving particles...')

    pid_particle_pairs = list(simulator.world)
    pid_particle_pairs.sort()

    for pp in pid_particle_pairs:

        pid = pp[0]
        particle = pp[1]

        species_id_int   = id_to_int(particle.sid)
        structure_id_int = id_to_int(pp[1].structure_id)

        pid_int = id_to_int(pid)

        sectionname = 'PARTICLE_' + str(pid_int)
        cp.add_section(sectionname)

        cp.set(sectionname, 'id', pid_int)
        cp.set(sectionname, 'species_id', species_id_int)
        cp.set(sectionname, 'radius', particle.radius)
        cp.set(sectionname, 'D', particle.D)
        cp.set(sectionname, 'v', particle.v)
        cp.set(sectionname, 'position', list(particle.position))
        cp.set(sectionname, 'structure_id', structure_id_int)

    #### PARTICLE ORDER IN SCHEDULER ####
    # Pop the events from the scheduler and remember
    # them to push them back again afterwards.
    # This seems awkward but it is convenient because
    # we will do the same procedure when rebuilding the
    # system at loading / read-in.
    # Event IDs will change but this does not matter
    # for the simulated trajectory as such.    
    if __debug__:
        log.info('Saving particle order in scheduler...')

    eventlist = []
    scheduler_order = []
    while not simulator.scheduler.size == 0 :
        
        id, event = simulator.scheduler.pop()
        domain = simulator.domains[event.data]

        if __debug__:
            log.info('  Removed domain %s, particle id = %s from scheduler.' % (domain, domain.pid_particle_pair[0]) )

        # Store the events in a list
        # inserting from the beginning; thus it
        # will already have the right order for
        # reinsertion
        eventlist.insert(0, (event, domain))

        pid, particle = domain.pid_particle_pair        
        pid_int = id_to_int(pid)

        # Also output the scheduler order in reverse
        # fashion already, so we don't have to bother later
        scheduler_order.insert(0, pid_int)
    
    # Put the events back into the scheduler (in reverse order)
    # Note that the eventlist already has been correctly inverted at this point
    assert simulator.scheduler.size == 0   # just to be sure
    for (event, domain) in list(eventlist):
        if __debug__:
            log.info('  Re-instering domain %s, particle id = %s into scheduler.' % (domain, domain.pid_particle_pair[0]) )
        event_id = simulator.scheduler.add(DomainEvent(event.time, domain))
        # The event_id is changed in this process, so don't forget to update the domain
        domain.event_id = event_id

    # To be sure, check that event IDs in reconstructed scheduler
    # correspond to event IDs of the domains
    simulator.check_domains()

    # Create a new section in the save file
    sectionname = 'SCHEDULER'
    cp.add_section(sectionname)
    cp.set(sectionname, 'particle_order', list(scheduler_order))
    cp.set(sectionname, 't', simulator.t)
    cp.set(sectionname, 'dt', simulator.dt)
    cp.set(sectionname, 'step_counter', simulator.step_counter)

    #### ADD COUNTERS TO MODEL SECTION ####
    if __debug__:
        log.info('Saving object counters...')

    cp.set('MODEL', 'N_structure_types', N_structure_types)
    cp.set('MODEL', 'N_species', N_species)

    #### WRITE FILE ####
    if __debug__:
        log.info('Writing file...')

    with open(filename, 'wb') as outfile:
        cp.write(outfile)



def load_state(filename):
    """ Loads a file previously generated with save_state() and
        constructs a model, world and info needed to construct the 
        EGFRDSimulator.

        Returns: world, seed, time_info

        time_info is a tuple containing (t, dt, step_counter)
        where t, dt and step_counter are the values of these
        quantities at the moment of output / saving.
    """

    if __debug__:
        log.info('Attempting to load simulator state from file %s.' % str(filename))

    #### DEFINE THE CONFIG PARSER OBJECT ####
    cp = CP.ConfigParser()
    cp.optionxform = str # for upper case option names

    # Load the saved file
    # TODO Check for right format
    cp.read(filename)
  
    if __debug__:
        log.info('Sections in file: ' + str(cp.sections()) )

    #### FIRST GET GLOBAL INFO ####
    if __debug__:
        log.info('Reading global info...')

    world_size  = cp.getfloat('WORLD', 'world_size')
    matrix_size = cp.getint('WORLD',   'matrix_size')
    N_structure_types = cp.getint('MODEL', 'N_structure_types')
    N_species   = cp.getint('MODEL', 'N_species')
    seed        = cp.getint('SEED',  'seed')

    #### CREATE THE WORLD AND THE MODEL ####
    if __debug__:
        log.info('Creating the model...')

    m = model.ParticleModel(world_size)


    #### STRUCTURE_TYPES ####
    if __debug__:
        log.info('Reconstructing structure types...')

    structure_types_dict = {}  # will map the old (read-in) ID to the StructureType

    # First get the default structure type already defined at model creation
    # and add it to the internal dictionary
    def_structure_type_id = id_to_int(m.get_def_structure_type_id())
    def_structure_type    = m.get_structure_type_by_id(m.get_def_structure_type_id())
    structure_types_dict[def_structure_type_id] = def_structure_type
    if __debug__:
        log.info('  Added %s' % str(def_structure_type) )

    structure_types_sections = filter_sections(cp.sections(), 'STRUCTURETYPE')
    # Now read in the the structure types and add them to the new model
    # in the right order.
    # Note that add_structure_type() will create IDs automatically,
    # so it is important that we sort the sectionnames list by the ids.
    for sectionname in sorted(structure_types_sections, key = lambda name : name_to_int(name)):
        
        id      = cp.getint(sectionname, 'id')
        name    = cp.get(sectionname, 'name')        

        # Assert that the read id corresponds to the one in the sectionname
        assert name_to_int(sectionname) == id

        # Create a new StructureType object and add it to the model
        structure_type = _gfrd.StructureType()
        structure_type['name'] = name

        if id != def_structure_type_id:  # Skip adding the default structure_type a second time
            m.add_structure_type(structure_type)
            assert id_to_int(structure_type.id) == id            
            structure_types_dict[id] = structure_type
            if __debug__:
                log.info('  Added %s' % str(structure_type) )

    if __debug__:
        log.info('  structure_types_dict = %s' % str(structure_types_dict) )


    #### SPECIES ####
    if __debug__:
        log.info('Reconstructing species...')

    # Now the same for the species
    # Again we have to make sure we add the species in the right order => sort sections list
    species_dict = {}  # will map the old (read-in) ID to the Species

    species_sections = filter_sections(cp.sections(), 'SPECIES')
    for sectionname in sorted(species_sections, key = lambda name : name_to_int(name)):
        
        id      = cp.getint(sectionname, 'id')
        #name    = cp.getfloat(sectionname, 'name') # FIXME
        radius  = cp.getfloat(sectionname, 'radius')
        D       = cp.getfloat(sectionname, 'D')
        v       = cp.getfloat(sectionname, 'v')
        structure_type_id = cp.getint(sectionname, 'structure_type_id')

        # Assert that the read id corresponds to the one in the sectionname
        assert name_to_int(sectionname) == id

        # FIXME name cannot be output yet, so we have to make up one here
        name = 'Species_'+str(id)

        # Get the structure type for this species from the structure_type_id
        structure_type = structure_types_dict[structure_type_id]
        assert id_to_int(structure_type.id) == structure_type_id

        # Create a new Species object and add it to the model
        # Make sure to choose the right argument signature for each case
        if structure_type_id == def_structure_type_id:
            species = model.Species(name, D, radius)
        elif v == 0:
            species = model.Species(name, D, radius, structure_type)
        else:
            species = model.Species(name, D, radius, structure_type, v)

        m.add_species_type(species)
        assert id_to_int(species.id) == id
        species_dict[id] = species
        if __debug__:
            log.info('  Added %s' % str(species) )

    if __debug__:
        log.info('  species_dict = %s' % str(species_dict) )


    #### RULES ####
    if __debug__:
        log.info('Reconstructing rules network...')

    rules_dict = {}  # will map the old (read-in) ID to the reaction rule

    rules_sections = filter_sections(cp.sections(), 'REACTIONRULE')
    for sectionname in sorted(rules_sections, key = lambda name : name_to_int(name)):

        rtype        = cp.get(sectionname,      'type')
        rate         = cp.getfloat(sectionname, 'rate')
        reactant_ids = eval(cp.get(sectionname, 'reactant_ids')) # eval makes sure we get a list
        product_ids  = eval(cp.get(sectionname, 'product_ids'))

        rule = None
        if rtype == 'unimolecular':

            assert len(reactant_ids) == 1 and reactant_ids[0] != None and \
                   len(product_ids)  == 1 and product_ids[0]  != None, \
                     'Saved reaction rule does not have the proper signature: rule = %s' % sectionname

            # Get the involved species
            reactant = species_dict[int(reactant_ids[0])]
            product  = species_dict[int(product_ids[0])]

            # Create the rule
            rule = model.create_unimolecular_reaction_rule(reactant, product, rate)

        elif rtype == 'unbinding':

            assert len(reactant_ids) == 1 and reactant_ids[0] != None and \
                   product_ids[0] != None and product_ids[1]  != None, \
                     'Saved reaction rule does not have the proper signature: rule = %s' % sectionname

            # Get the involved species
            reactant = species_dict[int(reactant_ids[0])]
            product0 = species_dict[int(product_ids[0])]
            product1 = species_dict[int(product_ids[1])]

            # Create the rule
            rule = model.create_unbinding_reaction_rule(reactant, product0, product1, rate)

        elif rtype == 'binding_particle':
            
            assert reactant_ids[0] != None and reactant_ids[1] != None and \
                   len(product_ids) == 1   and product_ids[0]  != None, \
                     'Saved reaction rule does not have the proper signature: rule = %s' % sectionname

            # Get the involved species
            reactant0 = species_dict[int(reactant_ids[0])]
            reactant1 = species_dict[int(reactant_ids[1])]
            product   = species_dict[int(product_ids[0])]

            # Create the rule
            rule = model.create_binding_reaction_rule(reactant0, reactant1, product, rate)

        elif rtype == 'binding_surface':
            
            assert reactant_ids[0] != None and reactant_ids[1] != None and \
                   len(product_ids) == 1   and product_ids[0]  != None, \
                     'Saved reaction rule does not have the proper signature: rule = %s' % sectionname

            # Get the involved species
            structure_type = structure_types_dict[int(reactant_ids[0])]
            reactant       = species_dict[int(reactant_ids[1])]
            product        = species_dict[int(product_ids[0])]

            # Create the rule
            rule = model.create_binding_reaction_rule(reactant, structure_type, product, rate)

        else:

            raise LoadSaveError('Cannot create reaction rule when loading state: unknown reaction rule type.')

        # If we have successfully constructed a reaction rule, add it to the model
        # and to the internal rules dictionary
        if rule:
            
            m.add_reaction_rule(rule)
            rules_dict[name_to_int(sectionname)] = rule
            if __debug__:
                log.info('  Added %s, type = %s' % (rule, rtype) )

    if __debug__:
        log.info('  rules_dict = %s' % str(rules_dict) )

    #### CREATE THE WORLD ####
    if __debug__:
        log.info('Creating the world...')

    w = create_world(m, matrix_size)
    # This will come in handy later
    def_structure_id = w.get_def_structure_id()

    #### STRUCTURES ####
    if __debug__:
        log.info('Reconstructing structures...')

    structures_dict = {}  # will map the old (read-in) ID to the structure

    # The default structure should be always ID 1 and has to be added
    # from the beginning because it is created automatically at the init
    # of world.
    # To be sure we check that the default ID is indeed 1.
    # FIXME There must be a more elegant solution...
    assert id_to_int(def_structure_id) == 1
    def_structure = w.get_structure(def_structure_id)
    structures_dict[id_to_int(def_structure_id)] = def_structure

    structures_sections = filter_sections(cp.sections(), 'STRUCTURE')
    for sectionname in sorted(structures_sections, key = lambda name : name_to_int(name)):

        # Get the geometry key
        # TODO This assumes that the section name has the right signature
        #      => Implement a check
        geometry_key = sectionname.split(get_default_separator())[-2]

        # Find the geometric type of this structure according to its key
        for object_type, key in structure_keywords:
            if geometry_key == key:
                structure_object_type = object_type

        # Read in the information common to all object types
        id        = cp.getint(sectionname, 'id')        
        name      = cp.get(sectionname,    'name')
        parent_id = cp.getint(sectionname, 'parent_id')
        st_id     = cp.getint(sectionname, 'structure_type_id')        
        st_name   = cp.get(sectionname,    'structure_type_name')
        position  = vectorize(cp.get(sectionname, 'position'))
                    # eval makes sure we import a list

        # Assert that the read id corresponds to the one in the sectionname
        # This is particularly important in this part of the loading routine
        # because structures may be defined as substructures of previously
        # defined structures. Thus we want to define them in the order in
        # which they were defined before they have been saved to the file.
        assert name_to_int(sectionname) == id

        # Get the structure_type of this structure
        structure_type = structure_types_dict[st_id]
        if st_id != def_structure_type_id:
            assert structure_type['name'] == st_name

        # Get the parent structure from the structures dictionary
        # For this to work it is of critical importance that the order
        # in which the structures are defined is the same as the order
        # in which they were defined in before the output
        parent_structure = structures_dict[parent_id]

        structure = None
        if structure_object_type == CuboidalRegion:
        
            unit_x      = vectorize(cp.get(sectionname, 'unit_x'))
            unit_y      = vectorize(cp.get(sectionname, 'unit_y'))
            unit_z      = vectorize(cp.get(sectionname, 'unit_z'))            
            half_extent = vectorize(cp.get(sectionname, 'half_extent'))

            assert id == id_to_int(def_structure_id)
            # TODO Add more than default cuboidal regions
            # Right now we assume that there is only one which
            # is created automatically via create_world()            

        elif structure_object_type == SphericalSurface:
            
            radius = cp.getfloat(sectionname, 'radius')

            # Create the structure
            structure = model.create_spherical_surface(structure_type.id, name, \
                                                       position, radius, parent_structure.id)

        elif structure_object_type == CylindricalSurface:
            
            unit_z      = vectorize(cp.get(sectionname, 'unit_z'))            
            half_length = cp.getfloat(sectionname, 'half_length')
            radius      = cp.getfloat(sectionname, 'radius')

            # The cylinder creation function does not take the midpoint of
            # the cylinder (which is stored in position) and its half length
            # but the edge point and the total length. Thus we need:
            edge_pos = position - half_length * unit_z

            # Create the structure
            structure = model.create_cylindrical_surface(structure_type.id, name, edge_pos, radius, \
                                                         unit_z, 2.0*half_length, parent_structure.id)

        elif structure_object_type == DiskSurface:
            
            unit_z      = vectorize(cp.get(sectionname, 'unit_z'))            
            radius      = cp.getfloat(sectionname, 'radius')

            structure_type = structure_types_dict[st_id]            

            # Create the structure
            structure = model.create_disk_surface(structure_type.id, name, position, radius, \
                                                  unit_z, parent_structure.id)

        elif structure_object_type == PlanarSurface:

            unit_x      = vectorize(cp.get(sectionname, 'unit_x'))
            unit_y      = vectorize(cp.get(sectionname, 'unit_y'))
            unit_z      = vectorize(cp.get(sectionname, 'unit_z'))            
            half_extent = vectorize(cp.get(sectionname, 'half_extent'))

            # The plane creation function does not take the midpoint of
            # the plane (which is stored in position) and its half extents
            # but the corner point and the total extensions. Thus we need:
            corner_pos = position - half_extent[0] * unit_x - half_extent[1] * unit_y

            # Create the structure
            structure = model.create_planar_surface(structure_type.id, name, corner_pos, unit_x, unit_y, \
                                                    2.0*half_extent[0], 2.0*half_extent[1], parent_structure.id)

        # Add it to the world
        if structure:
            w.add_structure(structure)
            assert id_to_int(structure.id) == id # TODO
            structures_dict[id] = structure
            if __debug__:
                log.info('  Added %s, type = %s, id = %s, parent = %s' \
                         % (structure, structure_type, structure.id, structure.structure_id) )

    if __debug__:
        log.info('  structures_dict = %s' % str(structures_dict) )

    #### STRUCTURE CONNECTIONS ####
    if __debug__:
        log.info('Reconstructing structure connections...')

    connections_dict = {}  # will map the old (read-in) ID to the structure connection

    # Here we first read in all the connections and then figure out which sides are
    # connected in order to actually connect them.
    sc_sections = filter_sections(cp.sections(), 'STRUCTURECONNECTION')
    for sectionname in sorted(sc_sections, key = lambda name : name_to_int(name)):

        to_connect_id = name_to_int(sectionname)

        neighbor_id_0 = cp.getint(sectionname, 'neighbor_id_0')
        neighbor_id_1 = cp.getint(sectionname, 'neighbor_id_1')
        neighbor_id_2 = cp.getint(sectionname, 'neighbor_id_2')
        neighbor_id_3 = cp.getint(sectionname, 'neighbor_id_3')

        connections_dict[to_connect_id] = [neighbor_id_0, neighbor_id_1, \
                                           neighbor_id_2, neighbor_id_3] 

    # Figure out which sides are connected
    connections = []
    for to_connect_id, neighbor_list in connections_dict.iteritems():

        # For each neighbor_id find the side with which the neighbor is
        # connected to the structure with con_id
        for side in range(0,4):

            n_id = neighbor_list[side]
            nn_list = connections_dict[n_id]

            for n_side in range(0,4):
            
                nn_id = nn_list[n_side]
                if nn_id == to_connect_id:
                # We have found the neighbor_id that links to con_id
                    connections.append( (to_connect_id, side, n_id, n_side) )
                    break

    if __debug__:
        log.info('  connections = %s' % str(connections) )

    # Establish the connections
    for struct_id0, side0, struct_id1, side1 in connections:

        struct0 = structures_dict[struct_id0]
        struct1 = structures_dict[struct_id1]

        w.connect_structures(struct0, side0, struct1, side1)
        if __debug__:
            log.info('  Connected side %s of structure %s (id=%s) with side %s of structure %s (id=%s)' \
                     % (side0, struct0, struct0.id, side1, struct1, struct1.id) )

    #### PARTICLES ####
    if __debug__:
        log.info('Reconstructing particles...')

    # Here we first read in the scheduler order, i.e. the
    # order with which particles were sitting in the scheduler
    # before the output, and the particles as such.
    # Then we create the particles in the right order (FILO principle).
    # Note that the particle_order list of the SCHEDULER section
    # already contains the particles in reverse order as compared
    # to their order in the scheduler before saving, i.e. they
    # will be added to the system in the right order here.

    # Get the particle order and time counter values at output
    if __debug__:
        log.info('  Obtaining particle order in scheduler...')
    assert cp.has_section('SCHEDULER')
    particle_order = eval(cp.get('SCHEDULER', 'particle_order'))
    t  = cp.getfloat('SCHEDULER', 't')
    dt = cp.getfloat('SCHEDULER', 'dt')
    step_counter = cp.getfloat('SCHEDULER', 'step_counter')

    time_info = (t, dt, step_counter)

    if __debug__:
        log.info('  particle_order = %s' % str(particle_order) )
    
    # Read in the particle information
    if __debug__:
        log.info('  Reading particles info...')
    particles_dict = {}  # will map the old (read-in) ID to the particle info

    particle_sections = filter_sections(cp.sections(), 'PARTICLE')
    for sectionname in particle_sections:
    
        pid = cp.getint(sectionname, 'id')
        sid = cp.getint(sectionname, 'species_id')        
        radius   = cp.getfloat(sectionname, 'radius')
        D        = cp.getfloat(sectionname, 'D')
        v        = cp.getfloat(sectionname, 'v')
        position = vectorize(cp.get(sectionname, 'position'))
        str_id   = cp.getint(sectionname, 'structure_id')

        # Store the info in an internal particle dictionary
        # and add it to the collection (particles_dict)
        particle = { 'id'  : pid,
                     'species_id' : sid,
                     'radius'     : radius,
                     'D'   : D,
                     'v'   : v,
                     'position'   : position,
                     'structure_id' : str_id
                   }
        particles_dict[pid] = particle

    if __debug__:
        log.info('  Placing particles...')
    for pid in particle_order:

        particle = particles_dict[pid]
        particle_species = species_dict[ particle['species_id'] ]
        
        placed = place_particle(w, particle_species, particle['position'])
        # place_particle generates log info when placing


    # Return the world and the seed that was used 
    # to re-seed the simulator at output
    return w, seed, time_info



##########################
#### HELPER FUNCTIONS ####
##########################
def id_to_int(ID):
    """ Strips an eGFRD-type ID-specifier off everything
        but the integer ID numbers.

        I.e., PID(0:2) for example will be converted to 2.

        TODO Is the first number entry not important, too?
    """

    p = str(ID).partition(':')

    if p[2] != '':
        q = p[2].partition(')')

        if q[0] != '':
            return int(q[0])

    else:
        raise LoadSaveError('Could not extract number, probably the argument is not a valid ID.')


def get_default_separator():

    return '_'


def separate(sectionname):

    separator = get_default_separator()
    part = sectionname.partition(separator)

    return (part[0], part[2])

    # TODO Check for right tag format
    # raise LoadSaveError('Section name in input file did not have the required format.')


def filter_sections(sectionlist, tagstring):

    assert isinstance(tagstring, str)

    return [section for section in sectionlist if separate(section)[0] == tagstring]


def name_to_int(sectionname):

    separator = get_default_separator()

    return int(sectionname.split(separator)[-1])


def vectorize(input_list):

    return numpy.array(eval(input_list))


structure_keywords = [
        (CuboidalRegion,     'CUBE'),
        (SphericalSurface,   'SPHERE'),
        (CylindricalSurface, 'CYLINDER'),
        (DiskSurface,        'DISK'),
        (PlanarSurface,      'PLANE')
        ]


class LoadSaveError(Exception):
    pass


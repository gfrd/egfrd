class Single: 
    """ A Single contains:
        - 1 Particle of a certain Species
        - 1 shell (Sphere or Cylinder)
        - 1 Domain (Cartesian or Radial, wrapper around appropriate GF) 
        - public methods determineNextEvent and drawNewPosition

    """


class Pair:
    """A pair contains:
        - 2 Particles
        - 1 shell (Sphere or Cylinder)
        - 1 domain for center of mass (Cartesian or Radial)
        - 1 domain for interparticle vector (Composite: r and theta)
        - public methods determineNextEvent, drawNewPositions, drawNewCoM

    """


class EGFRDSimulator:
    """An EGFRDSimulator stores information about:

    surfaces:
        - surfaceList

    types of particles that can exist on each surface:
        - speciesList

    possible reactions each species can undergo:
        - reactionTypeMap1 (for the monomolecular reactions)
        - reactionTypeMap2 (for the bimolecular reactions)
        - interactionTypeMap (for the particle-surface interactions)

    position of each particle:
        - particleMatrix

    position and size of each shell:
        - sphereMatrix (for the spherical shells)
        - cylinderMatrix (for the cylindrical shells)

    order of events:
        - scheduler (contains Singles, Pairs and Multis)

    """
    def initialize:
        for particle in all particles:
            createSingle( particle )


    def step:
        # The event scheduler calls fireSingle, firePair or fireMulti, 
        # depending on the type of the top event.



    # singles

    def createSingle( particle ):
        # Create single with dt=0. The surface the particle is on determines 
        # the type of the single (CylindricalSurfaceSingle, 
        # PlanarSurfaceSingle, SphericalSingle).


    def createInteraction( single, surface ):
        # Create interaction single including the particle and the surface.
        # The surface the particle is on determines the type of the 
        # interaction single (CylindricalSurfaceInteraction, 
        # PlanarSurfaceInteraction).


    def fireSingle( single ):
        if eventType == REACTION:
            propagateSingle( single )
            fireSingleReaction( single )
            return
        if eventType == INTERACTION:
            propagateSingle( single )
            # The reactionType of the single distinguishes this case from the 
            # above in fireSingleReaction.
            fireSingleReaction( single )
            return

        # ESCAPE
        propagateSingle( single )

        # Find neighbors within burstVolume, and closestObject outside of 
        # burstVolume. BurstVolume depends on the surface the particle is on.  
        # Spherical by default.
        ( neighbors, closestObject, distanceToShellOfClosestObject ) =
            getNeighbors( position, burstVolume )

        if neighbors:
            burstedSingles = burstNonMultis( neighbor )
            if tryInteractionOrPairOrMulti( single, bursted ):
                return
            else:
                # If nothing was formed, recheck closest and restore shells.
                ( closestObject, distanceToShellOfClosestObject ) = 
                getClosestObj( position )

        updateSingle( single, closestObject, distanceToShellOfClosestObject )

        # Avoid livelock.
        restoreSingleShells( burstedSingles )
            

    def tryInteractionOrPairOrMulti( single1, neighbors ):
        # Try to make an interaction or a pair with single1 and the closest 
        # neighbor. If single1 and closest neighbor are too far apart, do 
        # nothing. If more than 2 objects are close together, make a multi.


    def burstSingle( single ):
        # There is a technical difference, not mentioned now.
        propagateSingle( single )


    def propagateSingle( single ):
        drawNewPosition()
        # If the single is an interaction single, remove it and create a 
        # normal single at the same position using createSingle().


    def restoreSingleShells( singles ):
        for single in singles:
            ( closestObject, distanceToShellOfClosestObject ) = getClosestObj()
            updateSingle( single, closestObject,
                          distanceToShellOfClosestObject )


    def updateSingle( single, closestObject, distanceToShellOfClosestObject ): 
        if closestObject is Single:
            radiusOfSingle = 
                calculateSingleShellSize( single, closestObject, 
                                          distanceToShellOfClosestObject )
        else:
            # closestObject is Pair or Multi or Surface.
            # It is harder to determine the optimal shellSize. Make single as 
            # big as possible for now.
            radiusOfSingle = distanceToShellOfClosestObject

        ( dt, eventType ) = determineNextEvent( )


    def calculateSingleShellSize( single, single2, 
                                  distanceToShellOfSingle2 ):
        # Optimal shellSize is half of the distance between the 2 particles if 
        # the radii and the diffusion constants are the same.
        return min( radius1 + sqrtD1 / ( sqrtD1 + sqrtD2 ) *
                              ( distanceBetweenParticles -
                                ( radius1 + radius2 ) ),
                    distanceToShellOfSingle2 )


    def fireSingleReaction( single ):
        if number of products is 1:
            if reactionType is SurfaceUnbindingReactionType:
                newpos = randomUnbindingSite( oldpos )
            elif reactionType is SurfaceBindingReactionType:
                newpos = projectedPoint( oldpos )
            else:
                # A -> B, no surfaces involved.
                newpos = oldpos

            createSingle( productParticle ) # at newpos

        elif number of products is 2:
            # A -> B + C, no surfaces involved.
            # Draw a random vector of a fixed length. The orientation depends 
            # on the surface the particle is on.
            length = ( particleRadius1 + particleRadius2 ) * 
                     MINIMAL_SEPERATION_FACTOR )
            vector = randomVector( length )

            newpos1 = oldpos + vector * ( D1 / D12 )
            newpos2 = oldpos - vector * ( D2 / D12 )

            createSingle( particle1 ) # at newpos1
            createSingle( particle2 ) # at newpos2



    # pairs

    def createPair( single1, single2, shellSize )
        # Create pair around the CoM of particle1 and particle2 with shellSize 
        # as radius. The surface the particles are on determines the type of 
        # the pair (CylindricalSurfacePair, PlanarSurfacePair, SphericalPair).


    def firePair( pair ):
        if eventType == SINGLE_REACTION:
            burstPair( pair )
            fireSingleReaction( reactingsingle )
        elif eventType == PAIR_REACTION:
            newCenterOfMass = drawNewCoM()
            createSingle( productParticle ) # at newCenterOfMass
        else:
            # Escape through either r or R.
            propagatePair( pair )


    def burstPair( pair ):
        # There is a technical difference, not mentioned now.
        propagatePair( pair )


    def propagatePair( pair ):
        drawNewPositions() # for both particles.



    # multis

    def createMulti( single, neighbors ):
        # Create multi and add the single and all objects in the neighbors 
        # list recursively. A multi has it's own Brownian Dynamics simulator.    
        addToMulti( single, multi )
        for neighbor in neighbors:
            addToMultiRecursive( neighbor, multi )


    def addToMultiRecursive( thisObject, multi ):
        if thisObject is Single:
            addToMulti( thisObject, multi )
            # First burst neighboring objects that have a shell within 
            # MULTI_SHELL_FACTOR * radius of thisObject.
            neighbors = burstNonMulti( neighbors )
            # Then select only those neighbors who have a particle that is 
            # within MULTI_SHELL_FACTOR * radius of thisObject. 
            (no pseudocode here)
            # Finally, add those neighbors recursively.
            for neighbor in neighbors:
                addToMultiRecursive( neighbor, multi )
        else:
            mergeMultis( thisObject, multi )


    def addToMulti( single, multi ):
        addParticle( particle )
        # Add a static shell with a radius that is the radius of the particle 
        # multiplied by MULTI_SHELL_FACTOR. When the particle leaves this 
        # shell, the multi is bursted (see fireMulti).  
        shellSize = radiusOfParticle * MULTI_SHELL_FACTOR
        addShell( position, shellSize )


    def mergeMultis( multi1, multi2 ):
        #Merge multi1 into multi2.


    def fireMulti( multi ):
        # Do one execution step of the Brownian Dynamics simulator. Not 
        # further discussed here.

        if reaction occured or particle escaped:
            burstMulti( multi )


    def burstMulti( multi ):
        for particle in all multi particles:
            createSingle( particle )


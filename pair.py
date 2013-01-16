from _gfrd import *
from constants import EventType
from _greens_functions import *
from greens_function_wrapper import *

from domain import (
    Domain,
    ProtectiveDomain)
from single import (
    NonInteractionSingle,
    InteractionSingle)
from multi import (
    Multi)
from shells import *
from utils import *

__all__ = [
    'CylindricalSurfacePair',
    'PlanarSurfacePair',
    'PlanarSurfaceTransitionPair',
    'SphericalPair',
    'Pair',
    'SimplePair',
    'MixedPair2D3D',
    'MixedPair1D3D',
    'MixedPair1DCap',
    ]

import logging
log = logging.getLogger('ecell')


class Pair(ProtectiveDomain, Others):
    """There are 3 types of pairs:
        * SphericalPair
        * PlanarSurfacePair
        * CylindricalSurfacePair

    """
    # CUTOFF_FACTOR is a threshold to choose between the real and 
    # approximate Green's functions.
    # H = 4.0: ~3e-5, 4.26: ~1e-6, 5.0: ~3e-7, 5.2: ~1e-7,
    # 5.6: ~1e-8, 6.0: ~1e-9
    CUTOFF_FACTOR = 5.6

    def __init__(self, domain_id, shell_id, rrs):
        # Calls required parent inits and sets variables.
        #
        # Sets:     multiplicity, single1, single2, pid_particle_pair1, 
        #           pid_particle_pair2, structure1, structure2, 
        #           interparticle_rrs, interparticle_ktot, com, iv, r0,
        #           sigma
        #
        # Requires nothing to be set.

        # Call init from parent class(es)
        ProtectiveDomain.__init__(self, domain_id, shell_id)
            # Note: for the Others superclass nothing is to be initialized.

        assert self.testShell       # assume that the testShell exists
        assert isinstance(self.testShell, testPair)  
        
        # Definitions
        self.multiplicity = 2

        self.single1 = self.testShell.single1
        self.single2 = self.testShell.single2 
        self.pid_particle_pair1 = self.testShell.pid_particle_pair1
        self.pid_particle_pair2 = self.testShell.pid_particle_pair2

        self.structure1 = self.testShell.structure1
        self.structure2 = self.testShell.structure2

        self.interparticle_rrs = rrs
        self.interparticle_ktot = self.calc_ktot(rrs)
        # copy constants from the testShell so that we don't have to recalculate them.
        self.com   = self.testShell.com
        self.iv    = self.testShell.iv
        self.r0    = self.testShell.r0
        self.sigma = self.testShell.sigma              

    def __del__(self):
        if __debug__:
            log.debug('del %s' % str(self))

    def get_D_tot(self):
        return self.pid_particle_pair1[1].D + \
               self.pid_particle_pair2[1].D
    D_tot = property(get_D_tot)
    D_r   = property(get_D_tot)

    def get_D_R(self):
        return (self.pid_particle_pair1[1].D *
                self.pid_particle_pair2[1].D) / self.D_tot
    D_R = property(get_D_R)

    def get_shell_direction(self):
        return self.testShell.get_orientation_vector()

    def initialize(self, t):
        # re-initialize the time parameters of the Domain to reuse it.
        # THIS IS PROBABLY USELESS
        self.last_time = t
        self.dt = 0
        self.event_type = None

    # draws a first event time and the (not fully specified) type of event
    def draw_com_escape_or_iv_event_time_tuple(self, r0):
        """Returns a (event time, event type, reactingsingle=None) tuple.
        
        """
        dt_com = draw_time_wrapper(self.com_greens_function())
        dt_iv = draw_time_wrapper(self.iv_greens_function(r0))  # uses the full solution to draw IV first event time
        if dt_com < dt_iv:
            return dt_com, EventType.COM_ESCAPE, None
        else:
            # Note: we are not calling pair.draw_iv_event_type yet, but 
            # postpone it to the very last minute (when this event is 
            # executed in fire_pair). So IV_EVENT can still be an iv 
            # escape or an iv reaction.
            return dt_iv, EventType.IV_EVENT, None

    def draw_single_reaction_time_tuple(self):
        """Return a (reaction time, event type, reactingsingle)-tuple.

        """
        dt_reaction1, event_type1 = self.single1.draw_reaction_time_tuple()
        dt_reaction2, event_type2 = self.single2.draw_reaction_time_tuple()
        if dt_reaction1 < dt_reaction2:
            return dt_reaction1, event_type1, self.single1
        else:
            return dt_reaction2, event_type2, self.single2

    def determine_next_event(self, r0):
        """Return a (event time, event type, reactingsingle)-tuple.

        """
        return min(self.draw_com_escape_or_iv_event_time_tuple(r0), 
                   self.draw_single_reaction_time_tuple()) 

    def draw_iv_event_type(self, r0):

        gf = self.iv_greens_function(r0)

        if self.dt < TIME_TOLERANCE**2:
            log.warning('Ridiculously small pair next-event time: pair.dt = %.30g < TIME_TOLERANCE^2, setting pair.dt = TIME_TOLERANCE^2 = %.30g' \
                       % (self.dt, TIME_TOLERANCE**2) )
            dt = TIME_TOLERANCE**2
        else:
            dt = self.dt

        event_kind = draw_event_type_wrapper(gf, dt)
        if event_kind == PairEventKind.IV_REACTION:
            return EventType.IV_REACTION
        elif event_kind == PairEventKind.IV_ESCAPE:
            return EventType.IV_ESCAPE
        raise NotImplemented()

    def draw_new_positions(self, dt, r0, old_iv, event_type):
        """Calculate new positions of the pair particles using a new 
        center-of-mass, a new inter-particle vector, and an old 
        inter-particle vector.

        """
        new_com = self.draw_new_com(dt, event_type)
        new_iv = self.draw_new_iv(dt, r0, old_iv, event_type)

        unit_z = self.get_shell_direction()
        newpos1, newpos2, sid1, sid2 = self.do_back_transform(new_com, new_iv,
                                                  self.pid_particle_pair1[1].D,
                                                  self.pid_particle_pair2[1].D,
                                                  self.pid_particle_pair1[1].radius,
                                                  self.pid_particle_pair2[1].radius,
                                                  self.structure1,
                                                  self.structure2,
                                                  unit_z, self.world)
        # sid = structure_id

        return newpos1, newpos2, sid1, sid2

    def draw_new_com(self, dt, event_type):
        # draws a new coordinate for the CoM in world coordinates
        if event_type == EventType.COM_ESCAPE:
            r = self.a_R
        else:
            gf = self.com_greens_function()
            r = draw_r_wrapper(gf, dt, self.a_R)

        # Add displacement to old CoM. This assumes (correctly) that 
        # r0=0 for the CoM. Compare this to 1D singles, where r0 is not  
        # necesseraly 0.

        # note that we need to make sure that the com and com_vector are in the structure to prevent
        # the particle leaving the structure due to numerical errors
        return self.com + self.create_com_vector(r)

    def draw_new_iv(self, dt, r0, old_iv, event_type):
        gf = self.choose_pair_greens_function(r0, dt)
        if event_type == EventType.IV_ESCAPE:
            r = self.a_r
        elif event_type == EventType.IV_REACTION:
            r = self.sigma      # maybe this should be zero (product particle is at CoM)
        else:
            r = draw_r_wrapper(gf, dt, self.a_r, self.sigma)

        # note that we need to make sure that the interparticle vector is in the structure to prevent
        # the particle leaving the structure due to numerical errors
        return self.create_interparticle_vector(gf, r, dt, r0, old_iv)

    def get_particles(self):
        return [self.pid_particle_pair1, self.pid_particle_pair2]
    particles = property(get_particles)

    def check(self):
    # perform internal consistency check 
        pass

    def __str__(self):
        sid1 = self.pid_particle_pair1[1].sid
        sid2 = self.pid_particle_pair2[1].sid
        name1 = self.world.model.get_species_type_by_id(sid1)["name"]
        name2 = self.world.model.get_species_type_by_id(sid2)["name"]
        return 'Pair[%s: %s, %s, (%s, %s)] with shell: %s' % (
            self.domain_id,
            self.pid_particle_pair1[0],
            self.pid_particle_pair2[0],
            name1, name2,
            self.shell)

class SimplePair(Pair):

    def __init__(self, domain_id, shell_id, rrs):
        # Calls required parent inits and sets variables.
        #
        # Sets:     a_R, a_r, strucutre, shell_size
        # Requires: r0, LD_MAX

        # Inits from parent classes 
        Pair.__init__(self, domain_id, shell_id, rrs) 
        
        # Definitions 
        self.structure = self.testShell.structure # structure on which both particles live                    
        
        shell_size = self.get_shell_size()                  # FIXME                      
                                                            # Needed by pair.init()?? No?
        self.a_R, self.a_r = self.determine_radii(self.r0, shell_size)
            # set the radii of the inner domains as a function of the outer protective domain
            # (self.r0 is determined in Pair.__init__)
            # determine_radii requires LD_MAX to be set


    @ classmethod
    def do_back_transform(cls, com, iv, D1, D2, radius1, radius2, structure1, structure2, unit_z, world):
    # here we assume that the com and iv are really in the structure and no adjustments have to be
    # made

        # For simple pairs the particles live on the same structure
        assert(structure1.id == structure2.id)
        D_tot = D1 + D2
        pos1 = com - iv * (D1 / D_tot)
        pos2 = com + iv * (D2 / D_tot)

        return pos1, pos2, structure1.id, structure2.id

    def get_shell_size(self):
        return self.shell.shape.radius

    def determine_radii(self, r0, shell_size):
        """Determine a_r and a_R from the size of the protective domain.

        Todo. Make dimension (1D/2D/3D) specific someday. Optimization only.

        """
        radius1 = self.pid_particle_pair1[1].radius
        radius2 = self.pid_particle_pair2[1].radius

        D1 = self.pid_particle_pair1[1].D
        D2 = self.pid_particle_pair2[1].D

        LD_MAX = self.LD_MAX # temporary value
        a_r_max = LD_MAX * (r0 - self.sigma) + self.sigma
        # a_r_max is the maximum size of a_r for a given maximum ratio of l/delta

        # Make sure that D1 != 0 to avoid division by zero in the followings.
        if D1 == 0:
            D1, D2 = D2, D1     # FIXME shouldn't we also here swap the particle radii if we're swapping diffusion contants?
            if __debug__:
                log.info('Swapping D1 and D2 because D1==0 (to avoid division by zero).')

        shell_size /= SAFETY    #??

        D_tot = D1 + D2
        D_geom = math.sqrt(D1 * D2)

        assert r0 >= self.sigma, \
            '%s;  r0 %g < sigma %g' % (self, r0, self.sigma)

        # equalize expected mean t_r and t_R.
        if ((D_geom - D2) * r0) / D_tot + shell_size +\
                math.sqrt(D2 / D1) * (radius1 - shell_size) - radius2 >= 0:
            Da = D1
            Db = D2
            radiusa = radius1
            radiusb = radius2
        else:
            Da = D2
            Db = D1
            radiusa = radius2
            radiusb = radius1


        #aR
        a_R = (D_geom * (Db * (shell_size - radiusa) + \
                         Da * (shell_size - r0 - radiusa))) /\
              (Da * Da + Da * Db + D_geom * D_tot)

        #ar
        a_r = (D_geom * r0 + D_tot * (shell_size - radiusa)) / (Da + D_geom)


        # Now if the planned domainsize for r is too large for proper convergence of the Green's
        # functions, make it the maximum allowed size
        if (a_r > a_r_max):
            a_r = a_r_max
            a_R = shell_size - radiusa - a_r * Da / D_tot
            if __debug__:
                log.info('Domainsize changed for convergence: a_r = a_r_max = %s, a_R = %s' % (a_r, a_R))

        assert a_R + a_r * Da / D_tot + radiusa >= \
               a_R + a_r * Db / D_tot + radiusb

        assert abs(a_R + a_r * Da / D_tot + radiusa - shell_size) \
            < 1e-12 * shell_size                        # here the shell_size is the relevant scale


        if __debug__:
          log.debug('shell_size = %g, a_r = %g, a_R = %g r0 = %g' % 
                 (shell_size, a_r, a_R, r0))
        if __debug__:
            tr = ((a_r - r0)**2) / (6 * self.D_r)       # the expected escape time of the IV
            if self.D_R == 0:
                tR = numpy.inf 
            else:
                tR = (a_R**2) / (6 * self.D_R)          # the expected escape time of the CoM
            log.debug('tr = %g, tR = %g' % (tr, tR))


        assert a_r > 0
        assert a_r > r0, '%g %g' % (a_r, r0)
        assert a_R > 0 or (a_R == 0 and (D1 == 0 or D2 == 0))

        return a_R, a_r

    def __str__(self):
        pass

class SphericalPair(SimplePair, hasSphericalShell):
    """2 Particles inside a (spherical) shell not on any surface.

    """
    def __init__(self, domain_id, shell_id, testShell, rrs):
        # Calls required parent inits and sets variables.
        #
        # Sets:         LD_MAX
        # Requires:     nothing to be set       

        assert isinstance(testShell, SphericalPairtestShell)
        
        # Parent classes
        hasSphericalShell.__init__(self, testShell, domain_id) 

        self.LD_MAX = numpy.inf # required by simplepair, set_radii           
        SimplePair.__init__(self, domain_id, shell_id, rrs)     # Always initialize AFTER hasSphericalShell

    def com_greens_function(self):
        # Green's function for centre of mass inside absorbing sphere.
        return GreensFunction3DAbsSym(self.D_R, self.a_R)

    def iv_greens_function(self, r0):
        # Green's function for interparticle vector inside absorbing 
        # sphere.  This exact solution is used for drawing times.
        return GreensFunction3DRadAbs(self.D_r, self.interparticle_ktot, r0,
                                              self.sigma, self.a_r)

    # selects between the full solution or an approximation where one of
    # the boundaries is ignored
    def choose_pair_greens_function(self, r0, t):
        distance_from_sigma = r0 - self.sigma
        distance_from_shell = self.a_r - r0

        threshold_distance = Pair.CUTOFF_FACTOR * \
            math.sqrt(6.0 * self.D_r * t)

        # if sigma reachable
        if distance_from_sigma < threshold_distance:
        
            # if shell reachable
            if distance_from_shell < threshold_distance:
                # near both a and sigma;
                # use GreensFunction3DRadAbs
                if __debug__:
                    log.debug('GF: normal')
                return self.iv_greens_function(r0)
            else:
                # near sigma; use GreensFunction3DRadInf
                if __debug__:
                    log.debug('GF: only sigma')
                return GreensFunction3DRadInf(self.D_r, self.interparticle_ktot, r0,
                                               self.sigma)

        # sigma unreachable
        else:
            if distance_from_shell < threshold_distance:
                # near a;
                if __debug__:
                    log.debug('GF: only a')
                return GreensFunction3DAbs(self.D_r, r0, self.a_r)
                
            else:
                # distant from both a and sigma; 
                if __debug__:
                    log.debug('GF: free')
                return GreensFunction3D(self.D_r, r0)

    def create_com_vector(self, r):
        return random_vector(r)

    def create_interparticle_vector(self, gf, r, dt, r0, old_iv): 

        if __debug__:
            log.debug("create_interparticle_vector: r=%g, dt=%g", r, dt)
        theta = draw_theta_wrapper(gf, r, dt)

        if theta == 0.0:
            # no rotation is necessary -> new_iv = new_iv
            new_iv = old_iv
        elif theta % numpy.pi != 0.0:
            # alternative calculation
            # rotate the old_iv to the new theta
            rotation_axis = crossproduct_against_z_axis(old_iv)
            rotation_axis = normalize(rotation_axis)
            new_iv = rotate_vector(old_iv, rotation_axis, theta)

            # rotate the new_iv around the old_iv with angle phi
            phi = myrandom.uniform() * 2 * Pi
            rotation_axis = normalize(old_iv)
            new_iv = rotate_vector(new_iv, rotation_axis, phi)
        else:
            # theta == numpi.pi -> just mirror the old_iv 
            new_iv = -old_iv

        # adjust length of the vector
        # note that r0 = length (old_iv)
        new_iv = (r/r0) * new_iv

        return new_iv

    def __str__(self):
        return 'Spherical' + Pair.__str__(self)


class PlanarSurfacePair(SimplePair, hasCylindricalShell):
    """2 Particles inside a (cylindrical) shell on a PlanarSurface. 
    (Hockey pucks).

    """
    def __init__(self, domain_id, shell_id, testShell, rrs):
        # Calls required parent inits and sets variables.
        #
        # Sets:         LD_MAX
        # Requires:     Nothing to be set.       

        assert isinstance(testShell, PlanarSurfacePairtestShell)
        hasCylindricalShell.__init__(self, testShell, domain_id)        
        
        self.LD_MAX = 20 # Required by SimplePair.__init__
        SimplePair.__init__(self, domain_id, shell_id, rrs)

    def com_greens_function(self):
        return GreensFunction2DAbsSym(self.D_R, self.a_R)

    def iv_greens_function(self, r0):
        return GreensFunction2DRadAbs(self.D_r, self.interparticle_ktot, r0,
                                      self.sigma, self.a_r)

    def choose_pair_greens_function(self, r0, t):
        # selects between the full solution or an approximation where one of
        # the boundaries is ignored
        distance_from_sigma = r0 - self.sigma
        distance_from_shell = self.a_r - r0

        threshold_distance = Pair.CUTOFF_FACTOR * \
            math.sqrt(4.0 * self.D_r * t)

        if distance_from_sigma < threshold_distance:
        
            if distance_from_shell < threshold_distance:
                # near both a and sigma;
                # use GreensFunction2DRadAbs
                if __debug__:
                    log.debug('GF2D: normal')
                return self.iv_greens_function(r0)
            else:
                # near sigma; use GreensFunction2DRadInf
                if __debug__:
                    log.debug('GF2D: only sigma')
                return self.iv_greens_function(r0)
                # Todo
                # return GreensFunction2DRadInf(self.D_r, self.interparticle_ktot, r0, self.sigma)
        else:
            if distance_from_shell < threshold_distance:
                # near a;
                if __debug__:
                    log.debug('GF2D: only a')
                return self.iv_greens_function(r0)
                # Todo
                # return GreensFunction2DAbs(self.D_r, r0, self.a_r)
                
            else:
                # distant from both a and sigma; 
                if __debug__:
                    log.debug('GF2D: free')
                return self.iv_greens_function(r0)
                # Todo
                # return GreensFunction2D(self.D_r, r0)

    def create_com_vector(self, r):
        x, y = random_vector2D(r)
        # project the com onto the surface unit vectors to make sure that the coordinates are in the surface
        return x * self.structure.shape.unit_x + y * self.structure.shape.unit_y

    def create_interparticle_vector(self, gf, r, dt, r0, old_iv): 

        if __debug__:
            log.debug("create_interparticle_vector: r=%g, dt=%g", r, dt)
        theta = draw_theta_wrapper(gf, r, dt)

        # note that r0 = length (old_iv)
        new_iv = (r/r0) * rotate_vector(old_iv, self.structure.shape.unit_z, theta)
        # note that unit_z can point two ways rotating the vector clockwise or counterclockwise
        # Since theta is symmetric this doesn't matter.

        # project the new_iv down on the unit vectors of the surface to prevent the particle from
        # leaving the surface due to numerical problem
        unit_x = self.structure.shape.unit_x
        unit_y = self.structure.shape.unit_y

        new_iv_x = unit_x * numpy.dot(new_iv, unit_x)
        new_iv_y = unit_y * numpy.dot(new_iv, unit_y)

        return new_iv_x + new_iv_y

    def __str__(self):
        return 'PlanarSurface' + Pair.__str__(self)


class PlanarSurfaceTransitionPair(SimplePair, hasSphericalShell):
    """2 Particles on two different but adjacent and orthogonal planar surfaces

    """
    # The test shell already defines structure1 and structure2 and should set
    # structure=structure1 by default.
    # We do all calculations here in structure1 assuming that pos2 has been 
    # "deflected back" properly into structure1 at test shell construction.
    # For safety each occurance of structure is replaced by structure1.
    # Otherwise this class is equal to PlanarSurfacePair inheriting from
    # hasSphericalShell instead of hasCylindricalShell and a modified 
    # do_back_transform() method.
    def __init__(self, domain_id, shell_id, testShell, rrs):
        # Calls required parent inits and sets variables.
        #
        # Sets:         structure1, strucutre2, structure, LD_MAX
        # Requires:     nothing to be set

        assert isinstance(testShell, PlanarSurfaceTransitionPairtestShell)

        hasSphericalShell.__init__(self, testShell, domain_id)
        self.LD_MAX = numpy.inf # TODO was 20. Why? Because of convergence issues! Maybe not really an issue here though
                                # Required by SimplePair.__init__        
        SimplePair.__init__(self, domain_id, shell_id, rrs)     # Always initialize AFTER hasSphericalShell

        # extra definitions for the TransitionPair
        self.structure1 = self.testShell.structure1
        self.structure2 = self.testShell.structure2
        # self.structure is equal to self.structure1 and inherited from the testShell
        # and SimplePair

    def com_greens_function(self):
        return GreensFunction2DAbsSym(self.D_R, self.a_R)

    def iv_greens_function(self, r0):
        return GreensFunction2DRadAbs(self.D_r, self.interparticle_ktot, r0,
                                      self.sigma, self.a_r)

    def choose_pair_greens_function(self, r0, t):
        # selects between the full solution or an approximation where one of
        # the boundaries is ignored
        distance_from_sigma = r0 - self.sigma
        distance_from_shell = self.a_r - r0

        threshold_distance = Pair.CUTOFF_FACTOR * \
            math.sqrt(4.0 * self.D_r * t)

        if distance_from_sigma < threshold_distance:
        
            if distance_from_shell < threshold_distance:
                # near both a and sigma;
                # use GreensFunction2DRadAbs
                if __debug__:
                    log.debug('GF2D: normal')
                return self.iv_greens_function(r0)
            else:
                # near sigma; use GreensFunction2DRadInf
                if __debug__:
                    log.debug('GF2D: only sigma')
                return self.iv_greens_function(r0)
                # Todo
                # return GreensFunction2DRadInf(self.D_r, self.interparticle_ktot, r0, self.sigma)
        else:
            if distance_from_shell < threshold_distance:
                # near a;
                if __debug__:
                    log.debug('GF2D: only a')
                return self.iv_greens_function(r0)
                # Todo
                # return GreensFunction2DAbs(self.D_r, r0, self.a_r)
                
            else:
                # distant from both a and sigma; 
                if __debug__:
                    log.debug('GF2D: free')
                return self.iv_greens_function(r0)
                # Todo
                # return GreensFunction2D(self.D_r, r0)

    def create_com_vector(self, r):
        x, y = random_vector2D(r)
        # project the com onto the surface unit vectors to make sure that the coordinates are in the surface
        # note that we use structure1 here, assuming that position of particle2 has been projected corretly
        # into this plane at construction of the test shell
        return x * self.structure1.shape.unit_x + y * self.structure1.shape.unit_y

    def create_interparticle_vector(self, gf, r, dt, r0, old_iv): 

        if __debug__:
            log.debug("create_interparticle_vector: r=%g, dt=%g", r, dt)
        theta = draw_theta_wrapper(gf, r, dt)

        # note that r0 = length (old_iv)
        new_iv = (r/r0) * rotate_vector(old_iv, self.structure1.shape.unit_z, theta)
        # note that unit_z can point two ways rotating the vector clockwise or counterclockwise
        # Since theta is symmetric this doesn't matter.

        # project the new_iv down on the unit vectors of the surface to prevent the particle from
        # leaving the surface due to numerical problem
        unit_x = self.structure1.shape.unit_x
        unit_y = self.structure1.shape.unit_y

        new_iv_x = unit_x * numpy.dot(new_iv, unit_x)
        new_iv_y = unit_y * numpy.dot(new_iv, unit_y)

        return new_iv_x + new_iv_y

    @ classmethod
    def do_back_transform(cls, com, iv, D1, D2, radius1, radius2, structure1, structure2, unit_z, world):

        # structure should be = structure1, but for safety we use the structures
        # inherited from the test shell

        # Calculate the new positions (still in structure1)
        D_tot = D1 + D2
        pos1 = com - iv * (D1 / D_tot)
        pos2 = com + iv * (D2 / D_tot)

        # Tentatively apply boundaries; this will also set the new structure IDs.
        # In some cases the new positions may end up on adjacent planar surfaces.
        # This may lead to an overlap which has to be checked for and removed if present
        # in the next step.
        new_pos1, new_sid1 = world.apply_boundary((pos1, structure1.id))
        new_pos2, new_sid2 = world.apply_boundary((pos2, structure1.id))

        # If the new positions lead to an overlap we have to enlarge the IV by a safety factor
        if world.distance(new_pos1, new_pos2) <= (radius1+radius2) * MINIMAL_SEPARATION_FACTOR :

            log.warn('Removing overlap in PlanarSurfaceTransitionPair')

            new_structure1 = world.get_structure(new_sid1)
            new_structure2 = world.get_structure(new_sid2)

            # Calculate the distances from the two new positions to the edge between the planes
            # in which the particles temporarily ended up
            # This is easily done by projecting new_pos1 into new_structure2 + vice versa and
            # computing the orthogonal components of the projections:
            _, (pos1_orth, _) = new_structure2.project_point(new_pos1)
            _, (pos2_orth, _) = new_structure1.project_point(new_pos2)

            l_iv = length(iv)
            assert(2.0*abs(pos1_orth*pos2_orth) < (l_iv*l_iv))  # That should never fail!
            iv_safety = 1.0/( 1.0 - 2.0*abs(pos1_orth*pos2_orth)/(l_iv*l_iv) )
            # Recalculate the interparticle vector and new positions
            iv_rescaled = iv_safety * iv
            pos1 = com - iv_rescaled * (D1 / D_tot)
            pos2 = com + iv_rescaled * (D2 / D_tot)

            # Re-apply boundaries
            new_pos1, new_sid1 = world.apply_boundary((pos1, structure1.id))
            new_pos2, new_sid2 = world.apply_boundary((pos2, structure1.id))               

        # TODO: Check that the new positions are in their new planes and 
        # maybe also that pos1 and pos2 are in structure1 in the first place?

        return new_pos1, new_pos2, new_sid1, new_sid2

    def draw_new_com(self, dt, event_type):
        # Draws a new coordinate for the CoM in world coordinates.
        # Since fire_pair_reaction() will call this function we have to
        # change the default draw_new_com() to ensure that the new CoM
        # is deflected into the right plane
        if event_type == EventType.COM_ESCAPE:
            r = self.a_R
        else:
            gf = self.com_greens_function()
            r = draw_r_wrapper(gf, dt, self.a_R)

        new_com = self.com + self.create_com_vector(r)

        if event_type == EventType.IV_REACTION:
            # Express the new CoM in plane1 and deflect it to plane2 if necessary
            # TODO This should be taken care of by apply_boundary() in the main routine            
            s1_ctr = self.structure1.shape.position
            new_com, changeflag = self.structure2.deflect(s1_ctr, new_com - s1_ctr)
            if __debug__:
                log.debug('Deflected new_com on structure2, changeflag=%s' % str(changeflag))

        # TODO: Check whether the new CoM is in the right structure
        return new_com

    def get_shell_size(self):
        return self.shell.shape.radius

    def __str__(self):
        return 'PlanarSurfaceTransition' + Pair.__str__(self)


class CylindricalSurfacePair(SimplePair, hasCylindricalShell):
    """2 Particles inside a (cylindrical) shell on a CylindricalSurface.  
    (Rods).

    """
    def __init__(self, domain_id, shell_id, testShell, rrs):
        # Calls required parent inits and sets variables.
        #
        # Sets:         LD_MAX
        # Requires:     nothing to be set        

        assert isinstance(testShell, CylindricalSurfacePairtestShell)

        hasCylindricalShell.__init__(self, testShell, domain_id)

        self.LD_MAX = numpy.inf # Required by SimplePair.__init__
        SimplePair.__init__(self, domain_id, shell_id, rrs)

    def get_shell_size(self):
        return self.shell.shape.half_length

    def get_v_tot(self):
        # 'direction' signifies the direction of the IV with respect to the surface unit_z
        # (which defines the direction of positive drift)
        direction = cmp(numpy.dot (self.iv, self.structure.shape.unit_z), 0)
        return (self.pid_particle_pair2[1].v - \
                self.pid_particle_pair1[1].v) * direction

    # NOTE the v_tot and v_r are now positive in the direction of the IV!
    v_tot = property(get_v_tot)
    
    v_r = property(get_v_tot)

    def get_v_R(self):
        # The v_R is always in the direction of the structure, no adjustment required.
        return (self.pid_particle_pair1[1].v * 
                self.pid_particle_pair2[1].D +
                self.pid_particle_pair2[1].v *
                self.pid_particle_pair1[1].D) / self.D_tot
    v_R = property(get_v_R)

    def com_greens_function(self):
        # The domain is created around r0, so r0 corresponds to r=0 within the domain
        return GreensFunction1DAbsAbs(self.D_R, self.v_R, 0.0, -self.a_R, self.a_R)

    def iv_greens_function(self, r0):
        return GreensFunction1DRadAbs(self.D_r, self.v_r, self.interparticle_ktot, r0, self.sigma, self.a_r)

    def choose_pair_greens_function(self, r0, t):
        # Todo
        return self.iv_greens_function(r0)

    def create_com_vector(self, z):
        # 'z' can be interpreted in two different ways here, it may a coordinate in the z direction or it may
        # be a displacement from the origin. In the first case it will already contain the drift information.
        if self.v_R == 0.0:
            # if there is no drift then we regard 'z' as a displacement from the origin. Although
            # it actually represents a coordinate when drawR was called. It is a displacement when z=a.
            # We randomize the direction.
            z = myrandom.choice(-1, 1) * z

        elif self.v_R != 0.0 and feq(z, self.a_R):
            # When there is drift and z=a, the 'z' actually represent the displacement from the origin and a
            # boundary must be chosen.

            # The Escape can be either to the left or to the right.
            # The side of escape should be decided on by the flux through both boundaries at the escape time
            gf = self.com_greens_function()
            event_kind = draw_event_type_wrapper(gf, self.dt)
            if event_kind == PairEventKind.IV_REACTION:     # IV_REACTION -> ESCAPE through 'left' boundary
                z = -z      # -self.a_R
            elif event_kind == PairEventKind.IV_ESCAPE:     # IV_ESCAPE -> ESCAPE through 'right' boundary
                #z = z      # self.a_R
                pass
            else:
                raise NotImplemented()
        else:
            # When there was drift and the particle was not at the boundary.
            # -> In this case the 'z' actually signifies a coordinate and nothing has to be done.
            pass

        # project the com onto the surface unit vector to make sure that the coordinates are in the surface
        return z * self.structure.shape.unit_z

    def create_interparticle_vector(self, gf, r, dt, r0, old_iv): 
        if __debug__:
            log.debug("create_interparticle_vector: r=%g, dt=%g", r, dt)

        # note that r0 = length (old_iv)
        new_iv = (r/r0) * old_iv

        # project the new_iv down on the unit vectors of the surface to prevent the particle from
        # leaving the surface due to numerical problem
        unit_z = self.structure.shape.unit_z
        new_iv_z = unit_z * numpy.dot(new_iv, unit_z)

        return new_iv_z

    def __str__(self):
        return 'CylindricalSurface' + Pair.__str__(self)


class MixedPair2D3D(Pair, hasCylindricalShell):

    def __init__(self, domain_id, shell_id, testShell, rrs):
        # Calls required parent inits and sets variables.
        #
        # Sets:         structure2D, structure3D, z_scaling_factor, a_R, a_r
        # Requires:     nothing to be set

        assert isinstance(testShell, MixedPair2D3DtestShell)
        hasCylindricalShell.__init__(self, testShell, domain_id)
        Pair.__init__(self, domain_id, shell_id, rrs)

        assert isinstance(self.testShell, testMixedPair2D3D)
        
        # Some definitions 
        self.structure2D = self.testShell.structure2D   # the surface on which particle1 lives
        self.structure3D = self.testShell.structure3D   # structure on which particle2 lives
        self.particle2D  = self.testShell.particle2D
        self.particle3D  = self.testShell.particle3D

        # initialize some useful constants
        self.z_scaling_factor = self.testShell.get_scaling_factor()

        # set the radii of the inner domains as a function of the outer protective domain
        self.a_R, self.a_r = self.determine_radii()
        
        # The latter is defined for completeness and used by check_domain() in egfrd.py:
        # (Note that that origin and target end up in same list.)
        self.origin_structure = self.testShell.origin_structure
        self.target_structure = self.testShell.target_structure

    def determine_radii(self):
        # Determines the dimensions of the domains used for the Green's functions 
        # from the dimensions of the cylindrical shell.
        # Note that the function assumes that the cylinder is dimensioned properly.
        radius2D = self.particle2D.radius
        radius3D = self.particle3D.radius
        D_2D = self.particle2D.D
        D_3D = self.particle3D.D

        shell_radius      = self.shell.shape.radius
        shell_half_length = self.shell.shape.half_length
        z_left            = radius2D
        z_right           = 2.0*shell_half_length - z_left

        # Use class methods to check dimensions of the cylinder
        r_check  = self.testShell.r_right(z_right)
        hl_check = (self.testShell.z_right(shell_radius) + self.testShell.z_left(shell_radius)) / 2.0
        if __debug__:
            assert feq(r_check, shell_radius) and feq(hl_check, shell_half_length), \
                'MixedPair2D3D: Domain did not obey scaling relationship. ' \
                'shell_radius = %s, shell_half_length = %s, radius_check = %s, half_length_check = %s, ' \
                'z_left = %s, z_right = %s' % \
                 (FORMAT_DOUBLE % shell_radius, FORMAT_DOUBLE % shell_half_length,
                  FORMAT_DOUBLE % r_check,      FORMAT_DOUBLE % hl_check,
                  FORMAT_DOUBLE % z_left,       FORMAT_DOUBLE % z_right          )

        # Partition the space in the protective domain over the IV and the CoM domains
        # The outer bound of the interparticle vector is set by the particle diffusing in 3D via:
        a_r = (z_right - radius3D) * self.z_scaling_factor

        # Calculate the maximal displacement of a particle given an interparticle vector bound a_r
        # taking into account the radius for both the 2D and 3D particle. This sets the remaining
        # space for the CoM diffusion.
        space_for_iv = max( (D_2D/self.D_tot) * a_r + radius2D,
                            (D_3D/self.D_tot) * a_r + radius3D)

        # It then determines the CoM vector bound via
        a_R = shell_radius - space_for_iv

        # Print the domain sizes and their estimated first passage times.
        if __debug__:
            tr = ((a_r - self.r0)**2.0) / (6.0 * self.D_r)  # the expected escape time of the iv
            if self.D_R == 0:
                tR = numpy.inf 
            else:
                tR = (a_R**2.0) / (4.0 * self.D_R)          # the expected escape time of the CoM
            log.debug('determine_radii: a_r= %s, tr= %s, a_R= %s, tR= %s, delta_tRr= %s' % \
                      (FORMAT_DOUBLE % a_r, FORMAT_DOUBLE % tr,
                       FORMAT_DOUBLE % a_R, FORMAT_DOUBLE % tR,
                       FORMAT_DOUBLE % (tr-tR) ))

#        assert feq(tr, tR), 'estimate first passage times were not equal'  ## FIXME

        # Some checks that shall never fail        
        assert (self.sigma < a_r) and (a_r < 2.0*shell_half_length * self.z_scaling_factor)
        assert (0 < a_R) and (a_R < shell_radius)

        return a_R, a_r

    @ classmethod 
    def do_back_transform(cls, com, iv, D1, D2, radius1, radius2, structure2D, structure3D, unit_z, world):
    # here we assume that the com and iv are really in the structure and no adjustments have to be
    # made
    #
    # structure3D is not used here but has to be passed because of the classmethod property        

        D_tot = D1 + D2
        weight1 = D1 / D_tot
        weight2 = D2 / D_tot

        z_btf = math.sqrt(D2/D_tot)     # z back transform factor

        # get the coordinates of the iv relative to the system of the surface (or actually the shell)
        iv_x = structure2D.shape.unit_x * numpy.dot(iv, structure2D.shape.unit_x)
        iv_y = structure2D.shape.unit_y * numpy.dot(iv, structure2D.shape.unit_y)

        iv_x_length = length(iv_x)
        iv_y_length = length(iv_y)

        # reflect the coordinates in the unit_z direction back to the side of the membrane
        # where the domain is. Note that it's implied that the origin of the coordinate system lies in the
        # plane of the membrane
        iv_z_length = abs(numpy.dot(iv, unit_z))   # FIXME maybe first project the shell unit_z onto the 
                                                   # surface unit_z to prevent numerical problems?
        # do the reverse scaling in z-direction
        iv_z_length_btf = iv_z_length * z_btf

        # The following is to avoid overlaps in case that the z-component of the IV is scaled down
        # so much that it would cause a particle overlap.
        # This can happen because the inner boundary in the space of transformed coordinates is not a
        # prolate spheroid, as it should, but approximated by a sphere. This can lead to IV exit points
        # within the prolate spheroid, i.e. within the sphere with radius sigma = radius1 + radius2 in
        # the real space after inverse transform.  FIXME This works, but is still somewhat of a HACK!
        min_iv_length = (radius1 + radius2) * MINIMAL_SEPARATION_FACTOR

        if (iv_x_length*iv_x_length + iv_y_length*iv_y_length + iv_z_length_btf*iv_z_length_btf) < min_iv_length*min_iv_length:

            z_safety_factor_sq = (min_iv_length*min_iv_length - iv_x_length*iv_x_length - iv_y_length*iv_y_length) \
                                 / (iv_z_length_btf*iv_z_length_btf)

            assert(z_safety_factor_sq >= 0.0)
            log.warn('MixedPair2D3D: Applying z_safety_factor = %s to enlarge too small interparticle vector.' % str(math.sqrt(z_safety_factor_sq)) )

            iv_z_length_btf = math.sqrt( z_safety_factor_sq ) * iv_z_length_btf
            

        # TODO Also check for membrane overlapping?

        # Construct the z-vector-component and the new position vectors
        iv_z = unit_z * iv_z_length_btf

        pos1 = com - weight1 * (iv_x + iv_y)            # the new 2D particle position
        pos2 = com + weight2 * (iv_x + iv_y) + iv_z     # the new 3D particle position

        # Class method; so don't use self.(..)
        return pos1, pos2, structure2D.id, structure3D.id

    @ classmethod
    def calc_z_scaling_factor(cls, D_2D, D_3D):
        # calculates the scaling factor to make the anisotropic diffusion problem into an isotropic one
        return math.sqrt( (D_2D + D_3D) / D_3D)

    def choose_pair_greens_function(self, r0, t):
        return self.iv_greens_function(r0)

    def com_greens_function(self):
        # diffusion of the CoM of a mixed pair is only two dimensional
        return GreensFunction2DAbsSym(self.D_R, self.a_R)

    def iv_greens_function(self, r0):
        # Diffusion of the interparticle vector is three dimensional
        # TODO Fix this ugly HACK to prevent particle overlap problem
        return GreensFunction3DRadAbs(self.D_r, self.interparticle_ktot, max(r0, self.sigma),
                                      self.sigma, self.a_r)

    def create_com_vector(self, r):
        x, y = random_vector2D(r)
        # project the com onto the surface unit vectors to make sure that the coordinates are in the surface
        return x * self.structure2D.shape.unit_x + y * self.structure2D.shape.unit_y

    def create_interparticle_vector(self, gf, r, dt, r0, old_iv):
        # FIXME This is actually the same method as for SphericalPair

        if __debug__:
            log.debug("create_interparticle_vector: r=%g, dt=%g", r, dt)

        theta = draw_theta_wrapper(gf, r, dt)   # draw_theta return a value between -Pi and Pi

        if theta == 0.0:
            new_iv = old_iv
        elif theta % numpy.pi != 0.0:
            # rotate the old_iv to the new theta
            rotation_axis = crossproduct_against_z_axis(old_iv)
            rotation_axis = normalize(rotation_axis)
            new_iv = rotate_vector(old_iv, rotation_axis, theta)

            # rotate the new_iv around the old_iv with angle phi
            phi = myrandom.uniform() * 2 * Pi
            rotation_axis = normalize(old_iv)
            new_iv = rotate_vector(new_iv, rotation_axis, phi)
        else:
            # theta == Pi -> just mirror the old_iv
            new_iv = -old_iv

        # adjust length of the vector to new r
        # note that r0 = length (old_iv)
        new_iv = (r/r0) * new_iv

        return new_iv

    def check(self):
    # perform internal consistency check 
        radius2D = self.particle2D.radius
        radius3D = self.particle3D.radius
        D_2D = self.particle2D.D
        D_3D = self.particle3D.D

        radius      = self.shell.shape.radius
        half_length = self.shell.shape.half_length

        # check that the shell obeys the scaling rules dimensions of the cylinder
        r_check  = self.testShell.r_right(2.0*half_length - radius2D)
        hl_check = (self.testShell.z_right(radius) + self.testShell.z_left(radius))/2.0
        assert feq(r_check, radius) and feq(hl_check, half_length), \
            'MixedPair2D3D: Domain did not obey scaling relationship. '

        # check that the CoM is in the surface
        assert numpy.dot( (self.com-self.structure2D.shape.position), self.structure2D.shape.unit_z) == 0.0
        # check that the IV is NOT in the surface

    def __str__(self):
        return 'Mixed2D3D' + Pair.__str__(self)

class MixedPair1D3D(Pair):

    @ classmethod
    def do_back_transform(cls, com, iv, D1, D2, radius1, radius2, structure1D, structure3D, unit_z, world):
    # here we assume that the com and iv are really in the structure and no adjustments have to be
    # made
    # unit_z is a putatively different unit_z than that available in the surface object and should not be
    # of use here (however the method is defined with the same parameter signature as elsewhere)
    #
    # note that structure1D is not a 1D object but the surface that the 1D particle is bound to
    #
    # structure3D is not used here but has to be passed because of the classmethod property

        D_tot = D1 + D2
        weight1 = D1 / D_tot
        weight2 = D2 / D_tot

        min_iv_r_length = radius2 + structure1D.shape.radius

        # get the coordinates of the iv relative to the system of the surface (or actually the shell)
        iv_z = structure1D.shape.unit_z * numpy.dot(iv, structure1D.shape.unit_z)
        iv_r = iv - iv_z

        # reflect the coordinates in the unit_z direction back to the side of the membrane
        # where the domain is. Note that it's implied that the origin of the coordinate system lies in the
        # plane of the membrane
        iv_r_length = length(iv_r)
        # do the reverse scaling
        iv_r_length_new = iv_r_length / cls.calc_r_scaling_factor(D1, D2)

        # if the particle is overlapping with the cylinder, make sure it doesn't
        if iv_r_length_new < min_iv_r_length:
            iv_r_length_new = min_iv_r_length * MINIMAL_SEPARATION_FACTOR

        iv_r =  (iv_r_length_new/iv_r_length) * iv_r

        pos1 = com - weight1 * iv_z
        pos2 = com + weight2 * iv_z + iv_r

        # Class method; so don't use self.(..)
        return pos1, pos2, structure1D.id, structure1D.id

    @ classmethod
    def calc_r_scaling_factor(cls, D1, D2):
        return math.sqrt((D1 + D2)/D2)


class MixedPair1DCap(Pair, hasCylindricalShell):
    # This domain is for the interaction of a particle drifting and diffusing on a cylinder (1D particle)
    # and another (static) particle which is bound to a cap of the same cylinder.
    # Since the cap particle is immobile, we effectively have a one-particle problem
    # and the CoM is not used. The displacement of the 1D particle is sampled via iv_greens_function.
    # However, we have to make this a pair domain because we have two particles involved which
    # can decay (fire single reactions) independently.
    # TODO: Is the new pos. of the cap particle sampled correctly in case of a decay event?

    def __init__(self, domain_id, shell_id, testShell, rrs):
        # Calls required parent inits and sets variables.
        #
        # Sets:             particle1D, structure1D, cap_particle, cap_structure
        # Requires:         nothing to be set       

        assert isinstance(testShell, MixedPair1DCaptestShell)
        hasCylindricalShell.__init__(self, testShell, domain_id)
        Pair.__init__(self, domain_id, shell_id, rrs)

        assert isinstance(self.testShell, MixedPair1DCaptestShell)

        # Some useful definitions 
        self.particle1D    = self.testShell.particle1D
        self.structure1D   = self.testShell.structure1D
        self.cap_particle  = self.testShell.cap_particle
        self.cap_structure = self.testShell.cap_structure

        # The latter is defined for completeness and used by check_domain() in egfrd.py:
        self.origin_structure = self.structure1D
        self.target_structure = self.cap_structure

        # The cap_particle should be immobile; check to be sure:
        assert self.cap_particle.D == 0
        assert self.cap_particle.v == 0
        # D_r and D_R are inherited from Pair class;
        # This results in D_r = particle1D.D and D_R = 0

        self.LD_MAX = numpy.inf

    def get_shell_size(self):   # TODO Is that even used?
        return self.shell.shape.half_length

    def get_a_r(self):
      
        ref_pt = self.testShell.get_referencepoint()
        assert ( all(self.cap_structure.shape.position - ref_pt == 0) )

        return abs(self.testShell.dz_right) - self.particle1D.radius
    a_r = property(get_a_r)
    a_R = property(get_a_r)

    def get_v_tot(self):
        # 'direction' signifies the direction of the IV with respect to the surface unit_z
        # (which defines the direction of positive drift).
        # Assuming self.pid_particle_pair2[1].v == 0 here, i.e. relative drift = drift of 1D part.
        direction = cmp(numpy.dot (self.iv, self.structure1D.shape.unit_z), 0)
        return self.particle1D.v * direction

    # NOTE the v_tot and v_r are now positive in the direction of the IV!
    v_tot = property(get_v_tot)    
    v_r = property(get_v_tot)

    def get_v_R(self):
        # With our transform D_R and v_R are zero if one of the particles is immobile
        return 0.0
    v_R = property(get_v_R)

    def com_greens_function(self):
        # The Green's function used is largely irrelevant, because D_R=0 and v_R=0        
        return GreensFunction1DAbsAbs(self.D_R, self.v_R, 0.0, -self.a_R, self.a_R)

    def iv_greens_function(self, r0):
        if( numpy.isinf(self.interparticle_ktot) ):
            return GreensFunction1DAbsAbs(self.D_r, self.v_r, r0, self.sigma, self.a_r)
        else:
            return GreensFunction1DRadAbs(self.D_r, self.v_r, self.interparticle_ktot, r0, self.sigma, self.a_r)

    def choose_pair_greens_function(self, r0, t):
        # Todo
        return self.iv_greens_function(r0)

    @ classmethod
    def do_back_transform(cls, com, iv, D1, D2, radius1, radius2, structure1, structure2, unit_z, world):
    # here we assume that structure1 = structure1D and structure2 = cap_structure
    # and that com has been correctly initialized with the initial (and constant) position
    # of the cap particle.

        pos2 = com          # cap particle (no. 2) stays in place
        pos1 = pos2 - iv    # IV always points from particle1 to particle2      

        return pos1, pos2, structure1.id, structure2.id
        # TODO Does this give the correct new positions in case of a single reaction of cap_particle ?

    def create_com_vector(self, z):        
    # This returns the displacement of the com, which is zero here.
        return 0.0

    def create_interparticle_vector(self, gf, r, dt, r0, old_iv): 
        if __debug__:
            log.debug("create_interparticle_vector: r=%g, dt=%g", r, dt)

        # note that r0 = length (old_iv)
        new_iv = (r/r0) * old_iv

        # project the new_iv down on the unit vectors of the surface to prevent the particle from
        # leaving the surface due to numerical problem
        unit_z = self.structure1D.shape.unit_z
        new_iv_z = unit_z * numpy.dot(new_iv, unit_z)

        return new_iv_z

    def __str__(self):
        return 'MixedPair1DCap' + Pair.__str__(self)


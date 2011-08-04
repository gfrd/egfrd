# Module Cyclops version 0.9.4.
# by Tim Peters (tim_one@email.msn.com).
# Released to the public domain 18-Jul-1999,

# Provided as-is; use at your own risk; no warranty; no promises; enjoy!

# Some parts of the code, and many ideas, were stolen with gratitude from:
#     Lars Marius Garshol's Plumbo.py.
#     Guido van Rossum.
#     Python's standard library module repr.py.

"""Module Cyclops -- stare at cycles in Python data structures.

Cyclops started life as a faster implementation of Lars Marius
Garshol's Plumbo.py, but almost instantly diverged due to feature
bloat.

The only object of real interest is class CycleFinder.  First get an
instance of that:

import Cyclops
z = Cyclops.CycleFinder()

This creates a cycle-tracker with an empty "root set".  Then tell z
which objects you're curious about (more on that later).

After running your regular code, tell z to look for cycles:

z.find_cycles(purge_dead_roots=0)
    Look for cycles reachable from the root set.
    Return 1 if a cycle was found, else 0.
    If optional arg purge_dead_roots is true, first remove objects
    from the root set whose true refcount (exclusive of CycleFinder
    internal references) is 0.
    See also method install_cycle_filter.

Then various methods can be used to get useful reports.  Which are
*most* useful you'll just have to play with.  If you have a few small
cycles, show_cycles is very revealing; if you have many cycles,
show_sccs and show_arcs *may* cut to the heart of the problem quickly.

z.show_stats()
    Print some stats about the last run of find_cycles:  number of
    distinct structured objects reachable from the root set, number of
    distinct objects involved in a cycle, number of cycles found,
    number of strongly connected components, number of objects in the
    root set, and number of edges (graph arcs) examined.

z.show_cycles()
    Print each cycle found to stdout.

z.show_cycleobjs(compare=typename_address_cmp)
    Print a listing of all objects involved in some cycle.  The
    objects are first sorted via the optional "compare" function.  By
    default, sorts first by type name and then by address (id); among
    objects of instance type, sorts by the instances' class names;
    among objects of class type, sorts by the classes' names.

x.show_sccs(compare=typename_address_cmp)
    Print a listing of cycle objects partitioned into strongly
    connected components (that is, the largest groupings possible such
    that each object in an SCC is reachable from every other object in
    that SCC).  Within each SCC, objects are sorted as for
    show_cycleobjs.

z.show_arcs(compare=None)
    Print a listing of all arc types involved in some cycle.  An arc
    is characterized by:
    1) its source object type or class,
    2) the selector (e.g. attribute) from the source, and
    3) its destination object type or class.
    The output is a table with four columns:  count, source, selector,
    destination.  The table is sorted by default on columns 2, 3, 4
    in that order; the first column (count) is ignored.  Specify an
    alternative compare function to change the sort; the compare
    function sees pairs of the form
        ((source_type, selector, destination_type), count)
    where each element of the triple is a string and count is an int.

z.show_obj(obj)
    Print a 2-line description of obj, including its address, adjusted
    refcount, type name and a short representation of its contents.
    See the method docstring for details.

There are two ways to add objects to the root set:

z.run(func, args=(), kwargs={})
    Run func (with optional args and keyword arguments kwargs), adding
    every object initialized by an __init__ function to the root set.
    This is the closest Cyclops gets to plumbo's "lazy" mode.  It's
    factored out here so it can be intermixed with the next method.

z.register(obj)
    Add obj to the root set.

To see the current root set,

z.get_rootset()
    Return the root set, as a list of (rc, cyclic?, obj) tuples.
    See method docstring for details.  In short, rc is the true
    refcount (CycleFinder internal references are subtracted out);
    cyclic is true iff the object is in a known cycle; and obj is the
    object.

Finally, to empty the root set again:

z.clear()
    Empty the root set, and release all other internal references to
    register'ed objects.

Notes:

+ Run Cyclops.py directly to exercise its _test() function; _test()
  sets up some common kinds of cycles, and should be easy to follow.
  Stare at the _test() code and the output until their relationship is
  clear.

+ find_cycles is linear-time, in the number of objects reachable from
  the root set plus the number of arcs connecting them.  This makes
  Cyclops pleasant for "real life" apps with tens of thousands of
  reachable objects, and at least usable well beyond that.

+ A (at least one) reference to each root-set object is maintained
  internally, so roots cannot die before invoking .clear() (or the
  CycleFinder is finalized).  This can distort the truth of your
  program, if a __del__ method of some root object that's not involved
  in a cycle could have caused cycles to get broken (this is unusual,
  but possible!).

  If you suspect it's happening, run find_cycles again passing true
  (optional arg purge_dead_roots):  find_cycles then releases all
  internal refs to root set objects whose true refcount is 0, thus
  allowing their __del__ methods to get invoked.  After that,
  find_cycles chases the remaining root set objects again.  You may
  need several iterations of invoking find_cycles(1) before reaching a
  steady state!

+ By default, cycles are chased through these (& only these) objects:

  - Lists.
  - Tuples.
  - Dicts (both keys and values).
  - Class instances (through their attributes).
  - Classes (through their attributes).
  - Instance method objects (through .im_self and .im_class).

  In particular, modules are not chased.  Maybe they should be.  See
  the next section for a way to force modules to get chased.


CHANGING THE SET OF INTERESTING TYPES

Methods of a CycleFinder object can be used to alter/query its view of
the set of types find_cycles "chases":

z.chase_type(t, t_refs_func, t_tag_func)
    Add type t to the set of interesting types.
    t_refs_func is a function of one argument x, where type(x) is t;
    it should return a sequence (typically a list or tuple) of all
    objects directly reachable from x.
    t_tag_func is a function of two arguments x and i, where type(x)
    is t and i in range(len(t_refs_func(x))).  It should return a
    brief string describing how t_refs_func(x][i] was obtained from
    x.
    See the _XXX_refs and _XXX_tag functions at the start of the
    module for examples of all this.

z.dont_chase_type(t)
    Remove type t from the set of interesting types.  find_cycles
    will not attempt to chase the objects reachable from objects of
    type t.

z.get_chased_types()
    Return the set of interesting types, as a list.


CHANGING THE SET OF INTERESTING CYCLES

Sometimes there are "expected" cycles you would like to ignore; e.g.,
this can happen if you install a module-chaser, because there are
cycles in Python's implementation you typically don't care about (for
example, if your module imports sys, there's a cycle because
sys.modules points back to your module!).

z.install_cycle_filter(filter_func=None)
    filter_func=None -> a way to ignore "expected" cycles.

    filter_func is a function of one argument, a cycle.  Each time
    find_cycles finds a cycle, the cycle is passed to filter_func.  The
    cycle is ignored unless filter_func returns true.  Passing None
    for filter_func restores the default behavior (do not ignore any
    cycle).

    A cycle is a list of (object, index) pairs, where the first object
    in the list is the same as the last object in the list, and where
    the object at cycle[i][0] is the cycle[i][1]'th object obtained
    from object cycle[i-1][0]. cycle[0][1] should be ignored (it tells
    us how we got to the first item in the cycle to begin with, but
    that's irrelevant to the cycle).


CASE STUDY

Below is the driver I used to track cycles in IDLE; it's a replacement
for IDLE's idle.py.

At first it didn't install function or module chasers, or a cycle filter,
and printed everything.  This turned up a bunch of easy cases, and the
show_sccs output was surprisingly revealing (all the cycles fell into a
handful of SCCs, which corresponded to distinct cycle-creating IDLE
subsystems).  show_arcs was also very helpful in getting the big picture.
show_cycles output was too voluminous to be helpful.

After those cycles were broken, the job got harder.  A module chaser was
added, which turned up another class of cycles, and then a function
chaser turned up 100s more.  Most of these involved expected cycles due
to Python's implementation, so a cycle filter was installed to ignore
cycles that didn't contain at least one class instance.  The remaining
cycles were isolated special cases, and only show_cycles output was
of real use.

After all cycles were purged, IDLE was still leaking, so driver output was
added to display the root-set objects still alive at the end.  This turned
up many cases where objects were living only because registration in the
root set was keeping them alive.  So a loop was added to the driver that
repeatedly purges dead root-set objects and tries again.  The __del__
methods of the purged roots caused other root objects to become trash,
and after several iterations of this the output reaches a steady state.

IDLE is still leaking (as of 18-Jul-1999), but ever less so, and while
Cyclops is no longer finding cycles, the driver's "live at the end"
output is still the best clue I've got for guessing what to do next.

Interesting:  At the start of this, it turned out that almost all cycles
were reachable from places outside themselves.  That is, they would not
have been considered trash even if Python used a mark-&-sweep form of
garbage collection.  IDLE's problems, in large part inherited from Tkinter,
are simply that "everything points to everything else".  The good news is
that Guido was able to clean most of this up just by adding reference-
purging code to widgets' explicitly-called destroy() methods.

#! /usr/bin/env python
import PyShell
import Cyclops
import types

def mod_refs(x):
    return x.__dict__.values()

def mod_tag(x, i):
    return "." + x.__dict__.keys()[i]

def func_refs(x):
    return x.func_globals, x.func_defaults

def func_tag(x, i):
    return (".func_globals", ".func_defaults")[i]

def instance_filter(cycle):
    for obj, index in cycle:
        if type(obj) is types.InstanceType:
            return 1
    return 0

# Note: PyShell binds ModifiedInterpreter.locals  to __main__.__dict__,
# and __main__ is *us*.  So if we don't keep the CycleFinder instance
# out of the global namespace here, z starts chewing over its own
# instance attributes.  Nothing breaks, but the output is at best
# surprising.

def hidez():
    z = Cyclops.CycleFinder()
    z.chase_type(types.ModuleType, mod_refs, mod_tag)
    z.chase_type(types.FunctionType, func_refs, func_tag)
    z.install_cycle_filter(instance_filter)
    z.run(PyShell.main)
    z.find_cycles()
    z.show_stats()
    z.show_cycles()
    # z.show_cycleobjs()
    # z.show_sccs()
    z.show_arcs()
    while 1:
        print "*" * 70
        print "non-cyclic root set objects:"
        sawitalready = {}
        numsurvivors = numdead = 0
        for rc, cyclic, x in z.get_rootset():
            if not sawitalready.has_key(id(x)):
                sawitalready[id(x)] = 1
                if rc == 0:
                    print "DEAD",
                    numdead = numdead + 1
                    z.show_obj(x)
                elif not cyclic:
                    numsurvivors = numsurvivors + 1
                    z.show_obj(x)
        x = None
        print numdead, "dead;", numsurvivors, "non-cycle & alive"
        if numdead == 0:
            break
        print "releasing dead root-set objects and trying again"
        z.find_cycles(1)
        z.show_stats()

hidez()
"""

# 0,9,4    18-Jul-1999
#    added purge_dead_roots arg to find_cycles
#    rearranged module docstring; added IDLE driver sample
# 0,9,3    29-Jun-1999
#    renamed print_obj to show_obj, and advertised it
#    redid adjusted refcount computations to account for root-set
#        objects too
#    changed get_rootset to return (refcount, cyclic?, obj) triples
# 0,9,2    27-Jun-1999
#    freed __init_tracer from dependence on name "self"
#    rearranged __find_cycles' inner loop for a nice little speed gain
#    which was promptly way more than lost by new code to compute &
#        display adjusted refcounts
#    which was partly regained by rewriting all that
# 0,9,1    26-Jun-1999
#    added SCC computation/display
#    added show_cycles; find_cycles no longer prints anything
#    changed all showXXX names to show_XXX
#    added install_cycle_filter
# 0,9,0    24-Jun-1999
#    first version I put under source control

__version__ = 0, 9, 4

#########################################################################
# Type-specific reference revealers.
#
# _T_refs(obj) should return a sequence of all objects directly
# reachable from obj.
#
# _T_tag(obj, i) should return a string briefly describing how
# _T_refs(obj][i] was obtained from obj.
#
# Why the separation?  Speed and space:  string tags are never
# computed unless a cycle is found and so something needs to be
# printed.

def _dict_refs(x):
    return x.keys() + x.values()
def _dict_tag(x, i):
    n = len(x)
    if i < n:
        return ".keys()[%d]" % i
    else:
        return "[%s]" % _quickrepr(x.keys()[i-n])

def _list_refs(x):
    return x
def _list_tag(x, i):
    return "[%d]" % i

_tuple_refs = _list_refs
_tuple_tag = _list_tag

def _instance_refs(x):
    # the keys are strings, so not interesting to return
    return x.__dict__.values()
def _instance_tag(x, i):
    return "." + x.__dict__.keys()[i]

_class_refs = _instance_refs
_class_tag = _instance_tag

def _instance_method_refs(x):
    return (x.im_self, x.im_class)
def _instance_method_tag(x, i):
    return (".im_self", ".im_class")[i]

import types
_default_refs_dispatcher = {
    types.DictType: _dict_refs,
    types.ListType: _list_refs,
    types.TupleType: _tuple_refs,
    types.InstanceType: _instance_refs,
    types.ClassType: _class_refs,
    types.MethodType: _instance_method_refs,
}
_default_tag_dispatcher = {
    types.DictType: _dict_tag,
    types.ListType: _list_tag,
    types.TupleType: _tuple_tag,
    types.InstanceType: _instance_tag,
    types.ClassType: _class_tag,
    types.MethodType: _instance_method_tag,
}
_InstanceType = types.InstanceType
_ClassType = types.ClassType
del types

del _dict_refs, _list_refs, _tuple_refs, _instance_refs, \
    _class_refs, _instance_method_refs
del _dict_tag, _list_tag, _tuple_tag, _instance_tag, \
    _class_tag, _instance_method_tag

#########################################################################
# A class to make short string representations of objects, for speed
# and brevity.  The std repr.repr sorts dicts by keys, but we refer to
# the keys via numeric subscript in cycle reports, so first a derived
# class that leaves dicts in raw order.  Also, instances, instance
# methods and classes frequently appear in cycle reports, so avoid
# chopping their reprs at all.  We're only trying to prevent massive
# expense for massive lists, tuples & dicts.

import repr
_repr = repr
del repr

class _CyclopsRepr(_repr.Repr):

    def __init__(self):
        _repr.Repr.__init__(self)
        # Boost the limit on types we don't know about; else it's too
        # easy to get a useless repr string when adding new types via
        # CycleFinder.chase_type.
        # Perhaps better to expose this class too, but-- sheesh --
        # this module is complicated enough!
        self.maxstring = self.maxother = 40

    # override base dictionary formatter; the code is almost the same,
    # we just leave the dict order alone

    def repr_dictionary(self, x, level):
        n = len(x)
        if n == 0:
            return '{}'
        if level <= 0:
            return '{...}'
        s = ''
        for k, v in x.items()[:min(n, self.maxdict)]:
            if s:
                s = s + ', '
            s = s + self.repr1(k, level-1)
            s = s + ': ' + self.repr1(v, level-1)
        if n > self.maxdict:
            s = s + ', ...'
        return '{' + s + '}'

    # don't chop instance, class or instance method reprs

    def repr_instance(self, x, level):
        try:
            return `x`
            # Bugs in x.__repr__() can cause arbitrary
            # exceptions -- then make up something
        except:
            return '<' + x.__class__.__name__ + ' instance at ' + \
                   hex(id(x))[2:] + '>'

    def repr_class(self, x, level):
        return `x`

    repr_instance_method = repr_class

_quickrepr = _CyclopsRepr().repr

#########################################################################
# CycleFinder is the workhorse.

def typename_address_cmp(x, y):
    if isinstance(x, _InstanceType) and isinstance(y, _InstanceType):
        xname, yname = x.__class__.__name__, y.__class__.__name__
    elif isinstance(x, _ClassType) and isinstance(y, _ClassType):
        xname, yname = x.__name__, y.__name__
    else:
        xname, yname = type(x).__name__, type(y).__name__
    return cmp((xname, id(x)), (yname, id(y)))

class CycleFinder:
    """Class for finding cycles in Python data structures.

    See Cyclops module docstring for details.
    """

    def __init__(self):
        """Create a cycle finder with empty root set."""

        self.clear()
        self.refs_dispatcher = _default_refs_dispatcher.copy()
        self.tag_dispatcher = _default_tag_dispatcher.copy()
        self.cycle_filter = None

    def clear(self):
        """Remove all internal references to external objects.

        Empties the root set.
        Does not change the set of types this CycleFinder chases.
        Does not change the cycle filter in effect.
        """

        self.roots = []
        self.__reset()

    def register(self, obj):
        """obj -> add object obj to the root set."""

        self.roots.append(obj)

    def run(self, func, args=(), kwargs={}):
        """func, args=(), kwargs={} -> add objects to root set by magic.

        Function func is invoked with arguments args and keyword
        arguments kwargs.  For the duration of the call, each class
        instance initialized by an __init__ call is automatically
        added to the root set.  The result of invoking func is
        returned. """

        # This clever method of trapping __init__ invocations was
        # stolen from Lars' plumbo.py.
        import sys
        sys.setprofile(self.__init_tracer)
        try:
            return apply(func, args, kwargs)
        finally:
            sys.setprofile(None)

    def find_cycles(self, purge_dead_roots=0):
        """purge_dead_roots=0 -> look for cycles, return true if found.

        Identify all cycles among objects reachable from the root set.
        Return true iff at least one cycle is found.

        This should be called before any of the show_XXX methods.
        Note that it's OK to add more objects to the root set and
        call it again, or to change the set of chased types, etc.
        find_cycles starts over from scratch each time it's called.

        If optional arg purge_dead_roots is true (default false),
        before searching for cycles the root set is purged of all
        objects that the preceding run of find_cycles determined had a
        true refcount of 0 (that is, the root set objects that are
        still alive only because they appear in the root set).
        Purging these allows their finalizers to get invoked, which
        may allow a cascade of other objects (including cycles) to go
        away too.

        See also method install_cycle_filter.
        """

        if purge_dead_roots and self.roots:
            id2rc = self.id2rc.get
            survivors = []
            for x in self.roots:
                if id2rc(id(x), None) != 0:
                    survivors.append(x)
            self.roots = survivors
            del x, survivors, id2rc

        self.__reset()
        self.__find_cycles(self.roots, 0)
        del self.seenids[id(self.roots)]    # not a user-visible object

        # Compute true refcounts for objects in cycles.
        id2rc = self.id2rc
        from sys import getrefcount
        for x in self.cycleobjs.values():
            # From the system refcount, subtract one for each of:
            #    being an element in the loop temp cycleobjs.values()
            #    being bound to "x"
            #    being an argument to getrefcount
            #    being a value in the cycleobjs dict
            #    showing up exactly once somewhere in self.sccno2objs
            xid = id(x)
            id2rc[xid] = getrefcount(x) - 5

        # Need also to subtract refs due to appearances in
        # self.cycles.  Complication:  some pairs in self.cycles may
        # be shared.
        seenpairs = {}
        isknown = seenpairs.has_key
        for cycle in self.cycles:
            for pair in cycle:
                pairid = id(pair)
                if not isknown(pairid):
                    seenpairs[pairid] = 1
                    xid = id(pair[0])
                    id2rc[xid] = id2rc[xid] - 1
        del isknown, seenpairs

        # And need also to subtract refs for membership is self.roots.
        # While we're at it, also compute adjusted refcounts for other
        # root set objects.
        for x in self.roots:
            xid = id(x)
            if id2rc.has_key(xid):
                id2rc[xid] = id2rc[xid] - 1
            else:
                assert not self.cycleobjs.has_key(xid)
                # Subtract one for each of:
                #     being bound to "x"
                #     being an argument to getrefcount
                #     being in self.roots
                id2rc[xid] = getrefcount(x) - 3

        return len(self.cycles) > 0

    def install_cycle_filter(self, filter_func=None):
        """filter_func=None -> a way to ignore "expected" cycles.

        See module docstring for details.  This is a callback function
        invoked whenever find_cycles() finds a cycle; the cycle is
        ignored unless the callback returns true.
        """

        self.cycle_filter = filter_func

    def show_stats(self):
        """Print statistics for the last run of find_cycles."""

        self._print_separator()
        print "# objects in root set:", len(self.roots)
        print "# distinct structured objects reachable:", len(self.seenids)
        print "# distinct structured objects in cycles:", len(self.cycleobjs)
        print "# cycles found:", len(self.cycles)
        print "# cycles filtered out:", self.ncyclesignored
        print "# strongly-connected components:", len(self.sccno2objs)
        print "# arcs examined:", self.narcs

    def show_cycles(self):
        """Print all cycles to stdout."""

        self._print_separator()
        print "# all cycles:"
        n = len(self.cycles)
        for i in xrange(n):
            self._print_cycle(self.cycles[i])
            if i < n-1:
                print "-" * 20

    def show_cycleobjs(self, compare=typename_address_cmp):
        """compare=typename_address_cmp -> print all objects in cycles.

        Prints to stdout.  Each distinct object find_cycles found in a
        cycle is displayed.  The set of objects found in cycles is
        first sorted by the optional "compare" function.  By default,
        objects are sorted using their type name as the primary key
        and their storage address (id) as the secondary key; among
        objects of instance type, sorts by the instances' class names;
        among objects of class type, sorts by the classes' names.
        """

        self._print_separator()
        print "# objects involved in cycles:"
        objs = self.cycleobjs.values()
        objs.sort(compare)
        for obj in objs:
            self.show_obj(obj)

    def show_sccs(self, compare=typename_address_cmp):
        """compare=typename_address_cmp -> print SCCs.

        Prints to stdout.  Shows the objects in cycles partitioned into
        strongly connected components (that is, the largest groupings
        possible such that each object in an SCC is reachable from every
        other object in that SCC).  Within each SCC, objects are sorted
        as for show_cycleobjs.
        """

        self._print_separator()
        print "# cycle objects partitioned into maximal SCCs:"
        sccs = self.sccno2objs.values()
        n = len(sccs)
        for i in xrange(n):
            print "--- SCC", i+1, "of", n
            objs = sccs[i]
            objs.sort(compare)
            for obj in objs:
                self.show_obj(obj)

    def show_arcs(self, compare=None):
        """compare=None -> print unique arc types in cycles.

        See module docstring for details.  Briefly, each arc in a
        cycle is categorized by the type of the source node, the kind
        of arc (how we got from the source to the destination), and
        the type of the destination node.  Each line of output
        consists of those three pieces of info preceded by the count
        of arcs of that kind.  By default, the rows are sorted first
        by column 2 (source node type), then by columns 3 and 4.
        """

        self._print_separator()
        print "# arc types involved in cycles:"
        items = self.arctypes.items()
        if compare:
            items.sort(compare)
        else:
            items.sort()
        for triple, count in items:
            print "%6d %-20s %-20s -> %-20s" % ((count,) + triple)

    def get_rootset(self):
        """Return the root set, as a list of (rc, cyclic?, obj) tuples.

        Should be called after find_cycles.  For each object in the
        root set, returns a triple consisting of
        refcount
            number of outstanding references less those due to
            CycleFinder internals; see show_obj docstring for more
            details; this will be None if find_cycles hasn't been
            run, or not since the last clear()
        cyclic?
            true (1) iff obj is known to be in a cycle
        obj
            the object
        """

        result = []
        getrc = self.id2rc.get
        incycle = self.cycleobjs.has_key
        for x in self.roots:
            xid = id(x)
            result.append((getrc(xid, None), incycle(xid), x))
        return result

    def chase_type(self, t, t_refs_func, t_tag_func):
        """t, t_refs_func, t_tag_func -> chase type t.

        See module docstring for details.
        """

        self.refs_dispatcher[t] = t_refs_func
        self.tag_dispatcher[t] = t_tag_func

    def dont_chase_type(self, t):
        """t -> remove type t from the set of chased types.

        See module docstring for details.
        """

        try:
            del self.refs_dispatcher[t], \
                self.tag_dispatcher[t]
        except KeyError:
            pass

    def get_chased_types(self):
        """Return the set of chased types, as a list."""

        return self.refs_dispatcher.keys()

    def __init_tracer(self, frame, event, args):
        try:
            a = event == "call"
            b = frame.f_code.co_name == "__init__"

            if a and b:
                # We want to capture the first argument -- this works whether
                # or not it's named "self", and in case the function is like
                #     def __init__(*args):
                # it's still OK:  we'll pick up the name 'args', and add the
                # tuple it's bound to to the root set; the tuple's first
                # element is "self", so will be found by the tuple-chaser.
                locals = frame.f_code.co_varnames
                if locals:
                    # first argname is first element of locals
                    self.register(frame.f_locals[locals[0]])
        except:
            print `event`, frame, frame.f_code, frame.f_code.co_name


    def __reset(self):
        # Clear out everything except:
        #       the root set
        #       the refs_dispatcher
        #       the tag_dispatcher
        #       the cycle filter

        # stack exactly mirrors __find_cycles' recursive calls; it's a
        # list of (object, index) pairs.
        self.stack = []

        # Map id of active object to its index in self.stack.  Since
        # it's a depth-first search, there's a cycle iff we hit an
        # object that's already on the stack.
        self.id2stacki = {}

        # Set of (addresses of) all interesting objects seen.  Since
        # we use a depth-first search, there's never a reason to
        # revisit a node.
        self.seenids = {}

        # List of all cycles found; each element is a stack slice (a
        # list of (object, index) pairs).
        self.cycles = []

        # Set of objects found in cycles (maps id(obj) -> obj).
        self.cycleobjs = {}

        # Classifies arcs found in cycles, mapping
        # (source_type, selector, destination_type) triples to a count
        # of how many times that triple appears in a cycle.
        self.arctypes = {}

        # Support for computing strongly-connected components (SCC).
        # We do this by merging cycles into equivalence classes.
        # Could be done faster by e.g. Tarjan's algorithm in
        # __find_cycles, but would rather keep that lean since cycles
        # are expected to be unusual.
        self.nextsccno = 1      # monotonically increasing SCC id
        self.id2sccno = {}      # map id(obj) to obj's sccno
        self.sccno2objs = {}    # map sccno back to list of objects

        # For objects in cycles and root set, map address to true
        # reference count.
        self.id2rc = {}

        # Number of arcs examined.
        self.narcs = 0

        # Number of cycles ignored (filtered out).
        self.ncyclesignored = 0

    def __find_cycles(self, obj, i, id=id, type=type, len=len):
        # This can be called an enormous number of times, so speed
        # tricks are appropriate.

        stack = self.stack

        # Set of ids of objects being, or formerly, chased.
        seenids = self.seenids
        already_seen = seenids.has_key

        # Maps active object id to index in stack.
        id2stacki = self.id2stacki
        currently_on_stack = id2stacki.has_key

        refs_dispatcher = self.refs_dispatcher
        is_interesting_type = refs_dispatcher.has_key

        myid = id(obj)
        seenids[myid] = 1
        id2stacki[myid] = len(stack)
        stack.append((obj, i))
        refs = refs_dispatcher[type(obj)](obj)
        self.narcs = self.narcs + len(refs)
        for i in xrange(len(refs)):
            child = refs[i]
            if is_interesting_type(type(child)):
                childid = id(child)
                if not already_seen(childid):
                    self.__find_cycles(child, i)
                elif currently_on_stack(childid):
                    cycle = stack[id2stacki[childid]:]
                    cycle.append((child, i)) # complete the cycle
                    self.__study_cycle(cycle)

        del stack[-1], id2stacki[myid]

    # a helper for __study_cycle
    def __obj2arcname(self, obj):
        if isinstance(obj, _InstanceType):
            name = obj.__class__.__name__ + "()"
        elif isinstance(obj, _ClassType):
            name = obj.__name__
        else:
            name = type(obj).__name__
        return name

    def __study_cycle(self, slice):
        assert len(slice) >= 2

        if self.cycle_filter is not None and \
           not self.cycle_filter(slice):
            self.ncyclesignored = self.ncyclesignored + 1
            return

        self.cycles.append(slice)

        # Pick (or create) an SCC equivalence class for this cycle.
        sccnowinner = self.id2sccno.get(id(slice[0][0]), None)
        if sccnowinner is None:
            sccnowinner = self.nextsccno
            self.nextsccno = self.nextsccno + 1
            self.sccno2objs[sccnowinner] = []
        classwinner = self.sccno2objs[sccnowinner]

        for i in xrange(len(slice)-1):
            obj1 = slice[i][0]
            key1 = self.__obj2arcname(obj1)

            obj2, index = slice[i+1]
            key2 = self.tag_dispatcher[type(obj1)](obj1, index)

            key3 = self.__obj2arcname(obj2)

            key = (key1, key2, key3)
            self.arctypes[key] = self.arctypes.get(key, 0) + 1

            self.cycleobjs[id(obj1)] = obj1

            # Merge the equivalence class for obj1 into classwinner.
            thissccno = self.id2sccno.get(id(obj1), None)
            if thissccno is None:
                self.id2sccno[id(obj1)] = sccnowinner
                classwinner.append(obj1)
            elif thissccno != sccnowinner:
                # merge obj1's entire equivalence class
                tomerge = self.sccno2objs[thissccno]
                for obj2 in tomerge:
                    self.id2sccno[id(obj2)] = sccnowinner
                classwinner.extend(tomerge)
                del self.sccno2objs[thissccno]

    ################################################################
    # Various output routines.  If someone is motivated enough to
    # change one of these, they're motivated enough to subclass us!

    def show_obj(self, obj):
        """obj -> print short description of obj to sdtout.

        This is of the form

        <address> rc:<refcount> <typename>
            repr: <shortrepr>

        where
        <address>
            hex address of obj
        <refcount>
            If find_cycles() has been run and obj is in the root set
            or was found in a cycle, this is the number of references
            outstanding less the number held internally by
            CycleFinder.  In most cases, this is what the true
            refcount would be had you not used CycleFinder at all.
            You can screw that up, e.g. by installing a cycle filter
            that holds on to references to one or more cycle elements.
            If find_cycles() has not been run, or has but obj wasn't
            found in a cycle and isn't in the root set, <refcount> is
            "?".
        <typename>
            type(obj), as a string.  If obj.__class__ exists, also
            prints the class name.
        <shortrepr>
            repr(obj), but using a variant of the std module repr.py
            that limits the number of characters displayed.
        """

        objid = id(obj)
        rc = self.id2rc.get(objid, "?")
        print hex(objid), "rc:" + str(rc), type(obj).__name__,
        if hasattr(obj, "__class__"):
            print obj.__class__,
        print
        print "    repr:", _quickrepr(obj)

    def _print_separator(self):
        print "*" * 70

    def _print_cycle(self, slice):
        n = len(slice)
        assert n >= 2
        print "%d-element cycle" % (n-1)
        for i in xrange(n):
            obj = slice[i][0]
            self.show_obj(obj)
            if i < n-1:
                index = slice[i+1][1]
                print "    this" + \
                      self.tag_dispatcher[type(obj)](obj, index), \
                      "->"

def _test():
    class X:
        def __init__(me, name):
            me.name = name
        def __repr__(self):
            return "X(" + `self.name` + ")"

    a, b, c, d = X('a'), X('b'), X('c'), X('d')
    a.k = b
    b.k = c
    c.k = a
    a.y = b.y = c.y = d
    d.k = {'harrumph': (1, 2, d, 3)}
    e = X('e')
    e.selfref = e
    a.gotoe = e
    a.__repr__ = a.__repr__

    class Y: pass
    X.gotoy = Y
    Y.gotox = X

    lonely = X('lonely')

    z = CycleFinder()
    z.register(a)
    z.register(lonely)
    del a, b, c, d, e, X, Y, lonely
    z.find_cycles()
    z.show_stats()
    z.show_cycles()
    z.show_cycleobjs()
    z.show_sccs()
    z.show_arcs()
    print "dead root set objects:"
    for rc, cyclic, x in z.get_rootset():
        if rc == 0:
            z.show_obj(x)
    z.find_cycles(1)
    z.show_stats()

if __name__ == "__main__":
    _test()

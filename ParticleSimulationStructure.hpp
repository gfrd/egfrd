#ifndef PARTICLE_SIMULATION_STRUCTURE_HPP
#define PARTICLE_SIMULATION_STRUCTURE_HPP

#include "Structure.hpp"

template<typename Ttraits_>
struct ImmutativeStructureVisitor;

template<typename Ttraits_>     // Note that Ttraits_ refers to a the World traits
struct MutativeStructureVisitor;

template<typename Ttraits_>
struct ParticleSimulationStructure: public Structure<Ttraits_>
{
    typedef Ttraits_ traits_type;
    typedef Structure<traits_type> base_type;

    typedef typename base_type::structure_name_type     structure_name_type;
    typedef typename base_type::structure_type_id_type  structure_type_id_type;
    typedef typename base_type::structure_id_type       structure_id_type;

    virtual ~ParticleSimulationStructure() {}

    virtual void accept(ImmutativeStructureVisitor<traits_type> const&) const = 0;

    virtual void accept(MutativeStructureVisitor<traits_type> const&) = 0;

    // Constructor
    ParticleSimulationStructure(structure_name_type const& name, structure_type_id_type const& sid, structure_id_type const& parent_struct_id):
        base_type(name, sid, parent_struct_id) {}
};

#endif /* PARTICLE_SIMULATION_STRUCTURE_HPP */

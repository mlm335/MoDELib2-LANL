/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Polycrystal_cpp_
#define model_Polycrystal_cpp_

#include <filesystem>
#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>
#include <utility>
#include <tuple>
#include <map>
#include <vector>
#include <deque>
#include <tuple>
#include <chrono>
#include <random>
#include <memory>
#include <Eigen/Core>
#include <SimplicialMesh.h>

#include <Grain.h>
#include <GrainBoundary.h>
#include <LatticeVector.h>
#include <TextFileParser.h>
#include <Polycrystal.h>

namespace model
{

    template <int dim>
    Polycrystal<dim>::Polycrystal(const std::string& polyFile,
                                  const SimplicialMeshType& mesh_in) :
    /* init */ PolycrystallineMaterialBase(std::filesystem::path(polyFile).parent_path().string()+"/"+TextFileParser(polyFile).readString("materialFile",false),
                                           TextFileParser(polyFile).readScalar<double>("absoluteTemperature",true))
    /* init */,mesh(mesh_in)
    /* init */,grains(getGrains(polyFile))
    /* init */,grainBoundaries(getGrainBoundaries())
    /* init */,Omega(getAtomicVolume())
    {
        std::cout<<greenBoldColor<<"Created Polycrystal"<<defaultColor<<std::endl;
        StaticID<Lattice<dim>>::force_count(0);
        StaticID<PlanarMeshFace<dim>>::force_count(0);
    }


    template <int dim>
    typename Polycrystal<dim>::GrainContainerType Polycrystal<dim>::getGrains(const std::string& polyFile) const
    {
        GrainContainerType temp;
        for(const auto& rIter : mesh.regions())
        {
            std::cout<<greenBoldColor<<"Creating Grain "<<rIter.second->regionID<<defaultColor<<std::endl;
            StaticID<Lattice<dim>>::set_count(rIter.second->regionID);
            const auto success(temp.emplace(rIter.second->regionID,std::shared_ptr<Grain<dim>>(new Grain<dim>(*rIter.second,*this,polyFile))));
            if(!success.second)
            {
                throw std::runtime_error("Cound not insert grain "+std::to_string(rIter.second->regionID));
            }
        }
        
        for(const auto& rgnBnd : mesh.regionBoundaries())
        {// loop over region boundaries
            for(const auto& face : rgnBnd.second.faces())
            {// loop over faces of each region boundary
                
                
                std::shared_ptr<GrainBoundary<dim>> gb(new GrainBoundary<dim>(rgnBnd.second,
                                                                              face.second,
                                                                              *grain(rgnBnd.first.first),
                                                                              *grain(rgnBnd.first.second)));
                
                temp.at(rgnBnd.first.first )->grainBoundaries().emplace(rgnBnd.first,gb);
                temp.at(rgnBnd.first.second)->grainBoundaries().emplace(rgnBnd.first,gb);
            }
        }        
        return temp;
    }

    template <int dim>
    double Polycrystal<dim>::getAtomicVolume() const
    {
        const double omg0(grains.begin()->second->latticeBasis.determinant());
        for(const auto& grain : grains)
        {
            const double omg(grain.second->latticeBasis.determinant());
            if(std::fabs(omg-omg0)>FLT_EPSILON)
            {
                throw std::runtime_error("grains have different atomic volume");
            }
        }
        
        return omg0/this->atomsPerUnitCell;
    }

    template <int dim>
    std::map<std::pair<size_t,size_t>,const GrainBoundary<dim>* const> Polycrystal<dim>::getGrainBoundaries() const
    {
        std::map<std::pair<size_t,size_t>,const GrainBoundaryType* const> temp;
        
        for(const auto& grain : grains)
        {
            for(const auto& gb : grain.second->grainBoundaries())
            {
                temp.emplace(gb.first,gb.second.get());
            }
        }
        
        return temp;
    }

    template <int dim>
    const std::shared_ptr<typename Polycrystal<dim>::GrainType>& Polycrystal<dim>::grain(const size_t& k) const
    {
        return grains.at(k);
    }

    template <int dim>
    const typename Polycrystal<dim>::GrainBoundaryType& Polycrystal<dim>::grainBoundary(const size_t& i,
                                                                                        const size_t& j) const
    {
        return (i<j)? *grainBoundaries.at(std::make_pair(i,j)) : *grainBoundaries.at(std::make_pair(j,i));
    }

    template <int dim>
    typename Polycrystal<dim>::LatticeVectorType Polycrystal<dim>::latticeVectorFromPosition(const VectorDim& p,
                                                                                             const Simplex<dim,dim>* const guess) const
    {
        const std::pair<bool,const Simplex<dim,dim>*> temp(mesh.searchWithGuess(p,guess));
        assert(temp.first && "Position not found in mesh");
        return grain(temp.second->region->regionID)->latticeVector(p);
    }

    template <int dim>
    typename Polycrystal<dim>::LatticeVectorType Polycrystal<dim>::latticeVectorFromPosition(const VectorDim& p) const
    {
        return latticeVectorFromPosition(p,&(mesh.simplices().begin()->second));
    }

    template <int dim>
    typename Polycrystal<dim>::ReciprocalLatticeVectorType Polycrystal<dim>::reciprocalLatticeVectorFromPosition(const VectorDim& p,
                                                                                                                 const Simplex<dim,dim>* const guess) const
    {
        const std::pair<bool,const Simplex<dim,dim>*> temp(mesh.searchWithGuess(p,guess));
        assert(temp.first && "Position not found in mesh");
        return grain(temp.second->region->regionID)->reciprocalLatticeVector(p);
    }

    template <int dim>
    typename Polycrystal<dim>::ReciprocalLatticeVectorType Polycrystal<dim>::reciprocalLatticeVectorFromPosition(const VectorDim& p) const
    {
        return reciprocalLatticeVectorFromPosition(p,&(mesh.simplices().begin()->second));
    }

    template <int dim>
    typename Polycrystal<dim>::VectorDim Polycrystal<dim>::randomPoint() const
    {
        std::mt19937 generator(std::chrono::system_clock::now().time_since_epoch().count());
        std::uniform_real_distribution<double> distribution(0.0,1.0);
        
        VectorDim P0(VectorDim::Zero());
        
        P0 << mesh.xMin(0)+distribution(generator)*(mesh.xMax(0)-mesh.xMin(0)),
        /* */ mesh.xMin(1)+distribution(generator)*(mesh.xMax(1)-mesh.xMin(1)),
        /* */ mesh.xMin(2)+distribution(generator)*(mesh.xMax(2)-mesh.xMin(2));
        
        return P0;
    }

    template <int dim>
    std::pair<LatticeVector<dim>,int> Polycrystal<dim>::randomLatticePointInMesh() const
    {
        const VectorDim P0=randomPoint();
        auto searchResult=mesh.search(P0);
        if(searchResult.first)
        {// point inside
            const LatticeVector<dim> L0 = grain(searchResult.second->region->regionID)->snapToLattice(P0);
            searchResult=mesh.searchRegionWithGuess(L0.cartesian(),searchResult.second);
            if(searchResult.first)
            {// point inside
                return std::make_pair(L0,searchResult.second->region->regionID);
            }
            else
            {
                return randomLatticePointInMesh();
            }
        }
        else
        {
            return randomLatticePointInMesh();
        }
        
    }

    template class Polycrystal<3>;
}
#endif


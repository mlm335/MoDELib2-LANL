/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_IrradiationPrismaticLoopGenerator_cpp_
#define model_IrradiationPrismaticLoopGenerator_cpp_

#include <numbers>
#include <chrono>
#include <random>
#include <cmath>
#include <list>
#include <assert.h>
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <limits>
#include <cmath>
#include <numbers> // std::numbers

//#include <Simplex.h>
#include <SimplicialMesh.h>
#include <Polycrystal.h>
#include <PolycrystallineMaterialBase.h>
#include <LatticeModule.h>
//#include <PlaneMeshIntersection.h>
#include <DislocationNodeIO.h>
#include <DislocationLoopIO.h>
#include <DislocationLoopLinkIO.h>
#include <DislocationLoopNodeIO.h>
#include <DDconfigIO.h>
#include <DDauxIO.h>

#include <DislocationLinkingNumber.h>
#include <TextFileParser.h>
#include <DislocationInjector.h>
#include <MeshBoundarySegment.h>
//#include <ConfinedDislocationObject.h>
#include <GlidePlaneModule.h>
#include <MeshModule.h>
#include <Plane.h>
#include <MicrostructureGenerator.h>
#include <PlanesIntersection.h>
#include <LatticePlane.h>
#include <IrradiationPrismaticLoopGenerator.h>

namespace model
{

IrradiationPrismaticLoopGenerator::IrradiationPrismaticLoopGenerator(const IrradiationPrismaticLoopDensitySpecification& spec, MicrostructureGenerator& mg)
{
    std::cout<<magentaBoldColor<<"Generating HCP prismatic loop density"<<defaultColor<<std::endl;
    
    // fractionLoops represents the fracton loops for HCP
    double targetVacancyBasalLoopDensity=spec.targetIHCPLoopDensity*spec.fractionBasalVacancy;
    double targetInterstitialPrismLoopDensity=spec.targetIHCPLoopDensity*spec.fractionPrismaticInterstitial;
    double targetVacancyPrismLoopDensity=spec.targetIHCPLoopDensity*spec.fractionPrismaticVacancy;
    double targetInterstitialBasalLoopDensity=spec.targetIHCPLoopDensity*(1-spec.fractionBasalVacancy-spec.fractionPrismaticInterstitial-spec.fractionPrismaticVacancy);
    
    //double targetBasalPrismLoopDensity=spec.targetIHCPLoopDensity*(1-spec.fractionBasalVacancy-spec.fractionPrismaticInterstitial);
    size_t nIHCPDefects=0, nInterstitialPrism=0, nVacancyPrism=0, nVacancyBasal=0, nInterstitialBasal=0;
    double defectsIHCPDensity=nIHCPDefects/  mg.ddBase.mesh.volume()/std::pow(mg.ddBase.poly.b_SI,3);
    double prismInterstitialDensity=nInterstitialPrism/  mg.ddBase.mesh.volume()/std::pow(mg.ddBase.poly.b_SI,3);
    double prismVacancyDensity=nVacancyPrism/  mg.ddBase.mesh.volume()/std::pow(mg.ddBase.poly.b_SI,3);
    double basalVacancyDensity=nVacancyBasal/  mg.ddBase.mesh.volume()/std::pow(mg.ddBase.poly.b_SI,3);
    double basalInterstitialDensity=nInterstitialBasal/  mg.ddBase.mesh.volume()/std::pow(mg.ddBase.poly.b_SI,3);

    
    std::cout<<"HCP loop generator with target Interstitial <a> loop density:"<<targetInterstitialPrismLoopDensity<<", and target Vacancy <a> loop density:"<<targetVacancyPrismLoopDensity<<", and target Vacancy <c> loop density:"<<targetVacancyBasalLoopDensity<<", and target Interstitial <c> loop density:"<<targetInterstitialBasalLoopDensity<<std::endl;
    
    
    while(defectsIHCPDensity<spec.targetIHCPLoopDensity)
    {
        
        std::normal_distribution<double> interstitialPrismSizeDistribution(spec.PrismaticInterstitialLoopSizeMean,spec.PrismaticInterstitialLoopSizeStd);
        std::normal_distribution<double> vacancyPrismSizeDistribution(spec.PrismaticVacancyLoopSizeMean,spec.PrismaticVacancyLoopSizeStd);
        std::normal_distribution<double> vacancyBasalSizeDistribution(spec.BasalVacancyLoopSizeMean,spec.BasalVacancyLoopSizeStd);
        std::normal_distribution<double> interstitialBasalSizeDistribution(spec.BasalInterstitialLoopSizeMean,spec.BasalInterstitialLoopSizeStd);

        std::mt19937 generator(std::chrono::system_clock::now().time_since_epoch().count());
        
        const std::pair<LatticeVector<3>,int> rp(mg.ddBase.poly.randomLatticePointInMesh());
        const int& grainID=rp.second;   // random grain ID
        const auto& grain(mg.ddBase.poly.grain(grainID));
        const LatticeVector<3>& L0=rp.first; // random lattice position in the grain
        const VectorDimD P0(L0.cartesian());   // cartesian position of L0

        std::uniform_int_distribution<> ssDist(0,mg.ddBase.poly.grain(grainID)->slipSystems().size()-1);
        const int pID(ssDist(generator)); // a random SlipSystem ID
        const auto& slipSystem(*grain->slipSystems().at(pID));

        const auto& planeBase(*grain->planeNormals().at(pID));
        const double planeSpacing(slipSystem.n.planeSpacing());
        
        
        if(fabs(planeSpacing-sqrt(3.0)/2.0)<FLT_EPSILON && prismInterstitialDensity<targetInterstitialPrismLoopDensity)
        {// Insert Prismatic Interstitial Loops
            
//            // Number of atoms in loop
//            const double size2radius=2.33e-29/(mg.ddBase.poly.b_SI*4.);   // Omega/pi/ba in m^-2
//            double xsize=interstitialPrismSizeDistribution(generator);
//            if(xsize<0) continue;
//            const double radius=sqrt(xsize*size2radius)/mg.ddBase.poly.b_SI; // Radius conversion from n atoms
            
            // Radius of loop
            double xsize=interstitialPrismSizeDistribution(generator);
            const double radius=xsize/mg.ddBase.poly.b_SI;
            
            VectorDimD b(slipSystem.s.cartesian());   // Prism axis
            const VectorDimD loopNorm(slipSystem.s.cartesian().normalized());   // Prism axis
            const VectorDimD R(slipSystem.unitNormal*radius);
            const ReciprocalLatticeDirection<3> r(grain->reciprocalLatticeDirection(loopNorm));
            const long int planeIndex(r.closestPlaneIndexOfPoint(P0));
            GlidePlaneKey<3> glidePlaneKey(planeIndex, r);
            std::shared_ptr<PeriodicGlidePlane<3>> glidePlane(mg.ddBase.periodicGlidePlaneFactory.get(glidePlaneKey));

            std::vector<VectorDimD> loopNodePos;
            for(int k=0;k<spec.IHCPLoopsNumberOfSides;++k)
            {
                loopNodePos.push_back(P0+Eigen::AngleAxisd(k*2.0*M_PI/spec.IHCPLoopsNumberOfSides, loopNorm)*R);
            }
            if (!mg.allPointsInGrain(loopNodePos, grainID)) {
                std::cout << "Some loop nodes are outside the grain. Skipping loop insertion." << std::endl;
                continue; 
            }
            
            const bool isVL(false);
            const VectorDimD Cycle_plane=(loopNodePos[1]-loopNodePos[0]).cross(loopNodePos[2]-loopNodePos[1]);
            if ((b.dot(Cycle_plane) > 0 && !isVL) || (b.dot(Cycle_plane) < 0 && isVL))
            {
                b *= -1.0;
            }

            mg.insertJunctionLoop(loopNodePos,glidePlane,b,loopNorm,P0,grainID,DislocationLoopIO<3>::SESSILELOOP);
            
            nInterstitialPrism++;
            prismInterstitialDensity=nInterstitialPrism/  mg.ddBase.mesh.volume()/std::pow(mg.ddBase.poly.b_SI,3);
            std::cout<<"Interstitial <a> loops density="<<prismInterstitialDensity<<std::endl;

        }
        else if(fabs(planeSpacing-sqrt(3.0)/2.0)<FLT_EPSILON && prismVacancyDensity<targetVacancyPrismLoopDensity)
        {// Insert Prismatic Vacancy Loops
            
//            // Number of Atoms
//            const double size2radius=2.33e-29/(mg.ddBase.poly.b_SI*4.);   // Omega/pi/ba in m^-2
//            double xsize=vacancyPrismSizeDistribution(generator);
//            if(xsize<0) continue;
//            const double radius=sqrt(xsize*size2radius)/mg.ddBase.poly.b_SI; // Radius conversion from n atoms
            
            // Radius of loop
            double xsize=vacancyPrismSizeDistribution(generator);
            const double radius=xsize/mg.ddBase.poly.b_SI; // Radius conversion from n atoms
            
            VectorDimD b(slipSystem.s.cartesian());   // Prism axis
            const VectorDimD loopNorm(slipSystem.s.cartesian().normalized());   // Prism axis
            const VectorDimD R(slipSystem.unitNormal*radius);
            const ReciprocalLatticeDirection<3> r(grain->reciprocalLatticeDirection(loopNorm));
            const long int planeIndex(r.closestPlaneIndexOfPoint(P0));
            GlidePlaneKey<3> glidePlaneKey(planeIndex, r);
            std::shared_ptr<PeriodicGlidePlane<3>> glidePlane(mg.ddBase.periodicGlidePlaneFactory.get(glidePlaneKey));

            std::vector<VectorDimD> loopNodePos;
            for(int k=0;k<spec.IHCPLoopsNumberOfSides;++k)
            {
                loopNodePos.push_back(P0+Eigen::AngleAxisd(k*2.0*M_PI/spec.IHCPLoopsNumberOfSides, loopNorm)*R);
            }
            if (!mg.allPointsInGrain(loopNodePos, grainID)) {
                std::cout << "Some loop nodes are outside the grain. Skipping loop insertion." << std::endl;
                continue; // Skip this loop
            }
            
            const bool isVL(true);
            const VectorDimD Cycle_plane=(loopNodePos[1]-loopNodePos[0]).cross(loopNodePos[2]-loopNodePos[1]);
            if ((b.dot(Cycle_plane) > 0 && !isVL) || (b.dot(Cycle_plane) < 0 && isVL))
            {
                b *= -1.0;
            }

            mg.insertJunctionLoop(loopNodePos,glidePlane,b,loopNorm,P0,grainID,DislocationLoopIO<3>::SESSILELOOP);
            
            nVacancyPrism++;
            prismVacancyDensity=nVacancyPrism/mg.ddBase.mesh.volume()/std::pow(mg.ddBase.poly.b_SI,3);
            std::cout<<"Vacancy <a> loops density="<<prismVacancyDensity<<std::endl;
    
        }
        else if(fabs(planeSpacing-sqrt(8.0/3.0))<FLT_EPSILON && basalVacancyDensity<targetVacancyBasalLoopDensity)
        {// Insert Vacancy Basal Loops
            
//            // Number of Atoms in Loop
//            const double size2radius=2.33e-29/(0.5*mg.ddBase.poly.b_SI*sqrt(8.0/3.0)*M_PI); // Omega/pi/bc in m^-2
//            double xsize=vacancyBasalSizeDistribution(generator);
//            if(xsize<0) continue;
//            const double radius=sqrt(xsize*size2radius)/mg.ddBase.poly.b_SI;
            
            // Radius of loop
            double xsize=vacancyBasalSizeDistribution(generator);
            const double radius=xsize/mg.ddBase.poly.b_SI; // Radius conversion from n atoms

            VectorDimD b(0.5*slipSystem.n.planeSpacing()*slipSystem.n.cartesian().normalized()); // 1/2 c-type loop
            const VectorDimD loopNorm(slipSystem.unitNormal);
            const VectorDimD R(slipSystem.s.cartesian().normalized()*radius);
            const ReciprocalLatticeDirection<3> r(grain->reciprocalLatticeDirection(loopNorm));
            const long int planeIndex(r.closestPlaneIndexOfPoint(P0));
            GlidePlaneKey<3> glidePlaneKey(planeIndex, r);
            std::shared_ptr<PeriodicGlidePlane<3>> glidePlane(mg.ddBase.periodicGlidePlaneFactory.get(glidePlaneKey));

            std::vector<VectorDimD> loopNodePos;
            for(int k=0;k<spec.IHCPLoopsNumberOfSides;++k)
            {
                loopNodePos.push_back(P0+Eigen::AngleAxisd(k*2.0*M_PI/spec.IHCPLoopsNumberOfSides, loopNorm)*R);
            }
            if (!mg.allPointsInGrain(loopNodePos, grainID)) {
                std::cout << "Some loop nodes are outside the grain. Skipping loop insertion." << std::endl;
                continue; // Skip this loop
            }

            const bool isVL(true);
            VectorDimD Cycle_plane=(loopNodePos[1]-loopNodePos[0]).cross(loopNodePos[2]-loopNodePos[1]);
            if ((b.dot(Cycle_plane) > 0 && !isVL) || (b.dot(Cycle_plane) < 0 && isVL))
            {
                b *= -1.0;
            }
            
            mg.insertJunctionLoop(loopNodePos,glidePlane,b,loopNorm,P0,grainID,DislocationLoopIO<3>::SESSILELOOP);
            
            nVacancyBasal++;
            basalVacancyDensity=nVacancyBasal/mg.ddBase.mesh.volume()/std::pow(mg.ddBase.poly.b_SI,3);
            std::cout<<"Vacancy <c> loops density="<<basalVacancyDensity<<std::endl;
            
        }
        else if(fabs(planeSpacing-sqrt(8.0/3.0))<FLT_EPSILON && basalInterstitialDensity<targetInterstitialBasalLoopDensity)
        {// Insert Interstitial Basal Loop
            
//            // Number of atoms in loop
//            const double size2radius=2.33e-29/(0.5*mg.ddBase.poly.b_SI*sqrt(8.0/3.0)*M_PI); // Omega/pi/bc in m^-2
//            double xsize=interstitialBasalSizeDistribution(generator);
//            if(xsize<0) continue;
//            const double radius=sqrt(xsize*size2radius)/mg.ddBase.poly.b_SI;
            
            // Radius of loop
            double xsize=interstitialBasalSizeDistribution(generator);
            const double radius=xsize/mg.ddBase.poly.b_SI; // Radius conversion from n atoms
            
            VectorDimD b(0.5*slipSystem.n.planeSpacing()*slipSystem.n.cartesian().normalized()); // 1/2 c-type loop
            const VectorDimD loopNorm(slipSystem.unitNormal);
            const VectorDimD R(slipSystem.s.cartesian().normalized()*radius);
            const ReciprocalLatticeDirection<3> r(grain->reciprocalLatticeDirection(loopNorm));
            const long int planeIndex(r.closestPlaneIndexOfPoint(P0));
            GlidePlaneKey<3> glidePlaneKey(planeIndex, r);
            std::shared_ptr<PeriodicGlidePlane<3>> glidePlane(mg.ddBase.periodicGlidePlaneFactory.get(glidePlaneKey));

            std::vector<VectorDimD> loopNodePos;
            for(int k=0;k<spec.IHCPLoopsNumberOfSides;++k)
            {
                loopNodePos.push_back(P0+Eigen::AngleAxisd(k*2.0*M_PI/spec.IHCPLoopsNumberOfSides, loopNorm)*R);
            }
            if (!mg.allPointsInGrain(loopNodePos, grainID)) {
                std::cout << "Some loop nodes are outside the grain. Skipping loop insertion." << std::endl;
                continue; // Skip this loop
            }
            
            const bool isVL(false);
            const VectorDimD Cycle_plane=(loopNodePos[1]-loopNodePos[0]).cross(loopNodePos[2]-loopNodePos[1]);
            if ((b.dot(Cycle_plane) > 0 && !isVL) || (b.dot(Cycle_plane) < 0 && isVL))
            {
                b *= -1.0;
            }
            
            mg.insertJunctionLoop(loopNodePos,glidePlane,b,loopNorm,P0,grainID,DislocationLoopIO<3>::SESSILELOOP);
            
            nInterstitialBasal++;
            basalInterstitialDensity=nInterstitialBasal/mg.ddBase.mesh.volume()/std::pow(mg.ddBase.poly.b_SI,3);
            std::cout<<"Interstitial <c> loops density="<<basalInterstitialDensity<<std::endl;
        }
        else
        {
            std::cout<<"Skipping a loop entry; not a prismatic plane or a basal plane"<<std::endl;
        }
        nIHCPDefects = nInterstitialPrism + nVacancyPrism + nVacancyBasal + nInterstitialBasal;
        defectsIHCPDensity=nIHCPDefects/  mg.ddBase.mesh.volume()/std::pow(mg.ddBase.poly.b_SI,3);
        std::cout<<"Total loop density="<<defectsIHCPDensity<<std::endl;
    }
    
    std::cout<<"Vacancy <c> loops number="<<nVacancyBasal<<std::endl;
    std::cout<<"Interstitial <c> loops number="<<nInterstitialBasal<<std::endl;
    std::cout<<"Vacancy <a> loops number="<<nVacancyPrism<<std::endl;
    std::cout<<"Interstitial <a> loops number="<<nInterstitialPrism<<std::endl;
}


IrradiationPrismaticLoopGenerator::IrradiationPrismaticLoopGenerator(const IrradiationPrismaticLoopIndividualSpecification& spec, MicrostructureGenerator& mg)
{
    std::cout<<magentaBoldColor<<"Generating individual HCP prismatic loops"<<defaultColor<<std::endl;
    if(spec.planeIDs.size())
    {
//            const std::vector<double> loopRadii(this->parser.readArray<double>("loopRadii_SI",true));
//            const Eigen::Matrix<double,Eigen::Dynamic,3> loopCenters(this->parser.readMatrix<double>("loopCenters",spec.planeIDs.size(),dim,true));
//            const std::vector<int> numberOfSides(this->parser.readArray<int>("numberOfSides",true));
//            const std::vector<int> isVacancyLoop(this->parser.readArray<int>("isVacancyLoop",true));

        if(spec.planeIDs.size()!=spec.loopRadii.size())
        {
            throw std::runtime_error("spec.planeIDs.size()="+std::to_string(spec.planeIDs.size())+" NOT EQUAL TO spec.spec.spec.loopRadii.size()="+std::to_string(spec.loopRadii.size()));
        }
        if(int(spec.planeIDs.size())!=spec.loopCenters.rows())
        {
            throw std::runtime_error("spec.planeIDs.size()="+std::to_string(spec.planeIDs.size())+" NOT EQUAL TO spec.loopCenters.rows()="+std::to_string(spec.loopCenters.rows()));
        }
        if(spec.planeIDs.size()!=spec.loopSides.size())
        {
            throw std::runtime_error("spec.planeIDs.size()="+std::to_string(spec.planeIDs.size())+" NOT EQUAL TO spec.loopSides.size()="+std::to_string(spec.loopSides.size()));
        }
        if(spec.planeIDs.size()!=spec.isVacancyLoop.size())
        {
            throw std::runtime_error("spec.planeIDs.size()="+std::to_string(spec.planeIDs.size())+" NOT EQUAL TO spec.isVacancyLoop.size()="+std::to_string(spec.isVacancyLoop.size()));
        }
        for(size_t k=0;k<spec.planeIDs.size();++k)
        {
            
            const bool isVL(spec.isVacancyLoop[k]>0);
//            generateSingle(mg,spec.planeIDs[k],spec.loopCenters.row(k),spec.loopRadii[k]/mg.ddBase.poly.b_SI,spec.loopSides[k],isVL);
            generateSingle(mg,spec.planeIDs[k],spec.loopCenters.row(k),spec.loopRadii[k],spec.loopSides[k],isVL);

        }
    }
    
}


void IrradiationPrismaticLoopGenerator::generateSingle(MicrostructureGenerator& mg,const int& pID,const VectorDimD& center,const double& radius,const size_t& sides,const bool& isVacancyLoop)
{
    std::pair<bool,const Simplex<3,3>*> found(mg.ddBase.mesh.search(center));
    if(!found.first)
    {
        std::cout<<"Point "<<center.transpose()<<" is outside mesh. EXITING."<<std::endl;
        exit(EXIT_FAILURE);
    }
    
    const int grainID(found.second->region->regionID);
    assert(mg.ddBase.poly.grains.size()==1 && "Periodic dislocations only supported for single crystals");
    const auto& grain(mg.ddBase.poly.grain(grainID));
    const VectorDimD P0(grain->snapToLattice(center).cartesian()); // WARNING: this may shift the point compared to the input.
    const auto& slipSystem(*grain->slipSystems().at(pID));
    const double planeSpacing(slipSystem.n.planeSpacing());
    const auto& planeBase(*grain->planeNormals().at(pID));

        if(fabs(planeSpacing-sqrt(3.0)/2.0)<FLT_EPSILON)
        { // prismatic plane spacing
            
            const double loopRadius(radius/mg.ddBase.poly.b_SI);            
            VectorDimD b(slipSystem.s.cartesian());   // Prism axis
            const VectorDimD loopNorm(slipSystem.s.cartesian().normalized());   // Prism axis
            const VectorDimD R(slipSystem.unitNormal*loopRadius);
            const ReciprocalLatticeDirection<3> r(grain->reciprocalLatticeDirection(loopNorm));
            const long int planeIndex(r.closestPlaneIndexOfPoint(P0));
            GlidePlaneKey<3> glidePlaneKey(planeIndex, r);
            std::shared_ptr<PeriodicGlidePlane<3>> glidePlane(mg.ddBase.periodicGlidePlaneFactory.get(glidePlaneKey));

            std::vector<VectorDimD> loopNodePos;
            for(size_t k=0;k<sides;++k)
            {
                loopNodePos.push_back(P0+Eigen::AngleAxisd(k*2.0*M_PI/sides, loopNorm)*R);
            }
            
            const VectorDimD Cycle_plane=(loopNodePos[1]-loopNodePos[0]).cross(loopNodePos[2]-loopNodePos[1]);
            if (b.dot(Cycle_plane)>0)
            {
                if (!isVacancyLoop) // 0 for interstitial
                {
                    b*=-1.0;
                }
            }
            else
            {
                if (isVacancyLoop) // 1 for vacancy
                {
                    b*=-1.0;
                }
            }
            std::cout << "Creating Individual Prismatic Loop" << std::endl;
            mg.insertJunctionLoop(loopNodePos,glidePlane,b,loopNorm,P0,grainID,DislocationLoopIO<3>::SESSILELOOP);
        }
        else if(fabs(planeSpacing-sqrt(8.0/3.0))<FLT_EPSILON)
        {// basal plane spacing
            
            VectorDimD b(0.5*slipSystem.n.planeSpacing()*slipSystem.n.cartesian().normalized()); // 1/2 c-type loop
            const VectorDimD loopNorm(slipSystem.unitNormal); 
            const ReciprocalLatticeDirection<3> r(grain->reciprocalLatticeDirection(loopNorm));
            const long int planeIndex(r.closestPlaneIndexOfPoint(center));
            GlidePlaneKey<3> glidePlaneKey(planeIndex, r);
            std::shared_ptr<PeriodicGlidePlane<3>> glidePlane(mg.ddBase.periodicGlidePlaneFactory.get(glidePlaneKey));
            const double loopRadius(radius/mg.ddBase.poly.b_SI);
            const VectorDimD R(slipSystem.s.cartesian().normalized()*loopRadius);
            
            std::vector<VectorDimD> loopNodePos;
            for(size_t k=0;k<sides;++k)
            {
                loopNodePos.push_back(P0+Eigen::AngleAxisd(k*2.0*M_PI/sides, loopNorm)*R);
            }
            
            VectorDimD Cycle_plane=(loopNodePos[1]-loopNodePos[0]).cross(loopNodePos[2]-loopNodePos[1]);
            if (b.dot(Cycle_plane)<0)
            {
                if (isVacancyLoop)
                {
                    b*=-1.0;
                }
            }
            else
            {
                if (!isVacancyLoop)
                {
                    b*=-1.0;
                }
            }
            mg.insertJunctionLoop(loopNodePos,glidePlane,b,loopNorm,P0,grainID,DislocationLoopIO<3>::SESSILELOOP);
            std::cout << "Created Individual Basal Loop" << std::endl;
        }

}


}
#endif

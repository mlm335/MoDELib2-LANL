/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_IrradiationPrismaticLoopDensitySpecification_cpp_
#define model_IrradiationPrismaticLoopDensitySpecification_cpp_

#include <IrradiationPrismaticLoopDensitySpecification.h>

namespace model
{

    IrradiationPrismaticLoopDensitySpecification::IrradiationPrismaticLoopDensitySpecification():
    /* init */ MicrostructureSpecificationBase("IrradiationPrismaticLoop","Density")
    /* init */,targetIHCPLoopDensity(0.0)
    /* init */,fractionBasalVacancy(0.0)
    /* init */,fractionPrismaticInterstitial(0.0)
    /* init */,fractionPrismaticVacancy(0.0)
    /* init */,PrismaticInterstitialLoopSizeMean(0.0)
    /* init */,PrismaticInterstitialLoopSizeStd(0.0)
    /* init */,PrismaticVacancyLoopSizeMean(0.0)
    /* init */,PrismaticVacancyLoopSizeStd(0.0)
    /* init */,BasalVacancyLoopSizeMean(0.0)
    /* init */,BasalVacancyLoopSizeStd(0.0)
    /* init */,BasalInterstitialLoopSizeMean(0.0)
    /* init */,BasalInterstitialLoopSizeStd(0.0)
    /* init */,IHCPLoopsNumberOfSides(0)
    {
        
    }

    IrradiationPrismaticLoopDensitySpecification::IrradiationPrismaticLoopDensitySpecification(const std::string& fileName):
    /* init */ MicrostructureSpecificationBase("IrradiationPrismaticLoop","Density",fileName)
    /* init */,targetIHCPLoopDensity(this->parser->readScalar<double>("targetIHCPLoopDensity",true))
    /* init */,fractionBasalVacancy(targetIHCPLoopDensity>0.0? this->parser->readScalar<double>("fractionBasalVacancy",true) : 0.0)
    /* init */,fractionPrismaticInterstitial(targetIHCPLoopDensity>0.0? this->parser->readScalar<double>("fractionPrismaticInterstitial",true) : 0.0)
    /* init */,fractionPrismaticVacancy(targetIHCPLoopDensity>0.0? this->parser->readScalar<double>("fractionPrismaticVacancy",true) : 0.0)
    /* init */,PrismaticInterstitialLoopSizeMean(targetIHCPLoopDensity>0.0? this->parser->readScalar<double>("PrismaticInterstitialLoopSizeMean",true) : 0.0)
    /* init */,PrismaticInterstitialLoopSizeStd(targetIHCPLoopDensity>0.0? this->parser->readScalar<double>("PrismaticInterstitialLoopSizeStd",true) : 0.0)
    /* init */,PrismaticVacancyLoopSizeMean(targetIHCPLoopDensity>0.0? this->parser->readScalar<double>("PrismaticVacancyLoopSizeMean",true) : 0.0)
    /* init */,PrismaticVacancyLoopSizeStd(targetIHCPLoopDensity>0.0? this->parser->readScalar<double>("PrismaticVacancyLoopSizeStd",true) : 0.0)
    /* init */,BasalVacancyLoopSizeMean(targetIHCPLoopDensity>0.0? this->parser->readScalar<double>("BasalVacancyLoopSizeMean",true) : 0.0)
    /* init */,BasalVacancyLoopSizeStd(targetIHCPLoopDensity>0.0? this->parser->readScalar<double>("BasalVacancyLoopSizeStd",true) : 0.0)
    /* init */,BasalInterstitialLoopSizeMean(targetIHCPLoopDensity>0.0? this->parser->readScalar<double>("BasalInterstitialLoopSizeMean",true) : 0.0)
    /* init */,BasalInterstitialLoopSizeStd(targetIHCPLoopDensity>0.0? this->parser->readScalar<double>("BasalInterstitialLoopSizeStd",true) : 0.0)
    /* init */,IHCPLoopsNumberOfSides(targetIHCPLoopDensity>0.0? this->parser->readScalar<int>("IHCPLoopsNumberOfSides",true) : 0)
    {
        
    }

}
#endif

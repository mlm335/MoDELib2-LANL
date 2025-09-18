/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_IrradiationPrismaticLoopIndividualSpecification_cpp_
#define model_IrradiationPrismaticLoopIndividualSpecification_cpp_

#include <IrradiationPrismaticLoopIndividualSpecification.h>

namespace model
{
    IrradiationPrismaticLoopIndividualSpecification::IrradiationPrismaticLoopIndividualSpecification():
    /* init */ MicrostructureSpecificationBase("IrradiationPrismaticLoop","Individual")
    {
        
    }

    IrradiationPrismaticLoopIndividualSpecification::IrradiationPrismaticLoopIndividualSpecification(const std::string& fileName):
    /* init */ MicrostructureSpecificationBase("IrradiationPrismaticLoop","Individual",fileName)
    {
        planeIDs=this->parser->readArray<int>("planeIDs",true);
        if(planeIDs.size())
        {
            loopRadii=this->parser->readArray<double>("loopRadii_SI",true);
            loopSides=this->parser->readArray<int>("loopSides",true);
            loopCenters=this->parser->readMatrix<double>("loopCenters",planeIDs.size(),3,true);
            isVacancyLoop=this->parser->readArray<int>("isVacancyLoop",true);
        }
    }
}
#endif

/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4 has been created by                                          */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4 is owned by                                 */
/*                 Charles Audet               - Polytechnique Montreal            */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD 4 has been funded by Rio Tinto, Hydro-Québec, Huawei-Canada,             */
/*  NSERC (Natural Sciences and Engineering Research Council of Canada),           */
/*  InnovÉÉ (Innovation en Énergie Électrique) and IVADO (The Institute            */
/*  for Data Valorization)                                                         */
/*                                                                                 */
/*  NOMAD v3 was created and developed by Charles Audet, Sebastien Le Digabel,     */
/*  Christophe Tribes and Viviane Rochon Montplaisir and was funded by AFOSR       */
/*  and Exxon Mobil.                                                               */
/*                                                                                 */
/*  NOMAD v1 and v2 were created and developed by Mark Abramson, Charles Audet,    */
/*  Gilles Couture, and John E. Dennis Jr., and were funded by AFOSR and           */
/*  Exxon Mobil.                                                                   */
/*                                                                                 */
/*  Contact information:                                                           */
/*    Polytechnique Montreal - GERAD                                               */
/*    C.P. 6079, Succ. Centre-ville, Montreal (Quebec) H3C 3A7 Canada              */
/*    e-mail: nomad@gerad.ca                                                       */
/*                                                                                 */
/*  This program is free software: you can redistribute it and/or modify it        */
/*  under the terms of the GNU Lesser General Public License as published by       */
/*  the Free Software Foundation, either version 3 of the License, or (at your     */
/*  option) any later version.                                                     */
/*                                                                                 */
/*  This program is distributed in the hope that it will be useful, but WITHOUT    */
/*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or          */
/*  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License    */
/*  for more details.                                                              */
/*                                                                                 */
/*  You should have received a copy of the GNU Lesser General Public License       */
/*  along with this program. If not, see <http://www.gnu.org/licenses/>.           */
/*                                                                                 */
/*  You can find information on the NOMAD software at www.gerad.ca/nomad           */
/*---------------------------------------------------------------------------------*/

#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Algos/Mads/MadsIteration.hpp"
#include "../../Algos/Mads/MadsMegaIteration.hpp"
#include "../../Algos/Mads/VNSSearchMethod.hpp"
#include "../../Algos/VNSMads/VNS.hpp"
#include "../../Cache/CacheBase.hpp"
#include "../../Output/OutputQueue.hpp"
//
// Reference: File VNS_Search.cpp in NOMAD 3.9.1
// Author: Christophe Tribes

void NOMAD::VNSSearchMethod::init()
{
    setStepType(NOMAD::StepType::SEARCH_METHOD_VNS_MADS);
    verifyParentNotNull();

    const auto parentSearch = getParentStep()->getParentOfType<NOMAD::VNSSearchMethod*>(false);
    
    // Do not perform if EVAL_SURROGATE_OPTIMIZATION is true
    if (nullptr != NOMAD::EvcInterface::getEvaluatorControl())
    {
        bool bBEval = ( NOMAD::EvcInterface::getEvaluatorControl()->getEvalType() == EvalType::BB ) ;
        
        // For some testing, it is possible that _runParams is null
        setEnabled((nullptr == parentSearch) && nullptr != _runParams && _runParams->getAttributeValue<bool>("VNS_MADS_SEARCH") && bBEval);
    }
    else
    {
        setEnabled(false);
    }
    if (isEnabled())
    {
        _trigger = _runParams->getAttributeValue<NOMAD::Double>("VNS_MADS_SEARCH_TRIGGER").todouble();
        
        // At first the reference frame center is not defined.
        // We obtain the frame center from the EvaluatorControl. If
        _refFrameCenter = NOMAD::Point();
        
        // Create the VNS algorithm with its own stop reason
        _vnsStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::VNSStopType>>();
        _vnsAlgo = std::make_unique<NOMAD::VNS>(this,
                                                _vnsStopReasons ,
                                                _runParams,
                                                _pbParams);
    }

}

bool NOMAD::VNSSearchMethod::runImp()
{
    bool foundBetter = false;
    
    if (isEnabled())
    {
        auto evalType = NOMAD::EvcInterface::getEvaluatorControl()->getEvalType();
        
        // check the VNS_trigger criterion:
        size_t bbEval = NOMAD::EvcInterface::getEvaluatorControl()->getBbEval();
        if (bbEval == 0 || double(_trialPointStats.getNbEvalsDone(evalType))/bbEval < _trigger)
        {
            
            EvalPointPtr frameCenter = nullptr;
               
            // Barrier of parent Mads not the same as the VNS Mads suboptimization
            std::shared_ptr<NOMAD::Barrier> barrier = nullptr;
            
            // Check that mesh from upper MadsMegaIteration is finer than initial
            auto madsMegaIter = getParentOfType<NOMAD::MadsMegaIteration*>(false);
            auto frameSize = madsMegaIter->getMesh()->getDeltaFrameSize();
            auto initialFrameSize = madsMegaIter->getMesh()->getInitialFrameSize();
            
            // Continue only if current frame size is small than initial frame size
            if ( initialFrameSize < frameSize )
            {
                OUTPUT_INFO_START
                AddOutputInfo("Current frame size larger than initial one. Stop VNS Mads Search.");
                OUTPUT_INFO_END
                setSuccessType(NOMAD::SuccessType::UNSUCCESSFUL);
                return foundBetter;
            }
            
            // Get barrier from upper MadsMegaIteration, if available.
            if (nullptr != madsMegaIter)
            {
                barrier = madsMegaIter->getBarrier();
            }
            else
            {
                throw NOMAD::Exception(__FILE__,__LINE__,"VNS Mads needs a barrier");
            }
            
            // MegaIteration's barrier member is already in sub dimension.
            auto bestXFeas = barrier->getFirstXFeas();
            auto bestXInf  = barrier->getFirstXInf();
            
            // Get the frame center for VNS sub optimization
            auto computeType = NOMAD::EvcInterface::getEvaluatorControl()->getComputeType();
            if (nullptr != bestXFeas
                && bestXFeas->getF(evalType, computeType).isDefined()
                && bestXFeas->getF(evalType, computeType) < MODEL_MAX_OUTPUT)
            {
                frameCenter = bestXFeas;
            }
            else if (nullptr != bestXInf
                     && bestXInf->getF(evalType, computeType).isDefined()
                     && bestXInf->getF(evalType, computeType) < MODEL_MAX_OUTPUT
                     && bestXInf->getH(evalType, computeType).isDefined()
                     && bestXInf->getH(evalType, computeType) < MODEL_MAX_OUTPUT)
            {
                frameCenter = bestXInf;
            }
            
            
            if ( nullptr != frameCenter )
            {
                
                _vnsAlgo->setEndDisplay(false);

                // VNS algo needs a frame center used as initial point for sub-optimization
                _vnsAlgo->setFrameCenter(frameCenter);

                // VNS conduct sub-optimization
                _vnsAlgo->start();
                _vnsAlgo->run();
                _vnsAlgo->end();
                
                // Get the success type and update Mads barrier with VNS Mads barrier
                auto vnsBarrier = _vnsAlgo->getBarrier();
                
                if (nullptr != vnsBarrier)
                {
                    auto vnsBestFeas = vnsBarrier->getFirstXFeas();
                    auto vnsBestInf = vnsBarrier->getFirstXInf();
                    NOMAD::SuccessType success = barrier->getSuccessTypeOfPoints(vnsBestFeas,
                                                                                 vnsBestInf,
                                                                                 NOMAD::EvalType::BB,
                                                                                 NOMAD::ComputeType::STANDARD);
                    setSuccessType(success);
                    if (success >= NOMAD::SuccessType::PARTIAL_SUCCESS)
                    {
                        foundBetter = true;
                    }
                    
                    // Update the barrier
                    barrier->updateWithPoints(vnsBarrier->getAllPoints(),
                                                                    NOMAD::EvalType::BB,
                                                                    NOMAD::ComputeType::STANDARD,
                                                                    _runParams->getAttributeValue<bool>("FRAME_CENTER_USE_CACHE"));
                    
                }
            }
        }
        else
        {
            OUTPUT_INFO_START
            AddOutputInfo("VNS trigger criterion not met. Stop VNS Mads Search.");
            OUTPUT_INFO_END
        }
    }
    return foundBetter;
}


void NOMAD::VNSSearchMethod::generateTrialPointsFinal()
{
    std::string s;
    NOMAD::EvalPointSet trialPoints;

    throw NOMAD::Exception(__FILE__,__LINE__,"VNS Mads generateTrialPointsFinal() not yet implemented.");
    
    // The trial points of one iteration of VNS are generated (not evaluated).
    // The trial points are obtained by shuffle + mads poll

    // auto madsIteration = getParentOfType<MadsIteration*>();

    /*
    // Note: Use first point of barrier as simplex center.
    NOMAD::VNSSingle singleVNS(this,
                            std::make_shared<NOMAD::EvalPoint>(getMegaIterationBarrier()->getFirstPoint()),
                            madsIteration->getMesh());
    singleVNS.start();
    singleVNS.end();

    // Pass the generated trial pts to this
    auto trialPts = singleVNS.getTrialPoints();
    for (auto point : trialPts)
    {
        insertTrialPoint(point);
    }
     */

}   // end generateTrialPoints



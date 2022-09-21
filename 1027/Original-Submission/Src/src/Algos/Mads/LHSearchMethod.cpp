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

#include "../../Algos/Mads/LHSearchMethod.hpp"
#include "../../Algos/SubproblemManager.hpp"
#include "../../Math/LHS.hpp"
#include "../../Type/LHSearchType.hpp"

void NOMAD::LHSearchMethod::init()
{
    setStepType(NOMAD::StepType::SEARCH_METHOD_LH);

    // For some testing, it is possible that _runParams is null
    if ( nullptr != _runParams)
    {
        auto lhSearch = _runParams->getAttributeValue<NOMAD::LHSearchType>("LH_SEARCH");
        setEnabled(lhSearch.isEnabled());
    }
    else
    {
        setEnabled(false);
    }
}


void NOMAD::LHSearchMethod::generateTrialPointsFinal()
{
    if (nullptr == _iterAncestor)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"LHSearchMethod: must have an iteration ancestor");
    }
    auto mesh = _iterAncestor->getMesh();
    if (nullptr == mesh)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"LHSearchMethod: must have a mesh");
    }

    // The frame center is only used to compute bounds, if they are not defined.
    // Use the first available point.
    auto barrier = getMegaIterationBarrier();
    if (nullptr == barrier)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"LHSearchMethod: must have a MadsMegaIteration ancestor with a barrier");
    }
    auto frameCenter = barrier->getFirstPoint();

    auto lhSearch = _runParams->getAttributeValue<NOMAD::LHSearchType>("LH_SEARCH");
    size_t n = _pbParams->getAttributeValue<size_t>("DIMENSION");
    size_t p = (0 == _iterAncestor->getK()) ? lhSearch.getNbInitial() : lhSearch.getNbIteration();
    auto lowerBound = _pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("LOWER_BOUND");
    auto upperBound = _pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("UPPER_BOUND");

    NOMAD::ArrayOfDouble deltaFrameSize = mesh->getDeltaFrameSize();
    NOMAD::Double scaleFactor = sqrt(-log(NOMAD::DEFAULT_EPSILON));
    // Apply Latin Hypercube algorithm (provide frameCenter, deltaFrameSize, and scaleFactor for updating bounds)
    NOMAD::LHS lhs(n, p, lowerBound, upperBound, frameCenter, deltaFrameSize, scaleFactor);
    auto pointVector = lhs.Sample();

    // Insert the point. Projection on mesh and snap to bounds is done in SearchMethod
    for (auto point : pointVector)
    {
        // Insert point (if possible)
        NOMAD::EvalPoint evalPoint(point);
        evalPoint.setPointFrom(std::make_shared<NOMAD::EvalPoint>(frameCenter), NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this));
        evalPoint.addGenStep(getStepType());
        insertTrialPoint(evalPoint);
    }
}

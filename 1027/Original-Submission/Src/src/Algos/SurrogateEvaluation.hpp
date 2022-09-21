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
#ifndef __NOMAD_4_2_SURROGATE_EVALUATION__
#define __NOMAD_4_2_SURROGATE_EVALUATION__

#include "../Algos/IterationUtils.hpp"
#include "../Algos/QuadModel/QuadModelIteration.hpp"
#include "../Algos/Step.hpp"

#include "../nomad_nsbegin.hpp"

/// Class to evaluate trial points using static SURROGATE or MODEL
class SurrogateEvaluation : public Step
{
private:
    EvalType _evalType;
    EvalPointPtr _frameCenter;
    EvaluatorPtr _evaluator;
    std::unique_ptr<QuadModelIteration> _quadModelIteration;
    
    EvalPointSet & _trialPoints;
    
public:
    /// Constructor
    explicit SurrogateEvaluation(const Step* parentStep,
                                 EvalPointSet & trialPoints,
                                 EvalPointPtr frameCenter = nullptr,
                                 EvalType evalType = EvalType::SURROGATE)
      : Step(parentStep),
       _trialPoints(trialPoints),
       _frameCenter(frameCenter),
       _evalType(evalType),
       _evaluator(nullptr),
       _quadModelIteration(nullptr)
    {
        init();
    }

private:
    void init();

    virtual void startImp() override;   ///< Construct Model if evalType==MODEL
    virtual bool runImp() override;     ///< Evaluate points using static surrogate
    virtual void endImp() override;     ///<  Do nothing

};

#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_2_SURROGATE_EVALUATION__

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

#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/Mads/Mads.hpp"
#include "../../Algos/QuadModel/QuadModelAlgo.hpp"
#include "../../Algos/QuadModel/QuadModelEvaluator.hpp"
#include "../../Algos/QuadModel/QuadModelIteration.hpp"
#include "../../Algos/QuadModel/QuadModelOptimize.hpp"
#include "../../Algos/SubproblemManager.hpp"
#include "../../Cache/CacheBase.hpp"
#include "../../Eval/ComputeSuccessType.hpp"
#include "../../Output/OutputQueue.hpp"
#include "../../Type/DirectionType.hpp"
#include "../../Type/EvalSortType.hpp"

void NOMAD::QuadModelOptimize::init()
{
    setStepType(NOMAD::StepType::OPTIMIZE);
    verifyParentNotNull();

    if (nullptr == _iterAncestor )
    {
        throw NOMAD::Exception(__FILE__,__LINE__,getName() + " must have an Iteration ancestor.");
    }

}


void NOMAD::QuadModelOptimize::startImp()
{
    auto modelDisplay = _runParams->getAttributeValue<std::string>("QUAD_MODEL_DISPLAY");
    _displayLevel = (std::string::npos != modelDisplay.find("O"))
        ? NOMAD::OutputLevel::LEVEL_INFO
        : NOMAD::OutputLevel::LEVEL_DEBUGDEBUG;

    OUTPUT_INFO_START
    std::string s;
    auto evcParams = NOMAD::EvcInterface::getEvaluatorControl()->getEvaluatorControlGlobalParams();
    s = "QUAD_MODEL_MAX_EVAL: " + std::to_string(evcParams->getAttributeValue<size_t>("QUAD_MODEL_MAX_EVAL"));
    AddOutputInfo(s, _displayLevel);
    s = "BBOT: " + NOMAD::BBOutputTypeListToString(NOMAD::QuadModelAlgo::getBBOutputType());
    AddOutputInfo(s, _displayLevel);
    OUTPUT_INFO_END

    generateTrialPoints();
}


bool NOMAD::QuadModelOptimize::runImp()
{
    std::string s;
    bool foundBetter = false;
    if ( ! _stopReasons->checkTerminate() )
    {
        foundBetter = evalTrialPoints(this);

        if (_modelFixedVar.nbDefined() > 0)
        {
            NOMAD::EvalPointSet evalPointSet;
            for (auto trialPoint : _trialPoints)
            {
                evalPointSet.insert(trialPoint.makeFullSpacePointFromFixed(_modelFixedVar));
            }
            _trialPoints.clear();
            _trialPoints = evalPointSet;
        }
        // Update barrier
        postProcessing();

        // If the oracle point cannot be evaluated the optimization has failed.
        if (_success==NOMAD::SuccessType::NOT_EVALUATED)
        {
            auto qmsStopReason = NOMAD::AlgoStopReasons<NOMAD::ModelStopType>::get ( getAllStopReasons() );
            qmsStopReason->set( NOMAD::ModelStopType::NO_NEW_POINTS_FOUND);
        }

    }


    return foundBetter;
}


void NOMAD::QuadModelOptimize::endImp()
{
    // Clean up the cache of points having only EvalType::SGTE
    NOMAD::CacheBase::getInstance()->deleteModelEvalOnly(NOMAD::getThreadNum());
}


void NOMAD::QuadModelOptimize::setupRunParameters()
{
    _optRunParams = std::make_shared<NOMAD::RunParameters>(*_runParams);

    _optRunParams->setAttributeValue("MEGA_SEARCH_POLL", false);
    
    _optRunParams->setAttributeValue("MAX_ITERATIONS", INF_SIZE_T);

    // Ensure there is no model, no NM and no VNS used in model optimization.
    _optRunParams->setAttributeValue("QUAD_MODEL_SEARCH", false);
    _optRunParams->setAttributeValue("SGTELIB_MODEL_SEARCH", false);
    _optRunParams->setAttributeValue("NM_SEARCH", false);
    _optRunParams->setAttributeValue("SPECULATIVE_SEARCH", true);
    
    // IMPORTANT: if VNS_MADS_SEARCH is changed to yes, the static members of VNSSearchMethod must be managed correctly
    _optRunParams->setAttributeValue("VNS_MADS_SEARCH", false);
    
    _optRunParams->setAttributeValue("ANISOTROPIC_MESH", false);

    // Set direction type to Ortho 2n
    _optRunParams->setAttributeValue("DIRECTION_TYPE",NOMAD::DirectionType::ORTHO_2N);

    // No hMax in the context of QuadModel
    _optRunParams->setAttributeValue("H_MAX_0", NOMAD::Double(NOMAD::INF));

    // Disable user calls
    _optRunParams->setAttributeValue("USER_CALLS_ENABLED", false);
    
    auto evcParams = NOMAD::EvcInterface::getEvaluatorControl()->getEvaluatorControlGlobalParams();
    _optRunParams->checkAndComply(evcParams, _optPbParams);
}


void NOMAD::QuadModelOptimize::setupPbParameters()
{
    _optPbParams = std::make_shared<NOMAD::PbParameters>(*_refPbParams);
    _optPbParams->setAttributeValue("LOWER_BOUND", _modelLowerBound);
    _optPbParams->setAttributeValue("UPPER_BOUND", _modelUpperBound);
    _optPbParams->setAttributeValue("FIXED_VARIABLE",_modelFixedVar);

    // Reset initial mesh and frame sizes
    // The initial mesh and frame sizes will be calculated from bounds and X0
    _optPbParams->resetToDefaultValue("INITIAL_MESH_SIZE");
    _optPbParams->resetToDefaultValue("INITIAL_FRAME_SIZE");
  
    // Use default min mesh and frame sizes
    _optPbParams->resetToDefaultValue("MIN_MESH_SIZE");
    _optPbParams->resetToDefaultValue("MIN_FRAME_SIZE");
    

    // Granularity is set to 0 and bb_input_type is set to all continuous variables. Candidate points are projected on the mesh before evaluation.
    _optPbParams->resetToDefaultValue("GRANULARITY");
    _optPbParams->resetToDefaultValue("BB_INPUT_TYPE");

    // No variable groups are considered for suboptimization
    _optPbParams->resetToDefaultValue("VARIABLE_GROUP");

    NOMAD::ArrayOfPoint x0s{_modelCenter};
    _optPbParams->setAttributeValue("X0", x0s);


    // We do not want certain warnings appearing in sub-optimization.
    _optPbParams->doNotShowWarnings();

    _optPbParams->checkAndComply();

}


void NOMAD::QuadModelOptimize::generateTrialPointsImp()
{

    // Set the bounds and fixed variables from the model
    setModelBoundsAndFixedVar();


    if ( _modelFixedVar.nbDefined() == _modelFixedVar.size() )
    {
        OUTPUT_INFO_START
        std::ostringstream oss;
        oss << "Effective dimension is null. No QuadModelOptimize" << std::endl;
        AddOutputInfo(oss.str(), NOMAD::OutputLevel::LEVEL_DEBUGDEBUG);
        OUTPUT_INFO_END

        return;
    }

    // Set specific evaluator control
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();

    // Enforce no opportunism and use no cache.
    auto previousOpportunism = evc->getOpportunisticEval();
    auto previousUseCache = evc->getUseCache();
    
    evc->setOpportunisticEval(true);
    evc->setUseCache(false);
    
    // Enforce specific sort type
    auto evcParams = evc->getEvaluatorControlGlobalParams();
    auto previousEvalSortType = evcParams->getAttributeValue<NOMAD::EvalSortType>("EVAL_QUEUE_SORT");
    evcParams->setAttributeValue("EVAL_QUEUE_SORT", NOMAD::EvalSortType::DIR_LAST_SUCCESS);
    evcParams->checkAndComply();
    

    auto modelDisplay = _runParams->getAttributeValue<std::string>("QUAD_MODEL_DISPLAY");

    OUTPUT_INFO_START
    std::string s = "Create QuadModelEvaluator with fixed variable = ";
    s += _modelFixedVar.display();
    AddOutputInfo(s);
    OUTPUT_INFO_END

    // Evaluations are in the quad model (local) full space: fixed variables detected when creating the training set (lb==ub) are not modified by the optimizer but evaluator receives points in local full space. The "local full space" is used because fixed variables from the parent problem are not seen/considered.
    // We send an empty Point for fixed variables to indicate that the evaluation are done in local full space.
    auto ev = std::make_shared<NOMAD::QuadModelEvaluator>(evc->getEvalParams(),
                                                          _model,
                                                          modelDisplay,
                                                          NOMAD::Point());

    // Replace the EvaluatorControl's evaluator with this one
    // we just created
    auto previousEvaluator = evc->setEvaluator(ev);
    if (nullptr == previousEvaluator)
    {
        std::cerr << "Warning: QuadModelOptimize: Could not set MODEL Evaluator" << std::endl;
        return;
    }

    auto madsStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::MadsStopType>>();


    // Setup Pb parameters just before optimization.
    setupPbParameters();

    // Set and verify run parameter values
    setupRunParameters();

    OUTPUT_INFO_START
    std::ostringstream oss;
    oss << "Run Parameters for QuadModelOptimize:" << std::endl;
    _optRunParams->display(oss, false);
    AddOutputInfo(oss.str(), NOMAD::OutputLevel::LEVEL_DEBUGDEBUG);
    OUTPUT_INFO_END

    // Run only if effective dimension not null.
    bool runOk = false;

    // Create a Mads step
    // Parameters for mads (_optRunParams and _optPbParams) are already updated.
    // We force mads to use local fixed variable to make sure we only work on fixed variables identified during construction of the training set. The evaluator works on the quad model, we don't need to map to the global full space for evaluation (not BB eval).
    auto mads = std::make_shared<NOMAD::Mads>(this, madsStopReasons, _optRunParams, _optPbParams, true /* use local fixed variables */);
    
    evc->resetModelEval();
    mads->start();
    runOk = mads->run();
    mads->end();

    evc->resetModelEval();
    evc->setEvaluator(previousEvaluator);

    // Note: No need to check the Mads stop reason: It is not a stop reason
    // for QuadModel.

    // Reset opportunism to previous values.
    evc->setOpportunisticEval(previousOpportunism);
    evc->setUseCache(previousUseCache);
    
    // Reset eval sort type to previous values
    evcParams->setAttributeValue("EVAL_QUEUE_SORT", previousEvalSortType);
    evcParams->checkAndComply();

    if (!runOk)
    {
        auto modelStopReasons = NOMAD::AlgoStopReasons<NOMAD::ModelStopType>::get(_stopReasons);
        modelStopReasons->set(NOMAD::ModelStopType::MODEL_OPTIMIZATION_FAIL);
    }
    else
    {
        // Get the best points in their reference dimension
        _bestXFeas = std::make_shared<NOMAD::EvalPoint>(mads->getBestSolution(true));
        _bestXInf  = std::make_shared<NOMAD::EvalPoint>(mads->getBestSolution(false));
        if (_bestXFeas->isComplete())
        {
            // New EvalPoint to be evaluated.
            // Add it to the list (local or in Search method).
            bool inserted = insertTrialPoint(*_bestXFeas);

            OUTPUT_INFO_START
            std::string s = "xt:";
            s += (inserted) ? " " : " not inserted: ";
            s += _bestXFeas->display();
            AddOutputInfo(s);
            OUTPUT_INFO_END
        }
        else
        {
            // Reset shared_ptr to default
            _bestXFeas.reset();
        }
        if (_bestXInf->isComplete())
        {
            // New EvalPoint to be evaluated.
            // Add it to the lists (local or in Search method).
            bool inserted = insertTrialPoint(*_bestXInf);
            
            OUTPUT_INFO_START
            std::string s = "xt:";
            s += (inserted) ? " " : " not inserted: ";
            s += _bestXInf->display();
            AddOutputInfo(s);
            OUTPUT_INFO_END
        
        }
        else
        {
            // Reset shared_ptr to default
            _bestXInf.reset();
        }
    }

}


// Set the bounds and the extra fixed variables (when lb=ub) of the model.
void NOMAD::QuadModelOptimize::setModelBoundsAndFixedVar()
{
    // When optWithScaleBounds is true, the training set is scaled with some generating directions. Warning: points are not necessarilly in [0,1]^n
    const SGTELIB::Matrix & X = _trainingSet->get_matrix_X();
    
    size_t n = _pbParams->getAttributeValue<size_t>("DIMENSION");
    if (n != (size_t)X.get_nb_cols())
    {
        throw NOMAD::Exception(__FILE__, __LINE__,
                               "QuadModel::setModelBounds() dimensions do not match");
    }
    
    int nbDim = X.get_nb_cols();
    int nbPoints = X.get_nb_rows();
    
    // Build model bounds & detect fixed variables
    NOMAD::Double lb;
    NOMAD::Double ub;
    bool isFixed = false;
    for (int j = 0; j < nbDim; j++)
    {
        lb = _modelLowerBound[j];
        ub = _modelUpperBound[j];
        for (int p = 0; p < nbPoints; p++)
        {
            // Use a comparison on regular double for better precision
            
            auto xpj = NOMAD::Double(X.get(p,j));
            if (lb.isDefined())
            {
                if (xpj.todouble() < lb.todouble())
                    lb = xpj;
            }
            else
            {
                lb = xpj;
            }
            if (ub.isDefined())
            {
                if (xpj.todouble() > ub.todouble())
                    ub = xpj;
            }
            else
            {
                ub = xpj;
            }
        }
        isFixed = false;
        // Comparison of Double at epsilon
        if (lb == ub)
        {
            _modelFixedVar[j] = ub;
            lb = NOMAD::Double(); // undefined
            ub = NOMAD::Double();
            isFixed = true;
        }
        
        
        if (!_optWithScaledBounds)
        {
            _modelLowerBound[j] = lb;
            _modelUpperBound[j] = ub;
        }
        else
        {
            // When scaling we force the bounds to be [0,1] whatever the training set except if we have a fixed variable.
            // If optWithScaledBounds is true and we detect a fixed variable, the modelFixedVar is not necessarily within [0,1] but the model center and the fixed variable must be equal.
            if (isFixed)
            {
                _modelLowerBound[j] = _modelUpperBound[j] = NOMAD::Double();
                _modelCenter[j] = _modelFixedVar[j];
            }
            // If not fixed, the model bounds and model center are pretermined.
            else
            {
                _modelLowerBound[j] = 0;
                _modelUpperBound[j] = 1;
                _modelCenter[j] = 0.5;
            }
        }
    }
    
    if (!_optWithScaledBounds)
    {
        // Detect the model center of the bounds
        // Scale the bounds around the model center
        auto reduction_factor = _runParams->getAttributeValue<NOMAD::Double>("QUAD_MODEL_SEARCH_BOUND_REDUCTION_FACTOR");
        for (int j = 0; j < nbDim; j++)
        {
            lb = _modelLowerBound[j];
            ub = _modelUpperBound[j];
            if (lb.isDefined() && ub.isDefined())
            {
                // The model center is the bounds middle point
                _modelCenter[j] = (lb + ub)/2.0;
                
                // Scale the bounds with respect to the bounds
                lb = _modelCenter[j] + (lb-_modelCenter[j])/reduction_factor;
                ub = _modelCenter[j] + (ub-_modelCenter[j])/reduction_factor;
                
                // Comparison of Double at epsilon
                if (lb == ub)
                {
                    _modelFixedVar[j] = ub;
                    _modelCenter[j] = ub;
                    lb = NOMAD::Double(); // undefined
                    ub = NOMAD::Double();
                    
                }
            }
            else
            {
                _modelCenter[j] = _modelFixedVar[j];
            }
            
            _modelLowerBound[j] = lb;
            _modelUpperBound[j] = ub;
        }
    }
    
    OUTPUT_INFO_START
    std::string s = "model lower bound: " + _modelLowerBound.display();
    AddOutputInfo(s);
    s = "model upper bound: " + _modelUpperBound.display();
    AddOutputInfo(s);
    s = "model center: " + _modelCenter.display();
    AddOutputInfo(s);
    OUTPUT_INFO_END
    
} // end setModelBounds

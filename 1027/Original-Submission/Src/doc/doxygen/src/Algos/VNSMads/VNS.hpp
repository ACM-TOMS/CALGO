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
#ifndef __NOMAD_4_2_VNS__
#define __NOMAD_4_2_VNS__


#include "../../Algos/Algorithm.hpp"
#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Algos/Mads/Mads.hpp"
#include "../../nomad_nsbegin.hpp"

/// Class implementing VNS Mads algorithm for constrained problems.

class VNS: public Algorithm
{
private:
    
    std::shared_ptr<AlgoStopReasons<MadsStopType>>    _madsStopReasons;
    
    std::shared_ptr<Barrier>    _barrier;
    
    std::shared_ptr<RunParameters>      _optRunParams; ///< run parameters for Mads sub optimization
    std::shared_ptr<PbParameters>       _optPbParams; ///< pb parameters for mads sub optimization
    
    EvalPointPtr          _frameCenter; ///< frame center to start mads sub optimization
    
    Point                           _refFrameCenter; ///<  The reference frame center to test if frame center is modified
    
    /**
     The neighborhood parameter is used to multiply the shake direction. Explore further away when neighborhood parameter is increased.
     */
    double _neighParameter;
    
public:
    /// Constructor
    /**
     \param parentStep          The parent of this Step -- \b IN.
     \param stopReasons         The stop reasons for NM -- \b IN.
     \param runParams           The run parameters that control NM -- \b IN.
     \param pbParams            The problem parameters that control NM -- \b IN.
     */
    explicit VNS(const Step* parentStep,
                std::shared_ptr<AlgoStopReasons<VNSStopType>> stopReasons,
                const std::shared_ptr<RunParameters>& runParams,
                const std::shared_ptr<PbParameters>& pbParams )
      : Algorithm(parentStep, stopReasons, runParams, pbParams),
        _barrier (nullptr),
        _frameCenter (nullptr),
        _neighParameter (0.0)
    {
        init();
    }

    /// Destructor
    virtual ~VNS() {}

    virtual void readInformationForHotRestart() override {}

    std::shared_ptr<Barrier> getBarrier() {return _barrier; }
    
    /// The frame center is used as initial point for the sub-obptimization
    void setFrameCenter(const EvalPointPtr frameCenter);
    
private:
    
    /// Helper for constructor
    void init();

    /// Implementation for run tasks.
    /**
     - Algorithm execution.
     - Shake current incumbents.
     - Perform mads.
     - Update the succes type
     - Perform Termination tasks (start, run, end)
     - Update the SearchMethod success type with best success found.
     \return \c true
     */
    virtual bool runImp() override;

    /// Implementation for start tasks.
    /**
     - Set the stop reason to STARTED
     - Reset sub-algorithm counter
     - Perform Initialization tasks (start, run, end)
     */
    virtual void startImp() override;

    virtual void endImp() override;
    
    // Helpers
    void setupRunParameters();
    void setupPbParameters(const NOMAD::Point & center, const NOMAD::ArrayOfDouble & currentMadsFrameSize);

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_2_VNS__

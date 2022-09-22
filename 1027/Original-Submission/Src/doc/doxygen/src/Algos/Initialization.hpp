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
#ifndef __NOMAD_4_2_INITIALIZATION__
#define __NOMAD_4_2_INITIALIZATION__

#include "../Algos/Step.hpp"
#include "../Algos/TrialPointStats.hpp"

#include "../nomad_nsbegin.hpp"

/// Class for initialization (step 0) of an Algorithm
/**
 This an abstract class, each algorithm should probably implement an initialization.
 */
class Initialization: public Step
{
protected:
    TrialPointStats                        _trialPointStats;   ///< The trial point counters stats for initialization
    
    std::shared_ptr<Barrier> _barrier;   ///< Barrier constructed from evaluated X0s

    
public:
    /// Constructor
    /*
     \param parentStep      The parent of this step -- \b IN.
     */
    explicit Initialization(const Step* parentStep)
      : Step(parentStep),
        _trialPointStats(parentStep),
        _barrier(nullptr)
    {
        init();
    }

    /// Destructor
    /**
     Upon destruction, print all that is in the output queue.
     */
    virtual ~Initialization();

    std::string getName() const override;

    const std::shared_ptr<Barrier>& getBarrier() const { return _barrier; }

protected:
    /// Helper for constructor
    void init();

public:
    virtual void startImp()    override {}
    virtual bool runImp()      override = 0;
    virtual void endImp()      override;

private:
    void incrementCounters () override;
    
};

#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_2_INITIALIZATION__

/*
////////////////////////////////////////////////////////////////////////////////////////
//
//  This file is part of Sandpile Simulator, a cellular automaton for sandpile dynamics.
//  Copyright (C) 2021 M. Frohne
//
//  Sandpile Simulator is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Affero General Public License as published
//  by the Free Software Foundation, either version 3 of the License,
//  or (at your option) any later version.
//
//  Sandpile Simulator is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//  See the GNU Affero General Public License for more details.
//
//  You should have received a copy of the GNU Affero General Public License
//  along with Sandpile Simulator. If not, see <https://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////////////
*/

#ifndef SIMULATIONMODEL_FWM_H
#define SIMULATIONMODEL_FWM_H

#include "aux.h"
#include "simulationmodel.h"
#include "avalanchestatistics.h"

/*! \brief Implements a customized sandpile model of ours, the "FW model" (FWM).
 *
 * In this model the integers stored at the lattice sites represent real stacks of sand grains,
 * i.e. each value stands for the height of a sand stack at the position of the lattice site.
 * As an example, a real-world sandpile exists in a 3D world and hence (when mapped on a
 * discrete lattice) consists of stacks of sand grains distributed on a 2-dim. lattice,
 * which can thus be simulated on a 2-dim. sandbox (see also Sandbox).
 * So, generally speaking the values at each lattice site of an n-dim. sandbox
 * quantify the (n+1)-dim. sandpile's extension in the (n+1)-th dimension.
 *
 * The dynamics of the simulated sandpile of course result from the driving and relaxation procedures,
 * which for this model are defined as follows (let s(r) be the grain stack height at lattice site r):
 *
 * - \b Driving:
 *   The driving procedure simply consists of choosing a (uniformly distributed) random
 *   lattice site r and increasing s(r) about 1. This is also called \e conservative \e driving.
 *   NB: \e non-conservative \e driving is typically used in the BTW model, though.
 *
 *   (see also drive())
 *
 * - \b Relaxation:
 *   The relaxation process is governed by the local slopes of the sandpile.
 *   Since we do not store slopes at the lattice sites but pile heights,
 *   these slopes first need to be calculated after the driving step
 *   before performing the actual relaxation. This however allows for
 *   directional slopes which isn't a feature of for instance the BTW model
 *   and thus may represent an essential improvement over the BTW model.
 *
 *   Because we operate on an equidistant lattice the slope between two neighboring
 *   sites r_1 and r_2 is defined as the height difference s(r_1) − s(r_2).
 *   If the absolute value of this slope is larger or equal the critical value q_crit
 *   the site with the higher stack is called \e critical \e site. The value q_crit
 *   is called critical slope and is a free parameter in this custom model (see also setModelParameter()).
 *
 *   If now one site r_0 is critical with respect to at least one of its nearest neighbors
 *   the stack at r_0 is considered unstable and about to collapse. All excess sand grains
 *   E = s(r_0) − s(r_<) with respect to that neighbor that has the least overall height, r_<,
 *   will then be redistributed to \e all nearest neighbors to which the slope is \e currently
 *   at least critical. One grain at a time is dropped to one of the critical neighbors
 *   and then the next critical neighbor is selected. This happens until none of the slopes
 *   to the selected neighbors is larger 0 anymore (neighbors reaching negative or 0 slope
 *   are gradually excluded) or E grains have been dropped. The initial site for that procedure
 *   is always selected randomly in order to hinder build-up of unintended directionality.
 *
 *   This whole procedure represents one unit time step. The first point r_0 to which
 *   the procedure is applied is that point whose stack has been modified in the
 *   previous driving step (as there the local slopes are initially changed).
 *
 *   When grains have been redistributed and thus slopes have changed, another relaxation step
 *   is performed immediately (i.e. before returning from relax()). For all nearest neighbors
 *   of the critical site from the previous step, {r_j}, it is checked if they have become
 *   critical (which also includes neighbors to which no grain has been dropped because the
 *   stack of the critical site has been lowered). The procedure for critical sites as described
 *   in the beginning is then applied to all critical sites of the {r_j} simultaneously.
 *
 *   Here simultaneously also means that the redistribution from one critical neighbor r_j^1
 *   doesn't affect the way of redistribution from another critical neighbor r_j^2.
 *   This is achieved by \e first calculating all the slopes for \e all the critical neighbors
 *   \e without updating these slopes upon redistribution of grains, such that \e then the
 *   simultaneous redistribution for all critical neighbors can be performed \e independently
 *   in an iterative manner. The grain redistributions of these independent redistribution processes
 *   are directly applied to the sandbox as the slopes have been fixed as described above.
 *   Another unit time step has passed.
 *
 *   All this is repeated until no site is critical anymore and thus the avalanche stops.
 *   The accumulated number of unit time steps sums up to the \e avalanche \e duration
 *   (see AvalancheStatistics::Event::duration).
 *
 *   (see also relax())
 *
 * Exceptions from the relaxation procedure described above:
 *
 * Boundary conditions: The FW model respects the boundary condition settings of the Sandbox
 * (see Sandbox::setBoundaryConditions()), which can be set open or closed for each edge of
 * the lattice independently. In this model the closed boundaries basically model an infinitely
 * high wall that does not let any sand grains pass whereas the open boundaries can be thought
 * of as the edge of a table where any excess sand can just drop from the table.
 *
 * This means that the whole stack of sand grains of a lattice site at an open edge is set to 0
 * when the height of that stack is larger or equal the critical slope q_crit (no dropping to the
 * nearest neighbors is considered here). This also counts as one unit time step.
 * At closed boundaries the relaxation follows the normal procedure except for the slope in direction
 * of that (non-existing) nearest neighbor that is not on the lattice being set to 0.
 *
 * In order to characterize the avalanches we still need to precisely define the observables in
 * AvalancheStatistics::Event. As explained above, the \e duration of an avalanche is defined
 * as the number of repetitions of the relaxation procedure (within a single call to relax())
 * needed to obtain a completely stable sandpile again. The \e linear size (see AvalancheStatistics::Event::linSize)
 * is the maximum distance between points of the set of avalanche sites (see calculateAvalancheLinSize()).
 * All lattice sites where sand grains have dropped from \e or onto during an avalanche at least once
 * we call "avalanche sites". The \e area (see AvalancheStatistics::Event::area) is the number of distinct avalanche sites
 * (see calculateAvalancheArea()). The \e size (see AvalancheStatistics::Event::size) in case of this model is the
 * accumulated number of individual collapsing grains stacks. So for each lattice site that becomes a critical site
 * the avalanche size is increased by 1. If the same site is critical for multiple iterations of the relaxation procedure
 * of the same avalanche, the avalanche size is increased multiple times accordingly.
 *
 * Note: This definition of the avalanche size may be considered as some sort of measure of the total mass movement
 * during an avalanche, however we actually don't count the \e number of moved grains but only \e how \e many \e times
 * these collapse events happen. So this interpretation is probably a little off. To actually make the avalanche size more
 * realistically represent mass movement, set the model parameter "FWM_relax_size_movedGrains" to 1 (see setModelParameter()).
 * This will make relax() set the size to the total/accumulated number of grain movements during the avalanche instead.
 */
class SimulationModel_FWM : public SimulationModel
{
public:
    SimulationModel_FWM(Sandbox& pSandbox, short pCritSlope);       ///< Constructor.
    virtual ~SimulationModel_FWM() = default;                       ///< Default destructor.
    //
    virtual Aux::Model id() const { return Aux::Model::_FWM; }      ///< \brief Get an ID identifying the SimulationModel
                                                                    ///         as defined in Aux::Model.
                                                                    ///  \returns Aux::Model::_FWM.
    //
    virtual void seed(std::seed_seq& pSeedSeq);     ///< Seed the random generator for initial topple direction randomization.
    //
    virtual bool setModelParameter(const std::string& pKey, int pVal);                      ///< Set additional model parameters.
    virtual bool getModelParameter(const std::string& pKey, int& pVal) const;               ///< Get additional model parameters.
    virtual std::vector<std::pair<std::string, std::string>> listModelParameters() const;   ///< List available model parameters.
    //
    virtual void drive() const;                                 ///< Drive the sandpile.
    virtual void relax(AvalancheStatistics::Event& pEvent);     ///< Relax the sandpile.
    //
    virtual double calculateSandboxCriticality() const;         ///< Quantify the criticality of the sandbox.

private:
    virtual double calculateAvalancheLinSize(const Avalanche& pAvalanche) const;    ///< Calculate the linear size of an avalanche.
    virtual long long calculateAvalancheArea(const Avalanche& pAvalanche) const;    ///< Calculate the area of an avalanche.

private:
    /*!
     * \brief Groups the position of a lattice site and if it is critical/unstable.
     */
    struct LatticeSite {
        std::vector<short> position;    //Cartesian position of the lattice site
        bool critical = false;          //Site is unstable?
    };

private:
    short critSlope;                            //Critical slope, i.e. local slope at which a lattice site becomes unstable
    bool avalancheSizeMovedGrains;              //Avalanche size is differently set to total number of individual grain movements
    //
    const std::vector<short>& currentGrain;     //The random lattice site that was used in the most recent driving step
    //
    std::mt19937_64 rndGenerator;                       //Mersenne Twister pRNG engine
    std::uniform_real_distribution<double> rndUnif;     //Generates unif. dist. numbers for randomly selecting initial drop neighbor
};

#endif // SIMULATIONMODEL_FWM_H

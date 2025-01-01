/*
////////////////////////////////////////////////////////////////////////////////////////
//
//  This file is part of Sandpile Simulator, a cellular automaton for sandpile dynamics.
//  Copyright (C) 2021, 2025 M. Frohne
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

#ifndef SANDSIM_SIMULATIONMODEL_BTW_H
#define SANDSIM_SIMULATIONMODEL_BTW_H

#include "aux.h"
#include "avalanche.h"
#include "avalanchestatistics.h"
#include "sandbox.h"
#include "simulationmodel.h"

#include <string>
#include <utility>
#include <vector>

/*! \brief Implements the Bak-Tang-Wiesenfeld (BTW) sandpile model.
 *
 * In this model the integers stored at the lattice sites represent sandpile slopes,
 * i.e. each value stands for the local height \e gradient at the position of the
 * lattice site. Since the gradients are only integer numbers they cannot contain
 * information about the direction of the slope. Compared with the picture of a real
 * sandpile (a grain stack at every lattice site) the critical state of the BTW model
 * must thus be principally thought of as a uniform “diagonal hillside” spanned
 * between one corner of the sandbox lattice and the opposite corner.
 *
 * The dynamics of the simulated sandpile of course result from the driving and relaxation procedures,
 * which for this model are defined as follows (let s(r) be the slope at lattice site r):
 *
 * - \b Driving:
 *   There are two modes for the driving procedure. The first one is called \e non-conservative \e driving and
 *   simply consists of choosing a (uniformly distributed) random lattice site r and increasing s(r) about 1.
 *   In the picture of a real sandpile this corresponds to an addition of half a sand grain at site r and
 *   also at \b all of its uphill neighbors, which for instance in a 2-dim. sandbox would span a triangle.
 *   The other, more realistic mode is called \e conservative \e driving, which corresponds to an addition
 *   of exactly one grain at site r only. Here s(r) at a random site r is increased about dim[sandbox] and
 *   the site's nearest uphill neighbor slopes s(r - e_i) are \e decreased about 1 each.
 *
 *   (see also drive())
 *
 * - \b Relaxation:
 *   The relaxation process is driven by the occurrence of so-called \e critical \e sites.
 *   These are sites that have a slope value larger or equal the \e critical \e slope q_crit,
 *   which in this model is fixed to 2 * dim[sandbox] (this equals the number of nearest
 *   neighbors of a site!). Sites of course mainly become critical because of a previous
 *   driving step but can also temporarily become critical during the relaxation process.
 *
 *   The criticality of all sites is checked in the beginning of the relaxation process.
 *   If no site is critical the sandpile is stable and the relaxation procedure ends immediately
 *   (note: with a duration of 0). If there are critical sites (maximally one \e directly after driving),
 *   the sandpile is unstable and grains will be redistributed. A grain redistribution corresponds to
 *   a slope redistribution in this model. This slope redistribution is constant for all critical sites
 *   (does not depend on the sandbox state) and so all required relaxations simply take place serially.
 *   For each critical site r_0 the slope s(r_0) is reduced by 2 * dim[sandbox] and all neighboring slopes
 *   of the critical site, s(r_0 +- e_i), are increased by 1. This whole procedure for all critical sites
 *   represents one unit time step. Note that the redistribution conserves slope units (except at boundaries)
 *   but not necessarily sand grains when trying to translate to the real sandpile picture!
 *
 *   When slopes have been redistributed (and thus changed), another relaxation step is performed
 *   immediately (i.e. before returning from relax()). This means that the criticality of all sites
 *   is re-evaluated and the previously described procedure is repeated. This in turn is repeated
 *   until no site is critical anymore and thus the avalanche stops. The accumulated number of unit
 *   time steps sums up to the \e avalanche \e duration (see AvalancheStatistics::Event::duration).
 *
 *   (see also relax())
 *
 * Boundary conditions: The boundary conditions of the sandbox (Sandbox::getBoundaryConditions()) are ignored.
 * The model's boundary conditions for slope units are open boundaries at all edges. In the picture of
 * a real sandpile (grain stack heights instad of slopes) the interpretation of boundary conditions
 * is a bit unclear. The \e effective boundary conditions for sand grains here seem to be something
 * between i) closed boundaries everywhere except for an open boundary at the corner of the "downhill side"
 * and ii) half-open boundary conditions with closed boundaries for edges meeting at the "uphill side"
 * corner and open boundaries for edges meeting at the "downhill side" corner.
 *
 * Note that the direction of the "diagonal hillside" (i.e. which corner is "uphill" or "downhill")
 * is actually undefined unless conservative driving (see drive()) is used!
 *
 * In order to characterize the avalanches we still need to precisely define the observables in
 * AvalancheStatistics::Event. As explained above, the \e duration of an avalanche is defined as
 * the number of repetitions of the relaxation procedure (within a single call of relax()) needed
 * to obtain a completely stable sandpile again. The \e linear size (see AvalancheStatistics::Event::linSize)
 * is the maximum distance between points of the set of avalanche sites (see calculateAvalancheLinSize()).
 * All lattice sites, where slopes were changed during an avalanche at least once, we call "avalanche sites".
 * The \e area (see AvalancheStatistics::Event::area) is the number of distinct avalanche sites (see calculateAvalancheArea()).
 * The \e size (see AvalancheStatistics::Event::size) in case of this model is the accumulated number of individual
 * critical lattice sites. This means that each lattice site that becomes a critical site the avalanche size
 * is increased by 1. If the same site is critical for multiple iterations of the relaxation procedure
 * of the same avalanche, the avalanche size is increased multiple times accordingly.
 *
 * Note: This definition of the avalanche size may be considered as
 * some sort of measure of the total mass movement during an avalanche.
 */
class SimulationModel_BTW : public SimulationModel
{
public:
    explicit SimulationModel_BTW(Sandbox& pSandbox);    ///< Constructor.
    ~SimulationModel_BTW() override = default;          ///< Default destructor.
    //
    Aux::Model id() const override { return Aux::Model::BakTangWiesenfeld; }    ///< \brief Get an ID identifying the SimulationModel
                                                                                ///         as defined in Aux::Model.
                                                                                ///  \returns Aux::Model::BakTangWiesenfeld.
    //
    bool setModelParameter(const std::string& pKey, int pVal) override;                     ///< Set additional model parameters.
    bool getModelParameter(const std::string& pKey, int& pVal) const override;              ///< Get additional model parameters.
    std::vector<std::pair<std::string, std::string>> listModelParameters() const override;  ///< List available model parameters.
    //
    void drive() const override;                                    ///< Drive the sandpile.
    void relax(AvalancheStatistics::Event& pEvent) override;        ///< Relax the sandpile.
    //
    double calculateSandboxCriticality() const override;            ///< Quantify the criticality of the sandbox.

private:
    double calculateAvalancheLinSize(const Avalanche& pAvalanche) const override;   ///< Calculate the linear size of an avalanche.
    long long calculateAvalancheArea(const Avalanche& pAvalanche) const  override;  ///< Calculate the area of an avalanche.

private:
    bool conservativeDriving;   //Drive the sandpile "conservatively" (instead of the simpler non-conservative driving)?
    //
    const std::vector<short>& currentGrain;     //The random lattice site that was used in the most recent driving step
};

#endif // SANDSIM_SIMULATIONMODEL_BTW_H

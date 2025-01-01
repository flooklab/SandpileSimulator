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

#include "simulationmodel_fwm.h"

/*!
 * \brief Constructor.
 *
 * \copydetails SimulationModel::SimulationModel()
 *
 * \param pCritSlope The value of the critical slope parameter q_crit (see also this class' documentation).
 */
SimulationModel_FWM::SimulationModel_FWM(Sandbox& pSandbox, short pCritSlope) :
    SimulationModel(pSandbox),
    critSlope(pCritSlope),
    avalancheSizeMovedGrains(false),
    currentGrain(sandbox.getCurrentGrain()),
    rndGenerator(),
    rndUnif(0, 1)
{
}

//Public

//

/*!
 * \brief Seed the random generator for initial topple direction randomization.
 *
 * Seeds the random generator that is invoked for every critical lattice site in order
 * to randomize/determine the neighbor in whose direction the grain redistribution starts.
 *
 * \param pSeedSeq The seed sequence.
 */
void SimulationModel_FWM::seed(std::seed_seq& pSeedSeq)
{
    rndGenerator.seed(pSeedSeq);
}

//

/*!
 * \brief Set additional model parameters.
 *
 * For available parameters see listModelParameters().
 *
 * \param pKey Name of the parameter.
 * \param pVal New value of the parameter.
 * \return True if the requested parameter is available and False else.
 */
bool SimulationModel_FWM::setModelParameter(const std::string& pKey, int pVal)
{
    if (pKey == "FWM_criticalSlope")
    {
        critSlope = pVal;
        return true;
    }
    else if (pKey == "FWM_relax_size_movedGrains")
    {
        avalancheSizeMovedGrains = (pVal == 0 ? false : true);
        return true;
    }

    return false;
}

/*!
 * \brief Get additional model parameters.
 *
 * For available parameters see listModelParameters().
 *
 * \param pKey Name of the parameter.
 * \param pVal Value of the parameter.
 * \return True if the requested parameter is available and False else.
 */
bool SimulationModel_FWM::getModelParameter(const std::string& pKey, int& pVal) const
{
    if (pKey == "FWM_criticalSlope")
    {
        pVal = critSlope;
        return true;
    }
    else if (pKey == "FWM_relax_size_movedGrains")
    {
        pVal = (avalancheSizeMovedGrains ? 1 : 0);
        return true;
    }

    return false;
}

/*!
 * \brief List available model parameters.
 *
 * Available parameters are:
 *
 * - FWM_critical_slope: The critical slope (q_crit) defines the height difference between neighboring lattice sites,
                         at which the sandpile becomes locally unstable (see also this class' documentation).
 *
 * - FWM_relax_size_movedGrains: Set to a value different from 0 to switch AvalancheStatistics::Event::size from representing
                                 \e accumulated \e number \e of \e critical \e sites to \e total \e number \e of \e grain
                                 \e movements. The default is 0 (size = number of critical sites). See also relax().
 *
 * \return List of available model parameters as pairs of {"NAME", "DESCRIPTION"}.
 */
std::vector<std::pair<std::string, std::string>> SimulationModel_FWM::listModelParameters() const
{
    return {{"FWM_criticalSlope", "The critical slope (q_crit) defines the height difference between neighboring lattice sites, "
                                  "at which the sandpile becomes locally unstable (which triggers an avalanche)."},
            {"FWM_relax_size_movedGrains", "Set to a value different from 0 to switch avalanche size from representing "
                                           "\"accumulated number of critical sites\" to \"total number of grain movements\". "
                                           "The default is 0 (size = number of critical sites)."}};
}

//

/*!
 * \brief Drive the sandpile.
 *
 * Performs conservative driving of the sandpile according to what is written
 * in the detailed description of this class, i.e. increases the height
 * of the sandpile by 1 at a randomly selected lattice site.
 */
void SimulationModel_FWM::drive() const
{
    sandbox.dropRandom(1);
}

/*!
 * \brief Relax the sandpile.
 *
 * Takes the random lattice site that was selected most recently by the sandbox (a new random lattice site
 * is for example selected for each driving step, which should happen right before the relaxation) as the
 * starting point and from there relaxes the sandpile iteratively until all lattice sites are stable again.
 * For more details about this procedure see this class' detailed description.
 *
 * The statistical information of the single avalanche produced by this relaxation
 * is written to \p pEvent, also if the initial site is not critical and thus
 * no avalanche happens (then \e duration = \e size = \e linSize = \e area = 0).
 * The event ID stays unchanged.
 *
 * Note: By default the \e avalanche \e size is set to the accumulated number of critical sites.
 * This behavior can be changed by toggling the model parameter "FWM_relax_size_movedGrains",
 * see listModelParameters() and SimulationModel_FWM.
 *
 * \param pEvent Destination for the statistical information belonging to the avalanche produced by this function.
 */
void SimulationModel_FWM::relax(AvalancheStatistics::Event& pEvent)
{
    //Empty avalanche
    Avalanche avalanche(sandbox);

    //Keep track of avalanche duration and size (calculate linSize and area at the end)
    decltype(pEvent.duration) tDuration = -1;   //Note: one increment happens below in any case
    decltype(pEvent.size)     tSize     =  0;

    //Count number of grain movements (used as avalanche size, if "FWM_relax_size_movedGrains" is enabled)
    decltype(pEvent.size) numMovedGrains = 0;

    //Potentially unstable lattice sites; begins with just driven site; is updated at the end of every iteration of below while loop
    std::map<long long, struct LatticeSite> potentialPointsMap {{sandbox.cartToLin(currentGrain), {currentGrain, false}}};

    //Relax, re-evaluate, relax, ... sandpile until no new critical sites were created / can be found in 'potentialPointsMap' anymore

    while (true)
    {
        decltype(tSize) numCriticalSites = 0;

        //Record pending grain redistributions (height changes) for the current
        //iteration here and apply them to the sandpile only after this iteration
        std::map<long long, short> heightChanges;

        for (auto& siteIt : potentialPointsMap)
        {
            //Need copy of position vector
            std::vector<short> tPoint = siteIt.second.position;

            //Boundary information for current point of loop
            char cBound = sandbox.getBoundaryInfoAt(tPoint);

            //Not on lattice, continue with next site
            if (cBound == 'x')
                continue;

            //At open boundary, need different procedure
            if (cBound == 'o')
            {
                //Relax site, if critical
                if (sandbox.getHeightAt(tPoint) >= critSlope)
                {
                    //Increase number of critical sites of current iteration
                    ++numCriticalSites;

                    //Remember this for the determination of points to be checked during the next iteration
                    siteIt.second.critical = true;

                    //Remove all grains at critical site (they fall off the table)
                    long long tIdx = sandbox.cartToLin(tPoint);
                    heightChanges[tIdx] -= sandbox.getHeightAt(tIdx);       //Note: have zero initialization if key doesn't exist yet

                    //Add site to the list of avalanche points
                    avalanche.addSite(tPoint);
                }

                continue;
            }

            //Normal procedure

            //Is current site of potentially critical sites actually critical? (determined further below)
            bool isCritical = false;

            //Height of current point of loop
            short centerHeight = sandbox.getHeightAt(tPoint);

            //Go through all nearest neighbors to record their slopes w.r.t. the central site and to see, if that one is critical

            std::map<long long, short> critSlopes;
            short maxSlope = 0;

            for (size_t i = 0; i < sandbox.getShape().size(); ++i)
            {
                //Move left
                --tPoint[i];

                //Boundary information for current neighbor
                char nBound = sandbox.getBoundaryInfoAt(tPoint);

                //Only consider neighbors on the lattice, since open boundary at center point already handled above
                if (nBound != 'x')
                {
                    short slope = centerHeight - sandbox.getHeightAt(tPoint);
                    if (slope >= critSlope)
                    {
                        critSlopes[sandbox.cartToLin(tPoint)] = slope;

                        if (slope > maxSlope)
                            maxSlope = slope;

                        isCritical = true;
                    }
                }

                //Move back and move right
                ++tPoint[i];
                ++tPoint[i];

                //Boundary information for current neighbor
                nBound = sandbox.getBoundaryInfoAt(tPoint);

                //Only consider neighbors on the lattice, since open boundary at center point already handled above
                if (nBound != 'x')
                {
                    short slope = centerHeight - sandbox.getHeightAt(tPoint);
                    if (slope >= critSlope)
                    {
                        critSlopes[sandbox.cartToLin(tPoint)] = slope;

                        if (slope > maxSlope)
                            maxSlope = slope;

                        isCritical = true;
                    }
                }

                //Move back again
                --tPoint[i];
            }

            //Perform the grain redistribution from the central site to its neighbors, if it was determined to be critical;
            //temporarily record the redistributions in 'heightChanges' instead of directly changing the sandbox, in order
            //to achieve simultaneous relaxation for all iterations of this while loop (see also code after while block).

            if (isCritical)
            {
                //Increase number of critical sites of current iteration
                ++numCriticalSites;

                //Remember this for the determination of points to be checked during the next iteration
                siteIt.second.critical = true;

                //Add site to the list of avalanche points
                avalanche.addSite(tPoint);

                //Below:
                //Actual relaxation of current critical site

                //Create iterator for gradually redistributing grains to the neighbors
                auto dropNeighbor = critSlopes.begin();

                //Randomize the choice of the neighbor to which the first drop happens

                int initialDropNeighbor = 0;
                if (critSlopes.size() > 1)
                    initialDropNeighbor = static_cast<int>(critSlopes.size()*rndUnif(rndGenerator));

                for (int i = 0; i < initialDropNeighbor; ++i)
                    ++dropNeighbor;

                //Redistribute at max 'grainsToRedist' grains in total, with 'grainsToRedist' equal
                //the number of excess grains with respect to the neighbor with the steepest slope
                short grainsToRedist = maxSlope;

                /* Repeatedly loop over all neighbors with initially critical slope (starting with the randomly
                   chosen neighbor from above) until maximally 'grainsToRedist' grains have been moved, but drop
                   a grain from the critical site to a neighbor only if the neighbor's height is *less* than the
                   height of the critical site, i.e. the drop goes locally "downhill" (remove the neighbor
                   from the loop otherwise). If no "downhill" neighbor is left, even before 'grainsToRedist'
                   grains have been redistributed, stop the redistribution nevertheless. First only apply
                   the redistributions to 'heightChanges' and not to the sandbox itself, in order to
                   achieve simultaneous relaxation of all currently critical sites (see also further below).
                */

                short centerOffset = 0; // Moving one grain from the central site to a neighbor changes the local slope by two units
                                        // but the removed grain at the central site also changes slopes w.r.t. all other neighbors
                                        // by one unit, so keep track of this contribution separately (*)

                while (grainsToRedist > 0)
                {
                    //Slope larger zero? Redistribute one grain
                    if ((dropNeighbor->second + centerOffset) > 0)
                    {
                        //Accumulate pending changes in the auxiliary map of height changes
                        //for later simultaneous application to the actual sandpile

                        --(heightChanges[siteIt.first]);            //Take one grain from the central, critical site
                        ++(heightChanges[dropNeighbor->first]);     //Add this grain to the current neighbor

                        --(dropNeighbor->second);   //Remove one slope unit for current neighbor, account for the other slope unit
                        --centerOffset;             //by reducing centerOffset (*); note: this offset appears in the loop condition

                        //Moved one grain --> one less grain to redistribute
                        --grainsToRedist;
                        ++numMovedGrains;

                        //Add neighbor to the list of avalanche points
                        avalanche.addSite(tPoint);
                    }
                    else
                    {
                        //Remove current neighbor from the loop if the sandpile has become locally even here
                        dropNeighbor = critSlopes.erase(dropNeighbor);

                        if (critSlopes.empty())
                        {
                            //Redistribution stops if no neighbors of the central, initially critical, site are "downhill" anymore
                            break;
                        }
                        else
                        {
                            //Move to the first neighbor before continuing the loop, if 'erase(...)' returned iterator to end()
                            if (dropNeighbor == critSlopes.end())
                                dropNeighbor = critSlopes.begin();

                            continue;
                        }
                    }

                    //Select the next neighbor; move to the first neighbor again, when out of range
                    if (++dropNeighbor == critSlopes.end())
                        dropNeighbor = critSlopes.begin();
                }
            }
        }

        //Now apply simultaneous grain redistributions to the sandpile
        for (const auto& it : heightChanges)
        {
            sandbox.changeHeightAt(it.first,       it.second    );
            //....................(grain position, height change);
        }

        //Increase avalanche's time duration by 1
        ++tDuration;

        //Increase avalanche size about the number of additional relaxation events / grain stack collapses
        tSize += numCriticalSites;

        //If no relaxation actually happened the avalanche stops at this iteration
        if (numCriticalSites == 0)
            break;

        //Record all neighbors of relaxed sites (i.e. where redistributions happened) as potentially unstable sites for the
        //next while-iteration, since slopes have changed there at least by means of the changed height of the central site

        decltype(potentialPointsMap) newPotentialPointsMap;

        for (auto siteIt : potentialPointsMap)
        {
            //Don't need to check neighbors of non-relaxed site
            if (!siteIt.second.critical)
                continue;

            //Need writable position vector
            std::vector<short> tNeighbor = std::move(siteIt.second.position);

            for (size_t i = 0; i < sandbox.getShape().size(); ++i)
            {
                //Move left
                --tNeighbor[i];

                newPotentialPointsMap[sandbox.cartToLin(tNeighbor)] = {tNeighbor, false};

                //Move back and move right
                ++tNeighbor[i];
                ++tNeighbor[i];

                newPotentialPointsMap[sandbox.cartToLin(tNeighbor)] = {tNeighbor, false};

                //Move back again
                --tNeighbor[i];
            }
        }

        potentialPointsMap.swap(newPotentialPointsMap);

    } // closing brace of 'while (true)'


    //Complete the relaxation: write statistical information to pEvent

    pEvent.duration = tDuration;

    if (avalancheSizeMovedGrains)
        pEvent.size     = numMovedGrains;
    else
        pEvent.size     = tSize;

    if (doLinSize)
        pEvent.linSize = calculateAvalancheLinSize(avalanche);

    if (doArea)
        pEvent.area = calculateAvalancheArea(avalanche);
}

//

/*!
 * \brief Quantify the criticality of the sandbox.
 *
 * Calculates a "criticality parameter" that is designed to range from 0 for an empty sandbox to roughly O(1) for a fully
 * critical sandbox and can hence be used to determine, if the sandpile has reached the state of self-organized criticality.
 * The criticality parameter is approximately compensated for different dimensions, lattice sizes, critical slopes, etc.
 * but the actual value for a critical sandpile might still vary with those parameters and for instance also depends
 * on the sandbox boundary conditions. Note also that the exact value always fluctuates while a simulation is running.
 * This means that sandbox criticality detection using this parameter should generally be based on the
 * \e saturation of the \e averaged value \e within \e some \e epsilon \e neighborhood.
 *
 * See also SimulationModel::calculateSandboxCriticality().
 *
 * \return %Sandbox criticality parameter between 0 and O(1).
 */
double SimulationModel_FWM::calculateSandboxCriticality() const
{
    /*
     * Below: Lambda expression for calculating the slope of a lattice point w.r.t to another lattice point.
     *
     * Defines a lambda expression that takes a point \p currentSite on the sandbox \p box
     * and an offset \p offset and calculates the height difference between \p currentSite
     * and the point shifted by \p offset. The height difference is added to \p slopeSum.
     *
     * Note: If the shifted point lies beyond the sandbox boundaries, the function returns without any action.
     */
    auto lFunc = [](const std::vector<short>& currentSite, const std::vector<short>& offset,
                    const Sandbox& box, long long& slopeSum) -> void
    {
        using OpOverloads_STL_vector::operator+;

        std::vector<short> neighborSite = currentSite + offset;

        //Neighboring point off lattice? Skip.
        for (size_t d = 0; d < box.getShape().size(); ++d)
        {
            if (neighborSite.at(d) < 0 || neighborSite.at(d) >= box.getShape().at(d))
                return;
        }

        slopeSum += std::abs(box.getHeightAt(currentSite) - box.getHeightAt(neighborSite));
    };

    long long slopeSum = 0;

    const size_t sandboxDim = sandbox.getShape().size();

    const std::vector<short> lowerBounds(sandboxDim, 0);
    const std::vector<short> upperBounds(sandbox.getShape());

    //Loop over all dimensions to sum all slopes in all directions
    for (size_t dim = 0; dim < sandboxDim; ++dim)
    {
        std::vector<short> shift(sandboxDim, 0);
        shift.at(dim) = 1;

        //Temporary storage for call to recursiveForLoop below
        std::vector<short> tmpSite(sandboxDim, 0);

        //Loop over the sandbox lattice and sum all slopes (slopeSum) in the direction
        //given by shift by applying lFunc to each point in the nested loop.
        Aux::recursiveForLoop(lFunc, lowerBounds, upperBounds, tmpSite, sandboxDim,
                              shift, sandbox, slopeSum);

        shift.at(dim) = -1; //Opposite direction

        tmpSite = std::vector<short>(sandboxDim, 0);

        //Loop over the sandbox lattice and sum all slopes (slopeSum) in the direction
        //given by shift by applying lFunc to each point in the nested loop.
        Aux::recursiveForLoop(lFunc, lowerBounds, upperBounds, tmpSite, sandboxDim,
                              shift, sandbox, slopeSum);
    }

    long long numSites = 1;
    for (short len : sandbox.getShape())
        numSites *= len;

    return static_cast<double>(slopeSum) / numSites / (critSlope - 1) / sandboxDim / 2;
}

//

/*!
 * \brief Calculate the linear size of an avalanche.
 *
 * The FW model's avalanche linear size equals the largest Euclidean distance
 * between any of the points that were part of the avalanche. These are all sites
 * where at least one sand grain has toppled from or onto during the avalanche.
 *
 * \param pAvalanche The Avalanche containing all the points on the lattice that were part of the considered avalanche.
 *
 * \return The AvalancheStatistics::Event::linSize of the avalanche.
 */
double SimulationModel_FWM::calculateAvalancheLinSize(const Avalanche& pAvalanche) const
{
    //Need only surface points of avalanche, "surface" provided by Avalanche sufficient
    std::vector<std::vector<short>> surfPoints;
    pAvalanche.getSurface(surfPoints);

    //Squared distances
    long long tDistSq = 0, maxDistSq = 0;

    //Calculate squared distances between all pairs of surface points, determine maximum of those

    std::vector<short> point1, point2;

    for (size_t i = 0; i < surfPoints.size(); ++i)  //Avoid double counting --> sum_{i,j,j<i}
    {
        point1 = surfPoints[i];

        for (size_t j = 0; j < i; ++j)
        {
            point2 = surfPoints[j];

            tDistSq = Aux::vectorDistanceSquared(point1, point2);

            //Adjust maximum squared distance if current is larger
            if (tDistSq > maxDistSq)
                maxDistSq = tDistSq;
        }
    }

    //Take square root and return
    return std::sqrt(static_cast<double>(maxDistSq));
}

/*!
 * \brief Calculate the area of an avalanche.
 *
 * The FW model's avalanche area equals the number of distinct points that were part of the avalanche.
 * These are all sites where at least one sand grain has toppled from or onto during the avalanche.
 *
 * \param pAvalanche The Avalanche containing all the points on the lattice that were part of the considered avalanche.
 *
 * \return The AvalancheStatistics::Event::area of the avalanche.
 */
long long SimulationModel_FWM::calculateAvalancheArea(const Avalanche& pAvalanche) const
{
    return pAvalanche.countSites();
}

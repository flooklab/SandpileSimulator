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

#include "simulationmodel.h"

#include <cmath>

/*!
 * \brief Constructor.
 *
 * Constructs a simulation model for a simulation operating on \p pSandbox.
 *
 * Calculation of AvalancheStatistics::Event::linSize and AvalancheStatistics::Event::area
 * is initially disabled here (see also enableRecordLinSize() and enableRecordArea()).
 *
 * \param pSandbox The sandbox containing the sandpile for which the simulation shall be performed.
 *                 Also defines the initial state of the simulated sandpile (which can in principle
 *                 be any (sensible) state). This means, in particular, it can also be a sandpile
 *                 resultant from a previously interrupted simulation, which can be continued later on.
 */
SimulationModel::SimulationModel(Sandbox& pSandbox) :
    sandbox(pSandbox),
    doLinSize(false),
    doArea(false)
{
}

//

/*!
 * \brief Enable/disable the calculation of AvalancheStatistics::Event::linSize for each avalanche event.
 *
 * The linear size calculation might be computationally rather expensive (compared to duration and size).
 * Hence the option for enabling or disabling this calculation is provided here (via setting a
 * control flag that might then of course be checked inside the relax() function).
 *
 * By default the control flag is set to false in the constructor (no linear size calculation),
 * which might be overridden by a concrete implementation of this class.
 * Of course this control flag could also simply be ignored by the implementation.
 *
 * \param pEnable True for enabled linSize calculation, False otherwise.
 */
void SimulationModel::enableRecordLinSize(const bool pEnable)
{
    doLinSize = pEnable;
}

/*!
 * \brief Enable/disable the calculation of AvalancheStatistics::Event::area for each avalanche event.
 *
 * The area calculation might be computationally rather expensive (compared to duration and size).
 * Hence the option for enabling or disabling this calculation is provided here (via setting
 * a control flag that might then of course be checked inside the relax() function).
 *
 * By default the control flag is set to false in the constructor (no area calculation),
 * which might be overridden by a concrete implementation of this class.
 * Of course this control flag could also simply be ignored by the implementation.
 *
 * \param pEnable True for enabled area calculation, False otherwise.
 */
void SimulationModel::enableRecordArea(const bool pEnable)
{
    doArea = pEnable;
}

/*!
 * \brief Check whether the linear size calculation is enabled.
 *
 * See also enableRecordLinSize().
 *
 * \return True if linSize calculation is enabled, False otherwise.
 */
bool SimulationModel::recordingLinSize() const
{
    return doLinSize;
}

/*!
 * \brief Check whether the area calculation is enabled.
 *
 * See also enableRecordArea().
 *
 * \return True if area calculation is enabled, False otherwise.
 */
bool SimulationModel::recordingArea() const
{
    return doArea;
}

//

/*!
 * \brief Calculate the linear size of an avalanche.
 *
 * Generic default implementation:
 * Here the avalanche linear size equals the largest distance between any of the points that were part of the avalanche.
 *
 * Note:
 * Which points are considered part of the avalanche (i.e. are included in \p pAvalanche)
 * of course itself depends on the specific implementation of SimulationModel!
 *
 * \param pAvalanche The Avalanche containing all the points on the lattice that were part of the considered avalanche.
 *
 * \return The AvalancheStatistics::Event::linSize of the avalanche.
 */
double SimulationModel::calculateAvalancheLinSize(const Avalanche& pAvalanche) const
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
 * Generic default implementation:
 * Here the avalanche area equals the number of distinct points that were part of the avalanche.
 *
 * Note:
 * Which points are considered part of the avalanche (i.e. are included in \p pAvalanche)
 * of course itself depends on the specific implementation of SimulationModel!
 *
 * \param pAvalanche The Avalanche containing all the points on the lattice that were part of the considered avalanche.
 *
 * \return The AvalancheStatistics::Event::area of the avalanche.
 */
long long SimulationModel::calculateAvalancheArea(const Avalanche& pAvalanche) const
{
    return pAvalanche.countSites();
}

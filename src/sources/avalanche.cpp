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

#include "avalanche.h"

/*!
 * \brief Constructor.
 *
 * Constructs an empty avalanche suitable for processing avalanches happening in the provided sandbox.
 *
 * \param pSandbox The sandbox in which the avalanche happened or any sandbox with the same shape as that one.
 */
Avalanche::Avalanche(const Sandbox& pSandbox) :
    sandbox(pSandbox),
    sandboxDim(pSandbox.getShape().size()),
    additionCounter(0)
{
}

//Public

/*!
 * \brief Remove all previously added sites from the avalanche.
 *
 * Resets addition counter and removes all avalanche site entries and by this
 * makes the instance ready for processing a new avalanche on the same sandpile.
 */
void Avalanche::clear()
{
    siteInAvalanche.clear();
    additionCounter = 0;
}

/*!
 * \brief Add a lattice site to the avalanche or increase it's addition count.
 *
 * If the site has not been added before it is added and is assigned a count of 1.
 * This count is increased by 1 every time this function is called for the same site again.
 * The total addition counter is increased by 1 as well.
 *
 * \param pSite A lattice site from an avalanche.
 */
void Avalanche::addSite(const std::vector<short>& pSite)
{
    long long idx = sandbox.cartToLin(pSite);
    auto site = siteInAvalanche.find(idx);

    if (site != siteInAvalanche.end())
        ++(site->second.numAdditions);
    else
        siteInAvalanche[idx] = {pSite, 1};

    ++additionCounter;
}

//

/*!
 * \brief Get the number of distinct lattice sites part of the avalanche.
 *
 * Sites that have been added multiple times are only counted \e once by this function.
 *
 * \return The number of distinct lattice sites that are part of the avalanche.
 */
long long Avalanche::countSites() const
{
    return siteInAvalanche.size();
}

/*!
 * \brief Get the total number of additions of sites, i.e. the number of calls to %addSite().
 *
 * Also equals the sum of all individual counts assigned to each added lattice site.
 *
 * Remeber: other than countSites(), this function counts multiple additions of the same site multiple times.
 *
 * \return The total number additions of lattice sites.
 */
long long Avalanche::countAdditions() const
{
    return additionCounter;
}

//

/*!
 * \brief Obtain a surface-like subset of the avalanche, preserving at least the maximum distance between two avalanche sites.
 *
 * Writes a set of points {S} forming a surface-like object to \p pSurfacePoints.
 * You can think of this surface-like object as something between the real surface
 * in the mathematical sense and the convex hull of the avalanche.
 * The diameter of the set {S} (the maximum distance between two points) equals that of the full avalanche.
 *
 * More precisely, the set is obtained by subdividing the avalanche into columns along the axis of the 0-th dimension
 * (could in principle be any available dimension) and then taking only the two boundary points of each column
 * (or one or zero when the column contains only one or zero avalanche points).
 *
 * This function can for instance be used for calculating the maximum distance between any two of the
 * avalanche sites on a reduced set of points. This approach might increase overall speed for large avalanches.
 *
 * \param pSurfacePoints Destination for the set of points forming the avalanche "surface" as defined above.
 */
void Avalanche::getSurface(std::vector<std::vector<short>>& pSurfacePoints) const
{
    //Prepare empty storage for surface points
    pSurfacePoints.clear();

    //Already done if no sites in avalanche
    if (countSites() == 0)
        return;

    //Before locating the surface points, constrain "search region" first

    std::vector<short> lowerBounds(sandbox.getShape());
    std::vector<short> upperBounds(sandboxDim, -1);

    for (const auto& it : siteInAvalanche)
    {
        for (size_t i = 0; i < sandboxDim; ++i)
        {
            if (it.second.position[i] < lowerBounds[i])
                lowerBounds[i] = it.second.position[i];
            if (it.second.position[i] > upperBounds[i])
                upperBounds[i] = it.second.position[i]+1;
        }
    }

    /*
     * Below: Lambda expression for determining the surface points of the avalanche.
     *
     * Defines a lambda expression that takes a point \p currentSite and iterates over
     * all possible points pt with 0-th dim. components pt[0] within [\p lower, \p upper)
     * and components of the other dim.s fixed to those of \p currentSite.
     *
     * The two surface points (w.r.t. 0-th dim.) of this column of the avalanche
     * are determined and only those are appended to \p surfPoints.
     * This is achieved by checking if the points exist as key in \p pointsMap and thus
     * the two outermost points within the interval \p lower <= pt[0] < \p upper that also
     * exist in \p pointsMap are the desired surface points.
     * (Choose the "search region" large enough (like above) in order not to cut off parts of the avalanche.)
     *
     * If there is only a single point / no points of the avalanche in the column,
     * only one point / no points will be added to \p surfPoints.
     */
    auto lFunc = [](std::vector<short> currentSite, const Sandbox& box, const std::map<long long, SiteAdditions>& pointsMap,
                    std::vector<std::vector<short>>& surfPoints, int lower, int upper) -> void
    {
        //Search for "lower" surface point
        for (int i = lower; i < upper; ++i)
        {
            //Update point to new position
            currentSite[0] = i;

            //Add the first point that is found in the map of avalanche sites
            //to the vector of surface points, adjust lower bound and break.
            if (pointsMap.find(box.cartToLin(currentSite)) != pointsMap.end())
            {
                surfPoints.push_back(currentSite);
                lower = i;
                break;
            }
        }

        //Search for "upper" surface point
        for (int i = upper-1; i > lower; --i)   //Dont process lower surface point again
        {
            //Update point to new position
            currentSite[0] = i;

            //Add the first point that is found in the map of avalanche sites to the vector of surface points.
            if (pointsMap.find(box.cartToLin(currentSite)) != pointsMap.end())
            {
                surfPoints.push_back(std::move(currentSite));
                return;
            }
        }
    };

    /* Because the lambda expression already takes care of the 0-th dimension
       of the search region, make the recursive loop below *not* loop over the
       0-th dimension, but keep the corresponding search area boundaries
       such that they can be passed as parameters for the lambda expression. */

    short lb = lowerBounds[0];
    short ub = upperBounds[0];

    lowerBounds[0] = 0;
    upperBounds[0] = 1;

    //Expect less or equal two surface points per column --> reserve memory

    using OpOverloads_STL_vector::operator-;

    int num = 2;
    for (short it : upperBounds-lowerBounds)
        num *= it;

    pSurfacePoints.reserve(num);

    //Temporary storage for call to recursiveForLoop below
    std::vector<short> tmp(sandboxDim);

    //Loop over the search region and add the surface points to
    //pSurfacePoints by applying lFunc to each point in the nested loop.
    Aux::recursiveForLoop(lFunc, lowerBounds, upperBounds, tmp, sandboxDim,
                          sandbox, siteInAvalanche, pSurfacePoints, lb, ub);
}

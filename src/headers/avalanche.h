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

#ifndef AVALANCHE_H
#define AVALANCHE_H

#include <map>

#include "sandbox.h"

/*!
 * \brief Store an avalanche and characterize a few of its features.
 *
 * This is basically a wrapper class for a set of points of the sandbox
 * (lattice sites) that are part of a single avalanche, i.e. sites that
 * were affected by the avalanche at least once. These points can be added
 * to the avalanche using addSite(). Sites that were affected
 * by the avalanche multiple times can be added with addSite() multiple times.
 * The individual number of additions is recorded for each site as well as
 * the total number of calls to addSite().
 *
 * Use clear() when you want to process another avalanche after the current one using the same Avalanche instance.
 *
 * The other member functions may be used to investigate some of the characteristic features of the avalanche.
 * Note that a complete characterization of an avalanche probably requires additional code in the relaxation
 * algorithm of SimulationModel (consider for instance AvalancheStatistics::Event::duration).
 */
class Avalanche
{
public:
    Avalanche(const Sandbox& pSandbox);             ///< Constructor.
    //
    void clear();                                   ///< Remove all previously added sites from the avalanche.
    void addSite(const std::vector<short>& pSite);  ///< Add a lattice site to the avalanche or increase it's addition count.
    //
    long long countSites() const;       ///< Get the number of distinct lattice sites part of the avalanche.
    long long countAdditions() const;   ///< Get the total number of additions of sites, i.e. the number of calls to %addSite().
    //
    void getSurface(std::vector<std::vector<short>>& pSurfacePoints) const; ///< \brief Obtain a surface-like subset of the avalanche,
                                                                            ///         preserving at least the maximum distance
                                                                            ///         between two avalanche sites.

private:
    const Sandbox& sandbox;     //Sandbox containing the sandpile on which the avalanche happens (or different sandbox with same shape)
    const size_t sandboxDim;    //Dimension of that sandbox
    //
    long long additionCounter;  //Total number of calls of addSite()
    //
    /*!
     * \brief Groups the position of a site and the number of times it was added to the avalanche.
     */
    struct SiteAdditions {
        std::vector<short> position;    //Cartesian position of the avalanche site
        int numAdditions = 0;           //Number of additions for the same site
    };
    std::map<long long, SiteAdditions> siteInAvalanche; //A map specifying lattice sites that are part of the avalanche
                                                        //together with the count of additions for the same sites
};

#endif // AVALANCHE_H

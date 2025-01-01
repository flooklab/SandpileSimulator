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

#include "simulationmodel_btw.h"

#include <cmath>

/*!
 * \brief Constructor.
 *
 * \copydetails SimulationModel::SimulationModel()
 */
SimulationModel_BTW::SimulationModel_BTW(Sandbox& pSandbox) :
    SimulationModel(pSandbox),
    conservativeDriving(false),
    currentGrain(sandbox.getCurrentGrain())
{
}

//Public

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
bool SimulationModel_BTW::setModelParameter(const std::string& pKey, const int pVal)
{
    if (pKey == "BTW_drive_conservative")
    {
        conservativeDriving = (pVal == 0 ? false : true);
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
bool SimulationModel_BTW::getModelParameter(const std::string& pKey, int& pVal) const
{
    if (pKey == "BTW_drive_conservative")
    {
        pVal = (conservativeDriving ? 1 : 0);
        return true;
    }

    return false;
}

/*!
 * \brief List available model parameters.
 *
 * Available parameters are:
 *
 * - BTW_drive_conservative: Set to a value different from 0 to enable the conservative driving mode of the BTW model.
                             The default is 0 (non-conservative driving). See also drive().
 *
 * For setting / getting the parameters see setModelParameter() / getModelParameter().
 *
 * \return List of available model parameters as pairs of {"NAME", "DESCRIPTION"}.
 */
std::vector<std::pair<std::string, std::string>> SimulationModel_BTW::listModelParameters() const
{
    return {{"BTW_drive_conservative", "Set to a value different from 0 to enable the conservative driving mode of the BTW model. "
                                       "The default is 0 (non-conservative driving)."}};
}

//

/*!
 * \brief Drive the sandpile.
 *
 * Performs conservative or non-conservative driving of the sandpile according to
 * what is written in the detailed description of this class. For non-conservative
 * driving the slope is increased by 1 at a randomly selected lattice site.
 * For conservative driving the slope is increased by dim[sandbox] at a randomly
 * selected site and decreased by 1 at the uphill neighbors of that site.
 *
 * For toggling between conservative and non-conservative driving use the
 * model parameter "BTW_drive_conservative" (see listModelParameters()).
 *
 * Note: For conservative driving this function reduces the slope of the neighbors to the "left"
 * (i.e. positions with \e lower values on the axes), which in this case effectively defines
 * the "uphill side" of the pile to be on the left for each axis/dimension respectively.
 */
void SimulationModel_BTW::drive() const
{
    if (!conservativeDriving)   //Non-conservative driving
        sandbox.dropRandom(1);
    else                        //Conservative driving
    {
        sandbox.dropRandom(sandbox.getShape().size());

        std::vector<short> tPoint = currentGrain;

        for (size_t i = 0; i < sandbox.getShape().size(); ++i)
        {
            //Move left
            --tPoint[i];

            sandbox.changeHeightAt(tPoint, -1);

            //Move back
            ++tPoint[i];
        }
    }
}

/*!
 * \brief Relax the sandpile.
 *
 * Iteratively relaxes the sandpile until all lattice sites are stable (again).
 * For each iteration the currently/still unstable sites are determined in the beginning
 * and then relaxed serially. A site is unstable, if the slope is larger or equal to
 * the critical slope 2*dim[sandbox]. For each unstable site the site's slope
 * is reduced by 2*dim[sandbox] and the neighboring slopes are increased by 1.
 * For more details about the procedure see this class' detailed description.
 *
 * The statistical information of the single avalanche produced by this relaxation
 * is written to \p pEvent, also if initially no site was critical and thus no
 * avalanche happened (then \e duration = \e size = \e linSize = \e area = 0).
 * The event ID stays unchanged.
 *
 * \param pEvent Destination for the statistical information belonging to the avalanche produced by this function.
 */
void SimulationModel_BTW::relax(AvalancheStatistics::Event& pEvent)
{
    //Empty avalanche
    Avalanche avalanche(sandbox);

    //Keep track of avalanche duration and size (calculate linSize and area at the end)
    decltype(pEvent.duration) tDuration = -1;   //Note: one increment happens below in any case
    decltype(pEvent.size)     tSize     =  0;

    //Fixed critical slope for BTW model
    const short critSlope = 2*sandbox.getShape().size();

    while (true)
    {
        //Determine critical sites
        std::vector<long long> criticalSites = sandbox.findSitesLargerEqualHeight(critSlope);

        //Relax all critical sites in series (possible since amount of redistribution constant for all critical sites)
        for (long long criticalSite : criticalSites)
        {
            //Calculate cartesian position of site (to access the site's neighbors etc.)
            std::vector<short> sitePos(sandbox.getShape().size());
            sandbox.linToCart(criticalSite, sitePos);

            //Relax critical site
            sandbox.changeHeightAt(criticalSite, -critSlope);

            //Add critical site to the list of avalanche points
            avalanche.addSite(sitePos);

            //Loop over sandbox dimensions to find (and relax) all nearest neighbor sites
            for (size_t i = 0; i < sandbox.getShape().size(); ++i)
            {
                //Move left
                --sitePos[i];

                //Only handle neighbors on the lattice
                if (sandbox.getBoundaryInfoAt(sitePos) != 'x')
                {
                    //Relax neighbor site
                    sandbox.changeHeightAt(sitePos, 1);

                    //Add neighbor site to the list of avalanche points
                    avalanche.addSite(sitePos);
                }

                //Move back and move right
                ++sitePos[i];
                ++sitePos[i];

                //Only handle neighbors on the lattice
                if (sandbox.getBoundaryInfoAt(sitePos) != 'x')
                {
                    //Relax neighbor site
                    sandbox.changeHeightAt(sitePos, 1);

                    //Add neighbor site to the list of avalanche points
                    avalanche.addSite(sitePos);
                }

                //Move back again
                --sitePos[i];
            }
        }

        //Increase avalanche's time duration by 1
        ++tDuration;

        //Increase avalanche size about the number of additional relaxation events (total number of critical sites)
        tSize += criticalSites.size();

        //If no relaxation actually happened the avalanche stops at this iteration
        if (criticalSites.empty())
            break;

    } // closing brace of 'while (true)'

    //Complete the relaxation: write statistical information to pEvent

    pEvent.duration = tDuration;
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
 * The criticality parameter is approximately compensated for different dimensions, lattice sizes, etc. but the actual value
 * for a critical sandpile might still vary with those parameters. Note also that the exact value always fluctuates while
 * a simulation is running. This means that sandbox criticality detection using this parameter should generally
 * be based on the \e saturation of the \e averaged value \e within \e some \e epsilon \e neighborhood.
 *
 * See also SimulationModel::calculateSandboxCriticality().
 *
 * \return %Sandbox criticality parameter between 0 and O(1).
 */
double SimulationModel_BTW::calculateSandboxCriticality() const
{
    long long numSites = 1;
    for (short len : sandbox.getShape())
        numSites *= len;

    long long slopeSum = 0;
    for (long long i = 0; i < numSites; ++i)
        slopeSum += sandbox.getHeightAt(i);

    double critSlope = 2. * sandbox.getShape().size();

    return static_cast<double>(slopeSum) / numSites / (critSlope - 1);
}

//

/*!
 * \brief Calculate the linear size of an avalanche.
 *
 * The BTW model's avalanche linear size equals the largest Euclidean distance between
 * any of the points that were part of the avalanche. These are all sites where the slope
 * changed at least once during the avalanche (i.e. in any of the iterations of the
 * relaxation procedure; can be the same slope before and after the complete avalanche).
 *
 * \param pAvalanche The Avalanche containing all the points on the lattice that were part of the considered avalanche.
 *
 * \return The AvalancheStatistics::Event::linSize of the avalanche.
 */
double SimulationModel_BTW::calculateAvalancheLinSize(const Avalanche& pAvalanche) const
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
 * The BTW model's avalanche area equals the number of distinct points that were part of the avalanche.
 * These are all sites where the slope changed at least once during the avalanche (i.e. in any of the
 * iterations of the relaxation procedure; can be the same slope before and after the complete avalanche).
 *
 * \param pAvalanche The Avalanche containing all the points on the lattice that were part of the considered avalanche.
 *
 * \return The AvalancheStatistics::Event::area of the avalanche.
 */
long long SimulationModel_BTW::calculateAvalancheArea(const Avalanche& pAvalanche) const
{
    return pAvalanche.countSites();
}

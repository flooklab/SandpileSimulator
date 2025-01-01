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

#ifndef SIMULATIONMODEL_H
#define SIMULATIONMODEL_H

#include <random>
#include <vector>

#include "noncopyable.h"
#include "sandbox.h"
#include "avalanchestatistics.h"
#include "avalanche.h"

/*! \brief Characterizes/models the dynamics of a sandpile.
 *
 * For characterizing sandpile dynamics using cellular automata one typically
 * investigates avalanches on these sandpiles (see also AvalancheStatistics).
 * As in reality, triggering those avalanches should be achieved by constantly
 * adding sand to the same sandpile. Thus, the simulation of a sandpile cellular
 * automaton must consist of two distinct steps (which will be repeated many times):
 *
 * - 1. "Driving":
 *   Specifies the manner, how sand is externally added to the pile.
 *   The procedure shall not depend on the current state of the sandpile.
 *
 *   [Note:
 *   Generally, as also explained in the documentation of Sandbox,
 *   it specifies any kind of perturbation of the considered system,
 *   as the sandbox lattice sites could also represent something
 *   different from sand grains.]
 *
 * - 2. "Relaxation":
 *   Specifies the sandpile's response on the driving step, assuming
 *   fast relaxation processes compared to the repitition rate of driving.
 *   The procedure shall depend on the current state of the sandpile.
 *
 * The dynamics of a sandpile cellular automaton is fully set by the implementation of those two processes.
 *
 * This class provides an abstract interface for repeatedly performing the driving
 * and relaxation steps on a sandpile. Specific sandpile models can be implemented
 * by deriving from this class.
 */
class SimulationModel : private NonCopyable
{
public:
    SimulationModel(Sandbox& pSandbox);             ///< Constructor.
    virtual ~SimulationModel() = default;           ///< Default destructor.
    //
    virtual Aux::Model id() const = 0;              ///< \brief Get an ID identifying the SimulationModel
                                                    ///         as defined in Aux::Model.
                                                    ///  \returns %SimulationModel ID.
    //
    void enableRecordLinSize(bool pEnable = false); ///< \brief Enable/disable the calculation of AvalancheStatistics::Event::linSize
                                                    ///         for each avalanche event.
    void enableRecordArea(bool pEnable = false);    ///< \brief Enable/disable the calculation of AvalancheStatistics::Event::area
                                                    ///         for each avalanche event.
    bool recordingLinSize() const;                  ///< Check whether the linear size calculation is enabled.
    bool recordingArea() const;                     ///< Check whether the area calculation is enabled.
    //
    /*!
     * \brief Seed any random generators that may be part of a specific model implementation.
     * \param pSeedSeq The seed sequence used to seed everything (in a fixed order).
     */
    virtual void seed(std::seed_seq& pSeedSeq) {
        (void)pSeedSeq; //Suppress unused parameter warning
        return;
    }
    //
    /*!
     * \brief Set additional parameters specific to the respective model.
     *
     * If a specific simulation model needs additional free model parameters,
     * these may be made accessible via the polymorphic interface
     * of the %SimulationModel base class.
     *
     * Doing so is possible by reimplementing this general setter function
     * (for the getter function see getModelParameter()).
     *
     * A list of all available parameters shall be made accessible via listModelParameters().
     *
     * \param pKey Name of the parameter.
     * \param pVal New value of the parameter.
     * \return True if the model provides the requested parameter and False else.
     */
    virtual bool setModelParameter(const std::string& pKey, int pVal) {
        (void)pKey; //Suppress unused parameter warning
        (void)pVal; //Suppress unused parameter warning
        return false;
    }
    /*!
     * \brief Get additional parameters specific to the respective model.
     *
     * If a specific simulation model needs additional free model parameters,
     * these may be made accessible via the polymorphic interface
     * of the %SimulationModel base class.
     *
     * Doing so is possible by reimplementing this general getter function
     * (for the setter function see setModelParameter()).
     *
     * The value of the parameter shall be directly written to the provided reference, the parameter val.
     *
     * A list of all available parameters shall be made accessible via listModelParameters().
     *
     * \param pKey Name of the parameter.
     * \param pVal Value of the parameter.
     * \return True if the model provides the requested parameter and False else.
     */
    virtual bool getModelParameter(const std::string& pKey, int& pVal) const {
        (void)pKey; //Suppress unused parameter warning
        (void)pVal; //Suppress unused parameter warning
        return false;
    }
    /*!
     * \brief List available additional parameters specific to the respective model.
     *
     * If a specific simulation model needs additional free model parameters,
     * these may be made accessible via the polymorphic interface
     * of the %SimulationModel base class.
     *
     * Doing so is possible by reimplementing setModelParameter() and getModelParameter().
     *
     * This function should be reimplemented accordingly, such that it returns a list of available
     * parameters i.e. all parameters that can be set/retrieved by the mentioned two functions.
     *
     * \return List of available model parameters as pairs of {"NAME", "DESCRIPTION"}.
     */
    virtual std::vector<std::pair<std::string, std::string>> listModelParameters() const {
        return {};
    }
    //
    virtual void drive() const = 0;         ///< \brief Drive the sandpile.
                                            ///  \details In a narrower sense, this function adds sand
                                            ///           to one or more (maybe random) lattice site(s).
                                            ///           Generally, this function could perform any
                                            ///           kind of external perturbation of the system/sandpile
                                            ///           that does not depend on the current state of the system.
                                            ///
    virtual void relax(AvalancheStatistics::Event& pEvent) = 0;
                                            ///< \brief Relax the sandpile.
                                            ///  \details In a narrower sense, if any lattice sites have become unstable
                                            ///           due to the driving step, this function redistributes sand grains
                                            ///           over the lattice until all the lattice sites are stable again.
                                            ///           Generally, this function describes the response of the system/sandpile
                                            ///           to the perturbation from the driving step. Contrary to the driving,
                                            ///           the relaxation shall depend on the current state of the system/sandpile.
                                            ///  \param pEvent Destination for the statistical information belonging to
                                            ///                the avalanche produced by this function.
    //
    virtual double calculateSandboxCriticality() const = 0; ///< \brief Quantify the criticality of the sandbox.
                                                            ///  \details Calculates a scalar that quantitatively represents the
                                                            ///           criticality of the sandbox. Ideally it should range from 0
                                                            ///           (empty, not critical) to 1 (every site just below critical),
                                                            ///           which should (also ideally) not depend on sandbox dimension
                                                            ///           or size. This is easy to achieve for simple models such as
                                                            ///           the BTW model but might be less straightforward for other
                                                            ///           models with configurable parameters or boundary conditions.
                                                            ///           Hence a range from 0 to ~O(1) should be ok. The parameter
                                                            ///           should in any case grow with continuing drives and then
                                                            ///           saturate as the sandbox becomes overall critical
                                                            ///           (i.e. stationary) on average. Note that for a critical
                                                            ///           sandbox the parameter does not have to be at its maximum
                                                            ///    possible value (not every single site must be just below critical).
                                                            ///  \return The criticality parameter.

protected:
    virtual double calculateAvalancheLinSize(const Avalanche& pAvalanche) const = 0;   ///< Calculate the linear size of an avalanche.
    virtual long long calculateAvalancheArea(const Avalanche& pAvalanche) const = 0;    ///< Calculate the area of an avalanche.

protected:
    Sandbox& sandbox;   //The sandbox containing the sandpile for which the simulation is performed
    //
    bool doLinSize;     //Calculate linear size for every avalanche event? (default: false)
    bool doArea;        //Calculate area for every avalanche event? (default: false)
};

#endif // SIMULATIONMODEL_H

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

#include "simulator.h"

/*!
 * \brief Construct a sandpile simulation for a certain type of simulation model and a sandbox of given shape.
 *
 * Prepares the %Simulator class for automatically performing the simulation on a sandbox
 * of given shape and using a SimulationModel according to the specified model string.
 * For this purpose, an appropriate empty sandbox and a simulation model are created,
 * which can be externally accessed via the pointers provided by getSandbox() and getModel().
 * The gathered statistics are written to an instance of AvalancheStatistics, which is internally
 * created as well and can be externally accessed via the pointer provided by getStatistics().
 *
 * Moreover, since %Simulator can take care of logging and also arranging plotting
 * of the sandpile if requested later on, this will be prepared here as well.
 *
 * \param pModelId Specifies the simulation model to be used.
 * \param pShape Dimensions of the sandbox to be created.
 * \param pSandboxSeedSeq Seed sequence used to seed the sandbox.
 * \param pLogger The logger to be used for logging. A new logger is created if \e nullptr is passed.
 *
 * \throws std::invalid_argument Invalid simulation model \p pModelId.
 */
Simulator::Simulator(Aux::Model pModelId, std::vector<short> pShape, std::seed_seq& pSandboxSeedSeq,
                     std::shared_ptr<Logger> pLogger) :
    sandbox(std::make_shared<Sandbox>(std::move(pShape), pSandboxSeedSeq)),
    stats(std::make_shared<AvalancheStatistics>()),
    model(nullptr),
    logger(pLogger),
#ifdef ENABLE_PLOTTER
    plotter(*sandbox),
#endif // ENABLE_PLOTTER
    plottingEnabled(false)
{
    if (logger == nullptr)
        logger = std::make_shared<Logger>(Logger::LogLevel::_INFO);

    try
    {
        model = createSimulationModel(pModelId, *sandbox);
    }
    catch (const std::invalid_argument&)
    {
        logger->logCritical("The specified simulation model is not available!");
        throw;
    }
}

/*!
 * \brief Construct a sandpile simulation from existing sandbox, statistics and simulation model.
 *
 * Prepares the %Simulator class for automatically performing the simulation on the provided sandbox using the
 * provided simulation model and writing the gathered statistics to the provided statistics collector and,
 * moreover, taking care of logging and also arranging plotting of the sandpile if requested later on.
 *
 * \param pSandbox An existing instance of Sandbox. Defines the initial state of the simulated sandpile,
 *                 i.e. will not be modified by this class until you call runSimulation().
 * \param pStatistics An existing instance of AvalancheStatistics,
 *                    destination for all statistics gathered during the simulation.
 * \param pSimulationModel An existing instance of SimulationModel used to drive and relax the sandpile.
 * \param pLogger The Logger to be used for logging. A new logger is created if \e nullptr is passed.
 */
Simulator::Simulator(std::shared_ptr<Sandbox> pSandbox, std::shared_ptr<AvalancheStatistics> pStatistics,
                     std::shared_ptr<SimulationModel> pSimulationModel, std::shared_ptr<Logger> pLogger) :
    sandbox(std::move(pSandbox)),
    stats(std::move(pStatistics)),
    model(std::move(pSimulationModel)),
    logger(pLogger),
#ifdef ENABLE_PLOTTER
    plotter(*sandbox),
#endif // ENABLE_PLOTTER
    plottingEnabled(false)
{
    if (logger == nullptr)
        logger = std::make_shared<Logger>(Logger::LogLevel::_INFO);
}

//Public

/*!
 * \brief Dynamically create a specific SimulationModel by its model ID.
 *
 * Factory for SimulationModel implementations.
 *
 * Note: Initializes the critical slope of SimulationModel_FWM to 8.
 *
 * \param pModelId Unique integer model ID.
 * \param pSandbox See SimulationModel::SimulationModel().
 * \return Shared base pointer to a new SimulationModel.
 *
 * \throws std::invalid_argument Invalid simulation model \p pModelId.
 */
std::shared_ptr<SimulationModel> Simulator::createSimulationModel(Aux::Model pModelId, Sandbox& pSandbox)
{
    switch (pModelId)
    {
        case Aux::Model::_BTW:
        {
            return std::make_shared<SimulationModel_BTW>(pSandbox);
            break;
        }
        case Aux::Model::_FWM:
        {
            return std::make_shared<SimulationModel_FWM>(pSandbox, 8);
            break;
        }
        default:
        {
            throw std::invalid_argument("The specified simulation model is not available!");
            break;
        }
    }
}

//

/*!
 * \brief Get the sandbox on which the simulation is done.
 *
 * \return A shared pointer to the used sandbox.
 */
std::shared_ptr<Sandbox> Simulator::getSandbox()
{
    return sandbox;
}

/*!
 * \brief Get the model used by simulator to actually simulate the sandpile.
 *
 * \return A shared pointer to the used simulation model.
 */
std::shared_ptr<SimulationModel> Simulator::getModel()
{
    return model;
}

/*!
 * \brief Get the simulation statistics about avalanches on the sandpile.
 *
 * \return A shared pointer to the avalanche statistics.
 */
std::shared_ptr<AvalancheStatistics> Simulator::getStatistics()
{
    return stats;
}

/*!
 * \brief Get the logger used for logging.
 *
 * \return A shared pointer to the logger.
 */
std::shared_ptr<Logger> Simulator::getLogger()
{
    return logger;
}

//

/*!
 * \brief Enable live plotting of the sandbox content.
 *
 * Uses an internal Plotter instance in order to start a separate plotting thread, which will plot
 * the evolution of the sandpile while the simulation is running. See Plotter::init() for details.
 *
 * \param pFramerateLimit Maximal update frequency of displayed sandbox content.
 */
void Simulator::enablePlotting(unsigned int pFramerateLimit)
{
#ifdef ENABLE_PLOTTER
    if (plottingEnabled)
        return;

    plottingEnabled = true;
    plotter.init(pFramerateLimit);
#else
    (void)pFramerateLimit;  //Suppress unused parameter warning
    logger->logWarning("Cannot enable plotting: Plotting feature has not been included in this build.");
#endif // ENABLE_PLOTTER
}

/*!
 * \brief Disable live plotting of the sandbox content.
 *
 * If a plotting thread is running, it will be terminated by this function call.
 */
void Simulator::disablePlotting()
{
#ifdef ENABLE_PLOTTER
    if (!plottingEnabled)
        return;

    plottingEnabled = false;
    plotter.stopPlotting();
#else
    logger->logWarning("Cannot disable plotting: Plotting feature has not been included in this build.");
#endif // ENABLE_PLOTTER
}

//

/*!
 * \brief Start the simulation.
 *
 * A loop is entered, which alternatingly triggers \p pNumDrives driving steps and \p pNumDrives relaxation steps.
 * The avalanche event characterization (AvalancheStatistics::Event) resulting from each driving step is appended
 * to the AvalancheStatistics instance (see constructors). Returns after \p pNumDrives relaxation steps have been performed.
 *
 * \param pNumDrives Number of driving & relaxation steps to be performed.
 *
 * \throws std::out_of_range \p pNumDrives is negative.
 */
void Simulator::runSimulation(long long pNumDrives)
{
    if (pNumDrives < 0)
        throw std::out_of_range("Number of drives/relaxations must not be negative!");

    std::string shapeString("{");
    shapeString += std::to_string(sandbox->getShape()[0]);

    for (size_t i = 1; i < sandbox->getShape().size(); ++i)
        shapeString += ", " + std::to_string(sandbox->getShape()[i]);

    shapeString += "}";

    logger->logDebug("#############################################");
    logger->logDebug("");
    logger->logDebug("Simulation configuration:");
    logger->logDebug("\tSandpile: " + shapeString);
    logger->logDebug("\tSimulation model: " + Aux::getModelStringFromId(model->id()));
    logger->logDebug("\tNumber of drives: " + std::to_string(pNumDrives));
    logger->logDebug("");
    logger->logDebug("+--------------------------------------------");
    logger->logMore("Start simulation.");
    logger->logDebug("");

    long long logProgressTriggerCtr = 1;
    long long logProgressTriggerThr = pNumDrives / 10;

    auto startingTime = std::chrono::high_resolution_clock::now();

    for (long long i = 0; i < pNumDrives; ++i)
    {
        if (logProgressTriggerCtr++ >= logProgressTriggerThr)
        {
            logProgressTriggerCtr = 1;
            auto currentDuration = std::chrono::duration_cast<std::chrono::minutes>(
                                       std::chrono::high_resolution_clock::now()-startingTime
                                       ).count();

            logger->logDebugDebug("Drive number " + std::to_string(i+1) + "\tafter " + std::to_string(currentDuration) + "min.");
        }

        AvalancheStatistics::Event event;

        model->drive();
        model->relax(event);

        stats->insert(event);
    }

    auto totalDuration = std::chrono::duration_cast<std::chrono::seconds>(
                             std::chrono::high_resolution_clock::now()-startingTime
                             ).count();
    auto durationPerDrive = std::chrono::duration_cast<std::chrono::microseconds>(
                                std::chrono::high_resolution_clock::now()-startingTime
                                ).count() / static_cast<double>(pNumDrives);

    int minutes = static_cast<int>(static_cast<double>(totalDuration)/60.);
    int seconds = static_cast<int>(static_cast<double>(totalDuration)-60*minutes);

    logger->logMore("Finished simulation.");
    logger->logDebug("+--------------------------------------------");
    logger->logDebug("");
    logger->logVerbose("Total duration: " + std::to_string(minutes) + "min " + std::to_string(seconds) + "s.");
    logger->logVerbose("Duration per drive: " + std::to_string(durationPerDrive) + "us.");
    logger->logDebug("");
    logger->logDebug("#############################################");
    logger->logDebug("");
}

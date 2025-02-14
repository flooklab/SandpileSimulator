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

#include "simulationmanager.h"

#include "momentanalysis.h"
#include "simulator.h"

#include <omp.h>

#include <chrono>
#include <cmath>
#include <deque>
#include <sstream>
#include <stdexcept>

/*!
 * \brief Constructor.
 *
 * Creates a new logger, if \p pLogger is \e nullptr.
 *
 * \param pLogger The Logger to be used for logging.
 */
SimulationManager::SimulationManager(std::shared_ptr<Logger> pLogger) :
    logger(std::move(pLogger))
{
    if (logger == nullptr)
        logger = std::make_shared<Logger>(Logger::LogLevel::Info);
}

//Public

/*!
 * \brief Get the logger used for logging.
 *
 * \return A shared pointer to the logger.
 */
std::shared_ptr<Logger> SimulationManager::getLogger()
{
    return logger;
}

//

/*!
 * \brief Perform repeated simulation runs for different lattice sizes to obtain avalanche observable distribution moment samples.
 *
 * Runs \p pNumSamples simulations (see Simulator::runSimulation()) for each sandbox in \p pCritSandboxes.
 * The used simulation model is \p pSimModel and each simulation comprises \p pNumDrives drives.
 * From the AvalancheStatistics of each simulation, the first and second distribution moments are calculated
 * (see MomentAnalysis::calcMoments()) and assigned to \p pSamplesPerSize_moments1 and \p pSamplesPerSize_moments2,
 * respectively, for each lattice size (map key) and for each repitition out of \p pNumSamples (vector index).
 *
 * All sandboxes must be hypercubes of dimension \p pSandboxDim and have edge lengths as given as keys of
 * \p pCritSandboxes. The repeated simulations all use copies of the same initial sandbox from \p pCritSandboxes.
 * To obtain varying results the simulation model is seeded differently each time.
 * The simulations for a single lattice size run in \p pNumThreads parallel OMP threads.
 * The number of threads is automatically set by OMP, if \p pNumThreads is 0 or negative.
 *
 * To speed up the simulations, recording of avalanche area and avalanche linear size statistics (slower than size and
 * duration calculation but \e probably \e small \e effect) can be disabled via \p pRecordArea and \p pRecordLinSize.
 *
 * Additional simulation model parameters might be provided as a vector
 * via \p pModelParameters (see also SimulationModel::setModelParameter()).
 *
 * Plotting/visualization of the sandpiles during simulations can be enabled via \p pEnablePlotting (creates extra thread!).
 * See also Simulator::enablePlotting(). Note: Enabling plotting requires \p pNumThreads == 1.
 *
 * \param pCritSandboxes Critical sandboxes with dimension \p pSandboxDim (hypercubes!) and different edge lengths (map keys).
 * \param pSandboxDim Sandbox dimension of all \p pCritSandboxes.
 * \param pSimModel Simulation model type to use for all simulations.
 * \param pSeedSeq Seed sequence used to seed the random generator for seeding the sandboxes and simulation models.
 * \param pSamplesPerSize_moments1 Resulting first avalanche observable moments (\p pNumSamples samples per lattice size).
 * \param pSamplesPerSize_moments2 Resulting second avalanche observable moments (\p pNumSamples samples per lattice size).
 * \param pNumDrives Number of drives for each simulation.
 * \param pNumSamples Number of repeated simulations per lattice size.
 * \param pRecordArea Enable recording of avalanche area.
 * \param pRecordLinSize Enable recording of avalanche linear size.
 * \param pModelParameters Additional parameters for the SimulationModel.
 * \param pNumThreads Number of parallel simulations.
 * \param pEnablePlotting Enable sandpile visualization during simulations (only works for \p pNumThreads == 1).
 *
 * \throws std::invalid_argument Invalid simulation model \p pSimModel.
 */
void SimulationManager::runFullSimulation(const std::vector<std::pair<short, Sandbox>>& pCritSandboxes,
                                          const short pSandboxDim, const Aux::Model pSimModel, std::seed_seq& pSeedSeq,
                                          std::map<short, std::vector<AvalancheStatistics::Moments>>& pSamplesPerSize_moments1,
                                          std::map<short, std::vector<AvalancheStatistics::Moments>>& pSamplesPerSize_moments2,
                                          const long long pNumDrives, const int pNumSamples, const bool pRecordArea, const bool pRecordLinSize,
                                          const std::vector<std::pair<std::string, int>> pModelParameters,
                                          int pNumThreads, const bool pEnablePlotting) const
{
    if (pSimModel == Aux::Model::InvalidModel)
    {
        logger->logCritical("The specified simulation model is not available!");
        throw std::invalid_argument("The specified simulation model is not available!");
    }

    logger->logInfo("Simulation configuration:");
    logger->logInfo("\tDimension: " + std::to_string(pSandboxDim));
    logger->logInfo("\tSimulation model: " + Aux::getModelStringFromId(pSimModel));
    logger->logInfo("\tNumber of samples: " + std::to_string(pNumSamples));
    logger->logInfo("\tDrives per sample: " + std::to_string(pNumDrives));
    logger->logInfo("+--------------------------------------------");

    //Use default number of threads, if not set to sensible value
    if (pNumThreads < 1)
        pNumThreads = omp_get_max_threads();

    omp_set_dynamic(0);
    omp_set_num_threads(pNumThreads);

    std::mt19937_64 rndGenerator(pSeedSeq);

    for (const auto& sandboxIt : pCritSandboxes)
    {
        const short len = sandboxIt.first;
        const Sandbox& initialSandbox = sandboxIt.second;

        std::vector<AvalancheStatistics::Moments> sample_moments1, sample_moments2;

        omp_lock_t vectorLock;
        omp_init_lock(&vectorLock);

        logger->logLess("Starting simulations for lattice size " + std::to_string(len) + ".");

        #pragma omp parallel for ordered
        for (int i = 0; i < pNumSamples; ++i)
        {
            //First use empty seed sequence for sandbox automatically generated in the simulator class.
            //Sandbox is overwritten and seeded later.
            std::seed_seq tSeedSeq {};

            Simulator simulator(pSimModel, std::vector<short>(pSandboxDim, len), tSeedSeq, logger);

            //Enable plotting, if requested and if only single-threaded simulation (crash otherwise)
            if (pEnablePlotting && pNumThreads == 1)
                simulator.enablePlotting();

            Sandbox& sandbox = *simulator.getSandbox();
            SimulationModel& model = *simulator.getModel();
            AvalancheStatistics& statistics = *simulator.getStatistics();

            //Start all simulations with the same critical sandbox but with a different random generator seed.

            sandbox = initialSandbox;

            #pragma omp ordered
            {
                std::seed_seq tSeedSeq1 = Aux::generateSeedSeq(rndGenerator);
                sandbox.seed(tSeedSeq1);
            }

            //Configure and seed simulation model

            #pragma omp ordered
            {
                std::seed_seq tSeedSeq2 = Aux::generateSeedSeq(rndGenerator);
                model.seed(tSeedSeq2);
            }

            model.enableRecordArea(pRecordArea);
            model.enableRecordLinSize(pRecordLinSize);

            for (const auto& modelParameter : pModelParameters)
            {
                if (!model.setModelParameter(modelParameter.first, modelParameter.second))
                {
                    logger->logWarning("Model parameter \"" + modelParameter.first + "\" is not available for simulation model \"" +
                                       Aux::getModelStringFromId(pSimModel) + "\".");
                }
            }

            //Run actual simulation
            simulator.runSimulation(pNumDrives);

            if (pEnablePlotting)
                simulator.disablePlotting();

            //Get list of avalanche events
            std::vector<AvalancheStatistics::Event> events;
            statistics.detachEvents(events);

            //Calculate first and second moments of event observables

            AvalancheStatistics::Moments firstMoments, secondMoments;
            MomentAnalysis::calcMoments(events, firstMoments, secondMoments);

            omp_set_lock(&vectorLock);
            sample_moments1.push_back(firstMoments);
            sample_moments2.push_back(secondMoments);
            omp_unset_lock(&vectorLock);
        }
        omp_destroy_lock(&vectorLock);

        pSamplesPerSize_moments1.insert(std::make_pair(len, sample_moments1));
        pSamplesPerSize_moments2.insert(std::make_pair(len, sample_moments2));

        logger->logLess("Finished simulations for lattice size " + std::to_string(len) + ".");
        logger->logInfo("");
    }
}

//

/*!
 * \brief Try to automatically drive sandboxes to state of self-organized criticality. Experimental!
 *
 * Drives (empty) sandboxes \p pSandboxes of dimension \p pSandboxDim until a critical state is "detected".
 * Detection of criticality is based on the "criticality parameter" SimulationModel::calculateSandboxCriticality(),
 * which is assumed to range from 0 for an empty sandbox to roughly O(1) for a critical sandbox. The criticality
 * parameter should be designed such that it saturates when criticality is reached and it should optimally be
 * independent of sandbox dimension and lattice size.
 *
 * The drives are accumulated by simulations (see Simulator::runSimulation()) of \p pDrivesBunch drives each,
 * using SimulationModel \p pSimModel (set model parameters via \p pModelParameters). After each simulation
 * the criticality parameter is evaluated and fed into a moving average over 50 simulations. The change of
 * this moving average between two simulations is compared to \p pEpsilon in order to check criticality.
 * At least \p pMinNumDrives will be performed in total. Maximally \p pMaxNumDrives will be performed in total.
 *
 * After reaching either criticality or \p pMaxNumDrives drives an additional simulation
 * will be performed to double the total number of drives to "ensure" criticality.
 *
 * Note:
 * The actual behavior of SimulationModel::calculateSandboxCriticality() depends on the used SimulationModel and possibly
 * on sandbox boundary conditions, might not be as independent of sandbox dimension and lattice size as expected,
 * saturates much slower (in terms of number of drives) per drive for larger lattice sizes, can fluctuate (not monotonic), etc.
 * This means that tuning parameters \p pEpsilon, \p pDrivesBunch, \p pMinNumDrives and \p pMaxNumDrives to obtain
 * reliable results might not always be possible, especially if \p pSandboxes spans a large range of lattice sizes.
 * %Sandbox criticality (i.e. stationary on average) should hence always be checked manually, especially for larger lattice sizes.
 *
 * Note:
 * Plotting the sandboxes during simulation (enable with \p pEnablePlotting) may help with checking, if criticality is really reached.
 *
 * \param pSandboxes Empty sandboxes with dimension \p pSandboxDim (hypercubes!) and different edge lengths (map keys).
 * \param pSandboxDim %Sandbox dimension of all \p pCritSandboxes.
 * \param pMinNumDrives (Half of) minimum number of drives.
 * \param pMaxNumDrives (Half of) maximum number of drives.
 * \param pDrivesBunch Number of drives between successive evaluations of criticality parameter.
 * \param pEpsilon Threshold to detect saturation of criticality parameter.
 * \param pSimModel Simulation model type to use for simulations.
 * \param pSeedSeq Seed sequence used to seed the random generator for seeding the sandboxes and simulation models.
 * \param pModelParameters Additional parameters for the SimulationModel.
 * \param pEnablePlotting Enable sandpile visualization during simulations.
 *
 * \throws std::invalid_argument Invalid simulation model \p pSimModel.
 *
 * \attention This function is not reliable! Always check manually (by plotting or so),
 *            if the sandbox state is approximately stationary!
 */
void SimulationManager::makeSandboxesCritical(std::vector<std::pair<short, Sandbox>>& pSandboxes,
                                              const short pSandboxDim, const long long pMinNumDrives, const long long pMaxNumDrives,
                                              const long long pDrivesBunch, const double pEpsilon, const Aux::Model pSimModel,
                                              std::seed_seq& pSeedSeq, const std::vector<std::pair<std::string, int>> pModelParameters,
                                              const bool pEnablePlotting) const
{
    if (pSimModel == Aux::Model::InvalidModel)
    {
        logger->logCritical("The specified simulation model is not available!");
        throw std::invalid_argument("The specified simulation model is not available!");
    }

    logger->logInfo("+--------------------------------------------");
    logger->logInfo("Driving sandpiles to state of SOC.");
    logger->logInfo("");
    logger->logInfo("Configuration:");
    logger->logInfo("\tDimension: " + std::to_string(pSandboxDim));
    logger->logInfo("\tSimulation model: " + Aux::getModelStringFromId(pSimModel));
    logger->logInfo("\tMin. number of drives: " + std::to_string(pMinNumDrives));
    logger->logInfo("\tMax. number of drives: " + std::to_string(pMaxNumDrives));
    logger->logInfo("\tEpsilon: " + Aux::numberToStringScientific(pEpsilon));
    logger->logInfo("");
    logger->logInfo("+--------------------------------------------");

    std::mt19937_64 rndGenerator(pSeedSeq);

    for (auto& sandboxIt : pSandboxes)
    {
        const short len = sandboxIt.first;
        Sandbox& initialSandbox = sandboxIt.second;

        //First use empty seed sequence for sandbox automatically generated in the simulator class.
        //Sandbox is overwritten and seeded later.
        std::seed_seq tSeedSeq {};

        Simulator simulator(pSimModel, std::vector<short>(pSandboxDim, len), tSeedSeq, logger);

        //Enable plotting to verify criticality
        if (pEnablePlotting)
            simulator.enablePlotting();

        Sandbox& sandbox = *simulator.getSandbox();
        SimulationModel& model = *simulator.getModel();

        //Seed the sandbox

        sandbox = initialSandbox;

        std::seed_seq tSeedSeq1 = Aux::generateSeedSeq(rndGenerator);
        sandbox.seed(tSeedSeq1);

        //Configure and seed simulation model

        std::seed_seq tSeedSeq2 = Aux::generateSeedSeq(rndGenerator);
        model.seed(tSeedSeq2);

        //Do not need statistics
        model.enableRecordArea(false);
        model.enableRecordLinSize(false);

        for (const auto& modelParameter : pModelParameters)
        {
            if (!model.setModelParameter(modelParameter.first, modelParameter.second))
            {
                logger->logWarning("Model parameter \"" + modelParameter.first + "\" is not available for simulation model \"" +
                                   Aux::getModelStringFromId(pSimModel) + "\".");
            }
        }

        logger->logLess("Starting \"criticalization\" for lattice size " + std::to_string(len) + ".");

        long long totalNumDrives = 0;

        std::deque<double> criticalityDeque;

        //Measure time
        auto timeBeginBunch = std::chrono::high_resolution_clock::now();

        //Repeat simulation with small number of drives until sandbox is critical
        for (double criticalityMovAv = 0, prevCriticalityMovAv = 0;;)
        {
            timeBeginBunch = std::chrono::high_resolution_clock::now();

            //Run actual simulation
            simulator.runSimulation(pDrivesBunch);

            totalNumDrives += pDrivesBunch;

            //Calculate criticality factor of sandbox using algorithm specific to simulation model
            double criticality = model.calculateSandboxCriticality();

            criticalityDeque.push_back(criticality);
            if (criticalityDeque.size() > 50)
                criticalityDeque.pop_front();

            criticalityMovAv = 0;
            for (double cr : criticalityDeque)
                criticalityMovAv += cr;
            criticalityMovAv /= criticalityDeque.size();

            logger->logInfo("Make critical... Drives: " + std::to_string(totalNumDrives) +
                            "; Criticality parameter (moving average): " + std::to_string(criticalityMovAv));

            if (totalNumDrives >= pMaxNumDrives)
            {
                logger->logWarning("Reached maximum number of drives.");
                break;
            }

            //Assume that sandbox is critical, if criticality factor saturates; drive at least pMinNumDrives times
            if (totalNumDrives >= pMinNumDrives && std::abs(criticalityMovAv - prevCriticalityMovAv) < pEpsilon)
                break;

            prevCriticalityMovAv = criticalityMovAv;
        }

        std::ostringstream tTimeStrm;
        tTimeStrm.precision(1);
        tTimeStrm<<std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-timeBeginBunch
                                                                         ).count() * (1/60.e6) * totalNumDrives / pDrivesBunch;

        logger->logInfo("Double the number of drives (estimated time: ~" + tTimeStrm.str() + " min) ...");

        //Criticality saturation not reliable, so run another simulation to double the number of drives
        simulator.runSimulation(totalNumDrives);
        totalNumDrives *= 2;

        if (pEnablePlotting)
            simulator.disablePlotting();

        logger->logLess("Finished \"criticalization\" for lattice size " + std::to_string(len) + ".");
        logger->logInfo("Reached criticality factor of " + std::to_string(model.calculateSandboxCriticality()) + " after " +
                          std::to_string(totalNumDrives) + " drives.");
        logger->logInfo("");

        initialSandbox = sandbox;
    }
}

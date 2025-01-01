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

#include "argumentparser.h"
#include "aux.h"
#include "avalanchestatistics.h"
#include "logger.h"
#include "momentanalysis.h"
#include "simulationmanager.h"
#include "simulationmodel.h"
#include "simulationmodel_fwm.h"
#include "simulator.h"

#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <functional>
#include <ios>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

int main(int argc, const char** argv)
{
    //Create and seed main pRNG, which will be used to generate all other required seed sequences
    //(may seed again explicitly via command line argument below)
    std::random_device rndDevice;
    std::seed_seq tSeedSeq {rndDevice(), rndDevice(), rndDevice(), rndDevice(),
                            rndDevice(), rndDevice(), rndDevice(), rndDevice()};
    std::mt19937_64 mainRndGenerator(tSeedSeq);

    //Set default values before parsing command line arguments

    std::string modelString = Aux::getModelStringFromId(Aux::Model::BakTangWiesenfeld);

    std::string modelParams = "";   //Additional model-specific parameters

    bool listParams = false;        //Only list available model parameters and return

    short sandboxDim = 2;           //Dimension of the sandbox

    short lowerLength = 10;         //Smallest lattice size
    short upperLength = 120;        //Largest lattice size
    short lengthIncrement = 10;     //Step size for simulations of different lattice sizes (between lowerLength and upperLength)

    std::string boundaryConditions = "cccc";    //Boundary conditions for each edge of each sandbox dimension

    int numDrives = 100000;         //Number of drives for each simulation sample for each lattice size

    int numSamples = 30;            //Number of samples per lattice size

    int numThreads = 1;             //Number of simulation threads for parallel computation of multiple samples

    bool makeCritical = false;      //Fill sandboxes (before simulation) until they reach state of SOC
    bool skipSim = false;           //Skip the simulation and moment analysis (to only create critical sandboxes)

    std::string criticalizationParams = ""; //Parameters to tune automatic sandbox 'criticalization'

    bool doMomentAnalysis = false;  //Perform moment analysis after simulation

    int numBootstraps = 50*numSamples;          //Number of bootstraps for moment analysis
    int numInnerBootstraps = 5*numBootstraps;   //Number of double bootstraps

    std::string loadSandboxPrefix = "";     //Path and filename prefix to load sandboxes from (postfix defined by model and size)
    std::string saveSandboxPrefix = "";     //Path and filename prefix to save sandboxes to (postfix defined by model and size)
    std::string saveResultsFileName = "";   //File name to save simulation results and moment analysis results to.

    Logger::LogLevel logLevel = Logger::LogLevel::Info;     //Log level

    std::string logFileName = "";                           //Save log output to this file name

    bool enablePlotting = false;    //Visualize sandpiles during simulations

    long long rndSeed = -1;         //Seed for main pRNG (only use positive seeds; keep automatic seed from above, if negative)

    try
    {
        typedef ArgumentParser::Argument Argument;
        typedef ArgumentParser::Argument::Type Type;

        ArgumentParser argParser("SandpileSimulator", argc, argv,
                                 {Argument(Type::String, 'm', "simulation-model",
                                                              "The simulation model to be used.",
                                                              "MODEL-ID",
                                                              true),
                                  Argument(Type::String, 'p', "model-parameters",
                                                              "List of integer parameters specific to the used simulation model.",
                                                              "PAR1=NUM1,PAR2=NUM2,...",
                                                              true,
                                                              {"simulation-model"}),
                                  Argument(Type::None,   'P', "list-parameters",
                                                              "List all model-specific parameters "
                                                              "(for parameter \"--model-parameters\") and return.",
                                                              "",
                                                              true,
                                                              {"simulation-model"}),
                                  Argument(Type::Int,    'd', "sandbox-dimension",
                                                              "Dimension of sandboxes for simulating DIM+1-dim. sandpiles.",
                                                              "DIM",
                                                              true),
                                  Argument(Type::Int,    'l', "lower-length",
                                                              "Minimal sandbox size to simulate (edge length of hypercube).",
                                                              "LEN",
                                                              true,
                                                              {"upper-length"}),
                                  Argument(Type::Int,    'L', "upper-length",
                                                              "Maximal sandbox size to simulate (edge length of hypercube).",
                                                              "LEN",
                                                              true,
                                                              {"lower-length"}),
                                  Argument(Type::Int,    'i', "length-increment",
                                                              "Step size of sandbox sizes.",
                                                              "DELTA-LEN",
                                                              true),
                                  Argument(Type::String, 'c', "boundary-conditions",
                                                              "Sandbox boundary conditions, c=closed, o=open "
                                                              "(default is all closed). Example: "
                                                              "2-dim. sandbox, all edges open except 2nd dim. upper edge: oooc.",
                                                              "LULU...",
                                                              true),
                                  Argument(Type::Int,    'n', "drives",
                                                              "Number of drives for every single simulation run.",
                                                              "DRIVES",
                                                              true),
                                  Argument(Type::Int,    'N', "num-samples",
                                                              "Number of repeated simulation runs for every lattice size.",
                                                              "SAMPLES",
                                                              true),
                                  Argument(Type::Int,    't', "num-threads",
                                                              "Number of parallel threads for i) main simulation and ii) random "
                                                              "number generation for bootstrap sampling in moment analysis.",
                                                              "THREADS",
                                                              true),
                                  Argument(Type::None,   'C', "make-critical",
                                                              "Drive the sandboxes to state of SOC prior to simulation.",
                                                              "",
                                                              true),
                                  Argument(Type::String, 'Y', "criticalization-parameters",
                                                              "Parameters for automatic sandbox \"criticalization\", which "
                                                              "determine \"criticality parameter\" saturation detection behavior. "
                                                              "Default is '500000,1000*LEN^DIM,10000,1e-6'.",
                                                              "MIN_NUM_DRIVES,MAX_NUM_DRIVES,DRIVES_PER_BUNCH,EPSILON",
                                                              true,
                                                              {"make-critical"}),
                                  Argument(Type::None,   'x', "skip-simulation",
                                                              "Skip the simulation and moment analysis. "
                                                              "Use together with \"make-critical\".",
                                                              "",
                                                              true,
                                                              {"make-critical"}),
                                  Argument(Type::None,   'a', "do-moment-analysis",
                                                              "Do moment analysis after all simulations have finished.",
                                                              "",
                                                              true),
                                  Argument(Type::Int,    'B', "bootstrap-samples",
                                                              "Number of (outer) bootstrap samples used in moment analysis.",
                                                              "SAMPLES",
                                                              true,
                                                              {"do-moment-analysis", "inner-bootstraps"}),
                                  Argument(Type::Int,    'b', "inner-bootstraps",
                                                              "Number of inner bootstrap samples used for double bootstrap.",
                                                              "SAMPLES",
                                                              true,
                                                              {"do-moment-analysis", "bootstrap-samples"}),
                                  Argument(Type::String, 'O', "load-sandboxes",
                                                              "Load the sandboxes before 'criticalization' / simulation. "
                                                              "Expected filename postfix format: \"_size_LENxLENx..."
                                                              "_model_MODEL-ID_params_PAR1=NUM1,PAR2=NUM2,....sbx\".",
                                                              "FILENAME-PREFIX",
                                                              true),
                                  Argument(Type::String, 'S', "save-sandboxes",
                                                              "Save the sandboxes after 'criticalization' / simulation. "
                                                              "File name postfix format: \"_size_LENxLENx..."
                                                              "_model_MODEL-ID_params_PAR1=NUM1,PAR2=NUM2,....sbx\".",
                                                              "FILENAME-PREFIX",
                                                              true),
                                  Argument(Type::String, 's', "save-results",
                                                              "Save the moment samples from simulation and the final results "
                                                              "from moment analysis.",
                                                              "FILENAME",
                                                              true),
                                  Argument(Type::String, 'v', "log-level",
                                                              "Global log level (verbosity) of the log output.",
                                                              "NONE|CRIT|ERROR|WARNG|LESS|INFO|MORE|VERB|DEBUG|DDBUG",
                                                              true),
                                  Argument(Type::String, 'f', "log-file",
                                                              "Write the program output to a log file.",
                                                              "FILENAME",
                                                              true),
                                  Argument(Type::None,   'G', "enable-plotting",
                                                              "Enable visualization of sandpiles during all simulations. "
                                                              "Requires a maximum threads count of 1 (\"--num-threads=1\").",
                                                              "",
                                                              true),
                                  Argument(Type::Int,    'r', "seed",
                                                              "Seed for main random generator used to seed all random generators. "
                                                              "If the seed is not set or a negative seed is provided, "
                                                              "a random device is used for seeding.",
                                                              "UNSIGNED_SEED",
                                                              true)},
                                  true);

        if (!argParser.parseArgs())
        {
            std::cout<<"\n";

            if (argParser.helpRequested())
                argParser.printHelp();
            else
                argParser.printUsage();
            return EXIT_FAILURE;
        }
        else if (argParser.helpRequested())
        {
            argParser.printHelp();
            return EXIT_SUCCESS;
        }

        //Set variables

        modelString = argParser.getParamValue("simulation-model", modelString);

        modelParams = argParser.getParamValue("model-parameters", modelParams);

        listParams = argParser.getParamCount("list-parameters", listParams ? 1 : 0) > 0;

        sandboxDim = static_cast<short>(std::stoi(argParser.getParamValue("sandbox-dimension", std::to_string(sandboxDim))));

        lowerLength = static_cast<short>(std::stoi(argParser.getParamValue("lower-length", std::to_string(lowerLength))));
        upperLength = static_cast<short>(std::stoi(argParser.getParamValue("upper-length", std::to_string(upperLength))));
        lengthIncrement = static_cast<short>(std::stoi(argParser.getParamValue("length-increment",
                                                                               std::to_string(lengthIncrement))));

        std::string defaultBoundaryConditions = "";
        for (short i = 0; i < sandboxDim; ++i)
            defaultBoundaryConditions.append("cc");
        boundaryConditions = argParser.getParamValue("boundary-conditions", defaultBoundaryConditions);

        numDrives = std::stoi(argParser.getParamValue("drives", std::to_string(numDrives)));

        numSamples = std::stoi(argParser.getParamValue("num-samples", std::to_string(numSamples)));

        numThreads = std::stoi(argParser.getParamValue("num-threads", std::to_string(numThreads)));

        makeCritical = argParser.getParamCount("make-critical", makeCritical ? 1 : 0) > 0;
        skipSim = argParser.getParamCount("skip-simulation", skipSim ? 1 : 0) > 0;

        criticalizationParams = argParser.getParamValue("criticalization-parameters", criticalizationParams);

        doMomentAnalysis = argParser.getParamCount("do-moment-analysis", doMomentAnalysis ? 1 : 0) > 0;

        numBootstraps = std::stoi(argParser.getParamValue("bootstrap-samples", std::to_string(50*numSamples)));
        numInnerBootstraps = std::stoi(argParser.getParamValue("inner-bootstraps", std::to_string(5*numBootstraps)));

        loadSandboxPrefix = argParser.getParamValue("load-sandboxes", loadSandboxPrefix);
        saveSandboxPrefix = argParser.getParamValue("save-sandboxes", saveSandboxPrefix);

        saveResultsFileName = argParser.getParamValue("save-results", saveResultsFileName);

        logLevel = Logger::labelToLogLevel(argParser.getParamValue("log-level", Logger::logLevelToLabel(logLevel)));

        logFileName = argParser.getParamValue("log-file", logFileName);

        enablePlotting = argParser.getParamCount("enable-plotting", enablePlotting ? 1 : 0) > 0;

        rndSeed = std::stoll(argParser.getParamValue("seed", "-1"));    //Do not re-seed, if no (or negative) seed provided

        //Check variables

        if (Aux::getModelIdFromString(modelString) == Aux::Model::InvalidModel)
            throw std::invalid_argument("The specified simulation model is not available!");

        if (sandboxDim <= 0)
            throw std::invalid_argument("Sandbox dimension must be positive!");

        if (lowerLength <= 0 || upperLength <= 0)
            throw std::invalid_argument("Sandbox edge lengths must be positive!");

        if (upperLength < lowerLength)
            throw std::invalid_argument("Sandbox upper edge length must not be smaller than lower edge length!");

        if (lengthIncrement <= 0)
            throw std::invalid_argument("Sandbox edge length increment must be positive!");

        if (boundaryConditions == "")
            boundaryConditions = defaultBoundaryConditions;

        if (boundaryConditions.size() != static_cast<size_t>(2*sandboxDim))
            throw std::invalid_argument("Wrong number of boundary conditions!");

        for (auto it = boundaryConditions.begin(); it != boundaryConditions.end(); ++it)
            if (*it != 'c' && *it != 'o')
                throw std::invalid_argument("Unrecognized boundary condition specifier!");

        if (numDrives < 0)
            throw std::invalid_argument("Number of drives must be positive!");

        if (numSamples <= 0)
            throw std::invalid_argument("Number of samples must be positive!");

        if (numThreads <= 0)
            throw std::invalid_argument("Number of threads must be positive!");

        if (numBootstraps <= 0)
            throw std::invalid_argument("Number of bootstrap samples must be positive!");

        if (numInnerBootstraps <= 0)
            throw std::invalid_argument("Number of inner bootstrap samples must be positive!");

        if (saveResultsFileName != "")
        {
            std::ofstream stream;
            stream.exceptions(std::ios_base::badbit | std::ios_base::failbit);
            stream.open(saveResultsFileName);
            stream.close();
        }

        if (logFileName != "")
        {
            std::ofstream stream;
            stream.exceptions(std::ios_base::badbit | std::ios_base::failbit);
            stream.open(logFileName);
            stream.close();
        }

        if (enablePlotting && numThreads > 1)
            throw std::invalid_argument("Plotting only works with maximally one thread (\"--num-threads=1\")!");
    }
    catch (const std::invalid_argument& exc)
    {
        std::cerr<<exc.what()<<std::endl;
        return EXIT_FAILURE;
    }
    catch (const std::ios_base::failure& exc)
    {
        std::cerr<<"ERROR: Could not create file \""<<saveResultsFileName<<"\"!\n";
        std::cerr<<exc.what()<<std::endl;
        return EXIT_FAILURE;
    }

    //Re-seed main pRNG with explicitly provided seed
    if (rndSeed >= 0)
        mainRndGenerator.seed(rndSeed);

    std::shared_ptr<Logger> logger = std::make_shared<Logger>(logLevel, logFileName);

    Aux::Model modelId = Aux::getModelIdFromString(modelString);

    if (listParams)
    {
        std::cout<<"Available model parameters for model \""<<modelString<<"\":\n\n";

        std::seed_seq tmpSeedSeq {};
        auto tParams = Simulator(modelId, {1}, tmpSeedSeq, logger).getModel()->listModelParameters();

        for (const auto& it : tParams)
            std::cout<<it.first<<":\n\t"<<it.second<<"\n\n";

        if (tParams.empty())
            std::cout<<"NONE\n";

        return EXIT_SUCCESS;
    }

    std::vector<std::pair<short, Sandbox>> critSandboxes;
    std::map<short, std::string> loadSandboxFileNames, saveSandboxFileNames;

    //Check if sandbox files can be read/written before simulation starts
    if (loadSandboxPrefix != "" || saveSandboxPrefix != "")
    {
        for (short len = lowerLength; len <= upperLength; len += lengthIncrement)
        {
            std::string tFileSuffix = "_size_";
            for (short i = 0; i < sandboxDim; ++i)
            {
                tFileSuffix.append(std::to_string(len));
                if (i < sandboxDim - 1)
                    tFileSuffix.append("x");
            }
            tFileSuffix.append("_model_").append(modelString);
            tFileSuffix.append("_params_").append(modelParams);
            tFileSuffix.append(".sbx");

            std::string tLoadFileName = loadSandboxPrefix + tFileSuffix;
            std::string tSaveFileName = saveSandboxPrefix + tFileSuffix;

            loadSandboxFileNames.insert({len, tLoadFileName});
            saveSandboxFileNames.insert({len, tSaveFileName});

            if (loadSandboxPrefix != "")
            {
                try
                {
                    std::ifstream stream;
                    stream.exceptions(std::ios_base::badbit | std::ios_base::failbit);
                    stream.open(tLoadFileName);
                    stream.close();
                }
                catch (const std::ios_base::failure& exc)
                {
                    logger->logCritical("ERROR: Could not open sandbox file \"" + tLoadFileName + "\"!");
                    std::cerr<<exc.what()<<std::endl;
                    return EXIT_FAILURE;
                }
            }
            if (saveSandboxPrefix != "" && saveSandboxPrefix != loadSandboxPrefix)  //Do not overwrite input file by new output file
            {
                try
                {
                    std::ofstream stream;
                    stream.exceptions(std::ios_base::badbit | std::ios_base::failbit);
                    stream.open(tSaveFileName);
                    stream.close();
                }
                catch (const std::ios_base::failure& exc)
                {
                    logger->logCritical("ERROR: Could not create sandbox file \"" + tSaveFileName + "\"!");
                    std::cerr<<exc.what()<<std::endl;
                    return EXIT_FAILURE;
                }
            }
        }
    }

    for (short len = lowerLength; len <= upperLength; len += lengthIncrement)
    {
        //Create sandbox; use empty seed sequence since we need to set different seed for each sample later either way
        std::seed_seq tmpSeedSeq {};
        Sandbox tSandbox(std::vector<short>(sandboxDim, len), tmpSeedSeq);

        //Load (critical) sandbox
        if (loadSandboxPrefix != "")
        {
            if (!tSandbox.load(loadSandboxFileNames.at(len)))
            {
                logger->logCritical("ERROR: Could not load sandbox file \"" + loadSandboxFileNames.at(len) + "\"!");
                return EXIT_FAILURE;
            }
        }

        std::vector<std::pair<bool, bool>> tBounds;
        for (auto it = boundaryConditions.begin(); it != boundaryConditions.end(); ++it)
            tBounds.push_back({*it == 'o' ? true : false, *(++it) == 'o' ? true : false});

        tSandbox.setBoundaryConditions(tBounds);

        critSandboxes.push_back(std::make_pair(len, std::move(tSandbox)));
    }

    //Get simulation model parameters from command line argument

    std::vector<std::pair<std::string, int>> modelParameters {};

    std::istringstream tIsstream(modelParams);

    for (std::string paramString; std::getline(tIsstream, paramString, ',');)
    {
        std::vector<std::string> tSubstrs;
        std::istringstream tIsstream2(paramString);

        for (std::string paramSubstring; std::getline(tIsstream2, paramSubstring, '=');)
            tSubstrs.push_back(paramSubstring);

        if (tSubstrs.size() != 2 || tSubstrs.at(0) == "")
        {
            logger->logCritical("ERROR: Could not parse simulation model parameters! \"" + paramString +
                                "\" does not match format \"PAR=NUM\".");
            return EXIT_FAILURE;
        }
        else
        {
            std::string paramName = tSubstrs.at(0);
            std::string paramValueStr = tSubstrs.at(1);
            int paramValue = 0;

            try
            {
                paramValue = std::stoi(paramValueStr);
            }
            catch (const std::exception& exc)
            {
                logger->logCritical("ERROR: Could convert simulation model parameter \"" + paramName + "\"! "
                                    "Integer number required.");
                std::cerr<<exc.what()<<std::endl;
                return EXIT_FAILURE;
            }

            modelParameters.push_back({paramName, paramValue});
        }
    }

    //Get criticalization ("c13n") parameters from command line argument

    long long c13nMinNumDrives = 50*10000;
    long long c13nMaxNumDrives = 1000*std::pow(upperLength, sandboxDim);
    long long c13nDrivesPerBunch = 10000;
    double c13nEpsilon = 1e-6;

    if (criticalizationParams != "")
    {
        std::vector<std::string> c13nParams;

        std::istringstream tIsstream2(criticalizationParams);

        for (std::string paramString; std::getline(tIsstream2, paramString, ',');)
            c13nParams.push_back(paramString);

        if (c13nParams.size() != 4)
        {
            logger->logCritical("ERROR: Could not parse 'criticalization' parameters! \"" + criticalizationParams +
                                "\" does not match format \"MIN_NUM_DRIVES,MAX_NUM_DRIVES,DRIVES_PER_BUNCH,EPSILON\".");
            return EXIT_FAILURE;
        }

        try
        {
            c13nMinNumDrives = std::stoll(c13nParams.at(0));
            c13nMaxNumDrives = std::stoll(c13nParams.at(1));
            c13nDrivesPerBunch = std::stoll(c13nParams.at(2));
            c13nEpsilon = std::stod(c13nParams.at(3));
        }
        catch (const std::exception& exc)
        {
            logger->logCritical("ERROR: Could not convert 'criticalization' parameters! \"" + criticalizationParams +
                                "\" types do not match \"LONGLONG,LONGLONG,LONGLONG,DOUBLE\".");
            std::cerr<<exc.what()<<std::endl;
            return EXIT_FAILURE;
        }
    }

    SimulationManager simulationManager(logger);

    //Measure elapsed times during "sandbox criticalization", simulation and moment analysis
    std::chrono::high_resolution_clock::time_point timeBeginCrit, timeEndCrit,
                                                   timeBeginSim, timeEndSim,
                                                   timeBeginAna, timeEndAna;

    timeBeginCrit = std::chrono::high_resolution_clock::now();

    //Drive sandboxes into state of SOC
    if (makeCritical)
    {
        logger->logLess("Starting \"criticalization\" of sandboxes.");

        std::seed_seq tSeedSeqC13n = Aux::generateSeedSeq(mainRndGenerator);

        simulationManager.makeSandboxesCritical(critSandboxes, sandboxDim, c13nMinNumDrives, c13nMaxNumDrives, c13nDrivesPerBunch,
                                                c13nEpsilon, modelId, tSeedSeqC13n, modelParameters, enablePlotting);

        logger->logLess("Finished \"criticalization\" of sandboxes.");
        logger->logInfo("+--------------------------------------------");
        logger->logInfo("");
    }
    timeEndCrit = std::chrono::high_resolution_clock::now();

    if (!skipSim)
    {
        logger->logLess("Starting sandpile simulations.");
        logger->logInfo("+--------------------------------------------");

        timeBeginSim = std::chrono::high_resolution_clock::now();

        std::seed_seq tSeedSeqSim = Aux::generateSeedSeq(mainRndGenerator);

        std::map<short, std::vector<AvalancheStatistics::Moments>> samplesPerSize_moments1, samplesPerSize_moments2;

        //Actual simulation(s) for collecting avalanche statistics
        simulationManager.runFullSimulation(critSandboxes, sandboxDim, modelId, tSeedSeqSim,
                                            samplesPerSize_moments1, samplesPerSize_moments2,
                                            numDrives, numSamples, true, true, modelParameters, numThreads, enablePlotting);

        timeEndSim = std::chrono::high_resolution_clock::now();

        logger->logLess("Finished sandpile simulations.");
        logger->logInfo("+--------------------------------------------");

        struct MomentAnalysis::Result sizeResult, linSizeResult, durationResult, areaResult;

        std::array<std::vector<double>, 10> rawFitResultsSize, rawFitResultsLinSize, rawFitResultsDuration, rawFitResultsArea;

        //Do moment analysis
        if (doMomentAnalysis)
        {
            logger->logLess("Starting moment analysis.");
            logger->logInfo("");
            logger->logInfo("Analysis configuration:");
            logger->logInfo("\tBootstrap sample length: " + std::to_string(numSamples));
            logger->logInfo("\tNumber of outer samples: " + std::to_string(numBootstraps));
            logger->logInfo("\tNumber of inner samples: " + std::to_string(numInnerBootstraps));
            logger->logInfo("");
            logger->logInfo("+--------------------------------------------");

            //Do moment analysis using numSimulations samples of the distributions' moments

            std::map<short, std::vector<double>> samplesPerSize_size_moment1, samplesPerSize_size_moment2;
            std::map<short, std::vector<double>> samplesPerSize_linSize_moment1, samplesPerSize_linSize_moment2;
            std::map<short, std::vector<double>> samplesPerSize_area_moment1, samplesPerSize_area_moment2;
            std::map<short, std::vector<double>> samplesPerSize_duration_moment1, samplesPerSize_duration_moment2;

            MomentAnalysis::splitMoments(samplesPerSize_moments1, samplesPerSize_size_moment1, samplesPerSize_linSize_moment1,
                                                                  samplesPerSize_area_moment1, samplesPerSize_duration_moment1);
            MomentAnalysis::splitMoments(samplesPerSize_moments2, samplesPerSize_size_moment2, samplesPerSize_linSize_moment2,
                                                                  samplesPerSize_area_moment2, samplesPerSize_duration_moment2);

            //Seed sequences for bootstrap sampling in moment analysis
            std::seed_seq tSeedSeqMA1 = Aux::generateSeedSeq(mainRndGenerator);
            std::seed_seq tSeedSeqMA2 = Aux::generateSeedSeq(mainRndGenerator);
            std::seed_seq tSeedSeqMA3 = Aux::generateSeedSeq(mainRndGenerator);
            std::seed_seq tSeedSeqMA4 = Aux::generateSeedSeq(mainRndGenerator);
            std::seed_seq tSeedSeqMA5 = Aux::generateSeedSeq(mainRndGenerator);
            std::seed_seq tSeedSeqMA6 = Aux::generateSeedSeq(mainRndGenerator);
            std::seed_seq tSeedSeqMA7 = Aux::generateSeedSeq(mainRndGenerator);
            std::seed_seq tSeedSeqMA8 = Aux::generateSeedSeq(mainRndGenerator);

            timeBeginAna = std::chrono::high_resolution_clock::now();

            logger->logLess("Running moment analysis for avalanche size distribution:");
            logger->logInfo("");

            sizeResult = MomentAnalysis::doMomentAnalysis(samplesPerSize_size_moment1, samplesPerSize_size_moment2,
                                                          rawFitResultsSize, numBootstraps, numInnerBootstraps,
                                                          tSeedSeqMA1, tSeedSeqMA2, logger, numThreads);

            logger->logInfo("+--------------------------------------------");
            logger->logLess("Running moment analysis for avalanche linear size distribution:");
            logger->logInfo("");

            linSizeResult = MomentAnalysis::doMomentAnalysis(samplesPerSize_linSize_moment1, samplesPerSize_linSize_moment2,
                                                             rawFitResultsLinSize, numBootstraps, numInnerBootstraps,
                                                             tSeedSeqMA3, tSeedSeqMA4, logger, numThreads);

            logger->logInfo("+--------------------------------------------");
            logger->logLess("Running moment analysis for avalanche duration distribution:");
            logger->logInfo("");

            durationResult = MomentAnalysis::doMomentAnalysis(samplesPerSize_duration_moment1, samplesPerSize_duration_moment2,
                                                              rawFitResultsDuration, numBootstraps, numInnerBootstraps,
                                                              tSeedSeqMA5, tSeedSeqMA6, logger, numThreads);

            logger->logInfo("+--------------------------------------------");
            logger->logLess("Running moment analysis for avalanche area distribution:");
            logger->logInfo("");

            areaResult = MomentAnalysis::doMomentAnalysis(samplesPerSize_area_moment1, samplesPerSize_area_moment2,
                                                          rawFitResultsArea, numBootstraps, numInnerBootstraps,
                                                          tSeedSeqMA7, tSeedSeqMA8, logger, numThreads);

            timeEndAna = std::chrono::high_resolution_clock::now();

            logger->logInfo("");
            logger->logInfo("+--------------------------------------------");
            logger->logLess("Finished moment analysis.");
            logger->logInfo("+--------------------------------------------");
            logger->logLess("");
            logger->logLess("Moment analysis results (critical exponents and observable dimensions):");
            logger->logLess("");
            logger->logLess("Size:");
            logger->logLess("\ttau:\t" + std::to_string(   sizeResult.critExp_mean) + " +- " +
                                           std::to_string(   sizeResult.critExp_std));
            logger->logLess("\tD:\t" + std::to_string(     sizeResult.critDim_mean) + " +- " +
                                         std::to_string(     sizeResult.critDim_std));
            logger->logLess("");
            logger->logLess("Linear size:");
            logger->logLess("\tlambda:\t" + std::to_string(linSizeResult.critExp_mean) + " +- " +
                                              std::to_string(linSizeResult.critExp_std));
            logger->logLess("\tM:\t" + std::to_string(     linSizeResult.critDim_mean) + " +- " +
                                         std::to_string(     linSizeResult.critDim_std));
            logger->logLess("");
            logger->logLess("Duration:");
            logger->logLess("\talpha:\t" + std::to_string( durationResult.critExp_mean) + " +- " +
                                             std::to_string( durationResult.critExp_std));
            logger->logLess("\tZ:\t" + std::to_string(     durationResult.critDim_mean) + " +- " +
                                         std::to_string(     durationResult.critDim_std));
            logger->logLess("");
            logger->logLess("Area:");
            logger->logLess("\tkappa:\t" + std::to_string( areaResult.critExp_mean) + " +- " +
                                             std::to_string( areaResult.critExp_std));
            logger->logLess("\tT:\t" + std::to_string(     areaResult.critDim_mean) + " +- " +
                                         std::to_string(     areaResult.critDim_std));
            logger->logInfo("");
            logger->logInfo("+--------------------------------------------");
        }

        //Save simulation and analysis results
        if (saveResultsFileName != "")
        {
            logger->logLess("Writing simulation/analysis results to file \"" + saveResultsFileName + "\"...");

            try
            {
                //Open results file
                std::ofstream stream;
                stream.exceptions(std::ios_base::badbit | std::ios_base::failbit);
                stream.open(saveResultsFileName);

                stream<<std::fixed<<std::setprecision(8);

                //Write simulation information

                stream<<"Results of simulation with command line \"";
                for (int i = 0; i < argc; ++i)
                {
                    if (i > 0)
                        stream<<' ';
                    stream<<argv[i];
                }
                stream<<"\"\n";

                //Write moment analysis results
                if (doMomentAnalysis)
                {
                    stream<<"\n\n- Moment analysis results (critical exponents 'K' and observable dimensions 'rho'):\n\n";

                    stream<<"Size:\n\ttau:\t"<<sizeResult.critExp_mean<<" +- "<<sizeResult.critExp_std<<"\n\tD:\t"
                                           <<sizeResult.critDim_mean<<" +- "<<sizeResult.critDim_std<<'\n';
                    stream<<"Linear size:\n\tlambda:\t"<<linSizeResult.critExp_mean<<" +- "<<linSizeResult.critExp_std<<"\n\tM:\t"
                                                     <<linSizeResult.critDim_mean<<" +- "<<linSizeResult.critDim_std<<'\n';
                    stream<<"Duration:\n\talpha:\t"<<durationResult.critExp_mean<<" +- "<<durationResult.critExp_std<<"\n\tZ:\t"
                                                 <<durationResult.critDim_mean<<" +- "<<durationResult.critDim_std<<'\n';
                    stream<<"Area:\n\tkappa:\t"<<areaResult.critExp_mean<<" +- "<<areaResult.critExp_std<<"\n\tT:\t"
                                             <<areaResult.critDim_mean<<" +- "<<areaResult.critDim_std<<'\n';
                }

                //Write raw moment samples from simulation

                stream<<"\n\n- Raw moment samples:\n\n";
                stream<<"Lattice Size,FirstMom-size,FirstMom-linSize,FirstMom-duration,FirstMom-area,"
                                     "SecondMom-size,SecondMom-linSize,SecondMom-duration,SecondMom-area\n";

                for (const auto& it : samplesPerSize_moments1)
                {
                    decltype(it.second)::size_type i = 0;
                    for (const AvalancheStatistics::Moments& moments1 : it.second)
                    {
                        const AvalancheStatistics::Moments& moments2 = samplesPerSize_moments2.at(it.first).at(i);
                        stream<<it.first<<','<<moments1.size<<','<<moments1.linSize<<','<<moments1.duration<<','<<moments1.area<<','
                                             <<moments2.size<<','<<moments2.linSize<<','<<moments2.duration<<','<<moments2.area
                                             <<'\n';
                    }
                }

                //Write bootstrapped moment fit results

                typedef std::array<std::vector<double>, 10> RawFitResults;
                typedef std::reference_wrapper<RawFitResults> RawFitResultsRef;
                typedef std::array<std::pair<std::string, RawFitResultsRef>, 4> RawFitResultsRefs;

                //Iteratively write bootstrapped fit results for each observable (size/linSize/duration/area)
                for (const auto& it : RawFitResultsRefs{{{"size", std::ref(rawFitResultsSize)},
                                                         {"linSize", std::ref(rawFitResultsLinSize)},
                                                         {"duration", std::ref(rawFitResultsDuration)},
                                                         {"area", std::ref(rawFitResultsArea)}}})
                {
                    stream<<"\n\nBootstrapped moment fit results ("<<it.first<<"):\n\n";
                    stream<<"FirstMom-slope,FirstMom-slopeStd,FirstMom-offset,FirstMom-offsetStd,FirstMom-slopeOffsetCov,"
                            "SecondMom-slope,SecondMom-slopeStd,SecondMom-offset,SecondMom-offsetStd,SecondMom-slopeOffsetCov\n";

                    const RawFitResults& fitResults = it.second;

                    //Write fit results for current observable for each (outer) bootstrap
                    for (std::remove_reference<decltype(fitResults[0])>::type::size_type i = 0; i < fitResults[0].size(); ++i)
                    {
                        stream<<fitResults[0][i]<<','<<fitResults[1][i]<<','<<fitResults[2][i]<<','<<fitResults[3][i]<<','
                              <<fitResults[4][i]<<','<<fitResults[5][i]<<','<<fitResults[6][i]<<','<<fitResults[7][i]<<','
                              <<fitResults[8][i]<<','<<fitResults[9][i]<<'\n';
                    }
                }

                logger->logLess("Written results to file \"" + saveResultsFileName + "\".");
                logger->logInfo("+--------------------------------------------");
            }
            catch (const std::ios_base::failure&)
            {
                logger->logError("ERROR: Could not write results to file \"" + saveResultsFileName + "\"!");
            }
        }
    }

    //Save sandboxes
    if (saveSandboxPrefix != "")
    {
        for (auto& it : critSandboxes)
            if (!it.second.save(saveSandboxFileNames.at(it.first)))
                logger->logError("ERROR: Could not write sandbox file \"" + saveSandboxFileNames.at(it.first) + "\"!");
    }

    //Print timing information

    logger->logInfo("");
    logger->logInfo("Summary:");
    logger->logInfo("");
    if (makeCritical)
    {
        logger->logInfo("\tSandbox \"criticalization\" took " + std::to_string(
                            std::chrono::duration_cast<std::chrono::seconds>(timeEndCrit-timeBeginCrit).count() / 60.) + " min.");
    }
    if (!skipSim)
    {
        logger->logInfo("\tSimulation took " + std::to_string(
                            std::chrono::duration_cast<std::chrono::seconds>(timeEndSim-timeBeginSim).count() / 60.) + " min.");

        if (doMomentAnalysis)
        {
            logger->logInfo("\tMoment analysis took " + std::to_string(
                                std::chrono::duration_cast<std::chrono::seconds>(timeEndAna-timeBeginAna).count() / 60.) + " min.");
        }
    }
    logger->logInfo("+--------------------------------------------");

    return EXIT_SUCCESS;
}

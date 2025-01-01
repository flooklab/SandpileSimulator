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

#ifndef SANDSIM_SIMULATOR_H
#define SANDSIM_SIMULATOR_H

#cmakedefine SANDSIM_ENABLE_PLOTTER @SANDSIM_ENABLE_PLOTTER@

#include "avalanchestatistics.h"
#include "logger.h"
#include "noncopyable.h"
#include "sandbox.h"
#include "simulationmodel.h"

#ifdef SANDSIM_ENABLE_PLOTTER
#include "plotter.h"
#endif // SANDSIM_ENABLE_PLOTTER

#include <memory>
#include <random>
#include <vector>

/*!
 * \brief Simulation class that manages the sandpile simulation.
 *
 * Instead of writing your own simulation loop you can also use this simple class,
 * which does the work for you. A variable number of consecutive, alternating driving
 * and relaxation steps can be simply triggered by runSimulation() and the results
 * can be obtained after that directly from the sandbox and statistics objects.
 *
 * You can maintain the full control over sandbox/statistics/model if you directly
 * provide respective shared pointers to the constructor in the beginning and/or
 * obtain shared pointers using the corresponding member functions of this class,
 * even for those sandbox/statistics/model automatically created by the constructor
 * Simulator(Aux::Model, std::vector<short>, std::seed_seq&, std::shared_ptr<Logger>);
 *
 * The sandpile can be plotted (especially during the simulation) in a dedicated window
 * in form of a heat map, which you can toggle with enablePlotting() / disablePlotting().
 */
class Simulator : private NonCopyable
{
public:
    Simulator(Aux::Model pModelId, std::vector<short> pShape, std::seed_seq& pSandboxSeedSeq,
              std::shared_ptr<Logger> pLogger = nullptr);               ///< \brief Construct a sandpile simulation for a certain
                                                                        ///  type of simulation model and a sandbox of given shape.
    Simulator(std::shared_ptr<Sandbox> pSandbox, std::shared_ptr<AvalancheStatistics> pStatistics,
              std::shared_ptr<SimulationModel> pSimulationModel,
              std::shared_ptr<Logger> pLogger = nullptr);               ///< \brief Construct a sandpile simulation from existing
                                                                        ///  sandbox, statistics and simulation model.
    //
    static std::shared_ptr<SimulationModel> createSimulationModel(Aux::Model pModelId, Sandbox& pSandbox);
                                                                                            ///< \brief Dynamically create a specific
                                                                                            ///  SimulationModel by its model ID.
    //
    std::shared_ptr<Sandbox> getSandbox();                  ///< Get the sandbox on which the simulation is done.
    std::shared_ptr<SimulationModel> getModel();            ///< Get the model used by simulator to actually simulate the sandpile.
    std::shared_ptr<AvalancheStatistics> getStatistics();   ///< Get the simulation statistics about avalanches on the sandpile.
    std::shared_ptr<Logger> getLogger();                    ///< Get the logger used for logging.
    //
    void enablePlotting(unsigned int pFramerateLimit = 30); ///< Enable live plotting of the sandbox content.
    void disablePlotting();                                 ///< Disable live plotting of the sandbox content.
    //
    void runSimulation(long long pNumDrives);               ///< Start the simulation.

private:
    std::shared_ptr<Sandbox> sandbox;               //Sandbox containing the sandpile to be simulated
    std::shared_ptr<AvalancheStatistics> stats;     //Statistics collector capturing the statistics of avalanches
    std::shared_ptr<SimulationModel> model;         //Model according to which the simulation behaves
    std::shared_ptr<Logger> logger;                 //For logging the progress and stuff
#ifdef SANDSIM_ENABLE_PLOTTER
    Plotter plotter;                                //Plots the sandpile
#endif // SANDSIM_ENABLE_PLOTTER
    //
    bool plottingEnabled;                           //Plot the sandpile during the simulation using the plotter?
};

#endif // SANDSIM_SIMULATOR_H

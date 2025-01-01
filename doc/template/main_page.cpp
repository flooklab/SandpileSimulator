////////////////////////////////////////////////////////////
/// \mainpage
///
/// Sandpile %Simulator is a cellular automaton (library) for simulating and analyzing the dynamics of
/// sandpiles with special focus on self-organized criticality (SOC) and avalanche scaling behavior.
/// Sandpiles with arbitrary dimensions can be simulated using different models and from the
/// collected statistics scaling exponents can be extracted using a "moment analysis" approach.
///
/// Copyright (C) 2021, 2025 M. Frohne
///
/// Sandpile %Simulator is free software: you can redistribute it and/or modify
/// it under the terms of the GNU Affero General Public License as published
/// by the Free Software Foundation, either version 3 of the License,
/// or (at your option) any later version.
///
/// Sandpile %Simulator is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
/// See the GNU Affero General Public License for more details.
///
/// You should have received a copy of the GNU Affero General Public License
/// along with Sandpile %Simulator. If not, see <https://www.gnu.org/licenses/>.
///
/// \section example Example
/// Below you find a short example code for performing repeated simulations with the "FW model"
/// on 2-dim. sandboxes of different lattice sizes for extracting the \e critical \e exponent
/// and \e observable \e dimension for the avalanche observable \e linear \e size.
/// The sandboxes loaded from disk should already be in the state of SOC.
///
/// For a much more comprehensive "example" code please also look at this project's \p main.cpp.
///
/// \code
///
/// #include <iostream>
/// #include <array>
/// #include <vector>
/// #include <map>
/// #include <random>
///
/// #include "aux.h"
/// #include "sandbox.h"
/// #include "simulationmodel.h"
/// #include "avalanchestatistics.h"
/// #include "simulator.h"
/// #include "momentanalysis.h"
///
/// int main()
/// {
///     std::map<short, std::vector<AvalancheStatistics::Moments>> samplesPerSize_moments1;
///     std::map<short, std::vector<AvalancheStatistics::Moments>> samplesPerSize_moments2;
///
///     std::seed_seq seedSeq {};
///
///     //Run repeated simulations for different lattice sizes
///
///     for (short latticeSize = 10; latticeSize <= 70; latticeSize += 10)
///     {
///         Simulator simulator(Aux::Model::_FWM, {latticeSize, latticeSize}, seedSeq);
///
///         Sandbox& sandbox = *simulator.getSandbox();
///         SimulationModel& model = *simulator.getModel();
///         AvalancheStatistics& statistics = *simulator.getStatistics();
///
///         //Load a critical sandbox matching the current lattice size
///         sandbox.load("...criticalSandbox(latticeSize)....sbx");
///
///         sandbox.setBoundaryConditions({{false, true}, {false, true}});
///
///         model.setModelParameter("FWM_criticalSlope", 10);
///
///         model.enableRecordLinSize(true);
///
///         simulator.enablePlotting(30);
///
///         std::vector<AvalancheStatistics::Moments> sample_moments1, sample_moments2;
///
///         for (int i = 0; i < 20; ++i)
///         {
///             simulator.runSimulation(100000);
///
///             std::vector<AvalancheStatistics::Event> events;
///             statistics.detachEvents(events);
///
///             AvalancheStatistics::Moments firstMoments, secondMoments;
///             MomentAnalysis::calcMoments(events, firstMoments, secondMoments);
///
///             sample_moments1.push_back(firstMoments);
///             sample_moments2.push_back(secondMoments);
///         }
///
///         samplesPerSize_moments1.insert(std::make_pair(latticeSize, sample_moments1));
///         samplesPerSize_moments2.insert(std::make_pair(latticeSize, sample_moments2));
///     }
///
///     //Moment analysis
///
///     std::map<short, std::vector<double>> samplesPerSize_size_moment1, samplesPerSize_size_moment2;
///     std::map<short, std::vector<double>> samplesPerSize_linSize_moment1, samplesPerSize_linSize_moment2;
///     std::map<short, std::vector<double>> samplesPerSize_area_moment1, samplesPerSize_area_moment2;
///     std::map<short, std::vector<double>> samplesPerSize_duration_moment1, samplesPerSize_duration_moment2;
///
///     MomentAnalysis::splitMoments(samplesPerSize_moments1, samplesPerSize_size_moment1, samplesPerSize_linSize_moment1,
///                                                           samplesPerSize_area_moment1, samplesPerSize_duration_moment1);
///     MomentAnalysis::splitMoments(samplesPerSize_moments2, samplesPerSize_size_moment2, samplesPerSize_linSize_moment2,
///                                                           samplesPerSize_area_moment2, samplesPerSize_duration_moment2);
///
///     std::array<std::vector<double>, 10> rawFitResults;
///
///     struct MomentAnalysis::Result linSizeResult = MomentAnalysis::doMomentAnalysis(samplesPerSize_linSize_moment1,
///                                                                                    samplesPerSize_linSize_moment2,
///                                                                                    rawFitResults, 1000, 5000);
///
///     std::cout<<"Moment analysis results for linear size (critical exponents and observable dimensions):\n\n";
///
///     std::cout<<"\tlambda:\t"<<std::to_string(linSizeResult.critExp_mean)<<" +- "<<std::to_string(linSizeResult.critExp_std)<<"\n";
///     std::cout<<"\tM:\t"     <<std::to_string(linSizeResult.critDim_mean)<<" +- "<<std::to_string(linSizeResult.critDim_std)<<"\n";
///
///     return 0;
/// }
///
/// \endcode
////////////////////////////////////////////////////////////

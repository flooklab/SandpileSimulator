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

#ifndef SANDSIM_MOMENTANALYSIS_H
#define SANDSIM_MOMENTANALYSIS_H

#include "avalanchestatistics.h"
#include "logger.h"

#include <array>
#include <map>
#include <memory>
#include <random>
#include <tuple>
#include <utility>
#include <vector>

/*!
 * \brief Analyze power law scaling of an avalanche observable's dataset using moment analysis.
 *
 * For analysis of an observable Y whose distribution exhibits power law scaling and
 * is simulated on finite-size lattices a so-called "moment analysis" can be used.
 * We assume a so-called \e finite-size \e scaling of Y according to the following distribution function (see [1]):
 *
 * p(y) = a * y^{-rho} * G(y / (b * L^K))
 *
 * Here the power law scaling is given by the "critical exponent" 'rho' and the finite-size effects
 * (e.g. avalanche cut-off) are modeled by some function 'G' with the observable's value as argument
 * being scaled by a term dependent on the lattice size 'L' and the "observable dimension" 'K'.
 *
 * The most interesting quantities to extract here are the critical exponent 'rho' and the observable dimension 'K'.
 * As is shown in [1], this can be achieved by calculating the first and second \e moments of the measured/simulated
 * distribution p(y) for multiple different lattice sizes 'L' and then applying the actual "moment analysis".
 * This consists of doing straight line fits of the logarithmized moments, ln(moment[L_i]), vs. the logarithmized
 * lattice sizes, ln(L_i), for both the first and second moments. From the two resulting fit slopes one can
 * then calculate 'rho' and 'K' with simple equations (see calcCriticalExponent(), calcCriticalDimension()).
 * For more information see also doMomentAnalysis().
 *
 * The moments for the four different observables described by AvalancheStatistics::Event can be calculated
 * with calcMoments(). A moment analysis for a single observable that uses bootstrapping to determine the
 * uncertainties of 'rho' and 'K' can be done with doMomentAnalysis(). This requires a whole set of multiple
 * moment measurement samples for multiple lattice sizes for a \e single observable. Such sample sets
 * can be easily obtained from e.g. SimulationManager::runFullSimulation() for \e all four observables.
 * The function splitMoments() can be used to split these sets into individual sets for single observables.
 *
 * [1]: \e G. \e Pruessner, \e Self-Organized Criticality; \e Theory, \e Models \e and \e Characterisation \e (2012)
 */
class MomentAnalysis
{
private:
    typedef AvalancheStatistics::Event Event;
    typedef AvalancheStatistics::Moments Moments;

public:
    typedef std::vector<double> MomentSample;
    //
    /*!
     * \brief End result of a moment analysis for one observable.
     *
     * Contains the result in terms of means and uncertainties of the critical exponent and the observable dimension.
     */
    struct Result
    {
        double critExp_mean = 0;
        double critExp_std = 0;
        double critDim_mean = 0;
        double critDim_std = 0;
    };

public:
    MomentAnalysis() = delete;  ///< Deleted constructor.

public:
    static struct Result doMomentAnalysis(std::map<short, MomentSample> pSamplesPerSize_moment1,
                                          std::map<short, MomentSample> pSamplesPerSize_moment2,
                                          std::array<std::vector<double>, 10>& pRawFitResults,
                                          int pNumBootstraps, int pNumInnerBootstraps,
                                          std::seed_seq& pSeedSeq1, std::seed_seq& pSeedSeq2,
                                          std::shared_ptr<Logger> pLogger = nullptr, int pNumRNGThreads = 0);
                                                                                ///< \brief Perform a bootstrapped moment analysis
                                                                                ///  of moment samples for a single observable.
    //
    static void calcMoments(const std::vector<Event>& pEvents,
                            Moments& pFirstMoments, Moments& pSecondMoments);   ///< \brief Calculate moments of avalanche observable
                                                                                ///  distributions from (many) avalanche events.
    static void splitMoments(const std::map<short, std::vector<Moments>>& pSamplesPerSize,
                             std::map<short, std::vector<double>>& pSamplesPerSize_size,
                             std::map<short, std::vector<double>>& pSamplesPerSize_linSize,
                             std::map<short, std::vector<double>>& pSamplesPerSize_area,
                             std::map<short, std::vector<double>>& pSamplesPerSize_duration);   ///< \brief Split Moments samples into
                                                                                                /// samples for different observables.

private:
    static std::pair<double, double> calcMeanStdFromVector(const std::vector<double>& pVec,
                                                           const std::vector<double>& pUnc);    ///< \brief Calculate weighted
                                                                                                ///  sample mean and weighted
                                                                                                ///  sample standard deviation.
    //
    static double calcCriticalExponent(double pSlope_moment1, double pSlope_moment2);           ///< \brief Calculate critical
                                                                                                ///  exponent from moment fit slopes.
    static double calcCriticalExponentUncertainty(double pSlope_moment1, double pSlopeStd_moment1,
                                                  double pSlope_moment2, double pSlopeStd_moment2);
                                                                                            ///< \brief Calculate critical exponent
                                                                                            ///  uncertainty from moment fit slopes
                                                                                            ///  and uncertainties.
    static double calcCriticalDimension(double pSlope_moment1, double pSlope_moment2);  ///< \brief Calculate observable dimension
                                                                                        ///  from moment fit slopes.
    static double calcCriticalDimensionUncertainty(double pSlopeStd_moment1, double pSlopeStd_moment2);
                                                                                        ///< \brief Calculate observable dimension
                                                                                        ///  uncertainty from moment fit slopes
                                                                                        ///  and uncertainties.
    //
    static std::tuple<double, double, double, double, double> doFit(const std::vector<double>& pX, const std::vector<double>& pY,
                                                                    const std::vector<std::vector<double>>& pCovMat);
                                                                                        ///< \brief Perform a straight line fit
                                                                                        ///  for data with correlated y-values.
    //
    static void makeBootstrapSamples(const MomentSample& pOriginalSample,
                                     std::vector<MomentSample>& pBootstrapSamples,
                                     std::mt19937_64& pRndGenerator, int pNumThreads = 0);  ///< Generate bootstrap samples.
    static void generateBootstrappedFitSlopeSamplesFromMoments(std::map<short, MomentSample> pSamplesPerSize,
                                                               std::vector<double>& pFitSlopeSamples,
                                                               std::vector<double>& pFitSlopeStdSamples,
                                                               std::vector<double>& pFitOffsetSamples,
                                                               std::vector<double>& pFitOffsetStdSamples,
                                                               std::vector<double>& pFitSlopeOffsetCovSamples,
                                                               int pNumBootstraps, int pNumInnerBootstraps,
                                                               std::seed_seq& pSeedSeq,
                                                               std::shared_ptr<Logger> pLogger,
                                                               int pNumRNGThreads = 0); ///< \brief Generate/calculate straight line
                                                                                        /// fit slopes of bootstrapped, logarithmized
                                                                                        /// moment samples for a single observable.
};

#endif // SANDSIM_MOMENTANALYSIS_H

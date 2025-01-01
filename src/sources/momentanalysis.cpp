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

#include "momentanalysis.h"

//Public

/*!
 * \brief Perform a bootstrapped moment analysis of moment samples for a single observable.
 *
 * Determines the "critical exponent" and "observable dimension" of an avalanche observable.
 *
 * This is generally achieved by first calculating the first and second moments of the observable distribution
 * for different lattice sizes and for each moment doing a straight line fit of ln(moment[size]) vs. ln(size).
 * The critical exponent (and its uncertainty) then follow from the two resulting fit slopes via
 * calcCriticalExponent() (and calcCriticalExponentUncertainty()). Similarly, calcCriticalDimension()
 * (and calcCriticalDimensionUncertainty()) is used for obtaining the observable dimension.
 *
 * However, to obtain a more accurate estimate of the results this function does a bootstrapped version of
 * the above procedure. Therefore, sets of multiple measurements ("samples") of the first and second moments
 * must be provided for different lattice sizes via \p pSamplesPerSize_moment1 and \p pSamplesPerSize_moment2.
 * From every original sample a number of \p pNumBootstraps \e bootstrap \e samples will be derived and,
 * in turn, from every bootstrap sample \p pNumInnerBootstraps \e inner bootstrap samples will be derived.
 *
 * The inner bootstraps are used to determine the covariance matrix between moments of different lattice sizes
 * for every (outer) bootstrap. Using these covariance matrices, a fit as described before is then performed
 * for every (outer) bootstrap using the mean values of each sample. The bootstrapping and fitting is done
 * separately for the first and second moment via generateBootstrappedFitSlopeSamplesFromMoments().
 * This yields \p pNumBootstraps different fit slope values for the first and second moment each.
 *
 * From these \p pNumBootstraps pairs of slope values in turn \p pNumBootstraps values of
 * the critical exponent are derived as described before, resulting in a bootstrapped
 * critical exponent sample. Analogously bootstrapped samples for the critical exponent's
 * uncertainty, the observable dimension and the observable dimension's uncertainty are derived.
 *
 * To obtain the final result for the critical exponent, a weighted(!) bootstrap mean is then
 * calculated over the bootstrapped critical exponent sample (see calcMeanStdFromVector()).
 * The corresponding uncertainty sample is here only used for an inverse variance weighting of
 * the bootstrap mean and not to calculate the final uncertainty. That one is instead also obtained
 * from calcMeanStdFromVector() in terms of the unbiased, weighted sample standard deviation
 * (same weights as before). The same is done for the observavble dimension result as well.
 *
 * The results including the bootstrap-estimated uncertainties are returned as:
 *
 * {critExp, critExpStd, obsDim, obsDimStd}.
 *
 * The raw fit results from the moment fits for every (outer) bootstrap will be written
 * to \p pRawFitResults for both the first moment fits and the second moment fits.
 * The format is:
 *
 * {slopes_mom1, slopeStds_mom1, offsets_mom1, offsetStds_mom1, slopeOffsetCovs_mom1,
 *  slopes_mom2, slopeStds_mom2, offsets_mom2, offsetStds_mom2, slopeOffsetCovs_mom2}
 *
 * Every element in turn contains one value per bootstrap.
 *
 * Log output is generated via \p pLogger. If \p pLogger is \e nullptr, a new Logger is created!
 *
 * \param pSamplesPerSize_moment1 First moment samples for different lattice sizes.
 * \param pSamplesPerSize_moment2 Second moment samples for different lattice sizes.
 * \param pRawFitResults Raw moment fit results for every (outer) bootstrap for first and second moment; please see above.
 * \param pNumBootstraps Number of (outer) bootstraps.
 * \param pNumInnerBootstraps Number of inner/double bootstraps for each outer bootstrap.
 * \param pSeedSeq1 Seed sequence used to seed the random generator for bootstrap sampling (first moment).
 * \param pSeedSeq2 Seed sequence used to seed the random generator for bootstrap sampling (second moment).
 * \param pLogger The logger to be used for logging.
 * \param pNumRNGThreads Number of RNG threads, see makeBootstrapSamples().
 * \return Critical exponent and observable dimension plus uncertainties as {critExp, critExpStd, obsDim, obsDimStd}.
 */
struct MomentAnalysis::Result MomentAnalysis::doMomentAnalysis(std::map<short, MomentSample> pSamplesPerSize_moment1,
                                                               std::map<short, MomentSample> pSamplesPerSize_moment2,
                                                               std::array<std::vector<double>, 10>& pRawFitResults,
                                                               int pNumBootstraps, int pNumInnerBootstraps,
                                                               std::seed_seq& pSeedSeq1, std::seed_seq& pSeedSeq2,
                                                               std::shared_ptr<Logger> pLogger, int pNumRNGThreads)
{
    if (pLogger == nullptr)
        pLogger = std::make_shared<Logger>(Logger::LogLevel::_INFO);

    std::vector<double> fitSlopeSamples_moment1, fitSlopeSamples_moment2;
    std::vector<double> fitSlopeStdSamples_moment1, fitSlopeStdSamples_moment2;

    std::vector<double> fitOffsetSamples_moment1, fitOffsetSamples_moment2;
    std::vector<double> fitOffsetStdSamples_moment1, fitOffsetStdSamples_moment2;

    std::vector<double> fitSlopeOffsetCovSamples_moment1, fitSlopeOffsetCovSamples_moment2;

    pLogger->logInfo("[1/2] Bootstrapped fitting for first moment...");
    pLogger->logInfo("");

    generateBootstrappedFitSlopeSamplesFromMoments(pSamplesPerSize_moment1,
                                                   fitSlopeSamples_moment1, fitSlopeStdSamples_moment1,
                                                   fitOffsetSamples_moment1, fitOffsetStdSamples_moment1,
                                                   fitSlopeOffsetCovSamples_moment1,
                                                   pNumBootstraps, pNumInnerBootstraps, pSeedSeq1, pLogger, pNumRNGThreads);

    pLogger->logInfo("");
    pLogger->logInfo("[2/2] Bootstrapped fitting for second moment...");
    pLogger->logInfo("");

    generateBootstrappedFitSlopeSamplesFromMoments(pSamplesPerSize_moment2,
                                                   fitSlopeSamples_moment2, fitSlopeStdSamples_moment2,
                                                   fitOffsetSamples_moment2, fitOffsetStdSamples_moment2,
                                                   fitSlopeOffsetCovSamples_moment2,
                                                   pNumBootstraps, pNumInnerBootstraps, pSeedSeq2, pLogger, pNumRNGThreads);

    //Assign raw fit results to argument before calculating actual result
    pRawFitResults = {fitSlopeSamples_moment1,
                      fitSlopeStdSamples_moment1,
                      fitOffsetSamples_moment1,
                      fitOffsetStdSamples_moment1,
                      fitSlopeOffsetCovSamples_moment1,
                      fitSlopeSamples_moment2,
                      fitSlopeStdSamples_moment2,
                      fitOffsetSamples_moment2,
                      fitOffsetStdSamples_moment2,
                      fitSlopeOffsetCovSamples_moment2};

    //Calculate bootstrapped critical exponents and observable dimensions

    std::vector<double> critExponentSamples, critDimensionSamples;
    std::vector<double> critExponentStdSamples, critDimensionStdSamples;

    for (size_t i = 0; i < fitSlopeSamples_moment1.size(); ++i)
    {
        double slope_moment1 = fitSlopeSamples_moment1[i];
        double slopeStd_moment1 = fitSlopeStdSamples_moment1[i];

        double slope_moment2 = fitSlopeSamples_moment2[i];
        double slopeStd_moment2 = fitSlopeStdSamples_moment2[i];

        double critExponent = calcCriticalExponent(slope_moment1, slope_moment2);
        double critExponentStd = calcCriticalExponentUncertainty(slope_moment1, slopeStd_moment1, slope_moment2, slopeStd_moment2);

        double critDimension = calcCriticalDimension(slope_moment1, slope_moment2);
        double critDimensionStd = calcCriticalDimensionUncertainty(slopeStd_moment1, slopeStd_moment2);

        critExponentSamples.push_back(critExponent);
        critExponentStdSamples.push_back(critExponentStd);
        critDimensionSamples.push_back(critDimension);
        critDimensionStdSamples.push_back(critDimensionStd);
    }

    //Calculate end result by means of (weighted) bootstrap sample mean and (weighted) bootstrap sample standard deviation

    std::pair<double, double> critExponentResult = calcMeanStdFromVector(critExponentSamples, critExponentStdSamples);
    std::pair<double, double> critDimensionResult = calcMeanStdFromVector(critDimensionSamples, critDimensionStdSamples);

    return {critExponentResult.first, critExponentResult.second, critDimensionResult.first, critDimensionResult.second};
}

//

/*!
 * \brief Calculate moments of avalanche observable distributions from (many) avalanche events.
 *
 * Calculates the first and second Moments (\p pFirstMoments and \p pSecondMoments) of the avalanche
 * size/linSize/area/duration distributions described by a set of avalanche events \p pEvents.
 *
 * \param pEvents Avalanche event observables from a single simulation run.
 * \param pFirstMoments First distribution moments (one for each avalanche observable).
 * \param pSecondMoments Second distribution moments (one for each avalanche observable).
 */
void MomentAnalysis::calcMoments(const std::vector<Event>& pEvents, Moments& pFirstMoments, Moments& pSecondMoments)
{
    double tSize_first = 0, tSize_second = 0;
    double tLinSize_first = 0, tLinSize_second = 0;
    double tArea_first = 0, tArea_second = 0;
    double tDuration_first = 0, tDuration_second = 0;

    for (const Event& event : pEvents)
    {
        tSize_first += event.size;
        tSize_second += std::pow(event.size, 2.);
        tLinSize_first += event.linSize;
        tLinSize_second += std::pow(event.linSize, 2.);
        tArea_first += event.area;
        tArea_second += std::pow(event.area, 2.);
        tDuration_first += event.duration;
        tDuration_second += std::pow(event.duration, 2.);
    }

    pFirstMoments.size = tSize_first / pEvents.size();
    pFirstMoments.linSize = tLinSize_first / pEvents.size();
    pFirstMoments.area = tArea_first / pEvents.size();
    pFirstMoments.duration = tDuration_first / pEvents.size();

    pSecondMoments.size = tSize_second / pEvents.size();
    pSecondMoments.linSize = tLinSize_second / pEvents.size();
    pSecondMoments.area = tArea_second / pEvents.size();
    pSecondMoments.duration = tDuration_second / pEvents.size();
}

/*!
 * \brief Split Moments samples into samples for different observables.
 *
 * - \p pSamplesPerSize[len][i].size --> \p pSamplesPerSize_size[len][i]
 * - \p pSamplesPerSize[len][i].linSize --> \p pSamplesPerSize_linSize[len][i]
 * - \p pSamplesPerSize[len][i].area --> \p pSamplesPerSize_area[len][i]
 * - \p pSamplesPerSize[len][i].duration --> \p pSamplesPerSize_duration[len][i]
 *
 * \param pSamplesPerSize Avalanche observables' moment samples for different lattice sizes.
 * \param pSamplesPerSize_size Avalanche size's moment samples for different lattice sizes.
 * \param pSamplesPerSize_linSize Avalanche linear size's moment samples for different lattice sizes.
 * \param pSamplesPerSize_area Avalanche area's moment samples for different lattice sizes.
 * \param pSamplesPerSize_duration Avalanche duration's moment samples for different lattice sizes.
 */
void MomentAnalysis::splitMoments(const std::map<short, std::vector<Moments>>& pSamplesPerSize,
                                  std::map<short, std::vector<double>>& pSamplesPerSize_size,
                                  std::map<short, std::vector<double>>& pSamplesPerSize_linSize,
                                  std::map<short, std::vector<double>>& pSamplesPerSize_area,
                                  std::map<short, std::vector<double>>& pSamplesPerSize_duration)
{
    for (const auto& sampleIt : pSamplesPerSize)
    {
        int sampleLength = sampleIt.second.size();

        std::vector<double> tmpSample_size, tmpSample_linSize, tmpSample_area, tmpSample_duration;
        tmpSample_size.reserve(sampleLength);
        tmpSample_linSize.reserve(sampleLength);
        tmpSample_area.reserve(sampleLength);
        tmpSample_duration.reserve(sampleLength);

        for (const auto& moments : sampleIt.second)
        {
            tmpSample_size.push_back(moments.size);
            tmpSample_linSize.push_back(moments.linSize);
            tmpSample_area.push_back(moments.area);
            tmpSample_duration.push_back(moments.duration);
        }

        pSamplesPerSize_size.insert(std::make_pair(sampleIt.first, tmpSample_size));
        pSamplesPerSize_linSize.insert(std::make_pair(sampleIt.first, tmpSample_linSize));
        pSamplesPerSize_area.insert(std::make_pair(sampleIt.first, tmpSample_area));
        pSamplesPerSize_duration.insert(std::make_pair(sampleIt.first, tmpSample_duration));
    }
}

//Private

/*!
 * \brief Calculate weighted sample mean and weighted sample standard deviation.
 *
 * Calculates the weighted sample mean and the (unbiased) weighted sample standard deviation for a sample \p pVec
 * with corresponding uncertainties \p pUnc. The weights are inverse variance weights 1/(\p pUnc * \p pUnc).
 *
 * \param pVec The sample.
 * \param pUnc Uncertainties of the sample elements.
 * \return Result as {mean, std}.
 */
std::pair<double, double> MomentAnalysis::calcMeanStdFromVector(const std::vector<double>& pVec, const std::vector<double>& pUnc)
{
    double weightsSum = 0;
    for (double unc : pUnc)
        weightsSum += 1/(unc*unc);

    double squaredWeightsSum = 0;
    for (double unc : pUnc)
        squaredWeightsSum += 1/(unc*unc*unc*unc);

    //Calculate weighted sample mean
    double weightedSampleMean = 0;
    for (std::vector<double>::size_type i = 0; i < pVec.size(); ++i)
        weightedSampleMean += pVec.at(i) / (pUnc.at(i)*pUnc.at(i));
    weightedSampleMean /= weightsSum;

    //Calculate (unbiased) weighted sample variance
    double weightedSampleVar = 0;
    for (std::vector<double>::size_type i = 0; i < pVec.size(); ++i)
        weightedSampleVar += std::pow(pVec.at(i) - weightedSampleMean, 2.) / (pUnc.at(i)*pUnc.at(i));
    weightedSampleVar /= (weightsSum - (squaredWeightsSum / weightsSum));

    return {weightedSampleMean, std::sqrt(weightedSampleVar)};
}

//

/*!
 * \brief Calculate critical exponent from moment fit slopes.
 *
 * Assuming that \p pSlope_moment1 is the fit slope of a straight line fit of ln(distribution's first moment) vs. ln(lattice size)
 * and \p pSlope_moment2 the same for the second moment of the same distribution, the distribution's critical exponent is given by:
 *
 * 3.0 - \p pSlope_moment2 / (\p pSlope_moment2 - \p pSlope_moment1)
 *
 * \param pSlope_moment1 Slope of double-log_e fit for distribution's first moment vs. lattice size.
 * \param pSlope_moment2 Slope of double-log_e fit for distribution's second moment vs. lattice size.
 * \return Distribution's critical exponent.
 */
double MomentAnalysis::calcCriticalExponent(double pSlope_moment1, double pSlope_moment2)
{
    return 3.0 - pSlope_moment2 / (pSlope_moment2 - pSlope_moment1);
}

/*!
 * \brief Calculate critical exponent uncertainty from moment fit slopes and uncertainties.
 *
 * Calculates the uncertainty of calcCriticalExponent()'s result via (uncorrelated) propagation of uncertainty.
 *
 * \param pSlope_moment1 Slope of double-log_e fit for distribution's first moment vs. lattice size.
 * \param pSlopeStd_moment1 Slope uncertainty of double-log_e fit for distribution's first moment vs. lattice size.
 * \param pSlope_moment2 Slope of double-log_e fit for distribution's second moment vs. lattice size.
 * \param pSlopeStd_moment2 Slope uncertainty of double-log_e fit for distribution's second moment vs. lattice size.
 * \return Uncertainty of distribution's critical exponent.
 */
double MomentAnalysis::calcCriticalExponentUncertainty(double pSlope_moment1, double pSlopeStd_moment1,
                                                       double pSlope_moment2, double pSlopeStd_moment2)
{
    return std::sqrt(std::pow(pSlopeStd_moment1 * pSlope_moment2, 2.) +
                     std::pow(pSlopeStd_moment2 * pSlope_moment1, 2.)) / std::pow(pSlope_moment1 - pSlope_moment2, 2.);
}

/*!
 * \brief Calculate observable dimension from moment fit slopes.
 *
 * Assuming that \p pSlope_moment1 is the fit slope of a straight line fit of ln(distribution's first moment) vs. ln(lattice size)
 * and \p pSlope_moment2 the same for the second moment of the same distribution,
 * the distribution's "observable dimension" is given by:
 *
 * \p pSlope_moment2 - \p pSlope_moment1
 *
 * \param pSlope_moment1 Slope of double-log_e fit for distribution's first moment vs. lattice size.
 * \param pSlope_moment2 Slope of double-log_e fit for distribution's second moment vs. lattice size.
 * \return Distribution's observable dimension.
 */
double MomentAnalysis::calcCriticalDimension(double pSlope_moment1, double pSlope_moment2)
{
    return pSlope_moment2 - pSlope_moment1;
}

/*!
 * \brief Calculate observable dimension uncertainty from moment fit slopes and uncertainties.
 *
 * Calculates the uncertainty of calcCriticalDimension()'s result via (uncorrelated) propagation of uncertainty.
 *
 * \param pSlopeStd_moment1 Slope uncertainty of double-log_e fit for distribution's first moment vs. lattice size.
 * \param pSlopeStd_moment2 Slope uncertainty of double-log_e fit for distribution's second moment vs. lattice size.
 * \return Uncertainty of distribution's observable dimension.
 */
double MomentAnalysis::calcCriticalDimensionUncertainty(double pSlopeStd_moment1, double pSlopeStd_moment2)
{
    return std::sqrt(pSlopeStd_moment1 * pSlopeStd_moment1 + pSlopeStd_moment2 * pSlopeStd_moment2);
}

//

/*!
 * \brief Perform a straight line fit for data with correlated y-values.
 *
 * Fits a straight line to the data points (\p pX_i, \p pY_i)_i assuming
 * correlated y-values as described by the covariance matrix \p pCovMat.
 *
 * The resulting fit parameters and their uncertainties
 * are returned as {slope, slopeStd, offset, offsetStd, slopeOffsetCov}.
 *
 * The fit is based on a matrix equation being solved/calculated using GSL.
 * If a GSL operation or allocation fails, the function returns {0, 0, 0, 0, 0}.
 *
 * \param pX x-values of data points.
 * \param pY y-values of data points.
 * \param pCovMat Covariance matrix for \p pY.
 * \return The fit result as {slope, slopeStd, offset, offsetStd, slopeOffsetCov}.
 */
std::tuple<double, double, double, double, double> MomentAnalysis::doFit(const std::vector<double>& pX,
                                                                         const std::vector<double>& pY,
                                                                         const std::vector<std::vector<double>>& pCovMat)
{
    /*
     * Perform the fit for function y=ax+b by solving
     *
     *  pY = Z * beta
     *
     * for beta, with
     *
     *      /1 pX_1\         /b\
     *  Z = |1 ... |, beta = | |.
     *      \1 pX_n/         \a/
     *
     * One finds:
     *
     *  beta = (Z^T pCovMat^{-1} Z)^{-1} * Z^T * pCovMat^{-1} * pY
     *
     *  deltaBeta = (Z^T pCovMat^{-1} Z)^{-1}
     */

    //Check errors by return values
    gsl_set_error_handler_off();

    //Prepare matrix 'Z'

    gsl_matrix* zMat = gsl_matrix_alloc(pX.size(), 2);
    if (zMat == nullptr)
        return {0, 0, 0, 0, 0};

    for (size_t row = 0; row < zMat->size1; ++row)
    {
        gsl_matrix_set(zMat, row, 0, 1);
        gsl_matrix_set(zMat, row, 1, pX.at(row));
    }

    //Scale/normalize covariance matrix before inversion

    double meanCovDiag = 0;
    for (size_t row = 0; row < pCovMat.size(); ++row)
        meanCovDiag += pCovMat.at(row).at(row);
    meanCovDiag /= pCovMat.size();

    gsl_matrix* covMat = gsl_matrix_alloc(pY.size(), pY.size());
    if (covMat == nullptr)
        return {0, 0, 0, 0, 0};

    for (size_t row = 0; row < pCovMat.size(); ++row)
        for (size_t col = 0; col <= row; ++col) //Cholesky decomp. only needs lower triangular part
            gsl_matrix_set(covMat, row, col, pCovMat.at(row).at(col) / meanCovDiag);

    //Invert covariance matrix via Cholesky decomposition (pCovMat --> pCovMat^{-1})

    gsl_permutation* tmpGSLPermutation0 = gsl_permutation_alloc(covMat->size1);
    if (tmpGSLPermutation0 == nullptr)
        return {0, 0, 0, 0, 0};

    gsl_permutation_init(tmpGSLPermutation0);

    if (gsl_linalg_pcholesky_decomp(covMat, tmpGSLPermutation0) != 0)
        return {0, 0, 0, 0, 0};

    gsl_matrix* invCovMat = gsl_matrix_alloc(covMat->size1, covMat->size2);
    if (invCovMat == nullptr)
        return {0, 0, 0, 0, 0};

    if (gsl_linalg_pcholesky_invert(covMat, tmpGSLPermutation0, invCovMat) != 0)
        return {0, 0, 0, 0, 0};

    gsl_permutation_free(tmpGSLPermutation0);

    //Re-scale inverse covariance matrix after inversion
    for (size_t row = 0; row < pCovMat.size(); ++row)
        for (size_t col = 0; col < pCovMat.size(); ++col)
            gsl_matrix_set(invCovMat, row, col, gsl_matrix_get(invCovMat, row, col) / meanCovDiag);

    //Calculate (Z^T pCovMat^{-1})

    gsl_matrix* tMat1 = gsl_matrix_alloc(zMat->size2, invCovMat->size2);
    if (tMat1 == nullptr)
        return {0, 0, 0, 0, 0};

    if (gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1, zMat, invCovMat, 0, tMat1) != 0)
        return {0, 0, 0, 0, 0};

    //Calculate (Z^T pCovMat^{-1} Z)

    gsl_matrix* tMat2 = gsl_matrix_alloc(tMat1->size1, zMat->size2);
    if (tMat2 == nullptr)
        return {0, 0, 0, 0, 0};

    if (gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1, tMat1, zMat, 0, tMat2) != 0)
        return {0, 0, 0, 0, 0};

    //Scale/normalize matrix (Z^T pCovMat^{-1} Z) before inversion

    double meanVal = 0;
    meanVal += gsl_matrix_get(tMat2, 0, 0) + gsl_matrix_get(tMat2, 0, 1) + gsl_matrix_get(tMat2, 1, 0) + gsl_matrix_get(tMat2, 1, 1);
    meanVal /= 4;

    for (size_t row = 0; row < 2; ++row)
        for (size_t col = 0; col < 2; ++col)
            gsl_matrix_set(tMat2, row, col, gsl_matrix_get(tMat2, row, col) / meanVal);

    //Invert matrix (Z^T pCovMat^{-1} Z) via LU decomposition

    gsl_matrix* tMat3 = gsl_matrix_alloc(tMat2->size1, tMat2->size2);
    if (tMat3 == nullptr)
        return {0, 0, 0, 0, 0};

    int tmpInt = 0;
    gsl_permutation* tmpGSLPermutation1 = gsl_permutation_alloc(tMat2->size1);
    if (tmpGSLPermutation1 == nullptr)
        return {0, 0, 0, 0, 0};

    gsl_permutation_init(tmpGSLPermutation1);

    if (gsl_linalg_LU_decomp(tMat2, tmpGSLPermutation1, &tmpInt) != 0)
        return {0, 0, 0, 0, 0};

    if (gsl_linalg_LU_invert(tMat2, tmpGSLPermutation1, tMat3) != 0)
        return {0, 0, 0, 0, 0};

    gsl_permutation_free(tmpGSLPermutation1);

    //Re-scale matrix after inversion
    for (size_t row = 0; row < 2; ++row)
        for (size_t col = 0; col < 2; ++col)
            gsl_matrix_set(tMat3, row, col, gsl_matrix_get(tMat3, row, col) / meanVal);

    //Calculate (Z^T pCovMat^{-1} Z)^{-1} * Z^T

    gsl_matrix* tMat4 = gsl_matrix_alloc(tMat3->size1, zMat->size1);
    if (tMat4 == nullptr)
        return {0, 0, 0, 0, 0};

    if (gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1, tMat3, zMat, 0, tMat4) != 0)
        return {0, 0, 0, 0, 0};

    //Calculate (Z^T pCovMat^{-1} Z)^{-1} * Z^T * pCovMat^{-1}

    gsl_matrix* tMat5 = gsl_matrix_alloc(tMat4->size1, invCovMat->size2);
    if (tMat5 == nullptr)
        return {0, 0, 0, 0, 0};

    if (gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1, tMat4, invCovMat, 0, tMat5) != 0)
        return {0, 0, 0, 0, 0};

    //Calculate fit result [(Z^T pCovMat^{-1} Z)^{-1} * Z^T * pCovMat^{-1} * pY]

    gsl_matrix* yVec = gsl_matrix_alloc(pY.size(), 1);
    if (yVec == nullptr)
        return {0, 0, 0, 0, 0};

    for (size_t row = 0; row < yVec->size1; ++row)
        gsl_matrix_set(yVec, row, 0, pY.at(row));

    gsl_matrix* result = gsl_matrix_alloc(tMat5->size1, yVec->size2);
    if (result == nullptr)
        return {0, 0, 0, 0, 0};

    if (gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1, tMat5, yVec, 0, result) != 0)
        return {0, 0, 0, 0, 0};

    double fitOffset = gsl_matrix_get(result, 0, 0);
    double fitSlope = gsl_matrix_get(result, 1, 0);

    double fitOffsetStd = std::sqrt(gsl_matrix_get(tMat3, 0, 0));
    double fitSlopeStd = std::sqrt(gsl_matrix_get(tMat3, 1, 1));

    double fitSlopeOffsetCov = gsl_matrix_get(tMat3, 1, 0);

    gsl_matrix_free(result);
    gsl_matrix_free(yVec);
    gsl_matrix_free(tMat5);
    gsl_matrix_free(tMat4);
    gsl_matrix_free(tMat3);
    gsl_matrix_free(tMat2);
    gsl_matrix_free(tMat1);
    gsl_matrix_free(invCovMat);
    gsl_matrix_free(covMat);
    gsl_matrix_free(zMat);

    return {fitSlope, fitSlopeStd, fitOffset, fitOffsetStd, fitSlopeOffsetCov};
}

//

/*!
 * \brief Generate bootstrap samples.
 *
 * Generates uniformly distributed bootstrap samples \p pBootstrapSamples from an
 * original sample \p pOriginalSample using Mersenne Twister engines. The required
 * random indexes are generated in \p pNumThreads parallel OMP threads, each with a
 * separate RNG engine. All these engines are seeded in the beginning by \p pRndGenerator.
 * The length of the pre-filled \p pBootstrapSamples defines the number of generated samples.
 * Furthermore \p pBootstrapSamples is assumed to contain pre-filled vectors that all
 * have the same length as \p pOriginalSample; the values will be simply overwritten.
 *
 * Note: The number of threads is automatically set, if \p pNumThreads is 0 or negative.
 *
 * \param pOriginalSample Original data to draw the bootstrap samples from.
 * \param pBootstrapSamples (Pre-filled) Destination for the generated samples.
 * \param pRndGenerator Used Mersenne Twister engine.
 * \param pNumThreads Number of parallel RNG threads.
 */
void MomentAnalysis::makeBootstrapSamples(const MomentSample& pOriginalSample, std::vector<MomentSample>& pBootstrapSamples,
                                          std::mt19937_64& pRndGenerator, int pNumThreads)
{
    auto sampleLength = pOriginalSample.size();
    auto numSamples = pBootstrapSamples.size();

    //Use default number of threads, if not set to sensible value
    if (pNumThreads < 1)
        pNumThreads = omp_get_max_threads();

    //One RNG per thread, all seeded by 'pRndGenerator'
    std::vector<std::uniform_int_distribution<int>> rndUnifs(pNumThreads, std::uniform_int_distribution<int>(0, sampleLength-1));
    std::vector<std::mt19937_64> rndGens(pNumThreads, std::mt19937_64());
    for (auto& rndGen : rndGens)
        rndGen.seed(pRndGenerator());

    //Temporary storage for generated random indexes
    std::vector<std::remove_reference<decltype(pOriginalSample)>::type::size_type> rndIdxs(numSamples*sampleLength);

    omp_set_dynamic(0);
    omp_set_num_threads(pNumThreads);

    //Generate random indexes
    #pragma omp parallel for shared(rndIdxs)
    for (decltype(sampleLength) i = 0; i < rndIdxs.size(); ++i)
    {
        auto tNum = omp_get_thread_num();
        rndIdxs[i] = rndUnifs[tNum](rndGens[tNum]);
    }

    //Assign random bootstrap samples via previously generated indexes
    decltype(rndIdxs)::size_type i = 0;
    for (auto& newBootstrapSample : pBootstrapSamples)
        for (double& element : newBootstrapSample)
            element = pOriginalSample.at(rndIdxs[i++]);
}

/*!
 * \brief Generate/calculate straight line fit slopes of bootstrapped, logarithmized moment samples for a single observable.
 *
 * First logarithmizes (natural logarithm) the moment samples for every lattice size from \p pSamplesPerSize.
 * Then generates \p pNumBootstraps bootstrap samples from these logarithmized original samples via
 * makeBootstrapSamples(). For each (outer) bootstrap the sample means for every lattice size are then
 * fitted vs. the correspondingly logarithmized lattice sizes (keys of \p pSamplesPerSize) using doFit().
 * The covariance matrix for each fit is estimated via a double bootstrap procedure, drawing \p pNumInnerBootstraps
 * inner bootstrap samples for every (outer) bootstrap sample. The resulting fit slopes and slope uncertainties
 * for every (outer) bootstrap are finally written \p pFitSlopeSamples and \p pFitSlopeStdSamples, respectively.
 * Analogously, fit offsets, offset uncertainties and slope-offset covariances are written to
 * \p pFitOffsetSamples, \p pFitOffsetStdSamples and \p pFitSlopeOffsetCovSamples.
 *
 * Log output is generated via \p pLogger. If \p pLogger is \e nullptr, a new Logger is created!
 *
 * \param pSamplesPerSize Moment samples for different lattice sizes.
 * \param pFitSlopeSamples Fit slope results for every (outer) bootstrap.
 * \param pFitSlopeStdSamples Fit slope uncertainty results for every (outer) bootstrap.
 * \param pFitOffsetSamples Fit offset results for every (outer) bootstrap.
 * \param pFitOffsetStdSamples Fit offset uncertainty results for every (outer) bootstrap.
 * \param pFitSlopeOffsetCovSamples Fit slope-offset covariances for every (outer) bootstrap.
 * \param pNumBootstraps Number of (outer) bootstraps.
 * \param pNumInnerBootstraps Number of inner/double bootstraps for each outer bootstrap.
 * \param pSeedSeq Seed sequence used to seed the random generator for bootstrap sampling.
 * \param pLogger The logger to be used for logging.
 * \param pNumRNGThreads Number of RNG threads, see makeBootstrapSamples().
 */
void MomentAnalysis::generateBootstrappedFitSlopeSamplesFromMoments(std::map<short, MomentSample> pSamplesPerSize,
                                                                    std::vector<double>& pFitSlopeSamples,
                                                                    std::vector<double>& pFitSlopeStdSamples,
                                                                    std::vector<double>& pFitOffsetSamples,
                                                                    std::vector<double>& pFitOffsetStdSamples,
                                                                    std::vector<double>& pFitSlopeOffsetCovSamples,
                                                                    int pNumBootstraps, int pNumInnerBootstraps,
                                                                    std::seed_seq& pSeedSeq,
                                                                    std::shared_ptr<Logger> pLogger, int pNumRNGThreads)
{
    //Logarithmize moments such that linear fit can be applied

    for (auto& sample_it : pSamplesPerSize)
        for (double& moment : sample_it.second)
            moment = std::log(moment);

    std::vector<double> logLatticeSizes;
    for (const auto& sample_it : pSamplesPerSize)
        logLatticeSizes.push_back(std::log(sample_it.first));

    std::vector<short> boxSizes;
    for (const auto& it : pSamplesPerSize)
        boxSizes.push_back(it.first);

    MomentSample::size_type sampleLength = pSamplesPerSize.at(boxSizes.front()).size();

    //Generate and seed Mersenne Twister pRNG engine for bootstrap sampling
    std::mt19937_64 rndGenerator(pSeedSeq);

    //Need to store (only) means of bootstrap samples for every lattice size
    std::map<short, std::vector<double>> outerBootstrapMeansPerSize;
    std::map<short, std::vector<std::vector<double>>> innerBootstrapMeansPerSize;
    std::map<short, std::vector<double             >> meanInnerBootstrapMeansPerSize;

    //Generate outer and inner bootstrap samples, calculate corresponding sample means, discard actual samples again

    std::vector<MomentSample> outerBootstraps = std::vector<MomentSample>(pNumBootstraps, std::vector<double>(sampleLength));

    pLogger->logInfo("- Generating (double-)bootstrap samples...");

    int tLatSizeCtr = 1;
    for (const auto& sample_it : pSamplesPerSize)
    {
        pLogger->logMore("   Lattice size " + std::to_string(tLatSizeCtr++) + "/" + std::to_string(boxSizes.size()));

        //Generate outer samples
        makeBootstrapSamples(sample_it.second, outerBootstraps, rndGenerator, pNumRNGThreads);

        //Calculate outer sample means

        std::vector<double> outerBootstrapMeans;
        outerBootstrapMeans.reserve(pNumBootstraps);

        for (const MomentSample& bootstrapSample : outerBootstraps)
            outerBootstrapMeans.push_back(Aux::vectorSum(bootstrapSample) / bootstrapSample.size());

        //Do double-bootstrap for each outer sample

        std::vector<std::vector<double>> innerBootstrapMeans;
        std::vector<double             > meanInnerBootstrapMeans;

        //Reserve memory
        innerBootstrapMeans = std::vector<std::vector<double>>(pNumBootstraps, std::vector<double>(pNumInnerBootstraps));
        meanInnerBootstrapMeans.reserve(pNumBootstraps);

        std::vector<MomentSample> innerBootstraps = std::vector<MomentSample>(pNumInnerBootstraps, std::vector<double>(sampleLength));

        int tNumOuterSample = 0;
        for (const MomentSample& bootstrapSample : outerBootstraps)
        {
            //Generate inner samples
            makeBootstrapSamples(bootstrapSample, innerBootstraps, rndGenerator, pNumRNGThreads);

            //Calculate inner sample means

            std::vector<double> tmpInnerBootstrapMeans;
            tmpInnerBootstrapMeans.reserve(pNumInnerBootstraps);

            for (const MomentSample& innerBootstrapSample : innerBootstraps)
                tmpInnerBootstrapMeans.push_back(Aux::vectorSum(innerBootstrapSample) / innerBootstrapSample.size());

            innerBootstrapMeans.at(tNumOuterSample) = tmpInnerBootstrapMeans;

            //Calculate mean of inner sample means (needed for covariance matrix)
            meanInnerBootstrapMeans.push_back(Aux::vectorSum(tmpInnerBootstrapMeans) / tmpInnerBootstrapMeans.size());

            ++tNumOuterSample;
        }

        //Store sample means for current lattice size
        outerBootstrapMeansPerSize.insert(std::make_pair(sample_it.first, outerBootstrapMeans));
        innerBootstrapMeansPerSize.insert(std::make_pair(sample_it.first, innerBootstrapMeans));
        meanInnerBootstrapMeansPerSize.insert(std::make_pair(sample_it.first, meanInnerBootstrapMeans));
    }

    //Calculate double-bootstrapped covariance matrix for each outer bootstrap sample from inner bootstrap means

    typedef std::vector<std::vector<double>> CovMat;

    std::vector<CovMat> covMatPerOuterSample;
    covMatPerOuterSample.reserve(pNumBootstraps);

    pLogger->logInfo("- Calculating covariance matrices...");

    for (int i = 0; i < pNumBootstraps; ++i)
    {
        std::vector<std::vector<double>> covMat = std::vector<std::vector<double>>(boxSizes.size(),
                                                                                   std::vector<double>(boxSizes.size()));

        for (size_t row = 0; row < boxSizes.size(); ++row)
        {
            short lI = boxSizes.at(row);

            std::vector<double> tMatRow = std::vector<double>(boxSizes.size());

            for (size_t col = 0; col <= row; ++col) //doFit() below only needs lower triangular part of covariance matrix
            {
                short lJ = boxSizes.at(col);

                double tMatElement = 0;

                for (size_t dotProdIt = 0; dotProdIt < innerBootstrapMeansPerSize.at(lI).at(i).size(); ++dotProdIt)
                {
                    tMatElement += (
                                       (innerBootstrapMeansPerSize.at(lI).at(i).at(dotProdIt) -
                                        meanInnerBootstrapMeansPerSize.at(lI).at(i)
                                        ) *
                                       (innerBootstrapMeansPerSize.at(lJ).at(i).at(dotProdIt) -
                                        meanInnerBootstrapMeansPerSize.at(lJ).at(i)
                                        )
                                    );
                }

                tMatElement /= (pNumInnerBootstraps - 1);

                tMatRow.at(col) = tMatElement;
            }

            covMat.at(row) = tMatRow;
        }

        covMatPerOuterSample.push_back(covMat);
    }

    //Free memory
    innerBootstrapMeansPerSize.clear();
    meanInnerBootstrapMeansPerSize.clear();

    //Rearrange bootstrap means data to have outer sample index first and lattice size second

    typedef std::vector<double> BootstrapMeansPerSize;

    std::vector<BootstrapMeansPerSize> meansPerSizePerOuterSample;
    meansPerSizePerOuterSample.reserve(pNumBootstraps);

    for (int i = 0; i < pNumBootstraps; ++i)
    {
        BootstrapMeansPerSize outerBootstrapMeanPerSize;
        for (const auto& outerBootstrapMeans_it : outerBootstrapMeansPerSize)
            outerBootstrapMeanPerSize.push_back(outerBootstrapMeans_it.second[i]);

        meansPerSizePerOuterSample.push_back(outerBootstrapMeanPerSize);
    }

    //Free memory
    outerBootstrapMeansPerSize.clear();

    pLogger->logInfo("- Performing fits...");

    //Do a fit for each outer bootstrap sample using the double-bootstrapped covariance matrices from above
    for (int i = 0; i < pNumBootstraps; ++i)
    {
        //Do fit
        std::tuple<double, double, double, double, double> fitResult = doFit(logLatticeSizes, meansPerSizePerOuterSample.at(i),
                                                                             covMatPerOuterSample.at(i));

        if (std::get<0>(fitResult) == 0 && std::get<1>(fitResult) == 0 && std::get<2>(fitResult) == 0 && std::get<3>(fitResult) == 0)
            pLogger->logError("Fit (some GSL operation) failed!");

        double fitSlope = std::get<0>(fitResult);
        double fitSlopeStd = std::get<1>(fitResult);

        double fitOffset = std::get<2>(fitResult);
        double fitOffsetStd = std::get<3>(fitResult);

        double fitSlopeOffsetCov = std::get<4>(fitResult);

        pFitSlopeSamples.push_back(fitSlope);
        pFitSlopeStdSamples.push_back(fitSlopeStd);

        pFitOffsetSamples.push_back(fitOffset);
        pFitOffsetStdSamples.push_back(fitOffsetStd);

        pFitSlopeOffsetCovSamples.push_back(fitSlopeOffsetCov);
    }
}

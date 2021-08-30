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

#ifndef AUX_H
#define AUX_H

#include <iostream>
#include <sstream>
#include <type_traits>
#include <unordered_map>
#include <vector>
#include <random>

/*!
 * \brief Auxiliary class.
 *
 * Provides static member functions for various purposes.
 */
class Aux
{
private:
    Aux() = delete;         ///< Deleted constructor.

public:
    enum class Model        /// Enumeration of IDs of available simulation models.
    {
        __INVALID_MODEL,    ///< Error code / invalid model.
        _BTW,               ///< Bak-Tang-Wiesenfeld model (see SimulationModel_BTW).
        _FWM                ///< A customized sandpile model, the FW model (see SimulationModel_FWM).
    };

public:
    static std::string getModelStringFromId(Model pModelId);            ///< String representation of simulation model ID.
    static Model getModelIdFromString(const std::string& pModelString); ///< Integer representation of simulation model string.
    //
    static std::seed_seq generateSeedSeq(std::mt19937_64& pRndGen);     ///< Generate a seed sequence.
    //
    static long long vectorDistanceSquared(const std::vector<short>& pLeft, const std::vector<short>& pRight);
                                                                                            ///< \brief Compute squared Euclidean
                                                                                            ///  distance between two STL vectors.
    /*!
     * \brief Compute the sum over all vector elements.
     *
     * \tparam T Type of the vector elements.
     *
     * \param pVec Vector with arithmetic value type.
     */
    template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
    static T vectorSum(const std::vector<T>& pVec)
    {
        T retVal = 0;
        for (T num : pVec)
            retVal += num;
        return retVal;
    }
    //
    /*!
     * \brief Convert an arithmetic type to std::string in scientific notation.
     *
     * \tparam T Type of the number to be converted.
     *
     * \param pNumber Number to be converted.
     */
    template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
    static std::string numberToStringScientific(T pNumber)      //Convert an arithmetic type to std::string in scientific notation
    {
        char buf[15];
        std::snprintf(buf, 15, "%E", pNumber);
        return std::string(buf);
    }
    //
    /*!
     * \brief Loop over points of a hyperrectangle of variable dimension and execute a specified function for each point.
     *
     * Uses recursion to loop over all points of a hyperrectangle via nested for loops.
     * Within the innermost loop \p pFunction is called with first parameter the current
     * point of the loop and perhaps further parameters \p pArgs.
     * Hence a function call \p pFunction(pt, pArgs...) happens once for each point pt.
     * Non-void return values will be discarded. Use \p pArgs to communicate results of the function call.
     *
     * The hyperrectangle can be of variable dimension N = \p pRecursionDepth and it is bounded by the half-open interval
     * [\p pLowerBounds, \p pUpperBounds). This means the outermost loop ranges from pLowerBounds[N-1] to pUpperBounds[N-1]-1
     * and the innermost loop from pLowerBounds[0] to pUpperBounds[0]-1.
     *
     * \tparam FuncT Specific type of a FunctionObject (\p pFunction) that shall be called for each point.
     * \tparam ...Args Types of all additional parameters for \p pFunction if there are any.
     *
     * \param pFunction Any FunctionObject that accepts at least \p pCurrentSite as the first parameter
     *                  and optionally also further parameters \p pArgs. A lambda expression may be used.
     *                  The function's return value will be discarded.
     *                  Note: (some of) \p pArgs can be pass by reference in order to communicate the results.
     * \param pLowerBounds "Bottom-left" corner of the hyperrectangle that restricts the loop ranges;
     *                     the corner is included in the loop range. Bottom-left means that each component
     *                     of the point must be smaller than the corresponding component of \p pUpperBounds.
     *                     If any component is larger or equal then the recursion becomes ineffective
     *                     and no call of \p pFunction will happen.
     * \param pUpperBounds "Top-right" corner of the hyperrectangle that restricts the loop ranges;
     *                     neither the corner nor corresponding edges are included in the loop range.
     *                     Top-right means that each component of the point must be larger than the
     *                     corresponding component of \p pLowerBounds. If not, see \p pLowerBounds.
     * \param pCurrentSite Point resulting from the previous level of recursion. A point must be specified
     *                     for the initial call as well, while the actual coordinates of this initial point
     *                     don't matter as they are overwritten by the algorithm. But the point must be there
     *                     and of appropriate dimension in order to be usable for the algorithm as "temporary variable".
     * \param pRecursionDepth Current recursion depth. Set to N, initially, for looping over points of an N-dim. rectangle.
     * \param pArgs Optional arguments that will be passed to \p pFunction in addition to \p pCurrentSite.
     *
     *
     * \attention
     * - If the first parameter of \p pFunction (i.e. here: \p pCurrentSite) is not accepted by value
     *   but by reference and is then modified within \p pFunction, you might not obtain the expected
     *   results due to interference with the loops!
     * - \p pCurrentSite will be modified! Hence, pass a separate variable for \p pCurrentSite,
     *   not the same as for \p pLowerBounds or \p pUpperBounds!
     */
    template <typename FuncT, typename ... Args>
    static void recursiveForLoop(FuncT pFunction, const std::vector<short>& pLowerBounds, const std::vector<short>& pUpperBounds,
                                 std::vector<short>& pCurrentSite, size_t pRecursionDepth, Args& ... pArgs)
    {
        if (pRecursionDepth > 0)
        {
            for (int i = pLowerBounds[pRecursionDepth-1]; i < pUpperBounds[pRecursionDepth-1]; ++i)
            {
                pCurrentSite[pRecursionDepth-1] = i;
                recursiveForLoop<FuncT, Args ...>(pFunction, pLowerBounds, pUpperBounds, pCurrentSite, pRecursionDepth-1,
                                                  pArgs ...);
            }
        }
        else
        {
            pFunction(pCurrentSite, pArgs ...);
        }
    }

private:
    static const std::unordered_map<Model, std::string> mapIdToStringMap;   //Map simulation model integer IDs to model strings;
                                                                            //[definition in aux.cpp]
};

/*!
 * \brief Contains some useful operator overloads for STL vector arithmetics.
 */
namespace OpOverloads_STL_vector {
    /*! \relatesalso Aux */
    std::vector<short> operator+(const std::vector<short>& pLhs, const std::vector<short>& pRhs);   ///< \brief Overloaded '+'
                                                                                                /// operator: performs element-wise
                                                                                                /// addition of two STL vectors.
    /*! \relatesalso Aux */
    std::vector<double> operator+(const std::vector<double>& pLhs, const std::vector<double>& pRhs);///< \brief Overloaded '+'
                                                                                                /// operator: performs element-wise
                                                                                                /// addition of two STL vectors.
    /*! \relatesalso Aux */
    std::vector<short> operator-(const std::vector<short>& pLhs, const std::vector<short>& pRhs);   ///< \brief Overloaded '-'
                                                                                                /// operator: performs element-wise
                                                                                                /// subtraction of two STL vectors.
    /*! \relatesalso Aux */
    std::vector<double> operator-(const std::vector<double>& pLhs, const std::vector<double>& pRhs);///< \brief Overloaded '-'
                                                                                                    /// operator: performs
                                                                                                    /// element-wise subtraction
                                                                                                    /// of two STL vectors.
}

#endif // AUX_H

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

#include "aux.h"

//

/*!
 * \brief Maps simulation model integer IDs to model strings
 *
 * Defines what is returned by getModelStringFromId() and getModelIdFromString().
 * See also Aux::MODEL.
 */
const std::unordered_map<Aux::Model, std::string> Aux::mapIdToStringMap = {
    {Aux::Model::_BTW, "BTW"},
    {Aux::Model::_FWM, "FWM"}
};

/*!
 * \brief String representation of simulation model ID.
 *
 * Converts the unique integer/enum ID of a specific implementation of SimulationModel to the corresponding unique model string.
 *
 * If \p pModelId is invalid, the model string for the BTW model is returned.
 *
 * \param pModelId Unique integer/enum model ID.
 * \return Unique model string.
 */
std::string Aux::getModelStringFromId(Model pModelId)
{
    switch (pModelId)
    {
        case Model::_BTW:
        {
            return mapIdToStringMap.at(Model::_BTW);
            break;
        }
        case Model::_FWM:
        {
            return mapIdToStringMap.at(Model::_FWM);
            break;
        }
        default:
        {
            return mapIdToStringMap.at(Model::_BTW);
            break;
        }
    }
}

/*!
 * \brief Integer representation of simulation model string.
 *
 * Converts the unique model string of a specific implementation of SimulationModel to the corresponding unique integer/enum ID.
 *
 * If the provided string does not match any of the available models then Model::__INVALID_MODEL is returned.
 *
 * \param pModelString Unique model string.
 * \return Unique integer model ID.
 */
Aux::Model Aux::getModelIdFromString(const std::string& pModelString)
{
    //Search the predefined map for the model
    for (auto pair : mapIdToStringMap)
        if (pModelString == pair.second)
            return pair.first;

    //Return error code if model not found
    return Model::__INVALID_MODEL;
}

//

/*!
 * \brief Generate a seed sequence.
 *
 * \param pRndGen Random generator used to initialize the seed sequence.
 * \return New seed sequence {\p pRndGen(), \p pRndGen(), \p pRndGen(), \p pRndGen(),
 *                            \p pRndGen(), \p pRndGen(), \p pRndGen(), \p pRndGen()}.
 */
std::seed_seq Aux::generateSeedSeq(std::mt19937_64& pRndGen)
{
    return {pRndGen(), pRndGen(), pRndGen(), pRndGen(), pRndGen(), pRndGen(), pRndGen(), pRndGen()};
}

//

/*!
 * \brief Compute squared Euclidean distance between two STL vectors.
 *
 * \param pLeft One Cartesian vector.
 * \param pRight Another Cartesian vector of same length.
 * \return Squared Euclidean distance || \p pLeft - \p pRight || (or 0 if \p pLeft is empty).
 *
 * \attention Uses \p pLeft .size() to determine the size of \b both vectors.
 *            Hence undefined behavior for \p pRight smaller than \p pLeft.
 *            If \p pRight is larger than \p pLeft the function uses the truncated vector
 *            {pRight_0, ..., pRight_k}, k=\p pLeft .size()-1 for the computation.
 */
long long Aux::vectorDistanceSquared(const std::vector<short>& pLeft, const std::vector<short>& pRight)
{
    long long result = 0;

    for (size_t i = 0; i < pLeft.size(); ++i)
        result += (pLeft[i]-pRight[i])*(pLeft[i]-pRight[i]);

    return result;
}

//Non-members

namespace OpOverloads_STL_vector
{

/*!
 * \relatesalso Aux
 *
 * \brief Overloaded '+' operator: performs element-wise addition of two STL vectors.
 *
 * \param pLhs Vector on the left-hand side of the + operator.
 * \param pRhs Vector on the right-hand side of the + operator.
 * \return Element-wise sum of \p pLhs and \p pRhs.
 *
 * \attention Uses \p pLhs .size() to determine the size of \b both vectors.
 *            Hence undefined behavior for \p pRhs smaller than \p pLhs.
 *            If \p pRhs is larger than \p pLhs the function uses the truncated vector
 *            {pRhs_0, ..., pRhs_k}, k=\p pLhs .size()-1 for the computation.
 */
std::vector<short> operator+(const std::vector<short>& pLhs, const std::vector<short>& pRhs)
{
    size_t size = pLhs.size();
    std::vector<short> ret(size);

    for (size_t i = 0; i < size; ++i)
        ret[i] = pLhs[i] + pRhs[i];

    return ret;
}

/*!
 * \relatesalso Aux
 *
 * \brief Overloaded '+' operator: performs element-wise addition of two STL vectors.
 *
 * \param pLhs Vector on the left-hand side of the + operator.
 * \param pRhs Vector on the right-hand side of the + operator.
 * \return Element-wise sum of \p pLhs and \p pRhs.
 *
 * \attention Uses \p pLhs .size() to determine the size of \b both vectors.
 *            Hence undefined behavior for \p pRhs smaller than \p pLhs.
 *            If \p pRhs is larger than \p pLhs the function uses the truncated vector
 *            {pRhs_0, ..., pRhs_k}, k=\p pLhs .size()-1 for the computation.
 */
std::vector<double> operator+(const std::vector<double>& pLhs, const std::vector<double>& pRhs)
{
    size_t size = pLhs.size();
    std::vector<double> ret(size);

    for (size_t i = 0; i < size; ++i)
        ret[i] = pLhs[i] + pRhs[i];

    return ret;
}

/*!
 * \relatesalso Aux
 *
 * \brief Overloaded '-' operator: performs element-wise subtraction of two STL vectors.
 *
 * \param pLhs Vector on the left-hand side of the - operator.
 * \param pRhs Vector on the right-hand side of the - operator.
 * \return Element-wise difference of \p pLhs and \p pRhs.
 *
 * \attention Uses \p pLhs .size() to determine the size of \b both vectors.
 *            Hence undefined behavior for \p pRhs smaller than \p pLhs.
 *            If \p pRhs is larger than \p pLhs the function uses the truncated vector
 *            {pRhs_0, ..., pRhs_k}, k=\p pLhs .size()-1 for the computation.
 */
std::vector<short> operator-(const std::vector<short>& pLhs, const std::vector<short>& pRhs)
{
    size_t size = pLhs.size();
    std::vector<short> ret(size);

    for (size_t i = 0; i < size; ++i)
        ret[i] = pLhs[i] - pRhs[i];

    return ret;
}

/*!
 * \relatesalso Aux
 *
 * \brief Overloaded '-' operator: performs element-wise subtraction of two STL vectors.
 *
 * \param pLhs Vector on the left-hand side of the - operator.
 * \param pRhs Vector on the right-hand side of the - operator.
 * \return Element-wise difference of \p pLhs and \p pRhs.
 *
 * \attention Uses \p pLhs .size() to determine the size of \b both vectors.
 *            Hence undefined behavior for \p pRhs smaller than \p pLhs.
 *            If \p pRhs is larger than \p pLhs the function uses the truncated vector
 *            {pRhs_0, ..., pRhs_k}, k=\p pLhs .size()-1 for the computation.
 */
std::vector<double> operator-(const std::vector<double>& pLhs, const std::vector<double>& pRhs)
{
    size_t size = pLhs.size();
    std::vector<double> ret(size);

    for (size_t i = 0; i < size; ++i)
        ret[i] = pLhs[i] - pRhs[i];

    return ret;
}

} // namespace OpOverloads_STL_vector

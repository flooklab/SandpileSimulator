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

#include "randomizer.h"

/*!
 * \brief Constructor.
 *
 * \param pSandboxShape Shape of the sandbox in units of lattice sites (see also Sandbox).
 * \param pSeedSeq Seed sequence used to seed the internal random generator (see seed()).
 *
 * \throws std::out_of_range %Sandbox dimension is smaller than 1.
 */
Randomizer::Randomizer(std::vector<short> pSandboxShape, std::seed_seq& pSeedSeq) :
    sandboxShape(std::move(pSandboxShape)),
    rndGenerator(pSeedSeq)
{
    if (sandboxShape.size() < 1)
        throw std::out_of_range("Dimension of the sandbox must be at least 1!");

    long long n = 1;
    for (short len : sandboxShape)
    {
        rndUnifVec.push_back(std::uniform_int_distribution<short>(0, len-1));
        n *= len;
    }

    rndUnif = std::uniform_int_distribution<long long>(0, n-1);
}

//

/*!
 * \brief Seed the random generator engine.
 *
 * Seeds the Mersenne Twister engine used for generating the random lattice sites (both randLin() and randVec()).
 *
 * \param pSeedSeq The seed sequence.
 */
void Randomizer::seed(std::seed_seq& pSeedSeq)
{
    rndGenerator.seed(pSeedSeq);
}

//

/*!
 * \brief Generate random lattice site (access via linear index).
 *
 * Generates a random index that can be used for accessing a random lattice site of a sandbox with linear access.
 *
 * \return Integer random number in [0, N) with N the total number of lattice sites of the sandbox.
 */
long long Randomizer::randLin()
{
    return rndUnif(rndGenerator);
}

/*!
 * \brief Generate random lattice site (access via cartesian vector).
 *
 * Generates a random vector that represents a random lattice site of a sandbox of shape as specified on construction.
 *
 * \return Random vector with integer components v_i in [0, n_i) with n_i the sandbox's size in i-th dimension.
 */
std::vector<short> Randomizer::randVec()
{
    std::vector<short> tmpRndVec;
    for (size_t i = 0; i < sandboxShape.size(); ++i)
        tmpRndVec.push_back(rndUnifVec[i](rndGenerator));

    return tmpRndVec;
}

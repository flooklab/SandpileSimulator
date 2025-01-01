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

#ifndef SANDSIM_RANDOMIZER_H
#define SANDSIM_RANDOMIZER_H

#include <random>
#include <vector>

/*!
 * \brief Generate random sandbox lattice sites.
 *
 * Generates discrete, uniformly distributed positions on a lattice of specified shape.
 * Hence this class can be used especially for choosing random sites in a sandbox.
 *
 * There are two different (equivalent*) ways to obtain a newly generated lattice site,
 * namely either in form of a "linear index" (integer) or as a vector (n-dim. integer-vector).
 * For detailed information see randLin() and randVec() respectively.
 *
 * *Be careful:
 * 'equivalent' means: both ways in principle provide you with equally well suited random results,
 * just in a different format; so you might choose the method that fits best to your needs.
 * Nevertheless randVec() would generally be slower than randLin().
 * Thus, depending on the dimension and the actual application, always generating linear indexes
 * with manual conversion to a vector afterwards may still be the favorable option!
 */
class Randomizer
{
public:
    Randomizer(std::vector<short> pSandboxShape, std::seed_seq& pSeedSequence);     ///< Constructor.
    //
    void seed(std::seed_seq& pSeedSeq);     ///< Seed the random generator engine.
    //
    long long randLin();                    ///< Generate random lattice site (access via linear index).
    std::vector<short> randVec();           ///< Generate random lattice site (access via cartesian vector).

private:
    std::vector<short> sandboxShape;    //Shape of the sandbox for which the random lattice sites are generated
    //
    std::mt19937_64 rndGenerator;       //Mersenne Twister pRNG engine
    //
    std::uniform_int_distribution<long long> rndUnif;               //Generates uniformly distributed, integer random numbers (scalars)
    std::vector<std::uniform_int_distribution<short>> rndUnifVec;   //Generates uniformly distributed, integer random numbers (vectors)
};

#endif // SANDSIM_RANDOMIZER_H

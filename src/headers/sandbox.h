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

#ifndef SANDBOX_H
#define SANDBOX_H

#include <vector>
#include <array>
#include <exception>
#include <sstream>
#include <fstream>
#include <random>

#include "aux.h"
#include "randomizer.h"

/*!
 * \brief Class providing an interface to an n-dimensional sandbox.
 *
 * A sandpile in n+1 dimensions is represented by an n-dimensional array
 * of lattice sites (the "sandbox"), where each lattice site consists of
 * an integer. These itegers stand for whatever the chosen SimulationModel
 * uses as the characteristic property of the sandpile regarding the remaining
 * (n+1)-th dimension (e.g. for a 3-dim. sandpile: number of sand grains stacked
 * on top of each other, the slope, etc...).
 *
 * Note: The sandbox dimension must be at least 1.
 */
class Sandbox
{
public:
    Sandbox(std::vector<short> pShape, std::seed_seq& pSeedSeq);    ///< Constructor.
    //
    void fill(short pLevel = 0);                ///< Fill the whole lattice up to the specified level.
    void clear();                               ///< Clear the whole lattice to level 0.
    bool load(const std::string& pFileName);    ///< Load a priorly saved sandpile into the sandbox.
    bool save(const std::string& pFileName);    ///< Save the current state of the sandpile to a file.
    void seed(std::seed_seq& pSeedSeq);         ///< Seed the internal random lattice site generator.
    //
    const std::vector<short>& getShape() const; ///< Get the sandbox shape.
    //
    std::vector<std::pair<bool, bool>> getBoundaryConditions() const;       ///< Get the boundary conditions.
    void setBoundaryConditions(std::vector<std::pair<bool, bool>> pBounds); ///< Set the boundary conditions.
    void setBoundaryConditions(bool pAllOpen);                              ///< Set all boundary conditions collectively.
    //
    const std::vector<short>& getCurrentGrain() const;  ///< Get the latest random lattice site.
    void dropRandom(short pNumGrains);                  ///< Add integer units (grains, slope, ...) to a random lattice site.
    //
    short getHeightAt(const std::vector<short>& pPos) const;    ///< Get the stored integer at the specified lattice site.
    short getHeightAt(long long pIdx) const;                    ///< Get the stored integer at the specified lattice site.
    short getHeightCurrent() const;                             ///< Get the stored integer at the latest random lattice site.
    short changeHeightAt(const std::vector<short>& pPos, short pDiff);  ///< Change the stored integer at the specified lattice site.
    short changeHeightAt(long long pIdx, short pDiff);                  ///< Change the stored integer at the specified lattice site.
    short changeHeightCurrent(short pOffset);                   ///< Change the stored integer at the latest random lattice site.
    short clearColumnAt(const std::vector<short>& pPos);        ///< Reset the stored integer to 0 at the specified lattice site.
    short clearColumnAt(long long pIdx);                        ///< Reset the stored integer to 0 at the specified lattice site.
    short clearColumnCurrent();                                 ///< Reset the stored integer to 0 at the latest random lattice site.
    //
    std::vector<long long> findSitesLargerEqualHeight(short pHeight) const; ///< Find all lattice sites with a minimum stored value.
    //
    char getBoundaryInfoAt(const std::vector<short>& pPos) const;   ///< \brief Check whether a point lies in the sandbox and, if so,
                                                                    ///  whether it is at an open boundary or at a closed boundary.
    char getBoundaryInfoCurrent() const;                            ///< \brief Check whether the latest random lattice site is at an
                                                                    ///  open boundary or at a closed boundary or at no boundary.
    //
    long long cartToLin(const std::vector<short>& pCartPosition) const;                 ///< \brief Convert Cartesian position vector
                                                                                        ///         into linear lattice site index.
    void linToCart(long long pLinPosition, std::vector<short>& pCartPosition) const;    ///< \brief Convert linear lattice site index
                                                                                        ///         into Cartesian position vector.
    //
    std::array<short, 2> get2DSliceShape() const;           ///< Get the shape of a (fixed) 2D subset ("2D slice") of the sandpile.
    void get2DSlice(std::vector<std::vector<short>>& p2DSlice) const;   ///< Get a current 2D subset ("2D slice") of the sandpile.

private:
    void init2DSlice();                 ///< Initialize the (fixed) position/placement of a 2D subset ("2D slice") of the sandpile.
    //
    short& columnCurrent();                                             ///< Access the latest random lattice site.
    const short& columnCurrent() const;                                 ///< Access the latest random lattice site as const reference.
    short& columnAt(long long pLinPos);                                 ///< Access the lattice site at specified linear index.
    const short& columnAt(long long pLinPos) const;         ///< Access the lattice site at specified linear index as const reference.
    short& columnAt(const std::vector<short>& pCartPos);                ///< Access the lattice site at specified Cartesian position.
    const short& columnAt(const std::vector<short>& pCartPos) const;    ///< \brief Access the lattice site at specified
                                                                        ///  Cartesian position as const reference.

private:
    std::vector<short> sBox;        //The actual sandbox, accessed via a linear index
    std::vector<short> boxShape;    //Shape of the sandbox
    //
    std::vector<std::pair<bool, bool>> boxBounds;   //Boundaries (upper/lower edge) for each dimension (true: open, false: closed)
    //
    long long boxVolume;                        //Number of box's lattice points
    std::vector<long long> partialBoxVolumes;   //For index conversion: store volumes of n-dim subsets of sandbox
    //
    std::array<short, 2> sliceShape;            //Shape of the 2D slice
    long long sliceProjection_dimOffset;        //Offset for > 2-dim. sandbox, defining where the 2D slice is placed
    //
    std::vector<short> currentGrainVec;         //Position referring to the "current" random lattice site (changed by dropRandom())
    long long currentGrainLinIdx;               //Linear index representing 'currentGrainVec'
    //
    Randomizer randomizer;                      //pRNG wrapper class for choosing random sandbox lattice sites
};

#endif // SANDBOX_H

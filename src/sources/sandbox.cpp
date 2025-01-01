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

#include "sandbox.h"

#include <exception>
#include <fstream>
#include <ios>
#include <new>
#include <sstream>
#include <stdexcept>

/*!
 * \brief Constructor.
 *
 * Creates an empty sandbox of the given shape.
 *
 * Memory is allocated according to the given shape and each lattice site is set to zero (which is meant by "empty" sandbox).
 *
 * An appropriate mapping for obtaining a 2D subset (the "2D slice") is
 * defined here --- mainly for plotting purposes --- (see init2DSlice() for details).
 * The 2D slice can be obtained by the get2DSlice() function, see also its documentation.
 * Note that the "2D slice" is also defined for a 1-dim. sandbox.
 *
 * Note: The sandbox dimension (i.e. \p pShape .size()) must be at least 1
 * and every element of \p pShape must be positive and non-zero.
 *
 * \param pShape Shape of the sandbox in units of lattice sites.
 * \param pSeedSeq Seed sequence used to seed the internal random lattice site generator.
 *
 * \throws std::out_of_range %Sandbox dimension is smaller than 1.
 * \throws std::out_of_range %Sandbox shape is smaller than 1 in some direction.
 * \throws std::runtime_error Not enough memory for requested %Sandbox shape.
 */
Sandbox::Sandbox(std::vector<short> pShape, std::seed_seq& pSeedSeq) :
    boxShape(std::move(pShape)),
    boxBounds(boxShape.size(), {false, false}),
    boxVolume(0),
    sliceShape({0, 0}),
    sliceProjection_dimOffset(0),
    currentGrainLinIdx(0),
    randomizer(boxShape, pSeedSeq)
{
    if (boxShape.size() < 1)
        throw std::out_of_range("Dimension of the sandbox must be at least 1!");

    for (short len : boxShape)
        if (len <= 0)
            throw std::out_of_range("Extension of the sandbox must be at least 1 in each direction!");

    //Calculate volumes of n-dim subsets of sandbox needed for index conversion

    long long n = 1;        //Nicer: C++20 init-statement, long long n at loop scope
    for (short len : boxShape)
    {
        partialBoxVolumes.push_back(n);
        n *= len;
    }
    boxVolume = partialBoxVolumes.back() * boxShape.back();

    //Set the random lattice site to the origin initially
    currentGrainVec.assign(boxShape.size(), 0);

    //Calculate index offset for the 2D slice
    init2DSlice();

    try
    {
        //Fill box with zeros
        fill(0);
    }
    catch (std::bad_alloc&)     //Allocation failed
    {
        std::ostringstream stream;
        stream.precision(2);

        stream<<"Sandbox of shape {"<<boxShape[0];

        for (size_t i = 1; i < boxShape.size(); ++i)
            stream<<", "<<boxShape[i];

        stream<<"} too large. Not enough memory! (Needs ";
        stream<<boxVolume*sizeof(short)/1024./1024./1024.;
        stream<<"GiB)";

        throw std::runtime_error(stream.str());
    }
}

//Public

/*!
 * \brief Fill the whole lattice up to the specified level.
 *
 * Assigns the specified value to every lattice site.
 *
 * \param pLevel New value for every lattice site.
 */
void Sandbox::fill(const short pLevel)
{
    sBox.assign(boxVolume, pLevel);
}

/*!
 * \brief Clear the whole lattice to level 0.
 *
 * Resets each lattice site to value of 0.
 */
void Sandbox::clear()
{
    fill(0);
}

/*!
 * \brief Load a priorly saved sandpile into the sandbox.
 *
 * Reads (and sets) the height values for every lattice site from the sandbox file \p pFileName.
 * The file format must be according to save(). The Sandbox dimension/shape must exactly match that of the saved sandbox.
 *
 * \param pFileName File name for reading the sandbox file.
 * \return If successful.
 */
bool Sandbox::load(const std::string& pFileName)
{
    try
    {
        std::ifstream stream;
        stream.exceptions(std::ios_base::badbit | std::ios_base::failbit);
        stream.open(pFileName);

        std::string signature;
        if (!std::getline(stream, signature, ',') || signature != "SandpileSimulatorSandboxFile")
            return false;

        short tDim;
        std::string dimStr;
        if (!std::getline(stream, dimStr, ','))
            return false;

        try
        {
            tDim = std::stoi(dimStr);
        }
        catch (const std::exception&)
        {
            return false;
        }

        if (static_cast<size_t>(tDim) != boxShape.size())
            return false;

        std::vector<short> tShape;
        std::string shapeStr;
        for (short d = 0; d < tDim; ++d)
        {
            if (!std::getline(stream, shapeStr, ','))
                return false;

            try
            {
                tShape.push_back(std::stoi(shapeStr));
            }
            catch (const std::exception&)
            {
                return false;
            }
        }

        for (short d = 0; d < tDim; ++d)
            if (tShape.at(d) != boxShape.at(d))
                return false;

        std::vector<short> tBox;
        tBox.reserve(sBox.size());

        for (std::string substr; stream.good() && std::getline(stream, substr, ',');)
        {
            try
            {
                tBox.push_back(std::stoi(substr));
            }
            catch (const std::exception&)
            {
                return false;
            }
        }

        if (tBox.size() != sBox.size())
            return false;

        sBox.swap(tBox);

        stream.close();
    }
    catch (const std::ios_base::failure&)
    {
        return false;
    }

    return true;
}

/*!
 * \brief Save the current state of the sandpile to a file.
 *
 * Writes the height values for every lattice site to the sandbox file \p pFileName.
 *
 * See also load().
 *
 * \param pFileName File name for writing the sandbox file.
 * \return If successful.
 */
bool Sandbox::save(const std::string& pFileName)
{
    try
    {
        std::ofstream stream;
        stream.exceptions(std::ios_base::badbit | std::ios_base::failbit);
        stream.open(pFileName);

        stream<<"SandpileSimulatorSandboxFile"<<",";

        stream<<boxShape.size();

        for (short len : boxShape)
            stream<<","<<len;

        for (short val : sBox)
            stream<<","<<val;

        stream.close();
    }
    catch (const std::ios_base::failure&)
    {
        return false;
    }

    return true;
}

/*!
 * \brief Seed the internal random lattice site generator.
 *
 * Seeds the Randomizer used for generating random lattice sites.
 *
 * \param pSeedSeq The seed sequence.
 */
void Sandbox::seed(std::seed_seq& pSeedSeq)
{
    randomizer.seed(pSeedSeq);
}

//

/*!
 * \brief Get the sandbox shape.
 *
 * \return Const reference to the shape of the sandbox in units of lattice sites.
 */
const std::vector<short>& Sandbox::getShape() const
{
    return boxShape;
}

/*!
 * \brief Get the boundary conditions.
 *
 * See setBoundaryConditions(std::vector<std::pair<bool, bool>> pBounds) for the definition of these boundary conditions.
 *
 * \return The current boundary conditions.
 */
std::vector<std::pair<bool, bool>> Sandbox::getBoundaryConditions() const
{
    return boxBounds;
}

/*!
 * \brief Set the boundary conditions.
 *
 * For each dimension of the sandbox there exist two edges, a "lower" edge at i=0 and an "upper" edge at i=length-1.
 * You can set boundary conditions for each edge of each dimension independently, either to \e open or to \e closed
 * boundary conditions. What effects are caused by either open or closed boundaries is up to the used SimulationModel.
 *
 * For an n-dim. sandbox you have to provide n pairs of bools, where the first bool stands for the lower edge and the
 * second bool stands for the upper edge. True represents open boundaries and False represents closed boundaries.
 * The i-th pair of \p pBounds corresponds to the i-th dimension of the sandbox shape (getShape()).
 * If the size of \p pBounds does not match the sandbox dimension the boundary conditions will remain unchanged!
 *
 * Note: By default the constructor sets all boundaries to closed.
 *
 * \param pBounds The desired boundary conditions.
 */
void Sandbox::setBoundaryConditions(std::vector<std::pair<bool, bool>> pBounds)
{
    if (pBounds.size() == boxShape.size())
        boxBounds.swap(pBounds);
}

/*!
 * \brief Set all boundary conditions collectively.
 *
 * Sets either open boundary conditions for every edge or closed boundary conditions for every edge.
 * See also setBoundaryConditions(std::vector<std::pair<bool, bool>> pBounds).
 *
 * \param pAllOpen True for all boundaries open, False for all boundaries closed.
 */
void Sandbox::setBoundaryConditions(const bool pAllOpen)
{
    setBoundaryConditions(std::vector<std::pair<bool, bool>>(boxShape.size(), {pAllOpen, pAllOpen}));
}

//

/*!
 * \brief Get the latest random lattice site.
 *
 * Note: initial value after construction of this class is the origin.
 *
 * \return Const reference to the random lattice site that was selected most recently by dropRandom().
 */
const std::vector<short>& Sandbox::getCurrentGrain() const
{
    return currentGrainVec;
}

/*!
 * \brief Add integer units (grains, slope, ...) to a random lattice site.
 *
 * A random lattice site is selected and the specified number is added to the existing integer that constitutes the lattice site.
 *
 * \param pNumGrains Amount to be added to the content of the random lattice site.
 */
void Sandbox::dropRandom(const short pNumGrains)
{
    currentGrainLinIdx = randomizer.randLin();
    linToCart(currentGrainLinIdx, currentGrainVec);

    columnCurrent() += pNumGrains;
}

//

/*!
 * \brief Get the stored integer at the specified lattice site.
 *
 * \param pPos Position of the lattice site as n-dim. vector.
 * \return The integer stored at the specified lattice site.
 */
short Sandbox::getHeightAt(const std::vector<short>& pPos) const
{
    return columnAt(pPos);
}

/*!
 * \brief Get the stored integer at the specified lattice site.
 *
 * \param pIdx Position of the lattice site as linear index as obtained from cartToLin().
 * \return The integer stored at the specified lattice site.
 */
short Sandbox::getHeightAt(const long long pIdx) const
{
    return columnAt(pIdx);
}

/*!
 * \brief Get the stored integer at the latest random lattice site.
 *
 * \return Content of the random lattice site that was selected most recently by dropRandom().
 */
short Sandbox::getHeightCurrent() const
{
    return columnCurrent();
}

/*!
 * \brief Change the stored integer at the specified lattice site.
 *
 * \param pPos Position of the lattice site as n-dim. vector.
 * \param pDiff Amount to be added to the specified lattice site.
 * \return The new value (after addition of \p pDiff).
 */
short Sandbox::changeHeightAt(const std::vector<short>& pPos, const short pDiff)
{
    return (columnAt(pPos) += pDiff);
}

/*!
 * \brief Change the stored integer at the specified lattice site.
 *
 * \param pIdx Position of the lattice site as linear index as obtained from cartToLin().
 * \param pDiff Amount to be added to the specified lattice site.
 * \return The new value (after addition of \p pDiff).
 */
short Sandbox::changeHeightAt(const long long pIdx, const short pDiff)
{
    return (columnAt(pIdx) += pDiff);
}

/*!
 * \brief Change the stored integer at the latest random lattice site.
 *
 * The latest random lattice site is the lattice site that was selected most recently by dropRandom().
 *
 * \param pOffset Amount to be added to the current lattice site.
 * \return The new value (after addition of offset).
 */
short Sandbox::changeHeightCurrent(const short pOffset)
{
    return (columnCurrent() += pOffset);
}

/*!
 * \brief Reset the stored integer to 0 at the specified lattice site.
 *
 * \param pPos Position of the lattice site as n-dim. vector.
 * \return The old value (before resetting).
 */
short Sandbox::clearColumnAt(const std::vector<short>& pPos)
{
    short tmp = columnAt(pPos);
    columnAt(pPos) = 0;
    return tmp;
}

/*!
 * \brief Reset the stored integer to 0 at the specified lattice site.
 *
 * \param pIdx Position of the lattice site as linear index as obtained from cartToLin().
 * \return The old value (before resetting).
 */
short Sandbox::clearColumnAt(const long long pIdx)
{
    short tmp = columnAt(pIdx);
    columnAt(pIdx) = 0;
    return tmp;
}

/*!
 * \brief Reset the stored integer to 0 at the latest random lattice site.
 *
 * The latest random lattice site is the lattice site that was selected most recently by dropRandom().
 *
 * \return The old value (before resetting).
 */
short Sandbox::clearColumnCurrent()
{
    short tmp = columnCurrent();
    columnCurrent() = 0;
    return tmp;
}

//

/*!
 * \brief Find all lattice sites with a minimum stored value.
 *
 * Searches the lattice for sites with stored integer value larger than or equal to \p pHeight.
 * A vector of the linear indexes of these sites is returned.
 *
 * \return Linear indexes of matching lattice sites.
 */
std::vector<long long> Sandbox::findSitesLargerEqualHeight(const short pHeight) const
{
    std::vector<long long> sites;

    for (long long siteIdx = 0; siteIdx < static_cast<long long>(sBox.size()); ++siteIdx)
        if (sBox[siteIdx] >= pHeight)
            sites.push_back(siteIdx);

    return sites;
}

//

/*!
 * \brief Check whether a point lies in the sandbox and, if so, whether it is at an open boundary or at a closed boundary.
 *
 * Returns 'n' for \p pPos inside the sandbox, if \p pPos is not directly at a boundary of the sandbox.
 * If \p pPos is inside the sandbox and at least at one of the closed boundaries but at \b none of the open boundaries,
 * this function returns 'c'. It returns 'o' instead, if \p pPos is inside the sandbox and at least at one of the open boundaries.
 * For \p pPos outside the sandbox, 'x' is returned.
 *
 * See also getBoundaryConditions() for more information about the boundary conditions.
 *
 * \param pPos Position vector that shall be checked.
 * \return 'n', except when at least at one closed boundary (then returns 'c'),
 *              except when at least at one open boundary (then returns 'o'),
 *              except when beyond any boundary (then returns 'x').
 */
char Sandbox::getBoundaryInfoAt(const std::vector<short>& pPos) const
{
    //Assume: "normal" lattice site, not at a boundary or even beyond a boundary
    char retVal = 'n';

    //Change return value, if at least at one boundary, but *directly* return, if beyond any boundary

    for (size_t i = 0; i < boxShape.size(); ++i)
    {
        if (pPos[i] < 0)
            return 'x';
        else if (pPos[i] == 0)
        {
            if (retVal != 'o' && !boxBounds[i].first)   //An open boundary always "overwrites" closed boundaries
                retVal = 'c';
            else
                retVal = 'o';
        }
        else if (pPos[i] >= boxShape[i])
            return 'x';
        else if (pPos[i] == boxShape[i]-1)
        {
            if (retVal != 'o' && !boxBounds[i].second)  //An open boundary always "overwrites" closed boundaries
                retVal = 'c';
            else
                retVal = 'o';
        }
    }

    return retVal;
}

/*!
 * \brief Check whether the latest random lattice site is at an open boundary or at a closed boundary or at no boundary.
 *
 * Returns 'n' if the latest random lattice site (LRLS) is not directly at a boundary of the sandbox.
 * If LRLS is at least at one of the closed boundaries but at \b none of the open boundaries,
 * this function returns 'c'. It returns 'o' instead, if LRLS is at least at one of the open boundaries.
 *
 * See also getBoundaryConditions() for more information about the boundary conditions.
 *
 * Note: The LRLS is the lattice site that was selected most recently by dropRandom().
 *
 * \return 'n', except when at least at one closed boundary (then returns 'c'),
 *              except when at least at one open boundary (then returns 'o').
 */
char Sandbox::getBoundaryInfoCurrent() const
{
    //Assume: "normal" lattice site, not at a boundary or even beyond a boundary
    char retVal = 'n';

    //Change return value, if at least at one boundary.
    //Current grain should never be beyond any boundary, so don't check that.

    for (size_t i = 0; i < boxShape.size(); ++i)
    {
        if (currentGrainVec[i] == 0)
        {
            if (retVal != 'o' && !boxBounds[i].first)   //An open boundary always "overwrites" closed boundaries
                retVal = 'c';
            else
                retVal = 'o';
        }
        else if (currentGrainVec[i] == boxShape[i]-1)
        {
            if (retVal != 'o' && !boxBounds[i].second)  //An open boundary always "overwrites" closed boundaries
                retVal = 'c';
            else
                retVal = 'o';
        }
    }

    return retVal;
}

//

/*!
 * \brief Convert Cartesian position vector into linear lattice site index.
 *
 * This is a unique mapping for all coordinates within the bounds (getShape()) of the sandbox.
 * For the inverse mapping see linToCart().
 *
 * \param pCartPosition Cartesian position vector of a site on the sandbox lattice.
 * \return The corresponding linear index of the site.
 */
long long Sandbox::cartToLin(const std::vector<short>& pCartPosition) const
{
    long long retVal = 0;
    for (size_t i = 0; i < boxShape.size(); ++i)
        retVal += pCartPosition[i]*partialBoxVolumes[i];

    return retVal;
}


/*!
 * \brief Convert linear lattice site index into Cartesian position vector.
 *
 * This is the inverse mapping of cartToLin().
 *
 * \param pLinPosition Linear index of a site on the sandbox lattice.
 * \param pCartPosition Destination for the converted position.
 *                      The specified vector is assumed to have (at least) one element per sandbox dimension available already!
 *
 * \attention \p pCartPosition is assumed to have one element per sandbox dimension available already!
 */
void Sandbox::linToCart(long long pLinPosition, std::vector<short>& pCartPosition) const
{
    std::vector<long long>::const_iterator tPBoxVol_it = partialBoxVolumes.begin();
    for (int i = boxShape.size()-1; i >= 0; --i)
    {
        pCartPosition[i] = (pLinPosition / *(tPBoxVol_it + i));
        pLinPosition     =  pLinPosition % *(tPBoxVol_it + i);
    }
}

//

/*!
 * \brief Get the shape of a (fixed) 2D subset ("2D slice") of the sandpile.
 *
 * With the construction of this class a fixed 2D slice, a subset of the whole sandbox,
 * is defined (see init2DSlice() for details about and the purpose of this 2D slice).
 * Use this function to obtain the shape of the 2D slice. The slice itself is provided by get2DSlice().
 * Note that such a 2D slice is also defined for a 1-dim. sandbox (see init2DSlice()).
 *
 * \return Dimensions of the 2D slice.
 */
std::array<short, 2> Sandbox::get2DSliceShape() const
{
    return sliceShape;
}

/*!
 * \brief Get a current 2D subset ("2D slice") of the sandpile.
 *
 * The given 2-dim. array is assigned a 2D subset of the current sandpile state.
 * The 2D subset corresponds to the 2D slice defined in init2DSlice()
 * and must thus have the same shape as returned by get2DSliceShape().
 * Note that such a "2D subset" is also defined for a 1-dim. sandbox (see init2DSlice()).
 *
 * Note: \p p2DSlice is expected to be pre-filled according to that shape.
 *
 * \param p2DSlice Reference to where the 2D slice shall be written to.
 */
void Sandbox::get2DSlice(std::vector<std::vector<short>>& p2DSlice) const
{
    if (boxShape.size() == 1)
        p2DSlice[0] = sBox;
    else
    {
        for (short i = 0; i < boxShape[1]; ++i)
        {
            p2DSlice[i].assign(sBox.begin()+i*boxShape[0]+sliceProjection_dimOffset,
                               sBox.begin()+(i+1)*boxShape[0]+sliceProjection_dimOffset);
        }
    }
}

//Private

/*!
 * \brief Initialize the (fixed) position/placement of a 2D subset ("2D slice") of the sandpile.
 *
 * If the sandbox dimension is larger than or equal to 2, an appropriate mapping for obtaining a 2D subset (the "2D slice")
 * is defined here. For a 2-dim. sandbox the mapping is trivial, it simply corresponds to the sandbox itself. For higher
 * dimensions n the 2D slice is a 2-dim. slice, which is aligned at the middle point of each of the other n-2 dimensions.
 *
 * If, however, the sandbox dimension is equal to 1, the mapping is defined
 * as if the dimension was 2 with a "virtual" slice shape of {1, shape[0]}.
 *
 * The 2D slice is useful for plotting of the sandpile and can be obtained
 * by the get2DSlice() function, which relies on the mapping chosen here.
 */
void Sandbox::init2DSlice()
{
    sliceShape = {0, 0};

    if (boxShape.size() == 1)
        sliceShape = {1, boxShape[0]};
    else
    {
        //Select some 2D slice "in the middle" for higher dimensions via calculating an offset of sBox's index
        sliceProjection_dimOffset = 0;
        if (boxShape.size() > 2)
        {
            std::vector<short> tSlicePos = {0, 0};
            for (size_t i = 2; i < boxShape.size(); ++i)
                tSlicePos.push_back(boxShape[i]/2.);
            sliceProjection_dimOffset = cartToLin(tSlicePos);
        }

        sliceShape = {boxShape[1], boxShape[0]};
    }
}

//

/*!
 * \brief Access the latest random lattice site.
 *
 * The latest random lattice site is the lattice site that was selected most recently by dropRandom().
 *
 * \return Reference to the latest random lattice site.
 */
short& Sandbox::columnCurrent()
{
    return columnAt(currentGrainLinIdx);
}

/*!
 * \brief Access the latest random lattice site as const reference.
 *
 * Const version of columnCurrent().
 *
 * \return Const reference to the latest random lattice site.
 */
const short& Sandbox::columnCurrent() const
{
    return columnAt(currentGrainLinIdx);
}

/*!
 * \brief Access the lattice site at specified linear index.
 *
 * Directly returns the lattice site content via the lattice site container's linear index (see also cartToLin() / linToCart()).
 *
 * \param pLinPos Linear index of a site on the sandbox lattice.
 * \return Reference to the lattice site at pLinPos.
 */
short& Sandbox::columnAt(const long long pLinPos)
{
    return sBox[pLinPos];
}

/*!
 * \brief Access the lattice site at specified linear index as const reference.
 *
 * Const version of columnAt(long long pLinPos).
 *
 * \param pLinPos Linear index of a site on the sandbox lattice.
 * \return Const reference to the lattice site at pLinPos.
 */
const short& Sandbox::columnAt(const long long pLinPos) const
{
    return sBox[pLinPos];
}

/*!
 * \brief Access the lattice site at specified Cartesian position.
 *
 * First performs conversion of Cartesian position into linear lattice site index.
 * (Slower than columnAt(long long pLinPos)! See cartToLin() for details.)
 *
 * \param pCartPos Cartesian position vector of a site on the sandbox lattice.
 * \return Reference to the lattice site at pCartPos.
 */
short& Sandbox::columnAt(const std::vector<short>& pCartPos)
{
    return sBox[cartToLin(pCartPos)];
}

/*!
 * \brief Access the lattice site at specified Cartesian position as const reference.
 *
 * Const version of columnAt(const std::vector<short>& pCartPos).
 *
 * \param pCartPos Cartesian position vector of a site on the sandbox lattice.
 * \return Const reference to the lattice site at pCartPos.
 */
const short& Sandbox::columnAt(const std::vector<short>& pCartPos) const
{
    return sBox[cartToLin(pCartPos)];
}

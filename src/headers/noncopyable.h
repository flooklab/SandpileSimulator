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

#ifndef NONCOPYABLE_H
#define NONCOPYABLE_H

/*!
 * \brief Generic base class for creating non-copyable classes.
 *
 * When a derived class inherits from this base class,
 * copy construction of and assigning to the derived
 * class will be disabled due to deleted copy constructor
 * and assignment operator of the base class.
 */
class NonCopyable
{
protected:
    NonCopyable() = default;                                ///< Default constructor.
    ~NonCopyable() = default;                               ///< Default destructor.
    //
    NonCopyable(const NonCopyable&) = delete;               ///< Deleted copy constructor.
    NonCopyable& operator=(const NonCopyable&) = delete;    ///< Deleted assignment operator.
};

#endif // NONCOPYABLE_H

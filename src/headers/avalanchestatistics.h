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

#ifndef AVALANCHESTATISTICS_H
#define AVALANCHESTATISTICS_H

#include <vector>

#include "aux.h"

/*!
 * \brief Collect statistics of sandpile avalanches.
 *
 * In order to investigate and characterize sandpile dynamics using cellular automata
 * one typically triggers avalanches on those sandpiles and studies the avalanches'
 * characteristic properties (e.g. how wide does the avalanche spread, etc.).
 * More precisely, for each of the properties, the corresponding distribution function
 * (i.e. how often does property X take the value of x,
 * averaged over multiple avalanches) must be measured.
 *
 * This class provides an interface for collecting all the needed information.
 * For this purpose the Event struct defines a common classification of avalanche properties,
 * such that each avalanche can be represented by an instance of Event.
 * All the events are then gathered using this class.
 */
class AvalancheStatistics
{
public:
    struct Event;   //Forward declaration

public:
    AvalancheStatistics();          ///< Default constructor.
    //
    Event& insert();                ///< Create a new event and include it in the statistics.
    void insert(Event pEvent);      ///< Include the event in the statistics.
    //
    void detachEvents(std::vector<Event>& pDestination);    ///< Detach the collection of inserted events from this class's instance.

private:
    std::vector<Event> events;  //Vector of all inserted events

public:
    /*!
     * \brief Characterize an avalanche event.
     *
     * Note: Depending on the used simulation model, the exact definition of the observables
     * listed here may be slightly different. The descriptions provided here are only rough guidelines.
     */
    struct Event
    {
        long long _id = -1;     ///< (Consecutive) %Event ID.
        long long size = -1;    ///< \brief Abstract "size" of the avalanche.
                                ///  \details Number of critical lattice sites accumulated throughout
                                ///           all iterations of the relaxation procedure.
        double linSize = -1;    ///< \brief Diameter of the avalanche.
                                ///  \details Maximum distance between any of the lattice sites that took part in the avalanche
                                ///           (i.e. were modified at least once during the avalanche).
        long long area = -1;    ///< \brief Lattice area occupied by avalanche.
                                ///  \details The dim[sandbox]-dimensional volume that contains the avalanche,
                                ///           i.e. the number of distinct lattice sites that took part in the avalanche
                                ///           (i.e. were modified at least once during the avalanche).
        int duration = -1;      ///< \brief Time duration of the avalanche.
                                ///  \details Number of unit time steps (relaxation iterations) from the
                                ///           triggering of the avalanche until the avalanche stops.
    };
    /*!
     * \brief Characterize moments of distributions of event observables.
     *
     * After collecting many events, the distributions of the observables listed in Event
     * may be characterized by their \e moments. Use this struct to store a specific moment
     * (e.g. first moment) for each of the available observables. Different from Event,
     * every member of this struct is of floating point type such that one can also
     * store the logarithmized values of the moments, which is especially needed
     * for the moment analysis method (see MomentAnalysis).
     */
    struct Moments
    {
        double size = -1;       ///< \copybrief Event::size
                                ///  \copydetails Event::size
        double linSize = -1;    ///< \copybrief Event::linSize
                                ///  \copydetails Event::linSize
        double area = -1;       ///< \copybrief Event::area
                                ///  \copydetails Event::area
        double duration = -1;   ///< \copybrief Event::duration
                                ///  \copydetails Event::duration
    };
};

#endif // AVALANCHESTATISTICS_H

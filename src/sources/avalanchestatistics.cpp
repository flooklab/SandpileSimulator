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

#include "avalanchestatistics.h"

/*!
 * \brief Default constructor.
 */
AvalancheStatistics::AvalancheStatistics()
{
}

//

/*!
 * \brief Create a new event and include it in the statistics.
 *
 * Creates an %Event instance and adds it to the internal vector
 * that stores all the events included in the avalanche statistics
 * (see also insert(Event&)).
 *
 * The event's properties are modified using the returned reference.
 * The Event::_id of the new event is initially set to the ID of the lastly added event + 1.
 * If no event has been added before, the ID is set to 0.
 *
 * \return Reference to the newly created event.
 */
AvalancheStatistics::Event& AvalancheStatistics::insert()
{
    Event tEvent {};

    if (events.empty())
        tEvent._id = 0;
    else
        tEvent._id = events.back()._id + 1;

    events.push_back(tEvent);

    return events.back();
}

/*!
 * \brief Include the event in the statistics.
 *
 * Adds \p pEvent to the vector that stores all the events included in the avalanche statistics (see also insert()).
 *
 * The Event::_id of the added event is overwritten with the ID of the lastly added event + 1.
 * If no event has been added before, the ID is set to 0.
 *
 * \param pEvent The event to be inserted.
 */
void AvalancheStatistics::insert(Event pEvent)
{
    if (events.empty())
        pEvent._id = 0;
    else
        pEvent._id = events.back()._id + 1;

    events.push_back(pEvent);
}

/*!
 * \brief Detach the collection of inserted events from this class's instance.
 *
 * Moves all the previously inserted events to the given destination.
 * The internal list of events is then empty again.
 *
 * \param pDestination Destination for the detached events.
 */
void AvalancheStatistics::detachEvents(std::vector<Event>& pDestination)
{
    events.swap(pDestination);
    events.clear();
}

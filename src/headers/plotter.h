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

#ifndef PLOTTER_H
#define PLOTTER_H

#include <atomic>
#include <thread>
#include <functional>
#include <sstream>

#include <SFML/Graphics.hpp>

#include "noncopyable.h"
#include "aux.h"
#include "sandbox.h"

#define FONT_DEJAVU_DIR "/usr/share/fonts/TTF"

/*!
 * \brief Visualize a sandpile during simulation.
 *
 * While the simulation of a sandpile is running, its evolution can
 * be visualized for instance by plotting a heat map of the
 * lattice sites of the sandbox (or of a 2D subset of it).
 *
 * You can use this class to start another thread for plotting an n-dim. sandbox
 * simultaneously with the simulation. Here, a 2D heat map is plotted and for doing
 * so the 2D field of plot data is obtained using the function Sandbox::get2DSlice()
 * (see also updatePlotData()). Note: A 1-dim. sandbox will be plotted from the side instead.
 *
 * The plotting only starts after you've called init().
 * The plotting thread can also be stopped again via stopPlotting().
 *
 * The updating/plotting is \b not synchronized with the simulation and hence
 * the achievable frame rate does not (directly) limit the simulation speed, which is good.
 * But if you for instance wanted to carefully watch the avalanche evolution
 * you would have to actively slow down the simulation itself because
 * reducing the plotting frame rate would not slow down the simulation.
 */
class Plotter : private NonCopyable
{
public:
    Plotter(const Sandbox& pSandbox);           ///< Constructor.
    ~Plotter();                                 ///< Destructor.
    //
    void init(unsigned int pFramerateLimit);    ///< Start a thread and setup a window for plotting the sandpile evolution.
    void stopPlotting();                        ///< Close the plotting window and terminate the plotting thread.

private:
    void updatePlotData();                      ///< Update the plotted data, the 2D sandbox slice, with the current sandbox content.
    void plotLoop();                            ///< Plotting thread function. Perpetually updates and plots the sandbox slice.

private:
    const Sandbox& sandbox;                     //Sandbox whose content shall be plotted
    //
    std::thread plotThread;                     //The thread that continuously renders the sandbox
    std::atomic<bool> plotThreadRunning;        //Switch to control (i.e. terminate) the plotting thread
    std::unique_ptr<sf::RenderWindow> window;   //SFML window, in which the rendered sandbox is displayed
    unsigned int framerateLimit;                //Maximal allowed frame rate
    bool oneDimPlot;                            //Show a 1-dim. sandbox from the side instead of a heatmap of the 2D slice from the top
    //
    std::vector<std::vector<short>> plotData2d; //2D field containing those lattice sites of the sandbox that shall be plotted
    //
    sf::RectangleShape shape;                   //Simple rectangle object used to render/display each lattice site
};

#endif // PLOTTER_H

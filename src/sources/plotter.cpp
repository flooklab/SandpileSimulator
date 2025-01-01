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

#include "plotter.h"

#include <array>
#include <iostream>
#include <new>
#include <sstream>
#include <stdexcept>
#include <string>

/*!
 * \brief Constructor.
 *
 * Constructs a new %Plotter instance. This does not yet plot anything.
 * To open a window for plotting of the sandpile and starting the corresponding thread call init().
 * The plotting thread can be stopped again by calling stopPlotting().
 *
 * Note: After having stopped plotting you can of course always reactivate the thread by calling init() once again.
 *
 * \param pSandbox The sandbox containing the sandpile to be plotted.
 */
Plotter::Plotter(const Sandbox& pSandbox) :
    sandbox(pSandbox),
    plotThreadRunning(false),
    framerateLimit(30),
    oneDimPlot(sandbox.getShape().size() == 1)
{
}

/*!
 * \brief Destructor.
 *
 * Closes the plotting window and terminates the plotting thread, if still running.
 * See also stopPlotting().
 */
Plotter::~Plotter()
{
    stopPlotting();
}

//Public

/*!
 * \brief Start a thread and setup a window for plotting the sandpile evolution.
 *
 * Starts a separate thread that plots the evolution of the sandpile in a newly opened
 * render window (using the sf::RenderWindow class from the SFML library).
 * Only one running plotting thread at a time is possible.
 *
 * Memory for storing a perpetually updated 2D slice of the sandbox is allocated.
 * See also updatePlotData().
 *
 * Maximal time resolution (i.e. if/how well individual driving steps can be resolved on the screen) is limited
 * by the achievable frame rate of the render window, which can be further limited by the additional parameter.
 *
 * If a plotting thread is already running this function will just update its framerate limit.
 *
 * \param pFramerateLimit Maximal update frequency of displayed sandbox content in fps. Set to 0 to disable this limit.
 *
 * \throws std::runtime_error Not enough memory for required 2D slice shape.
 */
void Plotter::init(const unsigned int pFramerateLimit)
{
    framerateLimit = pFramerateLimit;

    if (plotThread.joinable())
    {
        if (plotThreadRunning)
        {
            if (window->isOpen())
                window->setFramerateLimit(framerateLimit);

            return;
        }
        else
            plotThread.join();
    }

    if (!window)
        window = std::make_unique<sf::RenderWindow>();

    if (!window->isOpen())
        window->create(sf::VideoMode(800, 600), sf::String("Sandpile Simulator - Sandbox Plot"));

    window->setFramerateLimit(framerateLimit);
    window->setVerticalSyncEnabled(false);

    //Initialize 2D projection of Sandbox
    std::array<short, 2> sliceShape = sandbox.get2DSliceShape();

    //Allocate memory according to slice shape
    std::vector<short> tmpVec;
    tmpVec.assign(sliceShape[1], 0);

    try
    {
        plotData2d.assign(sliceShape[0], tmpVec);
    }
    catch (std::bad_alloc&)     //Allocation failed
    {
        std::ostringstream stream;
        stream.precision(2);

        stream<<"2D slice of shape {"<<sliceShape[0];

        for (size_t i = 1; i < sliceShape.size(); ++i)
            stream<<", "<<sliceShape[i];

        stream<<"} too large. Not enough memory! (Needs ";
        stream<<sliceShape[0]*sliceShape[1]*sizeof(short)/1024./1024./1024.;
        stream<<"GiB)";

        throw std::runtime_error(stream.str());
    }

    plotThreadRunning = true;
    plotThread = std::thread(&Plotter::plotLoop, this);
}

/*!
 * \brief Close the plotting window and terminate the plotting thread.
 *
 * If a plotting thread is currently active the thread is terminated and the corresponding render window will be closed afterwards.
 */
void Plotter::stopPlotting()
{
    plotThreadRunning = false;
    if (plotThread.joinable())
        plotThread.join();

    if (window != nullptr)
        window->close();
}

//Private

/*!
 * \brief Update the plotted data, the 2D sandbox slice, with the current sandbox content.
 *
 * Obtains a current 2D projection ("slice") of the sandpile, which is to be be displayed in the render window, from the sandbox.
 * See Sandbox::get2DSlice() and Sandbox::init2DSlice() for details about how this 2D projection is defined.
 */
void Plotter::updatePlotData()
{
    sandbox.get2DSlice(plotData2d);
}

/*!
 * \brief Plotting thread function. Perpetually updates and plots the sandbox slice.
 *
 * Enters a loop that manages the user interaction with the render window, updating and rendering of the sandbox slice.
 * The loop is exited when the render window is closed by user interaction or if the function stopPlotting() is called.
 *
 * Note: Assumes that the \e DejaVu \e Sans font ("DejaVuSans.ttf") is installed on the system.
 * The font directory is defined as a macro in the class's header file.
 */
void Plotter::plotLoop()
{
    short upperHeightLimit = 1;
    short lowerHeightLimit = 0;

    if (oneDimPlot)
        shape.setFillColor(sf::Color(255., 255., 255.));
    else
        shape.setSize(sf::Vector2f(2, 2));

    std::string fontFileName = std::string(SANDSIM_FONT_DEJAVU_DIR) + "/DejaVuSans.ttf";

    sf::Font fnt;
    if (!fnt.loadFromFile(fontFileName))
        std::cerr<<"WARNING: Could not load font from file \""<<fontFileName<<"\"!\n";
    sf::Text txt(sf::String(), fnt);

    sf::Event event;

    sf::Clock clk;
    double frRt = 0, frRt2 = 0, frRt3 = 0, frRt4 = 0;

    while (plotThreadRunning && window->isOpen())
    {
        while(window->pollEvent(event))
        {
            switch (event.type)
            {
                case sf::Event::Closed:
                {
                    window->close();
                    plotThreadRunning = false;
                    break;
                }
                default:
                {
                    break;
                }
            }
        }

        updatePlotData();

        window->clear();

        if (oneDimPlot)
        {
            short maxHeight = 0;
            for (size_t j = 0; j < plotData2d[0].size(); ++j)
                if (plotData2d[0][j] > maxHeight)
                    maxHeight = plotData2d[0][j];

            double scaleFactor = 1.;
            if (maxHeight > 450)
                scaleFactor = 450. / maxHeight;

            for (size_t j = 0; j < plotData2d[0].size(); ++j)
            {
                for (short height = 0; height < plotData2d[0][j]; ++height)
                {
                    shape.setPosition(50+j, 550-height*scaleFactor);
                    shape.setSize(sf::Vector2f(1, height*scaleFactor));
                    window->draw(shape);
                }
            }
        }
        else
        {
            short currentHeight = 0;
            short tUpperHeightLimit = upperHeightLimit;
            short tLowerHeightLimit = tUpperHeightLimit;

            for (size_t i = 0; i < plotData2d.size(); ++i)
            {
                for (size_t j = 0; j < plotData2d[0].size(); ++j)
                {
                    currentHeight = plotData2d[i][j];

                    if (currentHeight > tUpperHeightLimit)
                    {
                        tUpperHeightLimit = currentHeight;

                        if (currentHeight > upperHeightLimit)
                            currentHeight = upperHeightLimit;
                    }
                    if (currentHeight < tLowerHeightLimit)
                    {
                        tLowerHeightLimit = currentHeight;

                        if (currentHeight < lowerHeightLimit)
                            currentHeight = lowerHeightLimit;
                    }

                    unsigned char colorVal = static_cast<unsigned char>(
                                                 (currentHeight-lowerHeightLimit)*254./(upperHeightLimit-lowerHeightLimit)
                                                 );
                    unsigned char colorVal2 = static_cast<unsigned char>(
                                                  254.-((currentHeight-lowerHeightLimit)*254./(upperHeightLimit-lowerHeightLimit))
                                                  );

                    shape.setPosition(50+2*j, 50+2*i);
                    shape.setFillColor(sf::Color(colorVal, 0, colorVal2));
                    window->draw(shape);
                }
            }

            upperHeightLimit = tUpperHeightLimit;
            lowerHeightLimit = tLowerHeightLimit;
        }

        //Display stabilized frame rate

        double tFr = 1000./clk.restart().asMilliseconds();
        frRt = (frRt + (frRt4 + frRt3 + frRt2 + tFr)/4.0) / 2.0;
        frRt4 = frRt3;
        frRt3 = frRt2;
        frRt2 = tFr;
        std::ostringstream sst;
        sst<<"Frame Rate: "<<static_cast<int>(frRt)<<"fps\t"<<"Fill Status (plotted slice):\n\t\t\t\t\t\t\t\t\t\t\t\t\t\t"
                                                            <<upperHeightLimit<<"/"<<INT16_MAX;
        txt.setString(sst.str());
        window->draw(txt);

        window->display();
    }

    plotThreadRunning = false;
}

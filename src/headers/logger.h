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

#ifndef SANDSIM_LOGGER_H
#define SANDSIM_LOGGER_H

#include <omp.h>

#include <fstream>
#include <string>

/*!
 * \brief Print log messages.
 *
 * Log messages can be printed to standard output via log().
 * If a log file name is provided to the constructor (see Logger()),
 * the log will also be written to this log file.
 *
 * Each log message begins with a GMT timestamp with 1s precision.
 * The exact format is "YYYY-MM-DDThh:mm:ssGMT".
 *
 * The log messages can have different log levels (see Logger::LogLevel).
 * The log level can be specified as an argument of log().
 * For every log level there is also a shortcut function:
 * logCritical(), ..., logDebugDebug().
 *
 * The global log level defines, which messages will be printed
 * (see Logger::LogLevel). It can be set via Logger() or via setLogLevel().
 *
 * The printing of each log message is protected by an OpenMP lock,
 * such that Logger can be used with OpenMP threads.
 * The thread number will be printed next to the timestamp.
 */
class Logger
{
public:
    enum class LogLevel;    //Forward declaration

public:
    explicit Logger(LogLevel pLogLevel = LogLevel::Info, const std::string& pLogFileName = "");     ///< Constructor.
    ~Logger();                                                                                      ///< Destructor.

public:
    /*!
     * \brief Enumeration of available log levels.
     *
     * Only log messages with log level below or equal to the current log level are accepted/logged.
     */
    enum class LogLevel
    {
        None       =  0,    ///< Do not print any log messages.
        Critical   = 10,    ///< Log critical errors only.
        Error      = 20,    ///< Log all errors.
        Warning    = 30,    ///< Log also warnings.
        Less       = 40,    ///< Log also important general notifications.
        Info       = 50,    ///< Log also normal notifications. The default.
        More       = 60,    ///< Log also less important notifications.
        Verbose    = 70,    ///< Log even more notifications.
        Debug      = 80,    ///< Log also debug messages.
        DebugDebug = 90     ///< Log even more debug messages.
    };

public:
    LogLevel getLogLevel() const;                                   ///< Get the log level.
    void setLogLevel(LogLevel pLogLevel);                           ///< Set the log level.
    //
    static std::string logLevelToLabel(LogLevel pLogLevel);         ///< Get the label for a log level.
    static LogLevel labelToLogLevel(const std::string& pLogLevel);  ///< Get the log level from its label.
    //
    void logCritical(const std::string& pMessage);      ///< Print a log message (LogLevel::Critical).
    void logError(const std::string& pMessage);         ///< Print a log message (LogLevel::Error).
    void logWarning(const std::string& pMessage);       ///< Print a log message (LogLevel::Warning).
    void logLess(const std::string& pMessage);          ///< Print a log message (LogLevel::Less).
    void logInfo(const std::string& pMessage);          ///< Print a log message (LogLevel::Info).
    void logMore(const std::string& pMessage);          ///< Print a log message (LogLevel::More).
    void logVerbose(const std::string& pMessage);       ///< Print a log message (LogLevel::Verbose).
    void logDebug(const std::string& pMessage);         ///< Print a log message (LogLevel::Debug).
    void logDebugDebug(const std::string& pMessage);    ///< Print a log message (LogLevel::DebugDebug).
    //
    void log(const std::string& pMessage, LogLevel pLogLevel = LogLevel::Info);     ///< Print a log message.

private:
    void logMessage(const std::string& pMessage, const std::string& pLevelString);  ///< Format and print an arbitrary log message.

private:
    LogLevel logLevel;      //Defines which log messages are accepted
    //
    omp_lock_t logLock;     //OMP lock to allow logging from different threads
    //
    std::ofstream logFile;  //Log file
};

#endif // SANDSIM_LOGGER_H

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

#include "logger.h"

#include <chrono>
#include <ctime>
#include <iomanip>
#include <ios>
#include <iostream>
#include <sstream>

/*!
 * \brief Constructor.
 *
 * Sets the logger's global log level to \p pLogLevel.
 * Initializes an OMP lock used to protect logging from different OMP threads.
 *
 * Tries to create/open a log file, if \p pLogFileName is not empty.
 *
 * \param pLogLevel Log level to use for all log messages.
 * \param pLogFileName File name for the log file. Can be empty (no log file).
 */
Logger::Logger(const LogLevel pLogLevel, const std::string& pLogFileName) :
    logLevel(pLogLevel)
{
    omp_init_lock(&logLock);

    //Try to open log file
    if (pLogFileName != "")
    {
        try
        {
            logFile.exceptions(std::ios_base::badbit | std::ios_base::failbit);
            logFile.open(pLogFileName);
        }
        catch (const std::ios_base::failure&)
        {
            std::cerr<<"ERROR: Could not open log file \"" + pLogFileName + "\"!"<<std::endl;
        }
    }
}

/*!
 * \brief Destructor.
 *
 * Closes the log file, if opened.
 * Destroys the OMP lock.
 */
Logger::~Logger()
{
    if (logFile.is_open())
        logFile.close();

    omp_destroy_lock(&logLock);
}

//Public

/*!
 * \brief Get the log level.
 *
 * \return Current log level.
 */
Logger::LogLevel Logger::getLogLevel() const
{
    return logLevel;
}

/*!
 * \brief Set the log level.
 *
 * \param pLogLevel New log level.
 */
void Logger::setLogLevel(const LogLevel pLogLevel)
{
    logLevel = pLogLevel;
}

//

/*!
 * \brief Get the label for a log level.
 *
 * Returns a (unique) formatted label for \p pLogLevel.
 * Converting back is possible using labelToLogLevel().
 *
 * \param pLogLevel The log level to get a label for.
 * \return The corresponding label for \p pLogLevel.
 */
std::string Logger::logLevelToLabel(const LogLevel pLogLevel)
{
    switch (pLogLevel)
    {
        case LogLevel::None:
            return "NONE";
        case LogLevel::Critical:
            return "CRIT";
        case LogLevel::Error:
            return "ERROR";
        case LogLevel::Warning:
            return "WARNG";
        case LogLevel::Less:
            return "LESS";
        case LogLevel::Info:
            return "INFO";
        case LogLevel::More:
            return "MORE";
        case LogLevel::Verbose:
            return "VERB";
        case LogLevel::Debug:
            return "DEBUG";
        case LogLevel::DebugDebug:
            return "DDBUG";
        default:
            return "INFO";
    }
}

/*!
 * \brief Get the log level from its label.
 *
 * Get the log level corresponding to a (unique) label \p pLogLevel,
 * which can in turn be obtained from logLevelToLabel().
 *
 * \param pLogLevel The label representing the requested log level.
 * \return The log level corresponding to label \p pLogLevel.
 */
Logger::LogLevel Logger::labelToLogLevel(const std::string& pLogLevel)
{
    if (pLogLevel == "NONE")
        return LogLevel::None;
    else if (pLogLevel == "CRIT")
        return LogLevel::Critical;
    else if (pLogLevel == "ERROR")
        return LogLevel::Error;
    else if (pLogLevel == "WARNG")
        return LogLevel::Warning;
    else if (pLogLevel == "LESS")
        return LogLevel::Less;
    else if (pLogLevel == "INFO")
        return LogLevel::Info;
    else if (pLogLevel == "MORE")
        return LogLevel::More;
    else if (pLogLevel == "VERB")
        return LogLevel::Verbose;
    else if (pLogLevel == "DEBUG")
        return LogLevel::Debug;
    else if (pLogLevel == "DDBUG")
        return LogLevel::DebugDebug;
    else
        return LogLevel::Info;
}

//

/*!
 * \brief Print a log message (LogLevel::Critical).
 *
 * Prints log message \p pMessage, if the Logger's global log level is higher than or equal to LogLevel::Critical.
 *
 * \param pMessage The message to log.
 */
void Logger::logCritical(const std::string& pMessage)
{
    if (logLevel >= LogLevel::Critical)
        logMessage(pMessage, logLevelToLabel(LogLevel::Critical));
}

/*!
 * \brief Print a log message (LogLevel::Error).
 *
 * Prints log message \p pMessage, if the Logger's global log level is higher than or equal to LogLevel::Error.
 *
 * \param pMessage The message to log.
 */
void Logger::logError(const std::string& pMessage)
{
    if (logLevel >= LogLevel::Error)
        logMessage(pMessage, logLevelToLabel(LogLevel::Error));
}

/*!
 * \brief Print a log message (LogLevel::Warning).
 *
 * Prints log message \p pMessage, if the Logger's global log level is higher than or equal to LogLevel::Warning.
 *
 * \param pMessage The message to log.
 */
void Logger::logWarning(const std::string& pMessage)
{
    if (logLevel >= LogLevel::Warning)
        logMessage(pMessage, logLevelToLabel(LogLevel::Warning));
}

/*!
 * \brief Print a log message (LogLevel::Less).
 *
 * Prints log message \p pMessage, if the Logger's global log level is higher than or equal to LogLevel::Less.
 *
 * \param pMessage The message to log.
 */
void Logger::logLess(const std::string& pMessage)
{
    if (logLevel >= LogLevel::Less)
        logMessage(pMessage, logLevelToLabel(LogLevel::Less));
}

/*!
 * \brief Print a log message (LogLevel::Info).
 *
 * Prints log message \p pMessage, if the Logger's global log level is higher than or equal to LogLevel::Info.
 *
 * \param pMessage The message to log.
 */
void Logger::logInfo(const std::string& pMessage)
{
    if (logLevel >= LogLevel::Info)
        logMessage(pMessage, logLevelToLabel(LogLevel::Info));
}

/*!
 * \brief Print a log message (LogLevel::More).
 *
 * Prints log message \p pMessage, if the Logger's global log level is higher than or equal to LogLevel::More.
 *
 * \param pMessage The message to log.
 */
void Logger::logMore(const std::string& pMessage)
{
    if (logLevel >= LogLevel::More)
        logMessage(pMessage, logLevelToLabel(LogLevel::More));
}

/*!
 * \brief Print a log message (LogLevel::Verbose).
 *
 * Prints log message \p pMessage, if the Logger's global log level is higher than or equal to LogLevel::Verbose.
 *
 * \param pMessage The message to log.
 */
void Logger::logVerbose(const std::string& pMessage)
{
    if (logLevel >= LogLevel::Verbose)
        logMessage(pMessage, logLevelToLabel(LogLevel::Verbose));
}

/*!
 * \brief Print a log message (LogLevel::Debug).
 *
 * Prints log message \p pMessage, if the Logger's global log level is higher than or equal to LogLevel::Debug.
 *
 * \param pMessage The message to log.
 */
void Logger::logDebug(const std::string& pMessage)
{
    if (logLevel >= LogLevel::Debug)
        logMessage(pMessage, logLevelToLabel(LogLevel::Debug));
}

/*!
 * \brief Print a log message (LogLevel::Debug_DEBUG).
 *
 * Prints log message \p pMessage, if the Logger's global log level is higher than or equal to LogLevel::Debug_DEBUG.
 *
 * \param pMessage The message to log.
 */
void Logger::logDebugDebug(const std::string& pMessage)
{
    if (logLevel >= LogLevel::DebugDebug)
        logMessage(pMessage, logLevelToLabel(LogLevel::DebugDebug));
}

//

/*!
 * \brief Print a log message.
 *
 * Prints log message \p pMessage, if \p pLogLevel is lower than or equal to the Logger's global log level.
 * Message with \p pLogLevel LogLevel::None will never be printed.
 *
 * \param pMessage The message to log.
 * \param pLogLevel The log level of \p pMessage.
 */
void Logger::log(const std::string& pMessage, const LogLevel pLogLevel)
{
    if (pLogLevel > logLevel)
        return;

    switch(pLogLevel)
    {
        case LogLevel::Critical:
        {
            logCritical(pMessage);
            break;
        }
        case LogLevel::Error:
        {
            logError(pMessage);
            break;
        }
        case LogLevel::Warning:
        {
            logWarning(pMessage);
            break;
        }
        case LogLevel::Less:
        {
            logLess(pMessage);
            break;
        }
        case LogLevel::Info:
        {
            logInfo(pMessage);
            break;
        }
        case LogLevel::More:
        {
            logMore(pMessage);
            break;
        }
        case LogLevel::Verbose:
        {
            logVerbose(pMessage);
            break;
        }
        case LogLevel::Debug:
        {
            logDebug(pMessage);
            break;
        }
        case LogLevel::DebugDebug:
        {
            logDebugDebug(pMessage);
            break;
        }
        default:
        {
            logInfo(pMessage);
            break;
        }
    }
}

//Private

/*!
 * \brief Format and print an arbitrary log message.
 *
 * Prints a formatted log message with text \p pMessage ("MESSAGE") prepended by
 * the current timestamp, \p pLevelString ("LEVEL") and the thread number ("n"):
 *
 * "[YYYY-MM-DDThh:mm:ssGMT, LEVEL|n] MESSAGE"
 *
 * If a log file was opened, the log message is also written to the log file.
 *
 * The logging is protected by the OMP lock, i.e. log messages from different threads will be displayed correctly.
 *
 * \param pMessage The message to log.
 * \param pLevelString Max. 5 character long string describing the log level.
 */
void Logger::logMessage(const std::string& pMessage, const std::string& pLevelString)
{
    auto now = std::chrono::system_clock::now();
    std::time_t time = std::chrono::system_clock::to_time_t(now);

    omp_set_lock(&logLock);

    std::ostringstream osstr;
    osstr<<"["<<std::put_time(std::gmtime(&time), "%FT%T%Z")<<", "<<std::setw(5)<<pLevelString
         <<"|"<<omp_get_thread_num()<<"] "<<pMessage<<"\n";

    //Output log message to both standard output and log file

    std::cout<<osstr.str();

    if (logFile.is_open())
        logFile<<osstr.str()<<std::flush;

    omp_unset_lock(&logLock);
}

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
Logger::Logger(LogLevel pLogLevel, const std::string& pLogFileName) :
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
void Logger::setLogLevel(LogLevel pLogLevel)
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
std::string Logger::logLevelToLabel(LogLevel pLogLevel)
{
    switch (pLogLevel)
    {
        case LogLevel::_NONE:
            return "NONE";
        case LogLevel::_CRITICAL:
            return "CRIT";
        case LogLevel::_ERROR:
            return "ERROR";
        case LogLevel::_WARNING:
            return "WARNG";
        case LogLevel::_LESS:
            return "LESS";
        case LogLevel::_INFO:
            return "INFO";
        case LogLevel::_MORE:
            return "MORE";
        case LogLevel::_VERBOSE:
            return "VERB";
        case LogLevel::_DEBUG:
            return "DEBUG";
        case LogLevel::_DEBUG_DEBUG:
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
        return LogLevel::_NONE;
    else if (pLogLevel == "CRIT")
        return LogLevel::_CRITICAL;
    else if (pLogLevel == "ERROR")
        return LogLevel::_ERROR;
    else if (pLogLevel == "WARNG")
        return LogLevel::_WARNING;
    else if (pLogLevel == "LESS")
        return LogLevel::_LESS;
    else if (pLogLevel == "INFO")
        return LogLevel::_INFO;
    else if (pLogLevel == "MORE")
        return LogLevel::_MORE;
    else if (pLogLevel == "VERB")
        return LogLevel::_VERBOSE;
    else if (pLogLevel == "DEBUG")
        return LogLevel::_DEBUG;
    else if (pLogLevel == "DDBUG")
        return LogLevel::_DEBUG_DEBUG;
    else
        return LogLevel::_INFO;
}

//

/*!
 * \brief Print a log message (LogLevel::_CRITICAL).
 *
 * Prints log message \p pMessage, if the Logger's global log level is higher than or equal to LogLevel::_CRITICAL.
 *
 * \param pMessage The message to log.
 */
void Logger::logCritical(const std::string& pMessage)
{
    if (logLevel >= LogLevel::_CRITICAL)
        logMessage(pMessage, logLevelToLabel(LogLevel::_CRITICAL));
}

/*!
 * \brief Print a log message (LogLevel::_ERROR).
 *
 * Prints log message \p pMessage, if the Logger's global log level is higher than or equal to LogLevel::_ERROR.
 *
 * \param pMessage The message to log.
 */
void Logger::logError(const std::string& pMessage)
{
    if (logLevel >= LogLevel::_ERROR)
        logMessage(pMessage, logLevelToLabel(LogLevel::_ERROR));
}

/*!
 * \brief Print a log message (LogLevel::_WARNING).
 *
 * Prints log message \p pMessage, if the Logger's global log level is higher than or equal to LogLevel::_WARNING.
 *
 * \param pMessage The message to log.
 */
void Logger::logWarning(const std::string& pMessage)
{
    if (logLevel >= LogLevel::_WARNING)
        logMessage(pMessage, logLevelToLabel(LogLevel::_WARNING));
}

/*!
 * \brief Print a log message (LogLevel::_LESS).
 *
 * Prints log message \p pMessage, if the Logger's global log level is higher than or equal to LogLevel::_LESS.
 *
 * \param pMessage The message to log.
 */
void Logger::logLess(const std::string& pMessage)
{
    if (logLevel >= LogLevel::_LESS)
        logMessage(pMessage, logLevelToLabel(LogLevel::_LESS));
}

/*!
 * \brief Print a log message (LogLevel::_INFO).
 *
 * Prints log message \p pMessage, if the Logger's global log level is higher than or equal to LogLevel::_INFO.
 *
 * \param pMessage The message to log.
 */
void Logger::logInfo(const std::string& pMessage)
{
    if (logLevel >= LogLevel::_INFO)
        logMessage(pMessage, logLevelToLabel(LogLevel::_INFO));
}

/*!
 * \brief Print a log message (LogLevel::_MORE).
 *
 * Prints log message \p pMessage, if the Logger's global log level is higher than or equal to LogLevel::_MORE.
 *
 * \param pMessage The message to log.
 */
void Logger::logMore(const std::string& pMessage)
{
    if (logLevel >= LogLevel::_MORE)
        logMessage(pMessage, logLevelToLabel(LogLevel::_MORE));
}

/*!
 * \brief Print a log message (LogLevel::_VERBOSE).
 *
 * Prints log message \p pMessage, if the Logger's global log level is higher than or equal to LogLevel::_VERBOSE.
 *
 * \param pMessage The message to log.
 */
void Logger::logVerbose(const std::string& pMessage)
{
    if (logLevel >= LogLevel::_VERBOSE)
        logMessage(pMessage, logLevelToLabel(LogLevel::_VERBOSE));
}

/*!
 * \brief Print a log message (LogLevel::_DEBUG).
 *
 * Prints log message \p pMessage, if the Logger's global log level is higher than or equal to LogLevel::_DEBUG.
 *
 * \param pMessage The message to log.
 */
void Logger::logDebug(const std::string& pMessage)
{
    if (logLevel >= LogLevel::_DEBUG)
        logMessage(pMessage, logLevelToLabel(LogLevel::_DEBUG));
}

/*!
 * \brief Print a log message (LogLevel::_DEBUG_DEBUG).
 *
 * Prints log message \p pMessage, if the Logger's global log level is higher than or equal to LogLevel::_DEBUG_DEBUG.
 *
 * \param pMessage The message to log.
 */
void Logger::logDebugDebug(const std::string& pMessage)
{
    if (logLevel >= LogLevel::_DEBUG_DEBUG)
        logMessage(pMessage, logLevelToLabel(LogLevel::_DEBUG_DEBUG));
}

//

/*!
 * \brief Print a log message.
 *
 * Prints log message \p pMessage, if \p pLogLevel is lower than or equal to the Logger's global log level.
 * Message with \p pLogLevel LogLevel::_NONE will never be printed.
 *
 * \param pMessage The message to log.
 * \param pLogLevel The log level of \p pMessage.
 */
void Logger::log(const std::string& pMessage, LogLevel pLogLevel)
{
    if (pLogLevel > logLevel)
        return;

    switch(pLogLevel)
    {
        case LogLevel::_CRITICAL:
        {
            logCritical(pMessage);
            break;
        }
        case LogLevel::_ERROR:
        {
            logError(pMessage);
            break;
        }
        case LogLevel::_WARNING:
        {
            logWarning(pMessage);
            break;
        }
        case LogLevel::_LESS:
        {
            logLess(pMessage);
            break;
        }
        case LogLevel::_INFO:
        {
            logInfo(pMessage);
            break;
        }
        case LogLevel::_MORE:
        {
            logMore(pMessage);
            break;
        }
        case LogLevel::_VERBOSE:
        {
            logVerbose(pMessage);
            break;
        }
        case LogLevel::_DEBUG:
        {
            logDebug(pMessage);
            break;
        }
        case LogLevel::_DEBUG_DEBUG:
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

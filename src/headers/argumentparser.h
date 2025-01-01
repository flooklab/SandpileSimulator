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

#ifndef ARGUMENTPARSER_H
#define ARGUMENTPARSER_H

#include <iostream>
#include <string>
#include <vector>
#include <deque>
#include <map>
#include <set>
#include <cctype>
#include <stdexcept>
#include <algorithm>


/*!
 * \brief Parse command line arguments.
 *
 * Automates parsing of command line switches and arguments by means of
 * a pre-defined list of possible arguments that is passed to the constructor.
 *
 * Possible argument types are simple switches and values of type string, integer or floating-point (see Argument::Type).
 * Each argument can be specified via its short form (dash + single character) or its long form (double dash + string/name).
 * Required argument dependencies can be specified for each argument.
 * See also Argument for more detailed information.
 *
 * For the short form, values must be provided as separate following arguments, e.g. "-a value".
 * For the long form values must be provided in the same argument separated by an equal sign, e.g. "--param=value".
 *
 * Multiple short-form arguments can be grouped together, e.g. "-abcd".
 * If some of those require values, they must be provided as a group of following arguments,
 * e.g. "-aBCd bValue cValue".
 *
 * Parsing of the command line arguments must be triggered by parseArgs().
 * %Argument values can then be queried with getParamValue(). Switch-type arguments
 * (Argument::Type::_NONE) do not have values. They can be queried with getParamCount(),
 * which returns the number of occurrences of the argument instead.
 * For error handling see parseArgs().
 *
 * Help/usage information can be printed with printHelp() and printUsage().
 * A help switch "-h/--help" can be automatically added to the accepted
 * arguments list (see ArgumentParser()) and queried with helpRequested().
 */
class ArgumentParser
{
public:
    struct Argument;    //Forward declaration

public:
    ArgumentParser(const std::string& pProgName, int pArgc, const char *const *const pArgv,
                   std::vector<Argument> pExpectedArgs, bool pAddHelpParam = false);        ///< Constructor.
    //
    bool parseArgs();                                           ///< Parse command line arguments and remember values and counts.
    //
    void printUsage() const;                                                    ///< Print valid/possible command line calls.
    void printHelp(size_t pLeftColumnLength = 8) const;                         ///< Print descriptions for all available arguments.
    //
    bool helpRequested() const;                                                 ///< Check, if a help switch was activated.
    //
    int getParamCount(const std::string& pParamName, int pDefault = 0) const;   ///< Get switch-type argument activation count.
    std::string getParamValue(const std::string& pParamName, std::string pDefault = "") const;  ///< Get argument value (as string).
    //
    std::vector<std::string> getTrailingArgs() const;                           ///< Get a list of trailing arguments/values.

private:
    bool unravelBunchOfArgs(const std::vector<std::string>& pBunch, bool pLastBunch = false);   ///< Interpret a group of arguments.

private:
    std::string progName;                                   //Name of the program/executable (for printUsage())
    std::vector<std::string> args;                          //Actual command line arguments from main()
    std::vector<Argument> expectedArgs;                     //Supported command line arguments
    //
    std::map<char, std::string> shortLongExpMap;            //Map from short form to long form of expected arguments
    std::map<std::string, char> longShortExpMap;            //Map from long form to short form of expected arguments
    //
    std::map<char, std::string> parsedArgsWithValues;       //Values of parsed value-type arguments
    std::map<char, int> parsedArgsFlags;                    //Counts of parsed switch-type arguments
    std::vector<std::string> floatingArgs, trailingArgs;    //Unassignable values (floating: between params, trailing: at the end)

public:
    /*!
     * \brief Definition of a command line argument/parameter.
     *
     * Defines a command line argument with a short (e.g. "-a") and a long (e.g. "--param") parameter expression,
     * a description of the argument, a value type (see Argument::Type; can be Type::_NONE for a "switch parameter")
     * and dependencies on other arguments. An argument can be mandatory or optional.
     *
     * Short parameter expressions (single char!) must not be a digit or a dash ('-')!
     *
     * See Argument() for more information.
     */
    struct Argument
    {
        /*!
         * \brief Argument types.
         *
         * Determines, if an argument requires a value and which type of value.
         */
        enum class Type
        {
            _NONE,          ///< Switch argument, which does not take a value.
            _STRING,        ///< %Argument with string as value.
            _INT,           ///< %Argument with integer as value.
            _FLOAT          ///< %Argument with floating-point number as value.
        };
        //
        const Type type;                            //Value type of the argument
        const char shortExp;                        //Short single-character form of the parameter
        const std::string longExp;                  //Long string form of the parameter
        const std::string usage;                    //Short expression describing the argument value (e.g. "FILENAME")
        const std::string description;              //Description of the argument
        const bool withValue;                       //Argument needs value (is set to true, if type is not _NONE)
        bool optional;                              //Optional argument (not const, change externally to reflect arg. dependencies!)
        const std::set<std::string> dependencies;   //Other arguments required by this one
        const std::set<std::string> conflicts;      //Other arguments conflicting with this one
        //
        /*!
         * \brief Constructor.
         *
         * Sets all argument properties (see also Argument). Member 'withValue' is set to true, iff \p pType is not Type::_NONE.
         *
         * Note: For a switch-type argument (\p pType == Type::_NONE) \p pOptional must be
         * true (mandatory switch does not make sense (in general)). After construction it can
         * however be manually set to false again e.g. in order to reflect an argument dependency.
         *
         * \param pType Value type of the argument.
         * \param pShortExp Single char as unique short expression for the parameter (must not be digit or '-').
         * \param pLongExp String as unique long expression for the parameter.
         * \param pDescription Description of the argument.
         * \param pUsage Short value description (e.g. "FILENAME", leads to printUsage() printing e.g. "--input==FILENAME").
         * \param pOptional %Argument is optional?
         * \param pDependencies %Argument dependencies expressed via their long parameter expressions.
         * \param pConflicts %Argument conflicts expressed via their long parameter expressions.
         *
         * \throws std::invalid_argument \p pShortExp must not be a digit or '-'.
         * \throws std::invalid_argument "Value description" \p pUsage must not be empty,
         *                               if argument requires value (\p pType not Type::_NONE).
         * \throws std::invalid_argument Switch-type argument must be optional.
         */
        Argument(Type pType, char pShortExp, std::string pLongExp, std::string pDescription, std::string pUsage = "",
                 bool pOptional = true, std::set<std::string> pDependencies = {}, std::set<std::string> pConflicts = {}) :
            type(pType),
            shortExp(pShortExp),
            longExp(std::move(pLongExp)),
            usage(std::move(pUsage)),
            description(std::move(pDescription)),
            withValue(type != Type::_NONE),
            optional(pOptional),
            dependencies(pDependencies),
            conflicts(pConflicts)
        {
            if (shortExp == '-' || std::isdigit(shortExp))
                throw std::invalid_argument("Invalid character for short argument name!");
            if (withValue && usage == "")
                throw std::invalid_argument("Argument with value requires usage string!");
            if (!withValue && !optional)
                throw std::invalid_argument("Switch-type argument must be optional (except, if dependency for other argument)!");
        }
    };
};

#endif // ARGUMENTPARSER_H

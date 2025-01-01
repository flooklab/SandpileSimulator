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

#include "argumentparser.h"

/*!
 * \brief Constructor.
 *
 * Registers a list of possible/valid command line arguments \p pExpectedArgs
 * and the actual command line arguments from the program call \p pArgv[\p pArgc].
 *
 * Checks \p pExpectedArgs for duplicate arguments (must have unique short form (e.g. "-a") and unique long form (e.g. "--paramA")).
 * Checks, if all specified argument dependencies and conflicts actually exist as valid parameters.
 *
 * If \p pAddHelpParam is true, a help parameter ("-h/--help") will be automatically added
 * to the \p pExpectedArgs list (can of course also be added manually). See also helpRequested().
 *
 * \p pProgName will be used as executable name for printUsage().
 *
 * \param pProgName Executable name for printUsage().
 * \param pArgc Command line arguments count (argc) from main().
 * \param pArgv Command line arguments (argv) from main().
 * \param pExpectedArgs List of possible arguments.
 * \param pAddHelpParam Add help parameter?
 *
 * \throws std::invalid_argument \p pArgc == 0 / missing argv[0]. Expecting at least one argument from OS.
 * \throws std::invalid_argument Duplicate argument in \p pExpectedArguments.
 * \throws std::invalid_argument Invalid argument dependency or conflict.
 */
ArgumentParser::ArgumentParser(const std::string& pProgName, int pArgc, const char *const *const pArgv,
                               std::vector<Argument> pExpectedArgs, bool pAddHelpParam) :
    progName(pProgName),
    args(pArgv, pArgv+pArgc),
    expectedArgs(pExpectedArgs)
{
    if (pArgc < 1)
        throw std::invalid_argument("You are trying to execute the progam in some weird execution environment. Missing argv[0].");
    if (args[0].find('-') == 0)
        throw std::invalid_argument("You are trying to execute the progam in some weird execution environment. Missing argv[0].");

    args.erase(args.begin());

    if (pAddHelpParam)
    {
        expectedArgs.push_back(Argument(Argument::Type::_NONE, 'h', "help",
                                        "Print a description of the command line arguments of this program."));
    }

    for (const Argument& arg : expectedArgs)
    {
        if (shortLongExpMap.find(arg.shortExp) != shortLongExpMap.end())
            throw std::invalid_argument("Same argument identifier specified twice in list of possible arguments: -\"" +
                                        std::string{arg.shortExp} + "\"");
        if (longShortExpMap.find(arg.longExp) != longShortExpMap.end())
            throw std::invalid_argument("Same argument identifier specified twice in list of possible arguments: --\"" +
                                        arg.longExp + "\"");

        shortLongExpMap[arg.shortExp] = arg.longExp;
        longShortExpMap[arg.longExp] = arg.shortExp;
    }

    for (const Argument& arg : expectedArgs)
    {
        for (const std::string& dep : arg.dependencies)
            if (longShortExpMap.find(dep) == longShortExpMap.end())
                throw std::invalid_argument("The argument dependency \"--" + dep + "\" is not a valid parameter!");

        for (const std::string& dep : arg.conflicts)
            if (longShortExpMap.find(dep) == longShortExpMap.end())
                throw std::invalid_argument("The argument conflict \"--" + dep + "\" is not a valid parameter!");
    }
}

//Public

/*!
 * \brief Parse command line arguments and remember values and counts.
 *
 * Parses the command line arguments passed to the constructor assuming the expected/possible arguments also passed
 * to the constructor (see ArgumentParser() and ArgumentParser). Interprets all valid and correctly formatted
 * arguments and internally assigns their values/counts (query by getParamValue(), getParamCount()).
 *
 * Invalid arguments and formatting errors are detected and also argument value types are checked,
 * as well as required argument dependencies and present argument conflicts. The function returns
 * false in case of such a problem. This means that checks for remaining/following arguments
 * will be skipped. However, values of other arguments (if correctly formatted!)
 * are still interpreted and can in principle be queried normally.
 *
 * Note: Wrongly formatted arguments that do not start with a dash will be simply ignored.
 * Also ignored will be values that are floating around (anything without leading dash) i.e. are not
 * needed by or assignable to a previous parameter. However, floating values that are not followed
 * by another parameter (i.e. values at the end of the command line) are "trailing arguments"
 * and will be added to a list without further processing (access by getTrailingArgs()).
 *
 * See also unravelBunchOfArgs().
 *
 * \return If parsing was successful and no formatting errors were detected.
 */
bool ArgumentParser::parseArgs()
{
    //Identify groups of argv elements argv[i] that belong together (e.g. "-a", "-bcd cValue dValue", "--e-param=12"),
    //which are separated by following elements starting with a dash; separately analyze each group as a whole
    //by unravelBunchOfArgs(), which extracts argument values and increments switch counts

    bool allOk = true;
    std::vector<std::string> bunchOfArgs;
    for (std::string arg : args)
    {
        //Start next group
        if (arg.find('-') == 0)
        {
            //Do NOT start next group, if not a parameter but a negative number!
            bool tIsNumber = false;
            if (arg.find("--") != 0)
            {
                try
                {
                    std::stoi(arg);
                    tIsNumber = true;
                }
                catch (const std::exception&)
                {
                    ;
                }
            }

            if (!tIsNumber)
            {
                bool ok = unravelBunchOfArgs(bunchOfArgs);
                if (!ok)
                    allOk = false;
                bunchOfArgs.clear();
            }
        }
        bunchOfArgs.push_back(std::move(arg));
    }
    bool ok = unravelBunchOfArgs(bunchOfArgs, true);
    if (!ok)
        allOk = false;

    //Return here in case of invalid arguments, missing values etc.
    if (!allOk)
        return false;

    //Handle argument dependencies: toggle 'optional' property of optional arguments that are required by other, entered arguments

    for (const auto& it : parsedArgsFlags)
    {
        std::set<std::string> tDependencies;
        for (const Argument& arg : expectedArgs)
        {
            if (arg.shortExp == it.first)
                tDependencies = arg.dependencies;
        }

        for (const std::string& longExp : tDependencies)
        {
            for (Argument& arg : expectedArgs)
            {
                if (arg.longExp == longExp)
                    arg.optional = false;
            }
        }
    }
    for (const auto& it : parsedArgsWithValues)
    {
        std::set<std::string> tDependencies;
        for (const Argument& arg : expectedArgs)
        {
            if (arg.shortExp == it.first)
                tDependencies = arg.dependencies;
        }

        for (const std::string& longExp : tDependencies)
        {
            for (Argument& arg : expectedArgs)
            {
                if (arg.longExp == longExp)
                    arg.optional = false;
            }
        }
    }

    //Check formatting (check string to number conversion etc.) and if all
    //required arguments are present and check for conflicting arguments

    for (const Argument& arg : expectedArgs)
    {
        if (arg.withValue && parsedArgsWithValues.find(arg.shortExp) == parsedArgsWithValues.end() &&
                parsedArgsFlags.find(arg.shortExp) != parsedArgsFlags.end())
        {
            std::cerr<<"Argument \"--"<<arg.longExp<<"\" requires a value!"<<std::endl;
            return false;
        }

        if (!arg.optional && parsedArgsWithValues.find(arg.shortExp) == parsedArgsWithValues.end() &&
                parsedArgsFlags.find(arg.shortExp) == parsedArgsFlags.end())
        {
            std::cerr<<"Missing required argument \"--"<<arg.longExp<<"\"!"<<std::endl;
            return false;
        }

        if (arg.withValue && parsedArgsWithValues.find(arg.shortExp) != parsedArgsWithValues.end())
        {
            if (arg.type == Argument::Type::_INT)
            {
                try
                {
                    std::stoi(parsedArgsWithValues.at(arg.shortExp));
                }
                catch (const std::exception&)
                {
                    std::cerr<<"Wrong format! Integer number required."<<std::endl;
                    return false;
                }
            }
            else if (arg.type == Argument::Type::_FLOAT)
            {
                try
                {
                    std::stod(parsedArgsWithValues.at(arg.shortExp));
                }
                catch (const std::exception&)
                {
                    std::cerr<<"Wrong format! Real number required."<<std::endl;
                    return false;
                }
            }
        }

        //If argument present, check for present conflicts
        if (parsedArgsWithValues.find(arg.shortExp) != parsedArgsWithValues.end() ||
                parsedArgsFlags.find(arg.shortExp) != parsedArgsFlags.end())
        {
            for (const std::string& conflict : arg.conflicts)
            {
                char tShortExp = longShortExpMap[conflict];

                if (parsedArgsWithValues.find(tShortExp) != parsedArgsWithValues.end() ||
                        parsedArgsFlags.find(tShortExp) != parsedArgsFlags.end())
                {
                    std::cerr<<"Found conflicting argument for argument \"--"<<arg.longExp<<"\": \"--"<<conflict<<"\""<<std::endl;
                    return false;
                }
            }
        }
    }

    return true;
}

//

/*!
 * \brief Print valid/possible command line calls.
 *
 * Prints an example command line invocation of the program with all possible arguments.
 * Optional arguments are indicated by square brackets.
 * Prints the executable name as defined in ArgumentParser().
 */
void ArgumentParser::printUsage() const
{
    std::string usageString = "Usage:\n " + progName;
    for (const Argument& arg : expectedArgs)
    {
        usageString.append(" ");
        if (arg.optional)
            usageString.append("[");
        usageString.append("--").append(arg.longExp);
        if (arg.withValue)
            usageString.append("=").append(arg.usage);
        if (arg.optional)
            usageString.append("]");
    }
    std::cout<<usageString<<std::endl;
}

/*!
 * \brief Print descriptions for all available arguments.
 *
 * Prints command line usage information (calls printUsage()) and then lists all available
 * command line arguments (shows short and long forms) together with a description for each one.
 *
 * For nicer formatting, all descriptions start with the same total indentation \p pLeftColumnLength.
 * If the short/long parameter name string exceeds this length, a new line is inserted before.
 *
 * \param pLeftColumnLength Common horizontal start position for all description strings.
 */
void ArgumentParser::printHelp(size_t pLeftColumnLength) const
{
    printUsage();

    std::string helpString = "\nArguments:\n";
    for (const Argument& arg : expectedArgs)
    {
        std::string currentLine;
        std::string full = arg.longExp;
        if (arg.withValue)
            full.append("=").append(arg.usage);
        currentLine.append(" ").append("-").append(std::string{arg.shortExp}).append(", ").append("--").append(full);

        if (!arg.description.empty())
        {
            if (currentLine.length() >= pLeftColumnLength)
            {
                currentLine.append("\n");
                for (size_t i = 0; i < pLeftColumnLength; ++i)
                    currentLine.append(" ");
            }
            else
            {
                for (auto i = currentLine.length(); i < pLeftColumnLength; ++i)
                    currentLine.append(" ");
            }
            currentLine.append(arg.description);
        }
        currentLine.append("\n\n");
        helpString.append(currentLine);
    }
    std::cout<<helpString<<std::endl;
}

//

/*!
 * \brief Check, if a help switch was activated.
 *
 * Checks entered command line arguments for switch-type parameter "--help"
 * in case such argument was specified as expected argument and for any
 * other switch-type argument with short form "-h" else.
 *
 * \return If help switch was found in entered arguments.
 */
bool ArgumentParser::helpRequested() const
{
    unsigned char helpParam = (longShortExpMap.find("help") != longShortExpMap.end()) ? longShortExpMap.at("help") : 'h';
    return parsedArgsFlags.find(helpParam) != parsedArgsFlags.end();
}

//

/*!
 * \brief Get switch-type argument activation count.
 *
 * Returns the number of times the switch-type argument \p pParamName was entered.
 * If this number is zero or if no switch \p pParamName exists, \p pDefault is returned.
 *
 * \param pParamName Long name of the parameter.
 * \param pDefault Default value.
 * \return Number of times the argument was entered (or \p pDefault).
 */
int ArgumentParser::getParamCount(const std::string& pParamName, int pDefault) const
{
    if (longShortExpMap.find(pParamName) != longShortExpMap.end() &&
            parsedArgsFlags.find(longShortExpMap.at(pParamName)) != parsedArgsFlags.end())
    {
        return parsedArgsFlags.at(longShortExpMap.at(pParamName));
    }

    return pDefault;
}

/*!
 * \brief Get argument value (as string).
 *
 * Returns the entered value for argument \p pParamName as string
 * (independent of Argument::Type; type conversion is checked in parseArgs()).
 * If the argument was not specified or if no argument \p pParamName exists, \p pDefault is returned.
 *
 * \param pParamName Long name of the parameter.
 * \param pDefault Default value.
 * \return Argument value (or \p pDefault) as string.
 */
std::string ArgumentParser::getParamValue(const std::string& pParamName, std::string pDefault) const
{
    if (longShortExpMap.find(pParamName) != longShortExpMap.end() &&
            parsedArgsWithValues.find(longShortExpMap.at(pParamName)) != parsedArgsWithValues.end())
        return parsedArgsWithValues.at(longShortExpMap.at(pParamName));

    return pDefault;
}

//

/*!
 * \brief Get a list of trailing arguments/values.
 *
 * \return Trailing arguments that were not assignable to a previous parameter.
 */
std::vector<std::string> ArgumentParser::getTrailingArgs() const
{
    return trailingArgs;
}

//Private

/*!
 * \brief Interpret a group of arguments.
 *
 * Parses a group of argv items argv[i] that belong together, which means the first item must start with a dash
 * and at least one short parameter (char) or a double dash and one long parameter (string). The following items
 * of the group must be values for short parameters from the first item. Values will be assigned to their parameters
 * in the order they appear and are requested (see below, example 2). There can be more value items than needed,
 * which will be simply ignored. Values for long parameters must be assigned in the same item with an equal sign (example 3).
 * Assigning a value to a no-value long parameter switch does not lead to an error and is treated as if entered without assignment.
 *
 * Here are three command line examples that would each form a valid argument group:
 * 1. "-a" (one argv item; single switch w/o value))
 * 2. "-bcd cValue dValue" (three argv items; one switch w/o value, two parameters with values)
 * 3. "--e-param=12" (one argv item, one (long) parameter with value)
 *
 * Values for each value argument will be put in 'parsedArgsWithValues'.
 * No-value arguments will be counted and the counts added to 'parsedArgsFlags'.
 *
 * If a (formatting) error is found, the function stops interpretation of the group and returns false.
 * Value or count assignments for the whole group might be missing/incomplete!
 *
 * Possible error conditions are:
 * - Parameter requires a value (see ArgumentParser(), Argument) but no value is assigned or too few value items follow.
 * - Value assignment to a long parameter with missing parameter name.
 * - Invalid (not expected) parameter (see ArgumentParser()).
 *
 * If the first item does not start with a dash, all items of the group will be ignored.
 *
 * Note: Formatting/content of values is not checked. Empty assignments ("--empty-assign=") are possible.
 *
 * \param pBunch List of arguments of the argument group.
 * \param pLastBunch Last argument group of command line? Decides whether unassignable
 *                   values are "floating" (discard) or "trailing" (remember).
 * \return If parsing of the argument group was successful.
 */
bool ArgumentParser::unravelBunchOfArgs(const std::vector<std::string>& pBunch, bool pLastBunch)
{
    std::deque<char> foundArgsWithValues;
    bool firstArg = true;

    for (std::string str : pBunch)
    {
        if (firstArg)
        {
            firstArg = false;

            if (str.find('-') == 0)
            {
                if (str.find("--") == 0)    //Long parameter
                {
                    //Remove double dash
                    str.erase(0, 2);

                    auto pos = str.find('=');
                    if (pos == std::string::npos)   //No '=' sign: must be parameter without value
                    {
                        if (longShortExpMap.find(str) == longShortExpMap.end())
                        {
                            std::cerr<<"Invalid argument \"--"<<str<<"\"!"<<std::endl;
                            return false;
                        }

                        char c = longShortExpMap.at(str);

                        for (const struct Argument& arg : expectedArgs)
                        {
                            if (c == arg.shortExp)
                            {
                                if (arg.withValue)
                                {
                                    std::cerr<<"Argument \"--"<<arg.longExp<<"\" requires a value!"<<std::endl;
                                    return false;
                                }

                                ++(parsedArgsFlags[c]);
                                break;
                            }
                        }
                    }
                    else    //Found '=' sign: split value from parameter
                    {
                        std::string longExp = str.substr(0, pos);
                        std::string value = str.substr(pos+1);

                        if (longExp.empty())
                        {
                            std::cerr<<"Missing parameter name for \""<<str<<"\"!"<<std::endl;
                            return false;
                        }
                        if (longShortExpMap.find(longExp) == longShortExpMap.end())
                        {
                            std::cerr<<"Invalid argument \"--"<<longExp<<"\"!"<<std::endl;
                            return false;
                        }

                        char c = longShortExpMap.at(longExp);

                        for (const struct Argument& arg : expectedArgs)
                        {
                            if (c == arg.shortExp)
                            {
                                if (arg.withValue)
                                    parsedArgsWithValues[c] = std::move(value);
                                else
                                    ++(parsedArgsFlags[c]);

                                break;
                            }
                        }
                    }
                }
                else    //Short parameter(s)
                {
                    //Remove dash
                    str.erase(0, 1);

                    for (auto it = str.begin(); it != str.end(); ++it)
                    {
                        char c = *it;
                        if (shortLongExpMap.find(c) == shortLongExpMap.end())
                        {
                            std::cerr<<"Invalid argument \"-"<<c<<"\"!"<<std::endl;
                            return false;
                        }

                        for (const struct Argument& arg : expectedArgs)
                        {
                            if (c == arg.shortExp)
                            {
                                if (arg.withValue)
                                    foundArgsWithValues.push_back(c);
                                else
                                    ++(parsedArgsFlags[c]);

                                break;
                            }
                        }
                    }
                }
            }
            else    //Not a parameter, add to floating (or trailing, if last bunch) arguments
            {
                if (pLastBunch)
                    trailingArgs.push_back(std::move(str));
                else
                    floatingArgs.push_back(std::move(str));
            }
        }
        else
        {
            if (!foundArgsWithValues.empty())
            {
                char param = *foundArgsWithValues.begin();
                parsedArgsWithValues[param] = std::move(str);
                foundArgsWithValues.erase(foundArgsWithValues.begin());
            }
            else
            {
                if (pLastBunch)
                    trailingArgs.push_back(std::move(str));
                else
                    floatingArgs.push_back(std::move(str));
            }
        }
    }

    //Missing values for short form arguments?
    if (!foundArgsWithValues.empty())
    {
        for (char c : foundArgsWithValues)
            std::cerr<<"Argument \"-"<<c<<"\" requires a value!"<<std::endl;

        return false;
    }

    return true;
}

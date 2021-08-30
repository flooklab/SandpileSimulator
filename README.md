# Sandpile Simulator

Sandpile Simulator is a cellular automaton (library) for simulating and analyzing the dynamics of
sandpiles with special focus on self-organized criticality (SOC) and avalanche scaling behavior.
Sandpiles with arbitrary dimensions can be simulated using different models and from the
collected statistics scaling exponents can be extracted using a "moment analysis" approach.

Sandpile Simulator is a rewrite (and extension) of the original sandpile simulation
project below, which was written in Python, Copyright (C) 2018 M. Frohne, P. Wolf:  
https://github.com/leloup314/Sandpiles

For more information on the topic see also this report:  
https://github.com/leloup314/Sandpiles/blob/sandpiles/report/report.pdf

Note that the information about the cellular automaton and simulation model algorithms (and the
data analysis approach) that is discussed in the report essentially still applies to this project.
However, some implementation and model details might differ (see this project's code documentation!).

## Build

Building the project can be configured with [CMake](https://cmake.org/).

Notes:
- Building the documentation requires [Doxygen](https://github.com/doxygen/doxygen) (optional).
- Enabling the feature for visualization of sandpiles requires [SFML](https://github.com/SFML/SFML) (optional).

## License information

Copyright (C) 2021 M. Frohne

Sandpile Simulator is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License,
or (at your option) any later version.

Sandpile Simulator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with Sandpile Simulator. If not, see <https://www.gnu.org/licenses/>.

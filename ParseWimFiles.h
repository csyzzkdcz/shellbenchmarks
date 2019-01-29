#ifndef PARSEWIMFILES_H
#define PARSEWIMFILES_H

#include <string>

class SimulationSetup;
class SimulationSetupNormal;
class SecondFundamentalFormDiscretization;

bool parseWimFiles(const std::string &prefixRes, const std::string &prefixTar, SimulationSetup &parsedSetup, const SecondFundamentalFormDiscretization &sff);

#endif
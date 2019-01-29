#ifndef SIMULATIONSETUPNORMAL_H
#define SIMULATIONSETUPNORMAL_H

#include "SimulationSetup.h"

struct SimulationSetupNormal : public SimulationSetup
{
    public:
    virtual void buildRestFundamentalForms(const SecondFundamentalFormDiscretization &sff) override;
};



#endif
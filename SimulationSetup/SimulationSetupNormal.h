#ifndef SIMULATIONSETUPNORMAL_H
#define SIMULATIONSETUPNORMAL_H

#include "SimulationSetup.h"

class SimulationSetupNormal : public SimulationSetup
{
    public:
    virtual void buildRestFundamentalForms(const SecondFundamentalFormDiscretization &sff) override;
    
    bool loadAbars() override;
};



#endif

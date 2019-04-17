#include "SimulationSetupNormal.h"
#include "../GeometryDerivatives.h"
#include "../SecondFundamentalForm/SecondFundamentalFormDiscretization.h"

void SimulationSetupNormal::buildRestFundamentalForms(const SecondFundamentalFormDiscretization &sff)
{
    int nfaces = mesh.nFaces();
    abars.resize(nfaces);
    bbars.resize(nfaces);

    for (int i = 0; i < nfaces; i++)
    {
        abars[i] = firstFundamentalForm(mesh, initialPos, i, NULL, NULL);
        bbars[i] = sff.secondFundamentalForm(mesh, initialPos, initialEdgeDOFs, i, NULL, NULL);
    }
}

bool SimulationSetupNormal::loadParams()
{
    return true;
}

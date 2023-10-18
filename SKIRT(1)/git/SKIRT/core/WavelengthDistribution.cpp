/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "WavelengthDistribution.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

void WavelengthDistribution::setupSelfBefore()
{
    SimulationItem::setupSelfBefore();

    _random = find<Random>();
}

//////////////////////////////////////////////////////////////////////

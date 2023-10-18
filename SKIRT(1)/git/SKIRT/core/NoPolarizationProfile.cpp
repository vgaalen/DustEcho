/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "NoPolarizationProfile.hpp"

////////////////////////////////////////////////////////////////////

int NoPolarizationProfile::dimension() const
{
    return 1;
}

////////////////////////////////////////////////////////////////////

StokesVector NoPolarizationProfile::polarizationForDirection(Direction /*bfk*/) const
{
    return StokesVector();
}

////////////////////////////////////////////////////////////////////

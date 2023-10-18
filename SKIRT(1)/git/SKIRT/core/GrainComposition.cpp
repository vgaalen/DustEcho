/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "GrainComposition.hpp"

////////////////////////////////////////////////////////////////////

string GrainComposition::resourceNameForMuellerMatrix() const
{
    return string();
}

////////////////////////////////////////////////////////////////////

bool GrainComposition::resourcesForSpheroidalEmission(bool& /*resource*/, double& /*interpol*/, string& /*tableName1*/,
                                                      string& /*tableName2*/) const
{
    return false;
}

////////////////////////////////////////////////////////////////////

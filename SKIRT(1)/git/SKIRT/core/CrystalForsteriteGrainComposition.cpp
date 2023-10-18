/*//////////////////////////////////////////////////////////////////
////       SKIRT -- an advanced radiative transfer code         ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "CrystalForsteriteGrainComposition.hpp"

//////////////////////////////////////////////////////////////////////

string CrystalForsteriteGrainComposition::name() const
{
    return "Crystalline_Forsterite";
}

//////////////////////////////////////////////////////////////////////

double CrystalForsteriteGrainComposition::bulkDensity() const
{
    return 3.33e3;
}

//////////////////////////////////////////////////////////////////////

string CrystalForsteriteGrainComposition::resourceNameForOpticalProps() const
{
    return "MinForsteriteOpticalProps";
}

//////////////////////////////////////////////////////////////////////

string CrystalForsteriteGrainComposition::resourceNameForEnthalpies() const
{
    return "DustEM_aSil_Entalphies";
}

//////////////////////////////////////////////////////////////////////

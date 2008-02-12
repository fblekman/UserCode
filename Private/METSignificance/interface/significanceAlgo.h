#ifndef SIGMET_ASIGNIFICANCE_H
#define SIGMET_ASIGNIFICANCE_H
// -*- C++ -*-
//
// Class: Significance
//
// Original Author:  Kyle Story
//   Created:        11/09/2007
//
//  Subroutine: ASignificance(vector<SigInputObj*> EventVec);
//  
//  Purpose:
//  
//  This subroutine takes in a vector of type SigInputObj that includes
//  all of the physics objects in the event.  It then calculates
//  the log significance of the MET of the event.
//
//  The Significance (S) is defined as:
//    ln(S) = 1/2 Chisq_0
//  where Chisq_0 is the value of Chi squared at MET=0.
//
//
// $Id: Significance.h,v 1.3 2007/12/07 00:20:10 kstory Exp $
// 
// Revision history
// 
// $Log: Significance.h,v $
// Revision 1.3  2007/12/07 00:20:10  kstory
// MET phi was changed to have a range of [-pi, pi]
//
// Revision 1.2  2007/11/30 22:02:33  kstory
// Changed the arguments to allow calculation of the total MET and the MET_phi.
//
// Revision 1.1  2007/11/28 22:51:02  xshi
// Copied from Kyle and put into CVS
//

#include <vector>
#include <string>
#include <iostream>
#include <sstream>

#include "Wenu/METSignificance/interface/SigInputObj.h"

namespace metsig{
  double chisq( double x,  double y,  double X0,  double Y0,  double rho,  double sigma_x,  double sigma_y);
  
  double ASignificance(const std::vector<metsig::SigInputObj*>& EventVec, double& met_r, double& met_phi, double& met_set);
}

#endif

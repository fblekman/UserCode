#include "Wenu/METSignificance/interface/significanceAlgo.h"
#include <iostream>
#include <string>
#include "math.h"
#include "TROOT.h"

//*** Chisq for a single object ******************************//
double 
metsig::chisq( double  x,  double  y,  double  X0,  double  Y0,  double  rho,  double  sigma_x,  double  sigma_y)
{
  
  double result = ( pow((sigma_y*(x-X0)),2) 
		    - 2.*rho*sigma_x*sigma_y*(x-X0)*(y-Y0) 
		    + pow((sigma_x*(y-Y0)),2) ) 
    / ( (1-pow(rho,2))*pow(sigma_x,2)*pow(sigma_y,2) );
  
  return result;
}
//*********************************************************//


//*************************************************//
double 
metsig::ASignificance(const std::vector<SigInputObj*>& EventVec, double &met_r, double &met_phi, double &met_set) 
{
  
  if(EventVec.size()<1) {
    std::cerr << "Event Vector is empty!  Return -1:\n";
    return(-1);
  }
  //=== Analytical Generation of Chisq Contours ===//
  std::vector<double> sigma_x;
  std::vector<double> sigma_y;
  std::vector<double> rho_xy;
  double set_worker;
  double sigma_xMET;
  double sigma_yMET;
  double rho_xyMET;
  double XMET=0, YMET=0;
  //--- Temporary variables ---//
  double sigma_u;
  double sigma_v;
  double phi_tmp;
  
  double sigma_xSum;
  double sigma_ySum;
  double rho_xySum;
  
  //--- Calculate the first SigInputObj: ---//

  XMET = EventVec.at(0)->get_energy() * cos(EventVec.at(0)->get_phi());
  YMET = EventVec.at(0)->get_energy() * sin(EventVec.at(0)->get_phi());
  set_worker = sqrt(XMET*XMET+YMET*YMET);
  sigma_u  = EventVec.at(0)->get_sigma_e();
  sigma_v  = EventVec.at(0)->get_sigma_tan();
  
  phi_tmp  = EventVec.at(0)->get_phi();
  sigma_x.push_back(sqrt( pow(cos(phi_tmp)*sigma_u,2) + pow(sin(phi_tmp)*sigma_v,2) ) );
  sigma_y.push_back( sqrt( pow(sin(phi_tmp)*sigma_u,2) + pow(cos(phi_tmp)*sigma_v,2) ) );
  rho_xy.push_back( ( cos(phi_tmp)*sin(phi_tmp)*(sigma_u*sigma_u - sigma_v*sigma_v) ) /
		    sqrt( (pow(cos(phi_tmp)*sigma_u,2) + pow(sin(phi_tmp)*sigma_v,2))*
			  (pow(sin(phi_tmp)*sigma_u,2) + pow(cos(phi_tmp)*sigma_v,2)) ) );
  
  sigma_xSum = sigma_x.at(0);
  sigma_ySum = sigma_y.at(0);
  rho_xySum = rho_xy.at(0);
  
  //--- Calculate sigma_x,y and rho of MET iteratively ---//
  for(unsigned int objnum=1; objnum < EventVec.size(); objnum++ ) {
    double xval = EventVec.at(objnum)->get_energy() * cos(EventVec.at(objnum)->get_phi());
    double yval =  EventVec.at(objnum)->get_energy() * sin(EventVec.at(objnum)->get_phi());
    XMET += xval;
    YMET += yval;
    set_worker+= sqrt(xval*xval+yval*yval);
    sigma_u = EventVec.at(objnum)->get_sigma_e();
    sigma_v = EventVec.at(objnum)->get_sigma_tan();
    phi_tmp = EventVec.at(objnum)->get_phi();
    sigma_x.push_back( sqrt( pow(cos(phi_tmp)*sigma_u,2) + pow(sin(phi_tmp)*sigma_v,2) ) );
    sigma_y.push_back( sqrt( pow(sin(phi_tmp)*sigma_u,2) + pow(cos(phi_tmp)*sigma_v,2) ) );
    rho_xy.push_back( ( cos(phi_tmp)*sin(phi_tmp)*(sigma_u*sigma_u - sigma_v*sigma_v) ) /
		      sqrt( (pow(cos(phi_tmp)*sigma_u,2) + pow(sin(phi_tmp)*sigma_v,2))*
			    (pow(sin(phi_tmp)*sigma_u,2) + pow(cos(phi_tmp)*sigma_v,2)) ) );
    
    sigma_xSum = sqrt(sigma_xSum*sigma_xSum + sigma_x.at(objnum)*sigma_x.at(objnum));
    sigma_ySum = sqrt(sigma_ySum*sigma_ySum + sigma_y.at(objnum)*sigma_y.at(objnum));
    rho_xySum  = (rho_xySum*sigma_xSum*sigma_ySum + rho_xy.at(objnum)*sigma_x.at(objnum)*sigma_y.at(objnum)) /
      sqrt( (sigma_xSum*sigma_xSum + sigma_x.at(objnum)*sigma_x.at(objnum))*
	    (sigma_ySum*sigma_ySum + sigma_y.at(objnum)*sigma_y.at(objnum)) );
  }
  
  sigma_xMET = sigma_xSum;
  sigma_yMET = sigma_ySum;
  rho_xyMET = rho_xySum;
  /*=============================================================*/
  
  //--- Calculate magnitude and angle of MET, store in returned variables ---//
  met_r = sqrt(XMET*XMET + YMET*YMET);
  met_set = set_worker;
  //--- Ensure met_phi is in [-pi, pi] ---//
  double tmp_met_phi;
  if( (XMET) >=0.0 ) {
    tmp_met_phi = TMath::ATan( (YMET) / (XMET) );
  }
  else {
    if( (YMET) >=0.0 ) {
      tmp_met_phi = TMath::ATan(YMET/XMET) + TMath::Pi();
    }
    else{ // => YMET<0
      tmp_met_phi = TMath::ATan(YMET/XMET) - TMath::Pi();
    }
  }
  met_phi = tmp_met_phi;
  
  //--- Calculate Significance ---//
  double chisq0 = chisq(0., 0., XMET, YMET, rho_xyMET, sigma_xMET, sigma_yMET);
  double chisqMET = chisq(XMET, YMET, XMET, YMET, rho_xyMET, sigma_xMET, sigma_yMET);
  double lnLikelihoodMax = -.5*chisqMET;
  double lnLikelihood0 = -.5*chisq0;
  double mySignificance = lnLikelihoodMax - lnLikelihood0;
  
  return mySignificance;
}
//*** End of ASignificance ********************************//

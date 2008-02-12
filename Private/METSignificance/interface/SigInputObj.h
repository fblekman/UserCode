#ifndef METSIG_SIGINPUTOBJ_H
#define METSIG_SIGINPUTOBJ_H
// -*- C++ -*-
//
// Class: SigInputObj
//
// Original Author:  Kyle Story
//   Created:        11/07/2007
//
// $Id: SigInputObj.h,v 1.1 2007/11/28 22:51:52 xshi Exp $
// 
// Revision history
// 
// $Log: SigInputObj.h,v $
// Revision 1.1  2007/11/28 22:51:52  xshi
// Copied from Kyle and put into CVS
//


#include <vector>
#include <string>
#include <iostream>
#include <sstream>

//=== Class SigInputObj ==============================//
namespace metsig{
  class SigInputObj{

  public:
    SigInputObj():
      type(""),energy(0.),phi(0.),sigma_e(0.),sigma_tan(0.)
      {;}// default constructor

    SigInputObj(std::string m_type, double m_energy, double m_phi,
	    double m_sigm_e, double m_sigma_phi);
    ~SigInputObj() {;}
    
    std::string get_type(){return(type);};
    double get_energy(){return(energy);};
    double get_phi(){return(phi);};
    double get_sigma_e(){return(sigma_e);};
    double get_sigma_tan(){return(sigma_tan);};
    
  private:
    std::string type;       //type of physics object
    /*** note: type = "jet", "uc" (un-clustered energy),
	 "electron", "muon", "hot-cell", "vertex" ***/
    double energy;    //magnitude of the energy
    double phi;       //azimuthal angle
    double sigma_e;   //gaus width in radial direction
    double sigma_tan; //gaus width in phi-hat direction (not in rad)
    
    void set_type(std::string m_type){type.clear(); type.append(m_type);};
    void set_energy(double m_energy){energy=m_energy;};
    void set_phi(double m_phi){phi=m_phi;};
    void set_sigma_e(double m_sigma_e){sigma_e=m_sigma_e;};
    void set_sigma_tan(double m_sigma_tan){sigma_tan=m_sigma_tan;};
  };
}

#endif

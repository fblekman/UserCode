#include "Wenu/METSignificance/interface/SigInputObj.h"

//=== Constructors ===============================//
metsig::SigInputObj::SigInputObj(std::string m_type, double m_energy, double m_phi,
		 double m_sigma_e, double m_sigma_tan) 
{
  type.clear(); 
  type.append(m_type);
  energy = m_energy;
  phi = m_phi;
  sigma_e = m_sigma_e;
  sigma_tan = m_sigma_tan;
}

//================================================//

//=== Methods ====================================//
// none yet...
//================================================//


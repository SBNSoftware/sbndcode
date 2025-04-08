//
// Build a dictionary.
//
// $Id: classes.h,v 1.8 2010/04/12 18:12:28  Exp $
// $Author:  $
// $Date: 2010/04/12 18:12:28 $
// 
// Original author Rob Kutschke, modified by wes
//

#include "canvas/Persistency/Common/Wrapper.h"

// data-products
// lardataobj
//#include "lardata/Utilities/AssociationUtil.h"
#include "canvas/Persistency/Common/Assns.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "sbndcode/BlipRecoSBND/Utils/DataTypes.h"
#include "lardataobj/RecoBase/SpacePoint.h"

//
// Only include objects that we would like to be able to put into the event.
// Do not include the objects they contain internally.
//

template class art::Assns<recob::Hit,blip::Blip,void>;
template class art::Wrapper<art::Assns<recob::Hit,blip::Blip,void> >;
template class std::vector<blip::Blip>;
template class art::Wrapper<std::vector<blip::Blip> >;
template class std::map<int, TVector3>;
template class art::Assns<blip::Blip,recob::Hit,void>;
template class art::Wrapper<art::Assns<blip::Blip,recob::Hit,void> >;
template class art::Assns<blip::Blip,recob::SpacePoint,void>;
template class art::Wrapper<art::Assns<blip::Blip,recob::SpacePoint,void> >;

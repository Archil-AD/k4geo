#include "DD4hep/DetFactoryHelper.h"
#include <DDRec/DetectorData.h>
#include "XML/Utilities.h"

namespace det {
/**
  Simple cylinder using Tube to be used to define cylinder composed of 1 single material
  @author Clement Helsens
**/
static dd4hep::Ref_t
createSimpleCylinder(dd4hep::Detector& lcdd, xml_h e, dd4hep::SensitiveDetector sensDet) {
  xml_det_t x_det = e;
  std::string name = x_det.nameStr();
  dd4hep::DetElement cylinderDet(name, x_det.id());

  dd4hep::Volume experimentalHall = lcdd.pickMotherVolume(cylinderDet);

  xml_comp_t cylinderDim(x_det.child(_U(dimensions)));

  dd4hep::Tube cylinder(
      cylinderDim.rmin(), cylinderDim.rmax(), cylinderDim.dz(), cylinderDim.phi0(), cylinderDim.deltaphi());

  dd4hep::Volume cylinderVol(
      x_det.nameStr() + "_SimpleCylinder", cylinder, lcdd.material(cylinderDim.materialStr()));

  if (x_det.isSensitive()) {
    dd4hep::xml::Dimension sdType(x_det.child(_U(sensitive)));
    cylinderVol.setSensitiveDetector(sensDet);
    sensDet.setType(sdType.typeStr());
  }

  dd4hep::PlacedVolume cylinderPhys;

  double zoff = cylinderDim.z_offset();
  if (fabs(zoff) > 0.000000000001) {
    dd4hep::Position trans(0., 0., zoff);
    cylinderPhys = experimentalHall.placeVolume(cylinderVol,
                                                dd4hep::Transform3D(dd4hep::RotationZ(0.), trans));
  } else
    cylinderPhys = experimentalHall.placeVolume(cylinderVol);

  cylinderPhys.addPhysVolID("system", x_det.id());

  cylinderDet.setPlacement(cylinderPhys);

  cylinderDet.setVisAttributes(lcdd, x_det.visStr(), cylinderVol);


  // Create caloData object
  auto caloData = new dd4hep::rec::LayeredCalorimeterData;
  caloData->extent[0] = cylinderDim.rmin();
  caloData->extent[1] = cylinderDim.rmax();
  if(name=="MuonTaggerBarrel")
  {
    caloData->extent[2] = 0;
    caloData->extent[3] = cylinderDim.dz();
    caloData->layoutType = dd4hep::rec::LayeredCalorimeterData::BarrelLayout;
  }
  else
  {
    caloData->extent[2] = zoff;
    caloData->extent[3] = zoff + cylinderDim.dz();
    caloData->layoutType = dd4hep::rec::LayeredCalorimeterData::EndcapLayout;
  }

  dd4hep::rec::LayeredCalorimeterData::Layer caloLayer;
  caloLayer.distance                  = cylinderDim.rmin();
  caloLayer.sensitive_thickness       = (cylinderDim.rmax() - cylinderDim.rmin());
  caloLayer.inner_thickness           = (cylinderDim.rmax() - cylinderDim.rmin()) / 2.0;
  caloLayer.outer_thickness           = (cylinderDim.rmax() - cylinderDim.rmin()) / 2.0;
  caloData->layers.push_back(caloLayer);

  cylinderDet.addExtension<dd4hep::rec::LayeredCalorimeterData>(caloData);

  // Set type flags
  dd4hep::xml::setDetectorTypeFlag(x_det, cylinderDet);

  return cylinderDet;
}
}
DECLARE_DETELEMENT(SimpleCylinder_o1_v01, det::createSimpleCylinder)


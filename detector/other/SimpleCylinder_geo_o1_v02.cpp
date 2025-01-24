#include "DD4hep/DetFactoryHelper.h"
#include <DDRec/DetectorData.h>
#include "XML/Utilities.h"
#include "DDRec/MaterialManager.h"
#include "DDRec/Vector3D.h"


namespace det {
/**
  Simple cylinder using Tube to be used to define cylinder composed of 1 single material
  Based on SimpleCylinder_geo_o1_v01.cpp prepared by Clement Helsens
  @author A. Durglishvili
**/
static dd4hep::Ref_t
createSimpleCylinder(dd4hep::Detector& lcdd, xml_h e, dd4hep::SensitiveDetector sensDet) {
  xml_det_t x_det = e;
  std::string name = x_det.nameStr();
  dd4hep::DetElement cylinderDet(name, x_det.id());
  dd4hep::Volume experimentalHall = lcdd.pickMotherVolume(cylinderDet);

  xml_comp_t cylinderDim(x_det.child(_U(dimensions)));

  // envelop shape used only for endcap
  dd4hep::Tube envelopShape(cylinderDim.rmin(), cylinderDim.rmax(),
                            (cylinderDim.z_offset()+cylinderDim.dz()),
                            cylinderDim.phi0(), cylinderDim.deltaphi());

  // envelop volume where we will place endcap cylinders
  dd4hep::Volume envelopeVolume(x_det.nameStr() + "_envelopeVolume", envelopShape, lcdd.air());

  // shape for barrel/endcap
  dd4hep::Tube cylinder(cylinderDim.rmin(), cylinderDim.rmax(),cylinderDim.dz(), cylinderDim.phi0(), cylinderDim.deltaphi());

  // barrel or endcap cylinder
  dd4hep::Volume cylinderVol(x_det.nameStr() + "_SimpleCylinder", cylinder, lcdd.material(cylinderDim.materialStr()));

  if (x_det.isSensitive()) {
    dd4hep::xml::Dimension sdType(x_det.child(_U(sensitive)));
    cylinderVol.setSensitiveDetector(sensDet);
    sensDet.setType(sdType.typeStr());
  }

  // Create caloData object
  auto caloData = new dd4hep::rec::LayeredCalorimeterData;
  dd4hep::rec::LayeredCalorimeterData::Layer caloLayer;
  dd4hep::rec::MaterialManager matMgr(experimentalHall);

  caloData->extent[0] = cylinderDim.rmin();
  caloData->extent[1] = cylinderDim.rmax();

  if(name=="MuonTaggerEndcap")
  {
    dd4hep::PlacedVolume cylinderPhys1; // negative endcap
    dd4hep::PlacedVolume cylinderPhys2; // positive endcap
    double zoff = cylinderDim.z_offset();
    dd4hep::Position trans1(0., 0., -zoff);
    dd4hep::Position trans2(0., 0., zoff);
    // place endcap cylinders into envelope volume
    cylinderPhys1 = envelopeVolume.placeVolume(cylinderVol,dd4hep::Transform3D(dd4hep::RotationZ(0.), trans1));
    cylinderPhys2 = envelopeVolume.placeVolume(cylinderVol,dd4hep::Transform3D(dd4hep::RotationZ(0.), trans2));

    //cylinderPhys1.addPhysVolID("system", x_det.id());
    cylinderPhys1.addPhysVolID("subsystem", 0); // negative endcap
    cylinderPhys1.addPhysVolID("layer", 0);

    //cylinderPhys2.addPhysVolID("system", x_det.id());
    cylinderPhys2.addPhysVolID("subsystem", 1); // positive endcap
    cylinderPhys2.addPhysVolID("layer", 0);

    // Place envelope volume into experimentalHall
    dd4hep::PlacedVolume placedEndcap = experimentalHall.placeVolume(envelopeVolume);
    placedEndcap.addPhysVolID("system", x_det.id());
    cylinderDet.setPlacement(placedEndcap);
    // FIXME! AD: maybe should be envelopeVolume instead of cylinderVol???
    cylinderDet.setVisAttributes(lcdd, x_det.visStr(), cylinderVol);

    caloData->extent[2] = zoff - cylinderDim.dz();
    caloData->extent[3] = zoff + cylinderDim.dz();
    caloData->layoutType = dd4hep::rec::LayeredCalorimeterData::EndcapLayout;

    dd4hep::rec::Vector3D ivr1 = dd4hep::rec::Vector3D(0., 0., zoff - cylinderDim.dz()); // defining starting vector points
    dd4hep::rec::Vector3D ivr2 = dd4hep::rec::Vector3D(0., 0., zoff + cylinderDim.dz());  // defining end vector

    const dd4hep::rec::MaterialVec &materials = matMgr.materialsBetween(ivr1, ivr2); // calling material manager to get material info between two points
    auto mat = matMgr.createAveragedMaterial(materials); // creating average of all the material between two points to calculate X0 and lambda of averaged material
    const double nRadiationLengths = cylinderDim.dz()*2. / mat.radiationLength();
    const double nInteractionLengths = cylinderDim.dz()*2. / mat.interactionLength();

    caloLayer.distance                  = zoff - cylinderDim.dz();
    caloLayer.sensitive_thickness       = cylinderDim.dz()*2.;
    caloLayer.inner_thickness           = cylinderDim.dz();
    caloLayer.outer_thickness           = cylinderDim.dz();

    caloLayer.inner_nRadiationLengths   = nRadiationLengths / 2.0;
    caloLayer.inner_nInteractionLengths = nInteractionLengths / 2.0;
    caloLayer.outer_nRadiationLengths   = nRadiationLengths / 2.0;
    caloLayer.outer_nInteractionLengths = nInteractionLengths / 2.0;
  }
  else
  {
    dd4hep::PlacedVolume cylinderPhys;
    cylinderPhys = experimentalHall.placeVolume(cylinderVol);

    cylinderPhys.addPhysVolID("system", x_det.id());
    cylinderPhys.addPhysVolID("layer", 0);
    cylinderDet.setPlacement(cylinderPhys);
    cylinderDet.setVisAttributes(lcdd, x_det.visStr(), cylinderVol);

    caloData->extent[2] = 0;
    caloData->extent[3] = cylinderDim.dz();
    caloData->layoutType = dd4hep::rec::LayeredCalorimeterData::BarrelLayout;

    dd4hep::rec::Vector3D ivr1 = dd4hep::rec::Vector3D(0., cylinderDim.rmin(), 0.); // defining starting vector points 
    dd4hep::rec::Vector3D ivr2 = dd4hep::rec::Vector3D(0., cylinderDim.rmax(), 0.);  // defining end vector

    const dd4hep::rec::MaterialVec &materials = matMgr.materialsBetween(ivr1, ivr2); // calling material manager to get material info between two points
    auto mat = matMgr.createAveragedMaterial(materials); // creating average of all the material between two points to calculate X0 and lambda of averaged material
    const double nRadiationLengths = (cylinderDim.rmax() - cylinderDim.rmin()) / mat.radiationLength();
    const double nInteractionLengths = (cylinderDim.rmax() - cylinderDim.rmin()) / mat.interactionLength();

    caloLayer.distance                  = cylinderDim.rmin();
    caloLayer.sensitive_thickness       = (cylinderDim.rmax() - cylinderDim.rmin());
    caloLayer.inner_thickness           = (cylinderDim.rmax() - cylinderDim.rmin()) / 2.0;
    caloLayer.outer_thickness           = (cylinderDim.rmax() - cylinderDim.rmin()) / 2.0;

    caloLayer.inner_nRadiationLengths   = nRadiationLengths / 2.0;
    caloLayer.inner_nInteractionLengths = nInteractionLengths / 2.0;
    caloLayer.outer_nRadiationLengths   = nRadiationLengths / 2.0;
    caloLayer.outer_nInteractionLengths = nInteractionLengths / 2.0;
  }

  caloLayer.cellSize0 = 20 * dd4hep::mm; // FIXME! AD: should be corrected from DDGeometryCreatorALLEGRO
  caloLayer.cellSize1 = 20 * dd4hep::mm; // FIXME! AD: should be corrected from DDGeometryCreatorALLEGRO

  // attach the layer to the caloData
  caloData->layers.push_back(caloLayer);

  // attach the layer to the cylinderDet
  cylinderDet.addExtension<dd4hep::rec::LayeredCalorimeterData>(caloData);

  // Set type flags
  dd4hep::xml::setDetectorTypeFlag(x_det, cylinderDet);

  return cylinderDet;
}
}
DECLARE_DETELEMENT(SimpleCylinder_o1_v02, det::createSimpleCylinder)


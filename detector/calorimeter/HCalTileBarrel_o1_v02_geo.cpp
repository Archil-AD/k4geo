// DD4hep
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DDRec/MaterialManager.h"
#include "DDRec/Vector3D.h"
#include "XML/Utilities.h"
#include <DDRec/DetectorData.h>

// k4geo
#include "detectorSegmentations/FCCSWHCalPhiRow_k4geo.h"
#include "detectorSegmentations/FCCSWHCalPhiTheta_k4geo.h"

using dd4hep::DetElement;
using dd4hep::PlacedVolume;
using dd4hep::Volume;
using dd4hep::xml::Dimension;

namespace det {

static dd4hep::Ref_t createHCal(dd4hep::Detector& lcdd, xml_det_t xmlDet, dd4hep::SensitiveDetector sensDet) {

  /////////////////// config parsing ///////////////////////////////////

  // Make volume that envelopes the whole barrel; set material to air
  Dimension xDimensions(xmlDet.dimensions());

  // sensitive detector type read from xml (for example "SimpleCalorimeterSD")
  Dimension xSensitive = xmlDet.child(_U(sensitive));
  sensDet.setType(xSensitive.typeStr());

  xml_comp_t xEndPlate = xmlDet.child(_Unicode(end_plate));
  double dZEndPlate = xEndPlate.thickness();
  xml_comp_t xFacePlate = xmlDet.child(_Unicode(face_plate));
  xml_comp_t xSpace = xmlDet.child(_Unicode(plate_space)); // to avoid overlaps
  double space = xSpace.thickness();
  xml_comp_t xSteelSupport = xmlDet.child(_Unicode(steel_support));
  double dSteelSupport = xSteelSupport.thickness();

  dd4hep::printout(dd4hep::DEBUG, "HCalTileBarrel_o1_v02", "steel support thickness (cm): %.2f", dSteelSupport);

  double sensitiveBarrelRmin = xDimensions.rmin() + xFacePlate.thickness() + space;

  // Hard-coded assumption that we have two different sequences for the modules
  std::vector<xml_comp_t> sequences = {xmlDet.child(_Unicode(sequence_a)), xmlDet.child(_Unicode(sequence_b))};
  // Check if both sequences are present
  if (!sequences[0] || !sequences[1]) {
    dd4hep::printout(dd4hep::ERROR, "HCalTileBarrel_o1_v02",
                     "The two sequences 'sequence_a' and 'sequence_b' must be present in the xml file.");
    throw std::runtime_error("Missing sequence_a or sequence_b in the xml file.");
  }
  // Check if both sequences have the same dimensions
  Dimension dimensionsA(sequences[0].dimensions());
  Dimension dimensionsB(sequences[1].dimensions());
  if (dimensionsA.dz() != dimensionsB.dz()) {
    dd4hep::printout(dd4hep::ERROR, "HCalTileBarrel_o1_v02",
                     "The dimensions of sequence_a and sequence_b do not match.");
    throw std::runtime_error("The dimensions of the sequence_a and sequence_b do not match.");
  }
  double dzSequence = dimensionsB.dz();
  dd4hep::printout(dd4hep::DEBUG, "HCalTileBarrel_o1_v02", "sequence thickness %.2f", dzSequence);

  // number of sequences fitting in Z
  unsigned int numSequencesZ = lcdd.constant<unsigned>("BarHCal_numSequencesZ");

  // number of modules along phi
  unsigned int numModules = lcdd.constant<unsigned>("BarHCal_n_phi_modules");

  // spacing between modules
  double moduleSpacing = lcdd.constant<double>("BarHCAL_module_spacing_phi");

  // get all 'layer' children of the 'layers' tag
  std::vector<xml_comp_t> Layers;
  for (xml_coll_t xCompColl(xmlDet.child(_Unicode(layers)), _Unicode(layer)); xCompColl; ++xCompColl) {
    Layers.push_back(xCompColl);
  }
  unsigned int numLayersR = 0;
  double moduleDepth = 0.;
  std::vector<double> layerDepths = std::vector<double>();
  std::vector<double> layerInnerRadii = std::vector<double>();
  for (std::vector<xml_comp_t>::iterator it = Layers.begin(); it != Layers.end(); ++it) {
    xml_comp_t layer = *it;
    Dimension layerDimension(layer.dimensions());
    numLayersR += layerDimension.nModules();
    for (int nLayer = 0; nLayer < layerDimension.nModules(); nLayer++) {
      moduleDepth += layerDimension.dr();
      layerDepths.push_back(layerDimension.dr());
    }
  }
  // Calculate correction along z based on the module size (can only have natural number of modules)
  double dzDetector = (numSequencesZ * dzSequence) / 2 + dZEndPlate + space;

  dd4hep::printout(dd4hep::INFO, "HCalTileBarrel_o1_v02", "dzDetector (cm): %.2f", dzDetector);
  dd4hep::printout(dd4hep::DEBUG, "HCalTileBarrel_o1_v02", "correction of dz in cm (negative = size reduced): %.2f",
                   dzDetector - xDimensions.dz());

  dd4hep::printout(dd4hep::DEBUG, "HCalTileBarrel_o1_v02",
                   "retrieved number of radial layers: %d , which end up to a full module depth in rho of %.2f cm",
                   numLayersR, moduleDepth);
  dd4hep::printout(dd4hep::DEBUG, "HCalTileBarrel_o1_v02", "retrieved number of radial layers: %d", layerDepths.size());
  dd4hep::printout(dd4hep::INFO, "HCalTileBarrel_o1_v02",
                   "constructing: %d sequences in Z, %d radial layers, in total %d tiles", numSequencesZ, numLayersR,
                   numLayersR * numSequencesZ);

  double rminSupport = sensitiveBarrelRmin + moduleDepth;
  double rmaxSupport = sensitiveBarrelRmin + moduleDepth + dSteelSupport;

  double sensitiveBarrelRmax = sensitiveBarrelRmin + moduleDepth;

  ////////////////////// detector building //////////////////////

  std::vector<dd4hep::PlacedVolume> layers;
  layers.reserve(layerDepths.size());
  std::vector<std::vector<dd4hep::PlacedVolume>> tileSeqTubesInLayers;
  tileSeqTubesInLayers.reserve(layerDepths.size());
  std::vector<std::vector<dd4hep::PlacedVolume>> seqInTileSeqTube;
  seqInTileSeqTube.reserve(layerDepths.size());
  std::vector<dd4hep::PlacedVolume> tilesPerLayer;
  tilesPerLayer.reserve(layerDepths.size());

  // top level det element representing whole hcal barrel
  DetElement caloDetElem(xmlDet.nameStr(), xmlDet.id());

  /// envelope shape
  dd4hep::Tube envelopeShape(xDimensions.rmin(), xDimensions.rmax(), xDimensions.dz());

  Volume envelopeVolume("HCalEnvelopeVolume", envelopeShape, lcdd.air());
  envelopeVolume.setVisAttributes(lcdd, xDimensions.visStr());

  // Add structural support made of steel inside of HCal
  dd4hep::Tube facePlateShape(xDimensions.rmin(), sensitiveBarrelRmin, (dzDetector - dZEndPlate - space));
  Volume facePlateVol("HCalFacePlateVol", facePlateShape, lcdd.material(xFacePlate.materialStr()));
  facePlateVol.setVisAttributes(lcdd, xFacePlate.visStr());
  PlacedVolume placedFacePlate = envelopeVolume.placeVolume(facePlateVol);
  DetElement facePlate_det(caloDetElem, "HCalFacePlate", 0);
  facePlate_det.setPlacement(placedFacePlate);

  // Add structural support made of steel at both ends of HCal
  dd4hep::Tube endPlateShape(xDimensions.rmin(), (xDimensions.rmax() - dSteelSupport), dZEndPlate / 2);
  Volume endPlateVol("HCalEndPlateVol", endPlateShape, lcdd.material(xEndPlate.materialStr()));
  endPlateVol.setVisAttributes(lcdd, xEndPlate.visStr());

  DetElement endPlatePos(caloDetElem, "HCalEndPlatePos", 0);
  dd4hep::Position posOffset(0, 0, dzDetector - (dZEndPlate / 2));
  PlacedVolume placedEndPlatePos = envelopeVolume.placeVolume(endPlateVol, posOffset);
  endPlatePos.setPlacement(placedEndPlatePos);

  DetElement endPlateNeg(caloDetElem, "HCalEndPlateNeg", 1);
  dd4hep::Position negOffset(0, 0, -dzDetector + (dZEndPlate / 2));
  PlacedVolume placedEndPlateNeg = envelopeVolume.placeVolume(endPlateVol, negOffset);
  endPlateNeg.setPlacement(placedEndPlateNeg);

  dd4hep::Tube supportShape(rminSupport, rmaxSupport, (dzDetector - dZEndPlate - space));
  Volume steelSupportVolume("HCalSteelSupportVol", supportShape, lcdd.material(xSteelSupport.materialStr()));
  steelSupportVolume.setVisAttributes(lcdd.invisible());
  PlacedVolume placedSupport = envelopeVolume.placeVolume(steelSupportVolume);
  DetElement support(caloDetElem, "HCalSteelSupport", 0);
  support.setPlacement(placedSupport);

  // loop over R ("layers")
  double layerR = 0.;
  for (unsigned int idxLayer = 0; idxLayer < layerDepths.size(); ++idxLayer) {
    std::string layerName = "HCalLayer" + std::to_string(idxLayer);

    // in Module rmin = 0  for first wedge, changed radius to the full radius starting at (0,0,0)
    double rminLayer = sensitiveBarrelRmin + layerR;
    double rmaxLayer = sensitiveBarrelRmin + layerR + layerDepths.at(idxLayer);
    layerR += layerDepths.at(idxLayer);
    layerInnerRadii.push_back(rminLayer);
    dd4hep::printout(dd4hep::INFO, "HCalTileBarrel_o1_v02", "layer %d (cm): %.2f - %.2f", idxLayer, rminLayer,
                     rmaxLayer);

    // alternate: even layers in first module consist of tile sequence b, odd layer of tile sequence a
    unsigned int sequenceIdx = idxLayer % 2;

    // create a tube volume for wedge shape tile sequences
    dd4hep::Tube tileSequenceShape(rminLayer, rmaxLayer, 0.5 * dzSequence);
    Volume tileSequenceVolume("HCalTileSequenceVol", tileSequenceShape, lcdd.air());

    // create a volume for a single wedge shape sequence of type a and b
    double phi_min = 0.0;
    double phi_max = 2. * M_PI / numModules;
    dd4hep::ConeSegment sequenceShape(0.5 * dzSequence, rminLayer, rmaxLayer, rminLayer, rmaxLayer,
                                            phi_min, phi_max);
    std::array<Volume,2> sequenceVol;
    for(unsigned int iseq=0; iseq < 2; iseq++){
      sequenceVol[iseq] = Volume("HCalSeqVol", sequenceShape, lcdd.air());
      sequenceVol[iseq].setVisAttributes(lcdd.invisible());
    }

    // create a layer volume
    dd4hep::Tube layerShape(rminLayer, rmaxLayer, dzDetector - dZEndPlate - space);
    Volume layerVolume("HCalLayerVol", layerShape, lcdd.air());
    layerVolume.setVisAttributes(lcdd.invisible());
    dd4hep::PlacedVolume placedLayerVolume = envelopeVolume.placeVolume(layerVolume);
    placedLayerVolume.addPhysVolID("layer", idxLayer);
    layers.push_back(placedLayerVolume);

    // place tiles in sequence_a and sequence_b volumes
    for(unsigned int iseq=0; iseq < 2; iseq++)
    {
      double tileZOffset = -0.5 * dzSequence;
      for (xml_coll_t xCompColl(sequences[iseq], _Unicode(module_component)); xCompColl; ++xCompColl) {
        xml_comp_t xComp = xCompColl;
        double phi1 = phi_min;
        double phi2 = phi_max;

        // reduce tile size by moduleSpacing size from both sides
        phi1 = asin(moduleSpacing/rminLayer);
        phi2 = phi2 - phi1;

        dd4hep::ConeSegment tileShape(0.5 * xComp.thickness(), rminLayer, rmaxLayer, rminLayer, rmaxLayer,
                                            phi1, phi2);
        Volume tileVol("HCalTileVol_" + xComp.nameStr(), tileShape, lcdd.material(xComp.materialStr()));
        tileVol.setVisAttributes(lcdd, xComp.visStr());
        if (xComp.isSensitive()) {
          tileVol.setSensitiveDetector(sensDet);
        }
        dd4hep::Position tileOffset( 0., 0., tileZOffset + 0.5 * xComp.thickness());
        PlacedVolume placedTileVol = sequenceVol[iseq].placeVolume(tileVol, tileOffset);
        tileZOffset += xComp.thickness();
      }
    }

    std::vector<dd4hep::PlacedVolume> sequence_vector;

    // place numModules=256 tile sequences into tube shape tileSequenceVolume
    for (unsigned int iWedge = 0; iWedge < numModules; ++iWedge) {
      double phi_angle = iWedge * (2. * M_PI / numModules);
      dd4hep::Transform3D transform(dd4hep::RotationZ(phi_angle), dd4hep::Position(0.,0.,0.));
      unsigned int iseq = 0;
      if(sequenceIdx==0) iseq = iWedge % 2;
      if(sequenceIdx==1) iseq = 1 - (iWedge % 2);
      PlacedVolume placedSequenceVol = tileSequenceVolume.placeVolume(sequenceVol[iseq], iWedge, transform);
      placedSequenceVol.addPhysVolID("phi", iWedge);
      sequence_vector.push_back(placedSequenceVol);
    }
    seqInTileSeqTube.push_back(sequence_vector);

    // second z loop (place sequences in layer)
    std::vector<dd4hep::PlacedVolume> tileSequence_vector;

    for (uint numSeq = 0; numSeq < numSequencesZ; numSeq++) {
      double zOffset = -dzDetector + numSeq * dzSequence + dzSequence / 2 + dZEndPlate + space;
      dd4hep::Position tileSequencePosition(0, 0, zOffset);
      dd4hep::PlacedVolume placedTileSequenceVolume = layerVolume.placeVolume(tileSequenceVolume, tileSequencePosition);
      placedTileSequenceVolume.addPhysVolID("row", numSeq);
      tileSequence_vector.push_back(placedTileSequenceVolume);
    }
    tileSeqTubesInLayers.push_back(tileSequence_vector);
  }

  // Place det elements within each other to recover volume positions later via cellID
  for (uint iLayer = 0; iLayer < numLayersR; iLayer++) {
    DetElement layerDet(caloDetElem, dd4hep::xml::_toString(iLayer, "layer%d"), iLayer);
    layerDet.setPlacement(layers[iLayer]);

    for (uint iSeqTube = 0; iSeqTube < tileSeqTubesInLayers[iLayer].size(); iSeqTube++) {
      DetElement tileSeqTubeDet(layerDet, dd4hep::xml::_toString(iSeqTube, "tileSeqTube%d"), iSeqTube);
      tileSeqTubeDet.setPlacement(tileSeqTubesInLayers[iLayer][iSeqTube]);

      for (uint iSeq = 0; iSeq < seqInTileSeqTube[iLayer].size(); iSeq++) {
        DetElement seqDet(tileSeqTubeDet, dd4hep::xml::_toString(iSeq, "tileSeq%d"), iSeq);
        seqDet.setPlacement(seqInTileSeqTube[iLayer][iSeq]);
      }
    }
  }

  // Place envelope (or barrel) volume
  Volume motherVol = lcdd.pickMotherVolume(caloDetElem);
  motherVol.setVisAttributes(lcdd.invisible());
  PlacedVolume envelopePhysVol = motherVol.placeVolume(envelopeVolume);
  envelopePhysVol.addPhysVolID("system", caloDetElem.id());
  caloDetElem.setPlacement(envelopePhysVol);

  // retrieve handle to segmentation, needed to get cell sizes
  dd4hep::Segmentation segHandle = sensDet.readout().segmentation();
  // try to retrieve segmentation itself
  std::string layerFieldName;
  dd4hep::DDSegmentation::FCCSWHCalPhiTheta_k4geo* seg_phitheta =
      dynamic_cast<dd4hep::DDSegmentation::FCCSWHCalPhiTheta_k4geo*>(segHandle.segmentation());
  dd4hep::DDSegmentation::FCCSWHCalPhiRow_k4geo* seg_phirow =
      dynamic_cast<dd4hep::DDSegmentation::FCCSWHCalPhiRow_k4geo*>(segHandle.segmentation());
  if (seg_phitheta) {
    dd4hep::printout(dd4hep::DEBUG, "HCalTileBarrel_o1_v02", "Segmentation is of type FCCSWHCalPhiTheta_k4geo");
    layerFieldName = seg_phitheta->fieldNameLayer();
  } else if (seg_phirow) {
    dd4hep::printout(dd4hep::DEBUG, "HCalTileBarrel_o1_v02", "Segmentation is of type FCCSWHCalPhiRow_k4geo");
    layerFieldName = seg_phirow->fieldNameLayer();
  } else {
    dd4hep::printout(dd4hep::ERROR, "HCalTileBarrel_o1_v02", "Unknown segmentation");
    throw std::runtime_error("Incorrect readout in calorimeter xml description!");
  }
  std::string cellIDEncoding = sensDet.readout().idSpec().fieldDescription();
  dd4hep::BitFieldCoder encoder(cellIDEncoding);

  // Create caloData object
  auto caloData = new dd4hep::rec::LayeredCalorimeterData;
  caloData->layoutType = dd4hep::rec::LayeredCalorimeterData::BarrelLayout;
  caloDetElem.addExtension<dd4hep::rec::LayeredCalorimeterData>(caloData);

  caloData->extent[0] = sensitiveBarrelRmin;
  caloData->extent[1] = sensitiveBarrelRmax;
  caloData->extent[2] = 0.; // NN: for barrel detectors this is 0
  caloData->extent[3] = dzDetector;

  dd4hep::rec::MaterialManager matMgr(envelopeVolume);
  dd4hep::rec::LayeredCalorimeterData::Layer caloLayer;

  dd4hep::printout(dd4hep::INFO, "HCalTileBarrel_o1_v02", "Layer structure information:");
  for (unsigned int idxLayer = 0; idxLayer < layerDepths.size(); ++idxLayer) {
    dd4hep::printout(dd4hep::INFO, "HCalTileBarrel_o1_v02", "  Layer %d", idxLayer);
    const double difference_bet_r1r2 = layerDepths.at(idxLayer);
    double thickness_sen = 0.;
    double absorberThickness = 0.;

    // Average material radiation length in a given layer depends on the polar angle.
    // Pandora uses it (but mainly for ECAL) to calculate expected amount of radiation from previous layer
    // However, it assumes that the material distribution is uniform vs z, and then corrects the radiation length
    // based on the direction of the particle.
    // In our case, this is not true. At the moment, we are going to calculate it at an angle of 60 degrees
    // so that a mixture of active and passive material is seen, but we'll probably have to do something different
    // in Pandora if we see that this has a significant impact on the results (probably it won't seen this information
    // does not seem to be relied upon a lot for the HCAL).
    const double angle = 60. * M_PI / 180.;
    dd4hep::rec::Vector3D ivr1 = dd4hep::rec::Vector3D(0., layerInnerRadii.at(idxLayer),
                                                       0); // defining starting vector points of the given layer
    dd4hep::rec::Vector3D ivr2 = dd4hep::rec::Vector3D(0., layerInnerRadii.at(idxLayer) + layerDepths.at(idxLayer),
                                                       layerDepths.at(idxLayer) * cos(angle) /
                                                           sin(angle)); // defining end vector points of the given layer

    const dd4hep::rec::MaterialVec& materials =
        matMgr.materialsBetween(ivr1, ivr2); // calling material manager to get material info between two points
    auto mat = matMgr.createAveragedMaterial(materials); // creating average of all the material between two points to
                                                         // calculate X0 and lambda of averaged material
    const double nRadiationLengths = layerDepths.at(idxLayer) / mat.radiationLength();
    const double nInteractionLengths = layerDepths.at(idxLayer) / mat.interactionLength();

    std::string str1("Polystyrene"); // sensitive material
    for (size_t imat = 0; imat < materials.size(); imat++) {
      std::string str2(materials.at(imat).first.name());
      if (str1.compare(str2) == 0) {
        thickness_sen += materials.at(imat).second;
      } else if (str2 != "Air") {
        absorberThickness += materials.at(imat).second;
      }
    }
    dd4hep::printout(dd4hep::INFO, "HCalTileBarrel_o1_v02", "    sensitive thickness is: %lf", thickness_sen);
    dd4hep::printout(dd4hep::INFO, "HCalTileBarrel_o1_v02", "    absorber thickness is: %lf", absorberThickness);
    dd4hep::printout(dd4hep::INFO, "HCalTileBarrel_o1_v02", "    number of radiation length is: %lf",
                     nRadiationLengths);
    dd4hep::printout(dd4hep::INFO, "HCalTileBarrel_o1_v02", "    number of interaction length is: %lf",
                     nInteractionLengths);

    caloLayer.distance = layerInnerRadii.at(idxLayer);   // radius of the current layer
    caloLayer.sensitive_thickness = difference_bet_r1r2; // radial dimension of the current layer
    // caloLayer.sensitive_thickness	= thickness_sen;
    caloLayer.absorberThickness = absorberThickness;

    caloLayer.inner_thickness = difference_bet_r1r2 / 2.0;
    caloLayer.inner_nRadiationLengths = nRadiationLengths / 2.0;
    caloLayer.inner_nInteractionLengths = nInteractionLengths / 2.0;
    caloLayer.outer_nRadiationLengths = nRadiationLengths / 2.0;
    caloLayer.outer_nInteractionLengths = nInteractionLengths / 2.0;
    caloLayer.outer_thickness = difference_bet_r1r2 / 2;

    if (seg_phitheta) {
      // cells have all the same phi-theta nominal sizes, so can pass dummy cell ID (it is ignored)
      std::vector<double> cellSizeVector = seg_phitheta->cellDimensions(0);
      double cellSizeTheta = cellSizeVector[1];
      double cellSizePhi = cellSizeVector[0];
      dd4hep::printout(dd4hep::INFO, "HCalTileBarrel_o1_v02", "    cell sizes in theta, phi: %lf rad , %lf rad",
                       cellSizeTheta, cellSizePhi);
      caloLayer.cellSize0 = cellSizeTheta;
      caloLayer.cellSize1 = cellSizePhi;
    } else if (seg_phirow) {
      // the merging of the rows into cells can differ layer by layer so need to pass
      // cellID with layer field properly filled in order to get good dimension.
      dd4hep::CellID cID = 0;
      encoder.set(cID, layerFieldName, idxLayer);
      std::vector<double> cellSizeVector = seg_phirow->cellDimensions(cID);
      double cellSizeZ = cellSizeVector[0];
      double cellSizePhi = cellSizeVector[1];
      dd4hep::printout(dd4hep::INFO, "HCalTileBarrel_o1_v02", "    cell sizes in z, phi: %.4f cm , %.4f cm", cellSizeZ,
                       cellSizePhi);
      caloLayer.cellSize0 = cellSizeZ;
      caloLayer.cellSize1 = cellSizePhi;
    }
    caloData->layers.push_back(caloLayer);
  }

  // Set type flags
  dd4hep::xml::setDetectorTypeFlag(xmlDet, caloDetElem);

  return caloDetElem;
}
} // namespace det

DECLARE_DETELEMENT(HCalTileBarrel_o1_v02, det::createHCal)

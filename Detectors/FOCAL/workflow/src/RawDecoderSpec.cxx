// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "Framework/CCDBParamSpec.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/ControlService.h"
#include "Framework/InputRecordWalker.h"
#include "Framework/DataRefUtils.h"
#include "Framework/Logger.h"
#include "Framework/WorkflowSpec.h"

#include "CommonConstants/Triggers.h"
#include "DetectorsRaw/RDHUtils.h"

#include "FOCALReconstruction/PadWord.h"
#include "FOCALReconstruction/HCALWord.h"          // <-- ADD
#include "FOCALWorkflow/RawDecoderSpec.h"

#include "ITSMFTReconstruction/GBTWord.h"

#include <algorithm>                               // <-- ADD (std::find, std::sort, std::min)
#include <iostream>
#include <sstream>
#include <set>

using namespace o2::focal::reco_workflow;
using namespace o2::focal;

void RawDecoderSpec::init(framework::InitContext& ctx)
{
  if (ctx.options().get<bool>("filterIncomplete")) {
    LOG(info) << "Enabling filtering of incomplete events in the pixel data";
    mFilterIncomplete = true;
  }
  if (ctx.options().get<bool>("displayInconsistent")) {
    LOG(info) << "Display additional information in case of inconsistency between pixel links";
    mDisplayInconsistent = true;
  }

  auto mappingfile = ctx.options().get<std::string>("pixelmapping");
  PixelMapper::MappingType_t mappingtype = PixelMapper::MappingType_t::MAPPING_UNKNOWN;
  auto chiptype = ctx.options().get<std::string>("pixeltype");
  if (chiptype == "IB") {
    LOG(info) << "Using mapping type: IB";
    mappingtype = PixelMapper::MappingType_t::MAPPING_IB;
  } else if (chiptype == "OB") {
    LOG(info) << "Using mapping type: OB";
    mappingtype = PixelMapper::MappingType_t::MAPPING_OB;
  } else {
    LOG(fatal) << "Unkown mapping type for pixels: " << chiptype;
  }

  if (mappingfile == "default") {
    LOG(info) << "Using default pixel mapping for pixel type " << chiptype;
    mPixelMapping = std::make_unique<PixelMapper>(mappingtype);
  } else {
    LOG(info) << "Using user-defined mapping: " << mappingfile;
    mPixelMapping = std::make_unique<PixelMapper>(PixelMapper::MappingType_t::MAPPING_UNKNOWN);
    mPixelMapping->setMappingFile(mappingfile, mappingtype);
  }
}

void RawDecoderSpec::run(framework::ProcessingContext& ctx)
{
  LOG(info) << "Running FOCAL decoding";
  resetContainers();

  mTimeframeHasPadData = false;
  mTimeframeHasPixelData = false;
  mTimeframeHasHcalData = false;                  // <-- ADD

  // FEE IDs (keep placeholders if you don't have final values yet)
  constexpr uint16_t FEE_PADS = 0xcafe;
  constexpr uint16_t FEE_HCAL = 0xbeef;           // <-- ADD (TODO: replace with real HCAL FEE ID)

  int inputs = 0;
  std::vector<char> rawbuffer;
  uint16_t currentfee = 0;
  o2::InteractionRecord currentIR;

  std::unordered_map<int, int> numHBFFEE, numEventsFEE;
  std::unordered_map<int, std::vector<int>> numEventsHBFFEE;

  int numHBFPadsTF = 0, numEventsPadsTF = 0;
  int numHBFHcalTF = 0, numEventsHcalTF = 0;      // <-- ADD

  std::vector<int> expectFEEs;

  for (const auto& rawData : framework::InputRecordWalker(ctx.inputs())) {
    if (rawData.header != nullptr && rawData.payload != nullptr) {
      const auto payloadSize = o2::framework::DataRefUtils::getPayloadSize(rawData);
      auto header = o2::framework::DataRefUtils::getHeader<o2::header::DataHeader*>(rawData);
      LOG(debug) << "Channel " << header->dataOrigin.str << "/" << header->dataDescription.str << "/" << header->subSpecification;

      gsl::span<const char> databuffer(rawData.payload, payloadSize);
      int currentpos = 0;
      bool firstHBF = true;

      while (currentpos < databuffer.size()) {
        auto rdh = reinterpret_cast<const o2::header::RDHAny*>(databuffer.data() + currentpos);
        if (mDebugMode) {
          o2::raw::RDHUtils::printRDH(rdh);
        }

        if (o2::raw::RDHUtils::getMemorySize(rdh) > o2::raw::RDHUtils::getHeaderSize(rdh)) {
          auto payloadsize = o2::raw::RDHUtils::getMemorySize(rdh) - o2::raw::RDHUtils::getHeaderSize(rdh);
          int endpoint = static_cast<int>(o2::raw::RDHUtils::getEndPointID(rdh));

          auto fee = o2::raw::RDHUtils::getFEEID(rdh);

          LOG(debug) << "Next RDH: ";
          LOG(debug) << "Found fee                   0x" << std::hex << fee << std::dec
                     << " (System " << (fee == FEE_PADS ? "Pads" : (fee == FEE_HCAL ? "HCAL" : "Pixels")) << ")";
          LOG(debug) << "Found trigger BC:           " << o2::raw::RDHUtils::getTriggerBC(rdh);
          LOG(debug) << "Found trigger Oribt:        " << o2::raw::RDHUtils::getTriggerOrbit(rdh);
          LOG(debug) << "Found payload size:         " << payloadsize;
          LOG(debug) << "Found offset to next:       " << o2::raw::RDHUtils::getOffsetToNext(rdh);
          LOG(debug) << "Stop bit:                   " << (o2::raw::RDHUtils::getStop(rdh) ? "yes" : "no");

          // ---- FIX: compute word size per FEE (pads/hcal/pixels)
          size_t wordSize = 0;
          if (fee == FEE_PADS) {
            wordSize = sizeof(o2::focal::PadGBTWord);
          } else if (fee == FEE_HCAL) {
            wordSize = sizeof(o2::focal::HCALGBTWord);
          } else {
            wordSize = sizeof(o2::itsmft::GBTWord);
          }
          LOG(debug) << "Number of GBT words:        " << (wordSize ? (payloadsize / wordSize) : 0);

          auto page_payload = databuffer.subspan(currentpos + o2::raw::RDHUtils::getHeaderSize(rdh), payloadsize);
          std::copy(page_payload.begin(), page_payload.end(), std::back_inserter(rawbuffer));
        }

        auto trigger = o2::raw::RDHUtils::getTriggerType(rdh);
        if (trigger & o2::trigger::HB) {
          if (o2::raw::RDHUtils::getStop(rdh)) {
            LOG(debug) << "Stop bit received - processing payload";

            if (rawbuffer.size()) {
              // PAD
              if (currentfee == FEE_PADS) {
                if (mUsePadData) {
                  LOG(debug) << "Processing PAD data";
                  auto nEventsPads = decodePadData(rawbuffer, currentIR);
                  mTimeframeHasPadData = true;

                  auto found = mNumEventsHBFPads.find(nEventsPads);
                  if (found != mNumEventsHBFPads.end()) {
                    found->second += 1;
                  } else {
                    mNumEventsHBFPads.insert({nEventsPads, 1});
                  }

                  numEventsPadsTF += nEventsPads;
                  numHBFPadsTF++;
                }
              }
              // HCAL
              else if (currentfee == FEE_HCAL) {
                if (mUseHcalData) {
                  LOG(debug) << "Processing HCAL data";
                  auto nEventsHcal = decodeHcalData(rawbuffer, currentIR);
                  mTimeframeHasHcalData = true;

                  auto found = mNumEventsHBFHcal.find(nEventsHcal);
                  if (found != mNumEventsHBFHcal.end()) {
                    found->second += 1;
                  } else {
                    mNumEventsHBFHcal.insert({nEventsHcal, 1});
                  }

                  numEventsHcalTF += nEventsHcal;
                  numHBFHcalTF++;
                }
              }
              // PIXELS (all other FEEs)
              else {
                if (mUsePixelData) {
                  auto feeID = o2::raw::RDHUtils::getFEEID(rdh);

                  if (firstHBF) {
                    if (std::find(expectFEEs.begin(), expectFEEs.end(), feeID) == expectFEEs.end()) {
                      expectFEEs.emplace_back(feeID);
                    }
                    firstHBF = false;
                  }

                  LOG(debug) << "Processing Pixel data from FEE " << feeID;
                  auto neventsPixels = decodePixelData(rawbuffer, currentIR, feeID);
                  mTimeframeHasPixelData = true;

                  auto found = numHBFFEE.find(feeID);
                  if (found != numHBFFEE.end()) {
                    found->second++;
                  } else {
                    numHBFFEE.insert({feeID, 1});
                  }

                  auto evFound = numEventsFEE.find(feeID);
                  if (evFound != numEventsFEE.end()) {
                    evFound->second += neventsPixels;
                  } else {
                    numEventsFEE.insert({feeID, neventsPixels});
                  }

                  auto evHBFFound = numEventsHBFFEE.find(feeID);
                  if (evHBFFound != numEventsHBFFEE.end()) {
                    evHBFFound->second.push_back(neventsPixels);
                  } else {
                    numEventsHBFFEE.insert({feeID, std::vector<int>{neventsPixels}});
                  }
                }
              }
            } else {
              LOG(debug) << "Payload size 0 - skip empty HBF";
            }

            rawbuffer.clear();
          } else {
            currentIR.bc = o2::raw::RDHUtils::getTriggerBC(rdh);
            currentIR.orbit = o2::raw::RDHUtils::getTriggerOrbit(rdh);
            currentfee = o2::raw::RDHUtils::getFEEID(rdh);
            LOG(debug) << "New HBF " << currentIR.orbit << " / " << currentIR.bc << ", FEE 0x" << std::hex << currentfee << std::dec;
          }
        }

        currentpos += o2::raw::RDHUtils::getOffsetToNext(rdh);
      }
    } else {
      LOG(error) << "Input " << inputs << ": Either header or payload is nullptr";
    }
    inputs++;
  }

  int numHBFPixelsTF = 0;
  if (mTimeframeHasPixelData) {
    if (!consistencyCheckPixelFEE(numHBFFEE)) {
      LOG(alarm) << "Mismatch in number of HBF / TF between pixel FEEs";
      if (mDisplayInconsistent) {
        printCounters(numHBFFEE);
      }
      mNumInconsistencyPixelHBF++;
    }

    numHBFPixelsTF = maxCounter(numHBFFEE);
    mNumHBFPixels += numHBFPixelsTF;

    if (!consistencyCheckPixelFEE(numEventsFEE)) {
      LOG(alarm) << "Mismatch in number of events / TF between pixel FEEs";
      if (mDisplayInconsistent) {
        printCounters(numEventsFEE);
      }
      mNumInconsistencyPixelEvent++;
    }

    mNumEventsPixels += maxCounter(numEventsFEE);

    if (!checkEventsHBFConsistency(numEventsHBFFEE)) {
      LOG(alarm) << "Mistmatch number of events / HBF between pixel FEEs";
      if (mDisplayInconsistent) {
        printEvents(numEventsHBFFEE);
      }
      mNumInconsistencyPixelEventHBF++;
    }

    fillPixelEventHBFCount(numEventsHBFFEE);

    if (mFilterIncomplete) {
      LOG(debug) << "Filtering incomplete pixel events";
      for (auto& hbf : mHBFs) {
        auto numErased = filterIncompletePixelsEventsHBF(hbf.second, expectFEEs);
        mNumEventsPixels -= numErased;
      }
    }
  }

  LOG(info) << "Found " << mHBFs.size() << " HBFs in timeframe";
  LOG(debug) << "EventBuilder: Pixels: " << (mTimeframeHasPixelData ? "yes" : "no");
  LOG(debug) << "EventBuilder: Pads:   " << (mTimeframeHasPadData ? "yes" : "no");
  LOG(debug) << "EventBuilder: HCAL:   " << (mTimeframeHasHcalData ? "yes" : "no"); // <-- ADD

  buildEvents();
  LOG(info) << "Found " << mOutputTriggerRecords.size() << " events in timeframe";

  sendOutput(ctx);

  mNumEventsPads += numEventsPadsTF;
  mNumHBFPads += numHBFPadsTF;

  mNumEventsHcal += numEventsHcalTF;              // <-- ADD
  mNumHBFHcal += numHBFHcalTF;                    // <-- ADD

  mNumTimeframes++;

  auto foundHBFperTFPads = mNumHBFperTFPads.find(numHBFPadsTF);
  if (foundHBFperTFPads != mNumHBFperTFPads.end()) {
    foundHBFperTFPads->second++;
  } else {
    mNumHBFperTFPads.insert({numHBFPadsTF, 1});
  }

  auto foundHBFperTFHcal = mNumHBFperTFHcal.find(numHBFHcalTF);  // <-- ADD
  if (foundHBFperTFHcal != mNumHBFperTFHcal.end()) {
    foundHBFperTFHcal->second++;
  } else {
    mNumHBFperTFHcal.insert({numHBFHcalTF, 1});
  }

  auto foundHBFperTFPixels = mNumHBFperTFPixels.find(numHBFPixelsTF);
  if (foundHBFperTFPixels != mNumHBFperTFPixels.end()) {
    foundHBFperTFPixels->second++;
  } else {
    mNumHBFperTFPixels.insert({numHBFPixelsTF, 1});
  }
}

void RawDecoderSpec::endOfStream(o2::framework::EndOfStreamContext& ec)
{
  std::cout << "Number of timeframes:             " << mNumTimeframes << std::endl;
  std::cout << "Number of pad HBFs:               " << mNumHBFPads << std::endl;
  std::cout << "Number of hcal HBFs:              " << mNumHBFHcal << std::endl; // <-- ADD
  std::cout << "Number of pixel HBFs:             " << mNumHBFPixels << std::endl;

  for (auto& [hbfs, tfs] : mNumHBFperTFPads) {
    std::cout << "Pads - Number of TFs with " << hbfs << " HBFs: " << tfs << std::endl;
  }
  for (auto& [hbfs, tfs] : mNumHBFperTFHcal) {                                  // <-- ADD
    std::cout << "HCAL - Number of TFs with " << hbfs << " HBFs: " << tfs << std::endl;
  }
  for (auto& [hbfs, tfs] : mNumHBFperTFPixels) {
    std::cout << "Pixels - Number of TFs with " << hbfs << " HBFs: " << tfs << std::endl;
  }

  std::cout << "Number of pad events:             " << mNumEventsPads << std::endl;
  std::cout << "Number of pixel events:           " << mNumEventsPixels << std::endl;
  std::cout << "Number of hcal events:            " << mNumEventsHcal << std::endl; // <-- ADD

  for (auto& [nevents, nHBF] : mNumEventsHBFPads) {
    std::cout << "Number of HBFs with " << nevents << " pad events:    " << nHBF << std::endl;
  }
  for (auto& [nevents, nHBF] : mNumEventsHBFHcal) {                              // <-- ADD
    std::cout << "Number of HBFs with " << nevents << " hcal events:    " << nHBF << std::endl;
  }
  for (auto& [nevents, nHBF] : mNumEventsHBFPixels) {
    std::cout << "Number of HBFs with " << nevents << " pixel events:    " << nHBF << std::endl;
  }

  std::cout << "Number of inconsistencies between pixel FEEs: "
            << mNumInconsistencyPixelHBF << " HBFs, "
            << mNumInconsistencyPixelEvent << " events, "
            << mNumInconsistencyPixelEventHBF << " events / HBF" << std::endl;
}

void RawDecoderSpec::sendOutput(framework::ProcessingContext& ctx)
{
  ctx.outputs().snapshot(framework::Output{o2::header::gDataOriginFOC, "PADLAYERS", mOutputSubspec}, mOutputPadLayers);
  ctx.outputs().snapshot(framework::Output{o2::header::gDataOriginFOC, "HCALPCBS",  mOutputSubspec}, mOutputHcalPCBs); // <-- ADD
  ctx.outputs().snapshot(framework::Output{o2::header::gDataOriginFOC, "PIXELHITS", mOutputSubspec}, mOutputPixelHits);
  ctx.outputs().snapshot(framework::Output{o2::header::gDataOriginFOC, "PIXELCHIPS", mOutputSubspec}, mOutputPixelChips);
  ctx.outputs().snapshot(framework::Output{o2::header::gDataOriginFOC, "TRIGGERS",  mOutputSubspec}, mOutputTriggerRecords);
}

void RawDecoderSpec::resetContainers()
{
  mHBFs.clear();
  mOutputPadLayers.clear();
  mOutputHcalPCBs.clear();                         // <-- ADD
  mOutputPixelChips.clear();
  mOutputPixelHits.clear();
  mOutputTriggerRecords.clear();
}

int RawDecoderSpec::decodePadData(const gsl::span<const char> padWords, o2::InteractionRecord& hbIR)
{
  LOG(debug) << "Decoding pad data for Orbit " << hbIR.orbit << ", BC " << hbIR.bc;
  constexpr std::size_t EVENTSIZEPADGBT = 1180,
                        EVENTSIZECHAR = EVENTSIZEPADGBT * sizeof(PadGBTWord) / sizeof(char);
  auto nevents = padWords.size() / (EVENTSIZECHAR);
  for (int ievent = 0; ievent < nevents; ievent++) {
    decodePadEvent(padWords.subspan(EVENTSIZECHAR * ievent, EVENTSIZECHAR), hbIR);
  }
  return nevents;
}

void RawDecoderSpec::decodePadEvent(const gsl::span<const char> padWords, o2::InteractionRecord& hbIR)
{
  gsl::span<const PadGBTWord> padWordsGBT(reinterpret_cast<const PadGBTWord*>(padWords.data()),
                                         padWords.size() / sizeof(PadGBTWord));
  mPadDecoder.reset();
  mPadDecoder.decodeEvent(padWordsGBT);

  std::map<o2::InteractionRecord, HBFData>::iterator foundHBF = mHBFs.find(hbIR);
  if (foundHBF == mHBFs.end()) {
    HBFData nexthbf;
    auto res = mHBFs.insert({hbIR, nexthbf});
    foundHBF = res.first;
  }
  foundHBF->second.mPadEvents.push_back(createPadLayerEvent(mPadDecoder.getData()));
}

// ========================
// HCAL decode + conversion
// ========================

int RawDecoderSpec::decodeHcalData(const gsl::span<const char> hcalWords, o2::InteractionRecord& hbIR)
{
  LOG(debug) << "Decoding hcal data for Orbit " << hbIR.orbit << ", BC " << hbIR.bc;

  constexpr std::size_t EVENTSIZEHCALGBT = 1180;
  constexpr std::size_t EVENTSIZECHAR = EVENTSIZEHCALGBT * sizeof(HCALGBTWord) / sizeof(char);

  auto nevents = hcalWords.size() / EVENTSIZECHAR;
  for (int ievent = 0; ievent < nevents; ++ievent) {
    decodeHcalEvent(hcalWords.subspan(EVENTSIZECHAR * ievent, EVENTSIZECHAR), hbIR);
  }
  return nevents;
}

void RawDecoderSpec::decodeHcalEvent(const gsl::span<const char> hcalWords, o2::InteractionRecord& hbIR)
{
  gsl::span<const HCALGBTWord> hcalWordsGBT(reinterpret_cast<const HCALGBTWord*>(hcalWords.data()),
                                           hcalWords.size() / sizeof(HCALGBTWord));

  mHcalDecoder.reset();
  mHcalDecoder.decodeEvent(hcalWordsGBT);

  auto foundHBF = mHBFs.find(hbIR);
  if (foundHBF == mHBFs.end()) {
    auto res = mHBFs.insert({hbIR, HBFData{}});
    foundHBF = res.first;
  }

  // IMPORTANT: HCALDecoder exposes getData() (HCALData), convert it to HCALPCBEvent array
  foundHBF->second.mHCALEvents.push_back(createHcalPCBEvent(mHcalDecoder.getData()));
}

std::array<HCALEvent, constants::HCAL_NPCBS>
RawDecoderSpec::createHcalPCBEvent(const o2::focal::HCALData& data) const
{
  std::array<HCALEvent, constants::HCAL_NPCBS> result{};
  std::array<uint8_t, 8> triggertimes{};

  // NOTE:
  // This assumes PCB 1..4 correspond to ASIC indices 0..3 in HCALData.
  // If a different mapping is needed  one couldchange `asicIndex` accordingly.
  for (std::size_t ipcb = 0; ipcb < constants::HCAL_NPCBS; ++ipcb) {
    const int asicIndex = static_cast<int>(ipcb); // <-- PCB=(ipcb+1)

    const auto& cont = data.getDataForASIC(asicIndex);
    const auto& asic = cont.getASIC();

    for (int ihalf = 0; ihalf < o2::focal::HCALASICData::NHALVES; ++ihalf) {
      const auto header = asic.getHeader(ihalf);
      const auto calib  = asic.getCalib(ihalf);
      const auto cmn    = asic.getCMN(ihalf);

      result[ipcb].setHeader(ihalf, header.getHeader(), header.getBCID(),
                             header.getWadd(), header.getFourbit(), header.getTrailer());
      result[ipcb].setCalib(ihalf, calib.getADC(), calib.getTOA(), calib.getTOT());
      result[ipcb].setCMN(ihalf, cmn.getADC(), cmn.getTOA(), cmn.getTOT());
    }

    for (int ich = 0; ich < o2::focal::HCALASICData::NCHANNELS; ++ich) {
      const auto ch = asic.getChannel(ich);
      result[ipcb].setChannel(ich, ch.getADC(), ch.getTOA(), ch.getTOT());
    }

    const auto triggers = cont.getTriggerWords();
    const auto nwin = std::min<std::size_t>(triggers.size(), constants::HCAL_WINDOW_LENGTH);

    for (std::size_t window = 0; window < nwin; ++window) {
      triggertimes.fill(0);
      triggertimes[0] = triggers[window].mTrigger0;
      triggertimes[1] = triggers[window].mTrigger1;
      triggertimes[2] = triggers[window].mTrigger2;
      triggertimes[3] = triggers[window].mTrigger3;
      triggertimes[4] = triggers[window].mTrigger4;
      triggertimes[5] = triggers[window].mTrigger5;
      triggertimes[6] = triggers[window].mTrigger6;
      triggertimes[7] = triggers[window].mTrigger7;

      result[ipcb].setTrigger(window, triggers[window].mHeader0, triggers[window].mHeader1, triggertimes);
    }
  }

  return result;
}

// ========================
// Pixel decode unchanged
// ========================

int RawDecoderSpec::decodePixelData(const gsl::span<const char> pixelWords, o2::InteractionRecord& hbIR, int feeID)
{
  // unchanged from your original
  LOG(debug) << "Decoding pixel data for Orbit " << hbIR.orbit << ", BC " << hbIR.bc;

  gsl::span<const o2::itsmft::GBTWord> pixelpayload(reinterpret_cast<const o2::itsmft::GBTWord*>(pixelWords.data()),
                                                   pixelWords.size() / sizeof(o2::itsmft::GBTWord));
  LOG(debug) << pixelWords.size() << " Bytes -> " << pixelpayload.size() << " GBT words";
  mPixelDecoder.reset();
  mPixelDecoder.decodeEvent(pixelpayload);

  std::map<o2::InteractionRecord, HBFData>::iterator foundHBF = mHBFs.end();

  int nevents = 0;
  for (auto& [trigger, chipdata] : mPixelDecoder.getChipData()) {
    LOG(debug) << "Found trigger orbit " << trigger.orbit << ", BC " << trigger.bc;
    if (trigger.orbit != hbIR.orbit) {
      LOG(debug) << "FEE 0x" << std::hex << feeID << std::dec
                 << ": Discarding spurious trigger with Orbit " << trigger.orbit
                 << " (HB " << hbIR.orbit << ")";
      continue;
    }

    if (foundHBF == mHBFs.end()) {
      foundHBF = mHBFs.find(hbIR);
      if (foundHBF == mHBFs.end()) {
        auto res = mHBFs.insert({hbIR, HBFData{}});
        foundHBF = res.first;
      }
    }

    auto triggerfound = std::find(foundHBF->second.mPixelTriggers.begin(),
                                  foundHBF->second.mPixelTriggers.end(), trigger);

    if (triggerfound != foundHBF->second.mPixelTriggers.end()) {
      auto index = triggerfound - foundHBF->second.mPixelTriggers.begin();
      for (const auto& chip : chipdata) {
        try {
          auto chipPosition = mPixelMapping->getPosition(feeID, chip);
          fillChipToLayer(foundHBF->second.mPixelEvent[index][chipPosition.mLayer], chip, feeID);
        } catch (PixelMapper::InvalidChipException& e) {
          LOG(warning) << e;
        }
      }
      foundHBF->second.mFEEs[index].push_back(feeID);
    } else {
      std::array<PixelLayerEvent, constants::PIXELS_NLAYERS> nextevent;
      foundHBF->second.mPixelEvent.push_back(nextevent);
      foundHBF->second.mPixelTriggers.push_back(trigger);

      auto& current = foundHBF->second.mPixelEvent.back();
      for (const auto& chip : chipdata) {
        try {
          auto chipPosition = mPixelMapping->getPosition(feeID, chip);
          fillChipToLayer(current[chipPosition.mLayer], chip, feeID);
        } catch (PixelMapper::InvalidChipException& e) {
          LOG(warning) << e;
        }
      }

      foundHBF->second.mFEEs.push_back({feeID});
    }
    nevents++;
  }
  return nevents;
}

std::array<o2::focal::PadLayerEvent, o2::focal::constants::PADS_NLAYERS>
RawDecoderSpec::createPadLayerEvent(const o2::focal::PadData& data) const
{
  // unchanged from your original
  std::array<PadLayerEvent, constants::PADS_NLAYERS> result;
  std::array<uint8_t, 8> triggertimes;
  for (std::size_t ilayer = 0; ilayer < constants::PADS_NLAYERS; ilayer++) {
    auto& asic = data.getDataForASIC(ilayer).getASIC();
    for (std::size_t ihalf = 0; ihalf < constants::PADLAYER_MODULE_NHALVES; ihalf++) {
      auto header = asic.getHeader(ihalf);
      auto calib = asic.getCalib(ihalf);
      auto cmn = asic.getCMN(ihalf);
      result[ilayer].setHeader(ihalf, header.getHeader(), header.getBCID(), header.getWadd(),
                               header.getFourbit(), header.getTrailer());
      result[ilayer].setCalib(ihalf, calib.getADC(), calib.getTOA(), calib.getTOT());
      result[ilayer].setCMN(ihalf, cmn.getADC(), cmn.getTOA(), cmn.getTOT());
    }
    for (std::size_t ichannel = 0; ichannel < constants::PADLAYER_MODULE_NCHANNELS; ichannel++) {
      auto channel = asic.getChannel(ichannel);
      result[ilayer].setChannel(ichannel, channel.getADC(), channel.getTOA(), channel.getTOT());
    }
    auto triggers = data.getDataForASIC(ilayer).getTriggerWords();
    for (std::size_t window = 0; window < constants::PADLAYER_WINDOW_LENGTH; window++) {
      std::fill(triggertimes.begin(), triggertimes.end(), 0);
      triggertimes[0] = triggers[window].mTrigger0;
      triggertimes[1] = triggers[window].mTrigger1;
      triggertimes[2] = triggers[window].mTrigger2;
      triggertimes[3] = triggers[window].mTrigger3;
      triggertimes[4] = triggers[window].mTrigger4;
      triggertimes[5] = triggers[window].mTrigger5;
      triggertimes[6] = triggers[window].mTrigger6;
      triggertimes[7] = triggers[window].mTrigger7;
      result[ilayer].setTrigger(window, triggers[window].mHeader0, triggers[window].mHeader1, triggertimes);
    }
  }
  return result;
}

void RawDecoderSpec::fillChipToLayer(o2::focal::PixelLayerEvent& pixellayer, const o2::focal::PixelChip& chipData, int feeID)
{
  // unchanged
  pixellayer.addChip(feeID, chipData.mLaneID, chipData.mChipID, chipData.mStatusCode, chipData.mHits);
}

void RawDecoderSpec::fillEventPixeHitContainer(std::vector<PixelHit>& eventHits,
                                               std::vector<PixelChipRecord>& eventChips,
                                               const PixelLayerEvent& pixelLayer, int layerIndex)
{
  // unchanged
  for (auto& chip : pixelLayer.getChips()) {
    auto starthits = eventHits.size();
    auto& chipHits = chip.mHits;
    std::copy(chipHits.begin(), chipHits.end(), std::back_inserter(eventHits));
    eventChips.emplace_back(layerIndex, chip.mFeeID, chip.mLaneID, chip.mChipID,
                            chip.mStatusCode, starthits, chipHits.size());
  }
}

void RawDecoderSpec::buildEvents()
{
  LOG(debug) << "Start building events" << std::endl;

  for (const auto& [hbir, hbf] : mHBFs) {

    // ---- CASE 1: PAD + PIXEL + HCAL
    if (mTimeframeHasPadData && mTimeframeHasPixelData && mTimeframeHasHcalData) {
      const auto n = hbf.mPixelTriggers.size();

      if (hbf.mPadEvents.size() != n ||
          hbf.mPixelEvent.size() != n ||
          hbf.mHCALEvents.size() != n) {
        LOG(error) << "Inconsistent number of events in HBF: "
                   << "pads=" << hbf.mPadEvents.size()
                   << ", pixels=" << hbf.mPixelEvent.size()
                   << ", hcal=" << hbf.mHCALEvents.size()
                   << ", triggers=" << n;
        continue;
      }

      for (std::size_t itrg = 0; itrg < n; ++itrg) {
        auto startPads  = mOutputPadLayers.size();
        auto startHCAL  = mOutputHcalPCBs.size();
        auto startHits  = mOutputPixelHits.size();
        auto startChips = mOutputPixelChips.size();

        // PAD: 18 layers (index 0..17 => layer 1..18 by convention)
        for (std::size_t ilayer = 0; ilayer < constants::PADS_NLAYERS; ++ilayer) {
          mOutputPadLayers.push_back(hbf.mPadEvents[itrg][ilayer]);
        }

        // HCAL: 4 pcbs (index 0..3 => pcb 1..4 by convention)
        for (std::size_t ipcb = 0; ipcb < constants::HCAL_NPCBS; ++ipcb) {
          mOutputHcalPCBs.push_back(hbf.mHCALEvents[itrg][ipcb]);
        }

        std::vector<PixelHit> eventHits;
        std::vector<PixelChipRecord> eventPixels;

        for (std::size_t ilayer = 0; ilayer < constants::PIXELS_NLAYERS; ++ilayer) {
          fillEventPixeHitContainer(eventHits, eventPixels, hbf.mPixelEvent[itrg][ilayer], ilayer);
        }

        std::copy(eventHits.begin(), eventHits.end(), std::back_inserter(mOutputPixelHits));
        std::copy(eventPixels.begin(), eventPixels.end(), std::back_inserter(mOutputPixelChips));

        // TriggerRecord ctor MUST be updated to include HCAL range:
        // (BC, firstPad, nPad, firstHcal, nHcal, firstChip, nChip, firstHit, nHit)
        mOutputTriggerRecords.emplace_back(hbf.mPixelTriggers[itrg],
                                           startPads,  constants::PADS_NLAYERS,
                                           startHCAL,  constants::HCAL_NPCBS,
                                           startChips, static_cast<int>(eventPixels.size()),
                                           startHits,  static_cast<int>(eventHits.size()));
      }
      continue;
    }

    // ---- CASE 2: PIXEL only
    if (mTimeframeHasPixelData) {
      if (hbf.mPixelEvent.size() != hbf.mPixelTriggers.size()) {
        LOG(error) << "Inconsistent number of pixel events (" << hbf.mPixelEvent.size()
                   << ") and triggers (" << hbf.mPixelTriggers.size() << ") in HBF";
        continue;
      }

      for (std::size_t itrg = 0; itrg < hbf.mPixelTriggers.size(); ++itrg) {
        auto startPads  = mOutputPadLayers.size();
        auto startHCAL  = mOutputHcalPCBs.size();
        auto startHits  = mOutputPixelHits.size();
        auto startChips = mOutputPixelChips.size();

        std::vector<PixelHit> eventHits;
        std::vector<PixelChipRecord> eventPixels;

        for (std::size_t ilayer = 0; ilayer < constants::PIXELS_NLAYERS; ++ilayer) {
          fillEventPixeHitContainer(eventHits, eventPixels, hbf.mPixelEvent[itrg][ilayer], ilayer);
        }

        std::copy(eventHits.begin(), eventHits.end(), std::back_inserter(mOutputPixelHits));
        std::copy(eventPixels.begin(), eventPixels.end(), std::back_inserter(mOutputPixelChips));

        mOutputTriggerRecords.emplace_back(hbf.mPixelTriggers[itrg],
                                           startPads,  0,
                                           startHCAL,  0,
                                           startChips, static_cast<int>(eventPixels.size()),
                                           startHits,  static_cast<int>(eventHits.size()));
      }
      continue;
    }

    // ---- CASE 3: PAD only
    if (mTimeframeHasPadData) {
      for (std::size_t itrg = 0; itrg < hbf.mPadEvents.size(); ++itrg) {
        auto startPads  = mOutputPadLayers.size();
        auto startHCAL  = mOutputHcalPCBs.size();
        auto startHits  = mOutputPixelHits.size();
        auto startChips = mOutputPixelChips.size();

        for (std::size_t ilayer = 0; ilayer < constants::PADS_NLAYERS; ++ilayer) {
          mOutputPadLayers.push_back(hbf.mPadEvents[itrg][ilayer]);
        }

        mOutputTriggerRecords.emplace_back(hbir,
                                           startPads, constants::PADS_NLAYERS,
                                           startHCAL, 0,
                                           startChips, 0,
                                           startHits, 0);
      }
      continue;
    }

    // ---- CASE 4: HCAL only
    if (mTimeframeHasHcalData) {
      for (std::size_t itrg = 0; itrg < hbf.mHCALEvents.size(); ++itrg) {
        auto startPads  = mOutputPadLayers.size();
        auto startHCAL  = mOutputHcalPCBs.size();
        auto startHits  = mOutputPixelHits.size();
        auto startChips = mOutputPixelChips.size();

        for (std::size_t ipcb = 0; ipcb < constants::HCAL_NPCBS; ++ipcb) {
          mOutputHcalPCBs.push_back(hbf.mHCALEvents[itrg][ipcb]);
        }

        mOutputTriggerRecords.emplace_back(hbir,
                                           startPads, 0,
                                           startHCAL, constants::HCAL_NPCBS,
                                           startChips, 0,
                                           startHits, 0);
      }
      continue;
    }
  }
}

int RawDecoderSpec::filterIncompletePixelsEventsHBF(HBFData& data, const std::vector<int>& expectFEEs)
{
  // unchanged from your original
  auto same = [](const std::vector<int>& lhs, const std::vector<int>& rhs) -> bool {
    bool missing = false;
    for (auto entry : lhs) {
      if (std::find(rhs.begin(), rhs.end(), entry) == rhs.end()) {
        missing = true;
        break;
      }
    }
    if (!missing) {
      for (auto entry : rhs) {
        if (std::find(lhs.begin(), lhs.end(), entry) == lhs.end()) {
          missing = true;
          break;
        }
      }
    }
    return missing;
  };

  std::vector<int> indexIncomplete;
  for (auto index = 0; index < data.mFEEs.size(); index++) {
    if (data.mFEEs[index].size() != expectFEEs.size()) {
      indexIncomplete.emplace_back(index);
      continue;
    }
    if (!same(data.mFEEs[index], expectFEEs)) {
      indexIncomplete.emplace_back(index);
    }
  }
  if (indexIncomplete.size()) {
    std::sort(indexIncomplete.begin(), indexIncomplete.end(), std::less<>());
    for (auto indexIter = indexIncomplete.rbegin(); indexIter != indexIncomplete.rend(); indexIter++) {
      data.mPixelEvent.erase(data.mPixelEvent.begin() + *indexIter);
      data.mPixelTriggers.erase(data.mPixelTriggers.begin() + *indexIter);
      data.mFEEs.erase(data.mFEEs.begin() + *indexIter);
    }
  }
  return indexIncomplete.size();
}

bool RawDecoderSpec::consistencyCheckPixelFEE(const std::unordered_map<int, int>& counters) const
{
  // unchanged
  bool initialized = false;
  bool discrepancy = false;
  int current = -1;
  for (auto& [fee, value] : counters) {
    if (!initialized) {
      current = value;
      initialized = true;
    }
    if (value != current) {
      discrepancy = true;
      break;
    }
  }
  return !discrepancy;
}

bool RawDecoderSpec::checkEventsHBFConsistency(const std::unordered_map<int, std::vector<int>>& counters) const
{
  // unchanged
  bool initialized = false;
  bool discrepancy = false;
  std::vector<int> current;
  for (auto& [fee, events] : counters) {
    if (!initialized) {
      current = events;
      initialized = true;
    }
    if (events != current) {
      discrepancy = true;
    }
  }
  return !discrepancy;
}

int RawDecoderSpec::maxCounter(const std::unordered_map<int, int>& counters) const
{
  // unchanged
  int maxCounter = 0;
  for (auto& [fee, counter] : counters) {
    if (counter > maxCounter) {
      maxCounter = counter;
    }
  }
  return maxCounter;
}

void RawDecoderSpec::printCounters(const std::unordered_map<int, int>& counters) const
{
  // unchanged
  for (auto& [fee, counter] : counters) {
    LOG(info) << "  FEE 0x" << std::hex << fee << std::dec << ": " << counter << " counts ...";
  }
}

void RawDecoderSpec::printEvents(const std::unordered_map<int, std::vector<int>>& counters) const
{
  // unchanged
  for (auto& [fee, events] : counters) {
    std::stringstream stringbuilder;
    bool first = true;
    for (auto ev : events) {
      if (first) {
        first = false;
      } else {
        stringbuilder << ", ";
      }
      stringbuilder << ev;
    }
    LOG(info) << "  FEE 0x" << std::hex << fee << std::dec << ": " << stringbuilder.str() << " events ...";
  }
}

void RawDecoderSpec::fillPixelEventHBFCount(const std::unordered_map<int, std::vector<int>>& counters)
{
  // unchanged
  int maxFEE = 0;
  int current = -1;
  for (auto& [fee, events] : counters) {
    int sum = 0;
    for (auto ev : events) {
      sum += ev;
    }
    if (sum > current) {
      maxFEE = fee;
      current = sum;
    }
  }
  auto en = counters.find(maxFEE);
  if (en != counters.end()) {
    for (auto nEventsPixels : en->second) {
      auto found = mNumEventsHBFPixels.find(nEventsPixels);
      if (found != mNumEventsHBFPixels.end()) {
        found->second += 1;
      } else {
        mNumEventsHBFPixels.insert({nEventsPixels, 1});
      }
    }
  }
}

// ----------------------------
// getRawDecoderSpec updated API
// ----------------------------
o2::framework::DataProcessorSpec o2::focal::reco_workflow::getRawDecoderSpec(bool askDISTSTF,
                                                                            uint32_t outputSubspec,
                                                                            bool usePadData,
                                                                            bool usePixelData,
                                                                            bool useHcalData,   // <-- ADD
                                                                            bool debugMode)
{
  constexpr auto originFOC = o2::header::gDataOriginFOC;
  std::vector<o2::framework::OutputSpec> outputs;

  outputs.emplace_back(originFOC, "PADLAYERS", outputSubspec, o2::framework::Lifetime::Timeframe);
  outputs.emplace_back(originFOC, "HCALPCBS",  outputSubspec, o2::framework::Lifetime::Timeframe); // <-- ADD
  outputs.emplace_back(originFOC, "PIXELHITS", outputSubspec, o2::framework::Lifetime::Timeframe);
  outputs.emplace_back(originFOC, "PIXELCHIPS", outputSubspec, o2::framework::Lifetime::Timeframe);
  outputs.emplace_back(originFOC, "TRIGGERS",  outputSubspec, o2::framework::Lifetime::Timeframe);

  std::vector<o2::framework::InputSpec> inputs{
    {"stf", o2::framework::ConcreteDataTypeMatcher{originFOC, o2::header::gDataDescriptionRawData},
     o2::framework::Lifetime::Timeframe}
  };
  if (askDISTSTF) {
    inputs.emplace_back("stdDist", "FLP", "DISTSUBTIMEFRAME", 0, o2::framework::Lifetime::Timeframe);
  }

  return o2::framework::DataProcessorSpec{
    "FOCALRawDecoderSpec",
    inputs,
    outputs,
    // IMPORTANT: pass useHcalData into task constructor
    o2::framework::adaptFromTask<o2::focal::reco_workflow::RawDecoderSpec>(
      outputSubspec, usePadData, usePixelData, useHcalData, debugMode
    ),
    o2::framework::Options{
      {"filterIncomplete",    o2::framework::VariantType::Bool,   false, {"Filter incomplete pixel events"}},
      {"displayInconsistent", o2::framework::VariantType::Bool,   false, {"Display information about inconsistent timeframes"}},
      {"pixeltype",           o2::framework::VariantType::String, "OB",  {"Pixel mapping type"}},
      {"pixelmapping",        o2::framework::VariantType::String, "default", {"File with pixel mapping"}}
    }
  };
}

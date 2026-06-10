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

#include <cstdint>
#include <cstdio>
#include <iostream>
#include <vector>

#include <gsl/span>

#include <fairlogger/Logger.h>
#include "FOCALReconstruction/HCALDecoder.h"
#include "FOCALReconstruction/HCALWord.h"

using namespace o2::focal;

void HCALDecoder::reset()
{
  mData.reset();
}

void HCALDecoder::decodeEvent(gsl::span<const HCALGBTWord> gbtdata)
{
  LOG(debug) << "decoding hcal data of size " << gbtdata.size() << "  GBT words - " << gbtdata.size() * sizeof(HCALGBTWord) / sizeof(uint64_t) << " 64 bit words";
  // first 35 GBT words : ASIC data
  // Other words: Trigger data

  // Note from Nina: 
  // The decoder was throwing an index error when run through the workflow because of a mismatch between the 64 HCAL channels and the 72 Pad channels
  // Quick fix: Changed hardcoded channel numbers in this file to be compatible with the number of HCAL channels, now the workflows execute and produce a root file
  // Thing to think about: Do we need to fundamentally change the data format to account for half a chip corresponding to 1 PCB in HCAL (move to 4x8 matrix of 32 channels with HCAL_NPCBS = 8 in the dataformat, instead of the current setup which is an 8x8 matrix of 64 channels with HCAL_NPCBS = 4 in the dataformat)

  std::size_t asicsize = 35 * HCALData::NASICS; 
  auto asicwords = gbtdata.subspan(0, asicsize);
  auto triggerwords = gbtdata.subspan(asicsize, gbtdata.size() - asicsize);
  for (int iasic = 0; iasic < HCALData::NASICS; iasic++) {
    // First part: ASIC words
    auto& asicdata = mData[iasic].getASIC();
    auto wordsthisAsic = asicwords.subspan(iasic * 35, 35);
    auto headerwords = wordsthisAsic[0].getASICData<HCALASICHeader>();
    asicdata.setFirstHeader(headerwords[0]);
    asicdata.setSecondHeader(headerwords[1]);
    int nchannels = 0;
    for (auto& datawords : wordsthisAsic.subspan(1, 32)) {
      for (auto& channelword : datawords.getASICData<HCALASICChannel>()) {
        asicdata.setChannel(channelword, nchannels);
        nchannels++;
      }
    }
    asicdata.setCMNs(wordsthisAsic[33].getASICData<HCALASICChannel>());
    asicdata.setCalibs(wordsthisAsic[34].getASICData<HCALASICChannel>());

    // Second part: Trigger words
    auto wordsTriggerThisAsic = triggerwords.subspan(iasic * mWin_dur, mWin_dur);
    auto& asiccont = mData[iasic];
    for (const auto trgword : wordsTriggerThisAsic) {
      asiccont.appendTriggerWord(trgword.getTriggerData());
    }
  }
}

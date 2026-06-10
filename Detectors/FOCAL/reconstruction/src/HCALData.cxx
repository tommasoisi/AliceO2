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
#include <algorithm>
#include "FOCALReconstruction/HCALData.h"

using namespace o2::focal;

HCALASICData::HCALASICData(HCALASICHeader firstheader, HCALASICHeader secondheader)
{
  setFirstHeader(firstheader);
  setSecondHeader(secondheader);
}

void HCALASICData::setHeader(HCALASICHeader header, int index)
{
  if (index >= NHALVES) {
    throw IndexException(index, NHALVES);
  }
  mHeaders[index] = header;
}

void HCALASICData::setChannel(HCALASICChannel data, int index)
{
  if (index >= NCHANNELS) {
    throw IndexException(index, NCHANNELS);
  }
  mChannels[index] = data;
}

void HCALASICData::setChannels(const gsl::span<const HCALASICChannel> channels)
{
  std::copy(channels.begin(), channels.end(), mChannels.begin());
}

void HCALASICData::setCMN(HCALASICChannel data, int index)
{
  if (index >= NHALVES) {
    throw IndexException(index, NHALVES);
  }
  mCMNChannels[index] = data;
}

void HCALASICData::setCMNs(const gsl::span<const HCALASICChannel> channels)
{
  std::copy(channels.begin(), channels.end(), mCMNChannels.begin());
}

void HCALASICData::setCalib(HCALASICChannel data, int index)
{
  if (index >= NHALVES) {
    throw IndexException(index, NHALVES);
  }
  mCalibChannels[index] = data;
}

void HCALASICData::setCalibs(const gsl::span<const HCALASICChannel> channels)
{
  std::copy(channels.begin(), channels.end(), mCalibChannels.begin());
}

HCALASICHeader HCALASICData::getHeader(int index) const
{
  if (index >= NHALVES) {
    throw IndexException(index, NHALVES);
  }
  return mHeaders[index];
}

gsl::span<const HCALASICHeader> HCALASICData::getHeaders() const
{
  return mHeaders;
}

gsl::span<const HCALASICChannel> HCALASICData::getChannels() const
{
  return mChannels;
}

HCALASICChannel HCALASICData::getChannel(int index) const
{
  if (index >= NCHANNELS) {
    throw IndexException(index, NCHANNELS);
  }
  return mChannels[index];
}

HCALASICChannel HCALASICData::getCalib(int index) const
{
  if (index >= NHALVES) {
    throw IndexException(index, NHALVES);
  }
  return mCalibChannels[index];
}

gsl::span<const HCALASICChannel> HCALASICData::getCalibs() const
{
  return mCalibChannels;
}

HCALASICChannel HCALASICData::getCMN(int index) const
{
  if (index >= NHALVES) {
    throw IndexException(index, NHALVES);
  }
  return mCMNChannels[index];
}

gsl::span<const HCALASICChannel> HCALASICData::getCMNs() const
{
  return mCMNChannels;
}

void HCALASICData::reset()
{
  /*
  for (auto& header : mHeaders) {
    header.mData = 0;
  }
  for (auto& channel : mChannels) {
    channel.mData = 0;
  }
  for (auto& channel : mCalibChannels) {
    channel.mData = 0;
  }
  for (auto& channel : mCMNChannels) {
    channel.mData = 0;
  }
  */
  std::fill(mHeaders.begin(), mHeaders.end(), HCALASICHeader(0));
  std::fill(mChannels.begin(), mChannels.end(), HCALASICChannel(0));
  std::fill(mCalibChannels.begin(), mCalibChannels.end(), HCALASICChannel(0));
  std::fill(mCMNChannels.begin(), mCMNChannels.end(), HCALASICChannel(0));
}

gsl::span<const HCALTriggerWord> HCALASICContainer::getTriggerWords() const
{
  return mTriggerData;
}

void HCALASICContainer::appendTriggerWords(gsl::span<const HCALTriggerWord> triggerwords)
{
  std::copy(triggerwords.begin(), triggerwords.end(), std::back_inserter(mTriggerData));
}

void HCALASICContainer::appendTriggerWord(HCALTriggerWord triggerword)
{
  mTriggerData.emplace_back(triggerword);
}

void HCALASICContainer::reset()
{
  mASIC.reset();
  mTriggerData.clear();
}

void HCALData::reset()
{
  for (auto& asic : mASICs) {
    asic.reset();
  }
}

const HCALASICContainer& HCALData::getDataForASIC(int index) const
{
  if (index >= NASICS) {
    throw IndexException(index, NASICS);
  }
  return mASICs[index];
}

HCALASICContainer& HCALData::getDataForASIC(int index)
{
  if (index >= NASICS) {
    throw IndexException(index, NASICS);
  }
  return mASICs[index];
}
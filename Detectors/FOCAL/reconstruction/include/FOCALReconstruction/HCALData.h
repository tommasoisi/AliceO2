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
#ifndef ALICEO2_FOCAL_HCALDATA_H
#define ALICEO2_FOCAL_HCALDATA_H

#include <array>
#include <exception>
#include <string>
#include <vector>

#include <gsl/span>

#include "Rtypes.h"

#include "FOCALReconstruction/HCALWord.h"

namespace o2::focal
{
class HCALASICData
{
 public:
  class IndexException : public std::exception
  {
   public:
    IndexException() = default;
    IndexException(int index, int maxindex) : std::exception(), mIndex(index), mMaxIndex(maxindex)
    {
      mMessage = "Invalid index " + std::to_string(mIndex) + ", max " + std::to_string(mMaxIndex);
    }
    ~IndexException() noexcept final = default;

    const char* what() const noexcept final
    {
      return mMessage.data();
    }
    int getIndex() const { return mIndex; }
    int getMaxIndex() const { return mMaxIndex; }

   private:
    int mIndex;
    int mMaxIndex;
    std::string mMessage;
  };

  static constexpr int NCHANNELS = 64;
  static constexpr int NHALVES = 2;

  HCALASICData() = default;
  HCALASICData(HCALASICHeader firstheader, HCALASICHeader secondheader);

  void setFirstHeader(HCALASICHeader header) { setHeader(header, 0); }
  void setSecondHeader(HCALASICHeader header) { setHeader(header, 1); }
  void setHeader(HCALASICHeader header, int index);

  void setChannel(HCALASICChannel data, int index);
  void setChannels(const gsl::span<const HCALASICChannel> channels);

  void setFirstCMN(HCALASICChannel data) { setCMN(data, 0); }
  void setSecondCMN(HCALASICChannel data) { setCMN(data, 1); }
  void setCMN(HCALASICChannel data, int index);
  void setCMNs(const gsl::span<const HCALASICChannel> channels);

  void setFirstCalib(HCALASICChannel data) { setCalib(data, 0); }
  void setSecondCalib(HCALASICChannel data) { setCalib(data, 1); }
  void setCalib(HCALASICChannel data, int index);
  void setCalibs(const gsl::span<const HCALASICChannel> channels);

  HCALASICHeader getFirstHeader() const { return getHeader(0); }
  HCALASICHeader getSecondHeader() const { return getHeader(1); }
  HCALASICHeader getHeader(int index) const;
  gsl::span<const HCALASICHeader> getHeaders() const;

  gsl::span<const HCALASICChannel> getChannels() const;
  HCALASICChannel getChannel(int index) const;

  HCALASICChannel getFirstCalib() const { return getCalib(0); }
  HCALASICChannel getSecondCalib() const { return getCalib(1); }
  HCALASICChannel getCalib(int index) const;
  gsl::span<const HCALASICChannel> getCalibs() const;

  HCALASICChannel getFirstCMN() const { return getCMN(0); }
  HCALASICChannel getSecondCMN() const { return getCMN(1); }
  HCALASICChannel getCMN(int index) const;
  gsl::span<const HCALASICChannel> getCMNs() const;

  void reset();

 private:
  std::array<HCALASICHeader, NHALVES> mHeaders;
  std::array<HCALASICChannel, NCHANNELS> mChannels;
  std::array<HCALASICChannel, NHALVES> mCalibChannels;
  std::array<HCALASICChannel, NHALVES> mCMNChannels;

  ClassDefNV(HCALASICData, 1);
};

class HCALASICContainer
{
 public:
  HCALASICContainer() = default;
  ~HCALASICContainer() = default;

  const HCALASICData& getASIC() const { return mASIC; }
  HCALASICData& getASIC() { return mASIC; }
  gsl::span<const HCALTriggerWord> getTriggerWords() const;
  void appendTriggerWords(gsl::span<const HCALTriggerWord> triggerwords);
  void appendTriggerWord(HCALTriggerWord triggerword);
  void reset();

 private:
  HCALASICData mASIC;
  std::vector<HCALTriggerWord> mTriggerData;

  ClassDefNV(HCALASICContainer, 1);
};

class HCALData
{
 public:
  class IndexException : public std::exception
  {
   public:
    IndexException() = default;
    IndexException(int index, int maxindex) : std::exception(), mIndex(index), mMaxIndex(maxindex)
    {
      mMessage = "Invalid index " + std::to_string(mIndex) + ", max " + std::to_string(mMaxIndex);
    }
    ~IndexException() noexcept final = default;
    const char* what() const noexcept final
    {
      return mMessage.data();
    }

    int getIndex() const { return mIndex; }
    int getMaxIndex() const { return mMaxIndex; }

   private:
    int mIndex;
    int mMaxIndex;
    std::string mMessage;
  };
  static constexpr int NASICS = 4;

  HCALData() = default;
  ~HCALData() = default;

  const HCALASICContainer& operator[](int index) const { return getDataForASIC(index); }
  HCALASICContainer& operator[](int index) { return getDataForASIC(index); }

  const HCALASICContainer& getDataForASIC(int index) const;
  HCALASICContainer& getDataForASIC(int index);
  void reset();

 private:
  std::array<HCALASICContainer, NASICS> mASICs;

  ClassDefNV(HCALData, 1);
};

} // namespace o2::focal
#endif // ALICEO2_FOCAL_HCALDATA_H
#pragma once

#include <chrono>
#include <nlohmann/json.hpp>
#include <fstream>
#include <NTL/ZZ_p.h>
#include <NTL/BasicThreadPool.h>

#include "netmp.h"

using namespace nlohmann;
using namespace NTL;
using namespace io;

struct TimePoint {
  using timepoint_t = std::chrono::high_resolution_clock::time_point;
  using timeunit_t = std::chrono::duration<double, std::milli>;

  TimePoint();
  double operator-(const TimePoint& rhs) const;

  timepoint_t time;
};

struct CommPoint {
  std::vector<uint64_t> stats;

  explicit CommPoint(NetIOMP& network);
  std::vector<uint64_t> operator-(const CommPoint& rhs) const;
};

class StatsPoint {
  TimePoint tpoint_;
  CommPoint cpoint_;

 public:
  explicit StatsPoint(io::NetIOMP& network);
  nlohmann::json operator-(const StatsPoint& rhs);
};

bool saveJson(const nlohmann::json& data, const std::string& fpath);
int64_t peakVirtualMemory();
int64_t peakResidentSetSize();
void initNTL(size_t num_threads);

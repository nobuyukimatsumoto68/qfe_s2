// rng.h

#pragma once

#include <cstdio>
#include <sstream>
#include <random>
#include <string>

class QfeRng {

public:
  QfeRng(int seed = 12345678);
  void WriteRng(FILE* file);
  void ReadRng(FILE* file);
  double RandReal(double min = 0.0, double max = 1.0);
  double RandNormal(double mean = 0.0, double stddev = 1.0);
  int RandInt(int min, int max);
  bool RandBool();

  std::mt19937 gen;
};

void QfeRng::WriteRng(FILE* file){

  std::stringstream rng_string;
  rng_string.setf(std::ios::dec | std::ios::left);
  rng_string.fill(' ');
  rng_string << gen;
  fputs(rng_string.str().c_str(), file);
}

void QfeRng::ReadRng(FILE* file){

  std::vector<char> rng_buf(0x2000);
  fgets(rng_buf.data(), 0x2000, file);
  std::stringstream rng_string(rng_buf.data());
  rng_string.setf(std::ios::dec);
  rng_string >> gen;
}

QfeRng::QfeRng(int seed) {
  gen = std::mt19937(seed);
}

double QfeRng::RandReal(double min, double max) {
  std::uniform_real_distribution<double> dist(min, max);
  return dist(gen);
}

double QfeRng::RandNormal(double mean, double stddev) {
  std::normal_distribution<double> dist(mean, stddev);
  return dist(gen);
}

int QfeRng::RandInt(int min, int max) {
  std::uniform_int_distribution<int> dist(min, max);
  return dist(gen);
}

bool QfeRng::RandBool() {
  std::uniform_int_distribution<int> dist(0, 1);
  return (dist(gen) == 1);
}

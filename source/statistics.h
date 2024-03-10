// statistics.h

#pragma once

#include <cmath>
#include <complex>
#include <deque>
#include <vector>

typedef std::complex<double> Complex;

double AutocorrTime(std::vector<double>& a);

class QfeMeasReal {

public:
  QfeMeasReal();
  void WriteMeasurement(FILE* file);
  void ReadMeasurement(FILE* file);
  void Reset();
  void Measure(double value, bool calc_autocorr = true);
  double Mean() const;
  double Error() const;
  double AutocorrFront();
  double AutocorrBack();

  double sum;  // sum of measurements
  double sum2;  // sum of squared measurements
  double last;  // most recent measurement
  int n;  // number of measurements
  std::vector<double> first_1000;  // first 1000 measurements
  std::deque<double> last_1000;  // last 1000 measurements
};

QfeMeasReal::QfeMeasReal() {
  Reset();
}

void QfeMeasReal::WriteMeasurement(FILE* file) {
  fprintf(file, "%.16e %.16e %d\n", Mean(), Error(), n);
}

void QfeMeasReal::ReadMeasurement(FILE* file) {
  Reset();
  double mean;
  double err;
  fscanf(file, "%lf %lf %d\n", &mean, &err, &n);
  double count = double(n);
  sum = mean * count;
  sum2 = (err * err * count + mean * mean) * count;
}

void QfeMeasReal::Reset() {
  sum = 0.0;
  sum2 = 0.0;
  last = 0.0;
  n = 0;
  first_1000.clear();
  last_1000.clear();
}

void QfeMeasReal::Measure(double value, bool calc_autocorr) {
  sum += value;
  sum2 += value * value;
  last = value;
  n++;

  if (calc_autocorr) {
    if (first_1000.size() < 1000) first_1000.push_back(value);
    last_1000.push_back(value);
    if (last_1000.size() > 1000) last_1000.pop_front();
  }
}

double QfeMeasReal::Mean() const {
  if (n == 0) return 0.0;
  return sum / double(n);
}

double QfeMeasReal::Error() const {
  if (n == 0) return 0.0;
  double mean = sum / double(n);
  double mean2 = sum2 / double(n);
  return sqrt((mean2 - mean * mean) / double(n));
}

double QfeMeasReal::AutocorrFront() {
  if (n == 0) return 0.0;
  return AutocorrTime(first_1000);
}

double QfeMeasReal::AutocorrBack() {
  if (n == 0) return 0.0;
  // convert last_1000 to a vector
  std::deque<double>::iterator it = last_1000.begin();
  std::vector<double> a;
  while (it != last_1000.end()) {
    a.push_back(*it++);
  }
  return AutocorrTime(a);
}

class QfeMeasComplex {

public:
  QfeMeasComplex();
  void Reset();
  void Measure(Complex value, bool calc_autocorr = true);
  Complex Mean();
  Complex Error();
  double AutocorrFront();
  double AutocorrBack();

  QfeMeasReal real_part;
  QfeMeasReal imag_part;
  Complex last;
};

QfeMeasComplex::QfeMeasComplex() {
  Reset();
}

void QfeMeasComplex::Reset() {
  real_part.Reset();
  imag_part.Reset();
  last = 0.0;
}

void QfeMeasComplex::Measure(Complex value, bool calc_autocorr) {
  real_part.Measure(real(value), calc_autocorr);
  imag_part.Measure(imag(value), calc_autocorr);
  last = value;
}

Complex QfeMeasComplex::Mean() {
  return Complex(real_part.Mean(), imag_part.Mean());
}

Complex QfeMeasComplex::Error() {
  return Complex(real_part.Error(), imag_part.Error());
}

double QfeMeasComplex::AutocorrFront() {
  return fmax(real_part.AutocorrFront(), imag_part.AutocorrFront());
}

double QfeMeasComplex::AutocorrBack() {
  return fmax(real_part.AutocorrBack(), imag_part.AutocorrBack());
}

// TODO: come up with a way to consolidate these so it's not just the same
//       function over and over again

double Mean(std::vector<double>& a) {
  double sum = 0.0;
  for (int i = 0; i < a.size(); i++) {
    sum += a[i];
  }
  return sum / double(a.size());
}

double LogMean(std::vector<double>& a) {
  return log(Mean(a));
}

double U4(std::vector<double>& m2, std::vector<double>& m4) {
  double m2_mean = Mean(m2);
  double m4_mean = Mean(m4);

  return 1.5 * (1.0 - m4_mean / (3.0 * m2_mean * m2_mean));
}

double Susceptibility(std::vector<double>& m2, std::vector<double>& m) {
  double m2_mean = Mean(m2);
  double m_mean = Mean(m);

  return m2_mean - m_mean * m_mean;
}

double JackknifeMean(std::vector<double>& a) {
  int n = a.size();
	double mean = Mean(a);
	double err = 0.0;

	for (int i = 0; i < n; i++) {
    std::vector<double> a_del = a;
    a_del.erase(a_del.begin() + i);
    double diff = Mean(a_del) - mean;
    err += diff * diff;
  }

	err = sqrt((double(n) - 1.0) / double(n) * err);
	return err;
}

double JackknifeLogMean(std::vector<double>& a) {
  int n = a.size();
  std::vector<double> bin_means(n);

  for (int i = 0; i < n; i++) {
    std::vector<double> a_del = a;
    a_del.erase(a_del.begin() + i);
    bin_means[i] = LogMean(a_del);
  }

  double mean = Mean(bin_means);
	double err = 0.0;

	for (int i = 0; i < n; i++) {
    double diff = bin_means[i] - mean;
    err += diff * diff;
  }

	err = sqrt((double(n) - 1.0) / double(n) * err);
	return err;
}

double JackknifeU4(std::vector<double>& m2, std::vector<double>& m4) {
  int n = m2.size();
	std::vector<double> bin_means(n);

  for (int i = 0; i < n; i++) {
    std::vector<double> m2_del = m2;
    std::vector<double> m4_del = m4;
    m2_del.erase(m2_del.begin() + i);
    m4_del.erase(m4_del.begin() + i);
    bin_means[i] = U4(m2_del, m4_del);
  }

  double mean = Mean(bin_means);
	double err = 0.0;

	for (int i = 0; i < n; i++) {
    double diff = bin_means[i] - mean;
    err += diff * diff;
  }

	err = sqrt((double(n) - 1.0) / double(n) * err);
	return err;
}

double JackknifeSusceptibility(std::vector<double>& m2, std::vector<double>& m) {
  int n = m2.size();
  std::vector<double> bin_means(n);

  for (int i = 0; i < n; i++) {
    std::vector<double> m2_del = m2;
    std::vector<double> m_del = m;
    m2_del.erase(m2_del.begin() + i);
    m_del.erase(m_del.begin() + i);
    bin_means[i] = Susceptibility(m2_del, m_del);
  }

	double mean = Mean(bin_means);
	double err = 0.0;

	for (int i = 0; i < n; i++) {
    double diff = bin_means[i] - mean;
    err += diff * diff;
  }

	err = sqrt((double(n) - 1.0) / double(n) * err);
	return err;
}

double AutocorrGamma(std::vector<double>& a, int n) {
  int N = a.size();
  double result = 0.0;
  double mean = Mean(a);
  int start = 0;
  int end = N - n;

  if (n < 0) {
    start = -n;
    end = N;
  }

  for (int i = start; i < end; i++) {
    result += (a[i] - mean) * (a[i + n] - mean);
  }

  return result / double(end - start);
}

double AutocorrTime(std::vector<double>& a) {
  double Gamma0 = AutocorrGamma(a, 0);
  double result = 0.5 * Gamma0;

  for (int n = 1; n < a.size(); n++) {
    double curGamma = AutocorrGamma(a, n);
    if (curGamma < 0.0) break;
    result += curGamma;
  }

  return result / Gamma0;
}

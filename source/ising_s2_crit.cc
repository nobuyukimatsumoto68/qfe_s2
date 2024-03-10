// ising_s2_crit.cc

#include <getopt.h>

#include <Eigen/Dense>
#include <cmath>
#include <complex>
#include <cstdio>
#include <string>
#include <vector>
#include <iostream>

#include "ising.h"
#include "s2.h"
#include "statistics.h"
#include "timer.h"
#include "util.h"

int main(int argc, char* argv[]) {
  // default parameters
  int q = 5;

  int n_refine = 1; // orig
  // int n_refine = 3;

  unsigned int seed = 1234u;
  bool cold_start = false;
  int l_max = 6;
  int n_therm = 2000;
  int n_traj = 20000;
  int n_skip = 10;
  int n_wolff = 5;
  int n_metropolis = 4;
  double wall_time = 0.0;
  std::string data_dir = "ising_s2_crit";
  std::string orbit_path = "";
  // std::string orbit_path = "/projectnb/qfe/nmatsumo/QFE_copy/grp/elem/o3q5.dat";
  // std::string orbit_path = "/projectnb/qfe/nmatsumo/QFE_copy/grp/site_g/o3q5k3.dat";
  // std::string orbit_path = "/projectnb/qfe/nmatsumo/QFE_copy/s2_refine/q5k3_orbit.dat";

  const struct option long_options[] = {
      {"q", required_argument, 0, 'q'},
      {"n_refine", required_argument, 0, 'N'},
      {"seed", required_argument, 0, 'S'},
      {"cold_start", no_argument, 0, 'C'},
      {"l_max", required_argument, 0, 'l'},
      {"n_therm", required_argument, 0, 'h'},
      {"n_traj", required_argument, 0, 't'},
      {"n_skip", required_argument, 0, 's'},
      {"n_wolff", required_argument, 0, 'w'},
      {"n_metropolis", required_argument, 0, 'e'},
      {"data_dir", required_argument, 0, 'd'},
      {"wall_time", required_argument, 0, 'W'},
      {"orbit_path", required_argument, 0, 'o'},
      {0, 0, 0, 0}};

  const char* short_options = "q:N:S:Cl:h:t:s:w:e:d:W:o:";

  while (true) {
    int o = 0;
    int c = getopt_long(argc, argv, short_options, long_options, &o);
    if (c == -1) break;

    switch (c) {
      case 'q':
        q = atoi(optarg);
        break;
      case 'N':
        n_refine = atoi(optarg);
        break;
      case 'S':
        seed = atol(optarg);
        break;
      case 'C':
        cold_start = true;
        break;
      case 'l':
        l_max = atoi(optarg);
        break;
      case 'h':
        n_therm = atoi(optarg);
        break;
      case 't':
        n_traj = atoi(optarg);
        break;
      case 's':
        n_skip = atoi(optarg);
        break;
      case 'w':
        n_wolff = atoi(optarg);
        break;
      case 'e':
        n_metropolis = atoi(optarg);
        break;
      case 'd':
        data_dir = optarg;
        break;
      case 'W':
        wall_time = std::stod(optarg);
        break;
      case 'o':
        orbit_path = optarg;
        break;
      default:
        break;
    }
  }

  std::string run_id = string_format("q%dk%d", q, n_refine);
  printf("run_id: %s\n", run_id.c_str());

  printf("n_therm: %d\n", n_therm);
  printf("n_traj: %d\n", n_traj);
  printf("n_skip: %d\n", n_skip);
  printf("n_wolff: %d\n", n_wolff);
  printf("n_metropolis: %d\n", n_metropolis);

  // number of spherical harmonics to measure
  int n_ylm = ((l_max + 1) * (l_max + 2)) / 2;
  printf("l_max: %d\n", l_max);
  printf("n_ylm: %d\n", n_ylm);

  // create a refined triangular lattice
  QfeLatticeS2 lattice(q, n_refine);

  if (!orbit_path.empty()) {
    std::cout << "orbit path = " << orbit_path << std::endl;
    FILE* orbit_file = fopen(orbit_path.c_str(), "r");
    assert(orbit_file != nullptr);
    lattice.ReadOrbits(orbit_file);
    fclose(orbit_file);
    std::cout << "orbit path opened" << std::endl;
  }
  else{
    std::cout << "orbit path=" << orbit_path << std::endl;
  }
  lattice.SeedRng(seed);
  lattice.UpdateWeights();

  printf("n_refine: %d\n", n_refine);
  printf("q: %d\n", q);
  printf("total sites: %d\n", lattice.n_sites);

  double vol = lattice.vol;
  double vol_sq = vol * vol;

  QfeIsing field(&lattice, 1.0);

  // check if there is an rng file
  std::string rng_path =
    string_format("%s/%s/%s_rng_%08X.dat", data_dir.c_str(), run_id.c_str(),
                  run_id.c_str(), seed);
  printf("opening file: %s\n", rng_path.c_str());
  FILE* rng_file = fopen(rng_path.c_str(), "r");
  if (rng_file != nullptr) {
    printf("loading rng state from file\n");
    lattice.rng.ReadRng(rng_file);
    fclose(rng_file);
  }
  std::cout << rng_path << std::endl;

  // check if there is a field checkpoint file
  std::string field_path =
      string_format("%s/%s/%s_field_%08X.dat", data_dir.c_str(), run_id.c_str(),
                    run_id.c_str(), seed);
  printf("opening file: %s\n", field_path.c_str());
  FILE* field_file = fopen(field_path.c_str(), "rb");
  if (field_file != nullptr) {
    printf("loading field from checkpoint file\n");
    field.ReadField(field_file);
    fclose(field_file);
  } else if (cold_start) {
    printf("cold start\n");
    field.ColdStart();
  } else {
    printf("hot start\n");
    field.HotStart();
  }

  printf("initial action: %.12f\n", field.Action());

  // numbering for adjacent sites, links, and faces
  //
  //                      s2
  //                    /    \
  //        f1a       /        \       f2a
  //               l1            l2
  //              /       f0       \
  //            /                    \
  //          s0 - - - -  l0  - - - - s1
  //            \                    /
  //              \       f1       /
  //               l3            l4
  //        f3b       \        /       f4b
  //                    \    /
  //                      s3
  //
  //

  // 3. this way uses dual links on the trianglular lattice with a kink
  for (int l0 = 0; l0 < lattice.n_links; l0++) {
    int f0 = lattice.links[l0].faces[0];
    int f1 = lattice.links[l0].faces[1];

    int s0 = lattice.links[l0].sites[0];
    int s1 = lattice.links[l0].sites[1];

    int s2;
    for (int e = 0; e < 3; e++) {
      s2 = lattice.faces[f0].sites[e];
      if (s2 != s0 && s2 != s1) break;
    }

    int s3;
    for (int e = 0; e < 3; e++) {
      s3 = lattice.faces[f1].sites[e];
      if (s3 != s0 && s3 != s1) break;
    }

    int l1 = lattice.FindLink(s0, s2);
    int l2 = lattice.FindLink(s1, s2);
    int l3 = lattice.FindLink(s0, s3);
    int l4 = lattice.FindLink(s1, s3);

    // calculate unit vectors from links to/from face circumcenters
    Vec3 v_f0 = lattice.FaceCircumcenter(f0);
    Vec3 v_f1 = lattice.FaceCircumcenter(f1);
    Vec3 v_l0 = 0.5 * (lattice.r[s0] + lattice.r[s1]);
    Vec3 v_l1 = 0.5 * (lattice.r[s0] + lattice.r[s2]);
    Vec3 v_l2 = 0.5 * (lattice.r[s1] + lattice.r[s2]);
    Vec3 v_l3 = 0.5 * (lattice.r[s0] + lattice.r[s3]);
    Vec3 v_l4 = 0.5 * (lattice.r[s1] + lattice.r[s3]);
    Vec3 v_l0_f0 = (v_f0 - v_l0).normalized();
    Vec3 v_l0_f1 = (v_f1 - v_l0).normalized();
    Vec3 v_f0_l1 = (v_l1 - v_f0).normalized();
    Vec3 v_f0_l2 = (v_l2 - v_f0).normalized();
    Vec3 v_f1_l3 = (v_l3 - v_f1).normalized();
    Vec3 v_f1_l4 = (v_l4 - v_f1).normalized();
    double cos1 = sqrt(0.5 * (1.0 + v_l0_f0.dot(v_f0_l1)));
    double cos2 = sqrt(0.5 * (1.0 + v_l0_f0.dot(v_f0_l2)));
    double cos3 = sqrt(0.5 * (1.0 + v_l0_f1.dot(v_f1_l3)));
    double cos4 = sqrt(0.5 * (1.0 + v_l0_f1.dot(v_f1_l4)));
    double cos12 = sqrt(0.5 * (1.0 - v_f0_l1.dot(v_f0_l2)));
    double cos34 = sqrt(0.5 * (1.0 - v_f1_l3.dot(v_f1_l4)));
    double cos_prod_num = cos1 * cos2 * cos3 * cos4;
    double cos_prod_den = cos12 * cos34;

    double len_l0_sq = lattice.EdgeSquared(l0);
    double len_l0 = sqrt(len_l0_sq);
    double len_l1 = lattice.EdgeLength(l1);
    double len_l2 = lattice.EdgeLength(l2);
    double len_l3 = lattice.EdgeLength(l3);
    double len_l4 = lattice.EdgeLength(l4);
    double len_f0 = len_l0 + len_l1 + len_l2;
    double len_f1 = len_l0 + len_l3 + len_l4;

    double tanh_sq_L_num = 4.0 * len_l0_sq * cos_prod_num;
    double tanh_sq_L_den = len_f0 * len_f1 * cos_prod_den;
    double L = atanh(sqrt(tanh_sq_L_num / tanh_sq_L_den));
    double K = 0.5 * asinh(1.0 / sinh(2.0 * L));

    // printf("%.12f %.12f\n", lattice.links[l0].wt, K);
    lattice.links[l0].wt = K;
  }

  // calculate spherical harmonics at each site. In all of these loops, (y_l,
  // y_m) are the integer eigenvalues of the spherical harmonics and y_i is a
  // sequential index over all of eigenvalues up to l_max (excluding negative
  // y_m).
  Eigen::MatrixXcd ylm(lattice.n_sites, n_ylm);
  for (int s = 0; s < lattice.n_sites; s++) {
    for (int y_i = 0, y_l = 0, y_m = 0; y_i < n_ylm; y_i++) {
      ylm(s, y_i) = lattice.CalcYlm(s, y_l, y_m);

      y_m++;
      if (y_m > y_l) {
        y_l++;
        y_m = 0;
      }
    }
  }

  // measurements
  std::vector<QfeMeasReal> legendre_2pt(l_max + 1);
  std::vector<QfeMeasReal> ylm_2pt(n_ylm);
  QfeMeasReal mag;     // magnetization
  QfeMeasReal mag_2;   // magnetization^2
  QfeMeasReal mag_4;   // magnetization^4
  QfeMeasReal mag_6;   // magnetization^6
  QfeMeasReal mag_8;   // magnetization^8
  QfeMeasReal mag_10;  // magnetization^10
  QfeMeasReal mag_12;  // magnetization^12
  QfeMeasReal action;
  QfeMeasReal cluster_size;
  QfeMeasReal accept_metropolis;

  // read bulk measurements
  std::string bulk_path =
      string_format("%s/%s/%s_bulk_%08X.dat", data_dir.c_str(), run_id.c_str(),
                    run_id.c_str(), seed);
  printf("opening file: %s\n", bulk_path.c_str());
  FILE* bulk_file = fopen(bulk_path.c_str(), "r");
  if (bulk_file != nullptr) {
    while (!feof(bulk_file)) {
      char meas_name[40];
      fscanf(bulk_file, "%s ", meas_name);
      if (strcmp(meas_name, "action") == 0) {
        action.ReadMeasurement(bulk_file);
      } else if (strcmp(meas_name, "mag") == 0) {
        mag.ReadMeasurement(bulk_file);
      } else if (strcmp(meas_name, "mag") == 0) {
        mag.ReadMeasurement(bulk_file);
      } else if (strcmp(meas_name, "mag^2") == 0) {
        mag_2.ReadMeasurement(bulk_file);
      } else if (strcmp(meas_name, "mag^4") == 0) {
        mag_4.ReadMeasurement(bulk_file);
      } else if (strcmp(meas_name, "mag^6") == 0) {
        mag_6.ReadMeasurement(bulk_file);
      } else if (strcmp(meas_name, "mag^8") == 0) {
        mag_8.ReadMeasurement(bulk_file);
      } else if (strcmp(meas_name, "mag^10") == 0) {
        mag_10.ReadMeasurement(bulk_file);
      } else if (strcmp(meas_name, "mag^12") == 0) {
        mag_12.ReadMeasurement(bulk_file);
      } else {
        printf("unknown measurement: %s\n", meas_name);
      }
    }
    fclose(bulk_file);
  }

  // read 2-point function legendre coefficients
  std::string legendre_2pt_path =
      string_format("%s/%s/%s_legendre_2pt_%08X.dat", data_dir.c_str(),
                    run_id.c_str(), run_id.c_str(), seed);
  printf("opening file: %s\n", legendre_2pt_path.c_str());
  FILE* legendre_2pt_file = fopen(legendre_2pt_path.c_str(), "r");
  if (legendre_2pt_file != nullptr) {
    int l;
    while (!feof(legendre_2pt_file)) {
      fscanf(legendre_2pt_file, "%d ", &l);
      if (l > l_max) {
        fscanf(legendre_2pt_file, "%*[^\n]");
      } else {
        legendre_2pt[l].ReadMeasurement(legendre_2pt_file);
      }
      fscanf(legendre_2pt_file, "\n");
    }
    fclose(legendre_2pt_file);
  }

  // read 2-point function ylm coefficients
  std::string ylm_2pt_path =
      string_format("%s/%s/%s_ylm_2pt_%08X.dat", data_dir.c_str(),
                    run_id.c_str(), run_id.c_str(), seed);
  printf("opening file: %s\n", ylm_2pt_path.c_str());
  FILE* ylm_2pt_file = fopen(ylm_2pt_path.c_str(), "r");
  if (ylm_2pt_file != nullptr) {
    int i_ylm, l, m;
    while (!feof(ylm_2pt_file)) {
      fscanf(ylm_2pt_file, "%d %d %d ", &i_ylm, &l, &m);
      if (i_ylm >= n_ylm) {
        fscanf(ylm_2pt_file, "%*[^\n]");
      } else {
        ylm_2pt[i_ylm].ReadMeasurement(ylm_2pt_file);
      }
      fscanf(ylm_2pt_file, "\n");
    }
    fclose(ylm_2pt_file);
  }

  Timer timer;

  for (int n = 0; n < (n_traj + n_therm); n++) {
    if (wall_time > 0.0 && timer.Duration() > wall_time) break;

    double metropolis_sum = 0.0;
    for (int j = 0; j < n_metropolis; j++) {
      metropolis_sum += field.Metropolis();
    }
    accept_metropolis.Measure(metropolis_sum);

    int cluster_size_sum = 0;
    for (int j = 0; j < n_wolff; j++) {
      cluster_size_sum += field.WolffUpdate();
    }
    cluster_size.Measure(double(cluster_size_sum) / vol);

    if (n % n_skip || n < n_therm) continue;

    // measure correlators
    std::vector<Complex> ylm_2pt_sum(n_ylm, 0.0);
    std::vector<Complex> ylm_4pt_sum(n_ylm, 0.0);
    double spin_sum = 0.0;

    for (int s = 0; s < lattice.n_sites; s++) {
      double wt_2pt = field.spin[s] * lattice.sites[s].wt;

      spin_sum += wt_2pt;

      for (int ylm_i = 0; ylm_i < n_ylm; ylm_i++) {
        Complex y = ylm(s, ylm_i);
        ylm_2pt_sum[ylm_i] += y * wt_2pt;
      }
    }

    double legendre_2pt_sum = 0.0;
    for (int ylm_i = 0, l = 0, m = 0; ylm_i < n_ylm; ylm_i++) {
      ylm_2pt[ylm_i].Measure(std::norm(ylm_2pt_sum[ylm_i]) / vol_sq);

      legendre_2pt_sum += ylm_2pt[ylm_i].last * (m == 0 ? 1.0 : 2.0);

      m++;
      if (m > l) {
        double coeff = 4.0 * M_PI / double(2 * l + 1);
        legendre_2pt[l].Measure(legendre_2pt_sum * coeff);
        legendre_2pt_sum = 0.0;
        l++;
        m = 0;
      }
    }

    double m = spin_sum / vol;
    double m_sq = m * m;
    mag.Measure(fabs(m));
    mag_2.Measure(m_sq);
    mag_4.Measure(m_sq * m_sq);
    mag_6.Measure(mag_4.last * m_sq);
    mag_8.Measure(mag_6.last * m_sq);
    mag_10.Measure(mag_8.last * m_sq);
    mag_12.Measure(mag_10.last * m_sq);
    action.Measure(field.Action());
    // printf("%06d %.12f %.4f %.4f\n", \
    //     n, action.last, \
    //     accept_metropolis.last, \
    //     cluster_size.last);
  }

  timer.Stop();
  printf("duration: %.6f\n", timer.Duration());

  printf("cluster_size/V: %.4f\n", cluster_size.Mean());
  printf("accept_metropolis: %.4f\n", accept_metropolis.Mean());

  // write rng state to file
  printf("writing rng state to file: %s\n", rng_path.c_str());
  rng_file = fopen(rng_path.c_str(), "w");
  // std::cout << "rng_path = " << rng_path << std::endl;
  assert(rng_file != nullptr);
  lattice.rng.WriteRng(rng_file);
  fclose(rng_file);

  // write field configuration to file
  printf("writing field to file: %s\n", field_path.c_str());
  field_file = fopen(field_path.c_str(), "wb");
  assert(field_file != nullptr);
  field.WriteField(field_file);
  fclose(field_file);

  double m_mean = mag.Mean();
  double m_err = mag.Error();
  double m2_mean = mag_2.Mean();
  double m2_err = mag_2.Error();
  double m4_mean = mag_4.Mean();
  double m4_err = mag_4.Error();

  // open an output file
  printf("opening file: %s\n", bulk_path.c_str());
  bulk_file = fopen(bulk_path.c_str(), "w");
  assert(bulk_file != nullptr);

  printf("action: %+.12e %.12e %.4f %.4f\n", action.Mean(), action.Error(),
         action.AutocorrFront(), action.AutocorrBack());
  fprintf(bulk_file, "action ");
  action.WriteMeasurement(bulk_file);

  printf("mag: %.12e %.12e %.4f %.4f\n", m_mean, m_err, mag.AutocorrFront(),
         mag.AutocorrBack());
  fprintf(bulk_file, "mag ");
  mag.WriteMeasurement(bulk_file);

  printf("m^2: %.12e %.12e %.4f %.4f\n", m2_mean, m2_err, mag_2.AutocorrFront(),
         mag_2.AutocorrBack());
  fprintf(bulk_file, "mag^2 ");
  mag_2.WriteMeasurement(bulk_file);

  printf("m^4: %.12e %.12e %.4f %.4f\n", m4_mean, m4_err, mag_4.AutocorrFront(),
         mag_4.AutocorrBack());
  fprintf(bulk_file, "mag^4 ");
  mag_4.WriteMeasurement(bulk_file);

  printf("m^6: %.12e %.12e %.4f %.4f\n", mag_6.Mean(), mag_6.Error(),
         mag_6.AutocorrFront(), mag_6.AutocorrBack());
  fprintf(bulk_file, "mag^6 ");
  mag_6.WriteMeasurement(bulk_file);

  printf("m^8: %.12e %.12e %.4f %.4f\n", mag_8.Mean(), mag_8.Error(),
         mag_8.AutocorrFront(), mag_8.AutocorrBack());
  fprintf(bulk_file, "mag^8 ");
  mag_8.WriteMeasurement(bulk_file);

  printf("m^10: %.12e %.12e %.4f %.4f\n", mag_10.Mean(), mag_10.Error(),
         mag_10.AutocorrFront(), mag_10.AutocorrBack());
  fprintf(bulk_file, "mag^10 ");
  mag_10.WriteMeasurement(bulk_file);

  printf("m^12: %.12e %.12e %.4f %.4f\n", mag_12.Mean(), mag_12.Error(),
         mag_12.AutocorrFront(), mag_12.AutocorrBack());
  fprintf(bulk_file, "mag^12 ");
  mag_12.WriteMeasurement(bulk_file);

  double U4_mean = 1.5 * (1.0 - m4_mean / (3.0 * m2_mean * m2_mean));
  double U4_err =
      0.5 * U4_mean *
      sqrt(pow(m4_err / m4_mean, 2.0) + pow(2.0 * m2_err / m2_mean, 2.0));
  printf("U4: %.12e %.12e\n", U4_mean, U4_err);

  double m_susc_mean = (m2_mean - m_mean * m_mean) * vol;
  double m_susc_err =
      sqrt(pow(m2_err, 2.0) + pow(2.0 * m_mean * m_err, 2.0)) * vol;
  printf("m_susc: %.12e %.12e\n", m_susc_mean, m_susc_err);

  // print 2-point function legendre coefficients
  printf("opening file: %s\n", legendre_2pt_path.c_str());
  legendre_2pt_file = fopen(legendre_2pt_path.c_str(), "w");
  assert(legendre_2pt_file != nullptr);
  for (int l = 0; l <= l_max; l++) {
    fprintf(legendre_2pt_file, "%02d ", l);
    legendre_2pt[l].WriteMeasurement(legendre_2pt_file);
  }
  fclose(legendre_2pt_file);

  // print 2-point function spherical harmonic coefficients
  printf("opening file: %s\n", ylm_2pt_path.c_str());
  ylm_2pt_file = fopen(ylm_2pt_path.c_str(), "w");
  assert(ylm_2pt_file != nullptr);
  for (int ylm_i = 0, l = 0, m = 0; ylm_i < n_ylm; ylm_i++) {
    fprintf(ylm_2pt_file, "%04d %02d %02d ", ylm_i, l, m);
    ylm_2pt[ylm_i].WriteMeasurement(ylm_2pt_file);
    m++;
    if (m > l) {
      l++;
      m = 0;
    }
  }
  fclose(ylm_2pt_file);

  return 0;
}

// s2.h

#pragma once

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <algorithm>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdio>
#include <string>
#include <unordered_map>
#include <vector>

#include "grp_o3.h"
#include "lattice.h"
#include "util.h"

// symmetry group data directory must be defined
#ifndef GRP_DIR
#error Error: GRP_DIR is not defined (it should be in the Makefile)
#endif

typedef Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic> ComplexMat;
typedef Eigen::Matrix<Complex, Eigen::Dynamic, 1> ComplexVec;

/// @brief Simplicial lattice discretization of a 2-sphere
class QfeLatticeS2 : public QfeLattice {
 public:
  QfeLatticeS2(int q = 5, int k = 1);
  void ReadBaseLattice(int q);
  void WriteSite(FILE* file, int s);
  void ReadSite(FILE* file, int s);
  int CreateOrbit(double xi1, double xi2);
  Vec3 CalcOrbitPos(int o);
  void ReadOrbits(FILE* file);
  void WriteOrbits(FILE* file);
  void UpdateOrbits();
  void ReadSymmetryData(int q, int k);
  void ResizeSites(int n_sites);
  void Inflate();
  void UpdateAntipodes();
  Vec3 FaceCircumcenter(int f);
  double EdgeSquared(int l);
  double EdgeLength(int l);
  double FlatArea(int f);
  void UpdateWeights();
  double CalcLatticeSpacing();
  void OptimizeIntegrator(int l_max);
  void UpdateYlm(int l_max);
  Complex GetYlm(int s, int l, int m);
  Complex CalcYlm(int s, int l, int m);
  double CosTheta(int s1, int s2);
  void PrintCoordinates();

  int q;                        // base polyhedron parameter
  std::vector<Vec3> r;          // vertex coordinates
  std::vector<int> antipode;    // antipode of each site (0 by default)
  std::vector<int> site_orbit;  // orbit id for each site
  std::vector<int> face_orbit;  // orbit id for each face
  std::vector<Vec3> orbit_xi;   // barycentric coordinates for each orbit
  Vec3 first_face_r[4];         // coordinates of first face vertices
  std::vector<GrpElemO3> G;     // symmetry group elements
  std::vector<int> site_g;      // group element for each site

  std::vector<std::vector<Complex>> ylm;  // spherical harmonics
};

/// @brief Create a simplicial discretization of a 2-sphere.
/// @param q Number of links meeting at each site. Valid values for @p q are 3,
/// 4, and 5 for a tetrahedron, octahedron, and icosahedron, respectively.
/// @param k Refinement level.
QfeLatticeS2::QfeLatticeS2(int q, int k) {
  // refinement level must be positive
  assert(k >= 1);

  // tetrahedron, octahedron, and icosahedran are the only valid base
  // polyhedrons
  assert(q >= 3 && q <= 5);

  this->q = q;

  // return base lattice if unrefined
  if (k == 1) {
    ReadBaseLattice(q);
    UpdateDistinct();
    ReadSymmetryData(q, k);
    return;
  }

  if (q == 3) {
    // k-refined tetrahedron has V = 2 k^2 + 2 vertices
    ResizeSites(2 * k * k + 2);
  } else if (q == 4) {
    // k-refined octahedron has V = 4 k^2 + 2 vertices
    ResizeSites(4 * k * k + 2);
  } else if (q == 5) {
    // k-refined icosahedron has V = 10 k^2 + 2 vertices
    ResizeSites(10 * k * k + 2);
  }

  // set the lattice volume
  vol = double(n_sites);

  // create an unrefined base lattice
  QfeLatticeS2 base_lattice(q);

  // maps to identify vertices, orbits, and faces
  std::unordered_map<std::string, int> coord_map;
  std::unordered_map<std::string, int> orbit_map;
  std::unordered_map<std::string, int> face_map;

  // index of next site and face to create
  int s_next = 0;
  int f_next = 0;

  // xy coordinates of each site
  std::vector<int> site_x(n_sites);
  std::vector<int> site_y(n_sites);

  // loop over faces of base polyhedron
  for (int f = 0; f < base_lattice.n_faces; f++) {
    // get the vertices of the base polyhedron face
    Vec3 face_r[3];
    for (int i = 0; i < 3; i++) {
      int b = base_lattice.faces[f].sites[i];
      face_r[i] = base_lattice.r[b];
      if (f == 0) {
        // set the coordinates of the first face's vertices
        first_face_r[i] = face_r[i];
      }
    }

    // unit vectors in the face in xy basis
    Vec3 n_x = (face_r[1] - face_r[0]) / double(k);
    Vec3 n_y = (face_r[2] - face_r[0]) / double(k);

    // list of sites in this face labeled by xy positions on the face
    int k1 = k + 1;
    int xy_max = k1 * k1;
    int xy_list[k1][k1];

    // loop over xy to find all sites and set their positions
    for (int xy = 0; xy <= xy_max; xy++) {
      int x = xy % k1;
      int y = xy / k1;

      // skip xy values outside the triangle
      if ((x + y) > k) continue;

      // calculate the coordinates of this vertex
      Vec3 v = (face_r[0] + x * n_x + y * n_y).normalized();

      // deal with negative zero
      std::string vec_name = Vec3ToString(v);

      // check if the site already exists
      if (coord_map.find(vec_name) == coord_map.end()) {
        // sorted barycentric coordinates define the site orbit
        int xi[3];
        xi[0] = x;
        xi[1] = y;
        xi[2] = k - x - y;
        std::sort(xi, xi + 3, std::greater<int>());
        std::string orbit_name = string_format("%d_%d_%d", xi[0], xi[1], xi[2]);

        // check if the orbit already exists
        if (orbit_map.find(orbit_name) == orbit_map.end()) {
          // create a new orbit
          double xi1 = double(xi[0]) / double(k);
          double xi2 = double(xi[1]) / double(k);
          orbit_map[orbit_name] = CreateOrbit(xi1, xi2);
        }

        // create a new site
        int orbit_id = orbit_map[orbit_name];
        coord_map[vec_name] = s_next;
        r[s_next] = v;
        sites[s_next].nn = 0;
        sites[s_next].wt = 1.0;
        sites[s_next].id = orbit_id;
        site_orbit[s_next] = orbit_id;
        s_next++;
        assert(s_next <= n_sites);
      }

      // get the site index
      int s = coord_map[vec_name];

      // save this site in the xy list
      xy_list[x][y] = s;
      site_x[s] = x;
      site_y[s] = y;
    }

    // add faces and links
    for (int xy = 0; xy <= xy_max; xy++) {
      int x = xy % k1;
      int y = xy / k1;
      if ((x + y) > k) continue;

      if ((x + y) != k) {
        // "forward" triangle
        int s1 = xy_list[x][y];
        int s2 = xy_list[x][y + 1];
        int s3 = xy_list[x + 1][y];
        AddFace(s1, s2, s3);
      }

      if ((x != 0) && (y != 0)) {
        // "backward" triangle
        int s1 = xy_list[x][y];
        int s2 = xy_list[x][y - 1];
        int s3 = xy_list[x - 1][y];
        AddFace(s1, s2, s3);
      }
    }

    // set the orbit id for each face
    int n_distinct_faces = 0;
    face_orbit.resize(n_faces);
    while (f_next != n_faces) {
      int xi[3] = {0, 0, 0};

      // compute the barycentric coordinates of this face
      for (int i = 0; i < 3; i++) {
        int s = faces[f_next].sites[i];
        int x = site_x[s];
        int y = site_y[s];

        xi[0] += x;
        xi[1] += y;
        xi[2] += k - x - y;
      }

      // sorted barycentric coordinates define the face orbit
      std::sort(xi, xi + 3, std::greater<int>());
      std::string face_name = string_format("%d_%d_%d", xi[0], xi[1], xi[2]);

      // check if the face already exists
      if (face_map.find(face_name) == face_map.end()) {
        face_map[face_name] = n_distinct_faces++;
      }
      face_orbit[f_next] = face_map[face_name];
      f_next++;
    }
  }

  // check that all of the sites and faces have been created
  assert(s_next == n_sites);
  assert(f_next == n_faces);

  // read the symmetry group data
  UpdateDistinct();
  ReadSymmetryData(q, k);
}

/// @brief Read base polyhedron data from grp directory
/// @param q polyhedron parameter
void QfeLatticeS2::ReadBaseLattice(int q) {
  // tetrahedron, octahedron, and icosahedron are the only valid base
  // polyhedrons
  assert(q >= 3 && q <= 5);
  this->q = q;

  // read the base lattice file in the symmetry group directory
  std::string lattice_path = string_format("%s/lattice/o3q%d.dat", GRP_DIR, q);
  FILE* file = fopen(lattice_path.c_str(), "r");
  assert(file != nullptr);
  ReadLattice(file);
  fclose(file);

  // set the lattice volume
  vol = double(n_sites);

  // create a single orbit
  CreateOrbit(0.0, 0.0);

  // initialize the first face vertex coordinates
  for (int i = 0; i < 3; i++) {
    int s = faces[0].sites[i];
    first_face_r[i] = r[s];
  }

  face_orbit.resize(n_faces);
  for (int f = 0; f < n_faces; f++) face_orbit[f] = 0;
}

/// @brief Write a site to a lattice file
/// @param file Lattice file
/// @param s Site index
void QfeLatticeS2::WriteSite(FILE* file, int s) {
  QfeLattice::WriteSite(file, s);
  double theta = acos(r[s].z());
  double phi = atan2(r[s].y(), r[s].x());
  fprintf(file, " %+.20f %+.20f", theta, phi);
}

/// @brief Read a site from a lattice file
/// @param file Lattice file
/// @param s Site index
void QfeLatticeS2::ReadSite(FILE* file, int s) {
  QfeLattice::ReadSite(file, s);
  double theta, phi;
  fscanf(file, " %lf %lf", &theta, &phi);

  r[s][0] = sin(theta) * cos(phi);
  r[s][1] = sin(theta) * sin(phi);
  r[s][2] = cos(theta);
  r[s].normalize();

  Vec3 north_pole(0.0, 0.0, 1.0);
  Vec3 south_pole(0.0, 0.0, -1.0);
  if (AlmostEq(r[s], north_pole)) r[s] = north_pole;
  if (AlmostEq(r[s], south_pole)) r[s] = south_pole;
}

/// @brief Create an orbit
/// @param xi1 1st barycentric coordinate
/// @param xi2 2nd barycentric coordinate
int QfeLatticeS2::CreateOrbit(double xi1, double xi2) {
  // barycentric coordinates, sorted to account for degeneracies
  double xi[3] = {xi1, xi2, 1.0 - xi1 - xi2};
  std::sort(xi, xi + 3, std::greater<double>());
  int o = orbit_xi.size();
  orbit_xi.push_back(Vec3(xi));
  return o;
}

/// @brief Calculate the coordinates of the first site in an orbit
/// @param o Orbit index
/// @return Normalized orbit coordinates
Vec3 QfeLatticeS2::CalcOrbitPos(int o) {
  Vec3 r = Vec3::Zero();
  Vec3 xi = orbit_xi[o];
  for (int i = 0; i < 3; i++) {
    r += xi(i) * first_face_r[i];
  }
  return r.normalized();
}

/// @brief Read an orbit file and convert to site coordinates
/// @param file Orbit file
void QfeLatticeS2::ReadOrbits(FILE* file) {
  // make sure we're at the beginning of the file
  fseek(file, 0L, SEEK_SET);

  // read orbit data
  orbit_xi.resize(n_distinct);
  for (int o = 0; o < n_distinct; o++) {
    // read barycentric coordinates
    int o_check;
    fscanf(file, "%d", &o_check);
    assert(o_check == o);

    // read dof values
    double xi_sum = 0.0;
    for (int i = 0; i < 2; i++) {
      double temp;
      fscanf(file, "%lf", &temp);
      xi_sum += temp;
      orbit_xi[o][i] = temp;
    }
    orbit_xi[o][2] = 1.0 - xi_sum;
    fscanf(file, "\n");
  }
  assert(feof(file));

  UpdateOrbits();
}

/// @brief Write orbit barycentric coordinates to a file that can be read
/// via ReadOrbits
/// @param file Orbit file
void QfeLatticeS2::WriteOrbits(FILE* file) {
  // read orbit data
  for (int o = 0; o < n_distinct; o++) {
    // read barycentric coordinates
    fprintf(file, "%d", o);

    // read dof values
    for (int i = 0; i < 2; i++) {
      fprintf(file, " %.16f", orbit_xi[o][i]);
    }
    fprintf(file, "\n");
  }
}

/// @brief Update positions of all sites using orbits and symmetry group data
void QfeLatticeS2::UpdateOrbits() {
  // calculate the orbit positions
  std::vector<Vec3> orbit_r(n_distinct);
  for (int o = 0; o < n_distinct; o++) {
    orbit_r[o] = CalcOrbitPos(o);
  }

  // use the symmetry group data to calculate site coordinates
  Vec3 north_pole(0.0, 0.0, 1.0);
  Vec3 south_pole(0.0, 0.0, -1.0);
  for (int s = 0; s < n_sites; s++) {
    int o = site_orbit[s];
    int g = site_g[s];
    r[s] = G[g] * orbit_r[o];
    r[s].normalize();
    if (AlmostEq(r[s], north_pole)) r[s] = north_pole;
    if (AlmostEq(r[s], south_pole)) r[s] = south_pole;
  }
}

/// @brief Read symmetry group data from pre-generated files in the grp
/// directory
/// @param q polyhedron parameter
/// @param k refinement level
void QfeLatticeS2::ReadSymmetryData(int q, int k) {
  // open the symmetry group data file
  std::string grp_path = string_format("%s/elem/o3q%d.dat", GRP_DIR, q);
  FILE* grp_file = fopen(grp_path.c_str(), "r");
  assert(grp_file != nullptr);

  // read group elements
  G.clear();
  while (!feof(grp_file)) {
    GrpElemO3 g;
    g.ReadGrpElem(grp_file);
    G.push_back(g);
  }
  fclose(grp_file);

  // calculate all of the orbit positions
  std::vector<Vec3> orbit_r(n_distinct);
  for (int o = 0; o < n_distinct; o++) {
    orbit_r[o] = CalcOrbitPos(o);
  }

  // open the site group element file
  std::string g_path = string_format("%s/site_g/o3q%dk%d.dat", GRP_DIR, q, k);
  FILE* g_file = fopen(g_path.c_str(), "r");
  bool site_g_success = true;
  if (g_file != nullptr) {
    // load pre-existing symmetry data
    site_g.resize(n_sites);
    std::vector<int>::iterator it = site_g.begin();
    while (!feof(g_file)) {
      assert(it != site_g.end());
      int g;
      fscanf(g_file, "%d\n", &g);
      *it++ = g;
    }
    fclose(g_file);

    // recalculate if the file was not long enough
    if (it != site_g.end()) {
      fprintf(stderr, "Rebuilding invalid data file: %s\n", g_path.c_str());
      site_g_success = false;
    } else {
      for (int s = 0; s < n_sites; s++) {
        int o = site_orbit[s];
        int g = site_g[s];
        Vec3 r_norm = r[s].normalized();
        Vec3 gr = G[g] * orbit_r[o];
        if (!AlmostEq(r_norm, gr, 1.0e-15)) {
          site_g_success = false;
          break;
        }
      }

      if (!site_g_success) {
        fprintf(stderr, "Rebuilding invalid data file: %s\n", g_path.c_str());
      }
    }

  } else {
    site_g_success = false;
  }

  if (!site_g_success) {
    // find the group element for each site
    site_g.resize(n_sites);
    for (int s = 0; s < n_sites; s++) {
      int o = site_orbit[s];
      Vec3 r_norm = r[s].normalized();

      // find the appropriate group element
      bool found_g = false;
      for (int g = 0; g < G.size(); g++) {
        Vec3 gr = G[g] * orbit_r[o];
        if (!AlmostEq(r_norm, gr, 1.0e-15)) continue;
        site_g[s] = g;
        found_g = true;
        break;
      }
      assert(found_g);
    }

    // write the site group elements to a file
    g_file = fopen(g_path.c_str(), "w");
    for (int s = 0; s < n_sites; s++) {
      fprintf(g_file, "%d\n", site_g[s]);
    }
    fclose(g_file);
  }

  UpdateOrbits();
}

/// @brief Change the number of sites.
/// @param n_sites New number of sites
void QfeLatticeS2::ResizeSites(int n_sites) {
  QfeLattice::ResizeSites(n_sites);
  r.resize(n_sites);
  ylm.resize(n_sites);
  antipode.resize(n_sites, 0);
  site_orbit.resize(n_sites);
}

/// @brief Project all site coordinates onto a unit sphere.
void QfeLatticeS2::Inflate() {
  for (int s = 0; s < n_sites; s++) {
    r[s].normalize();
  }
}

/// @brief Identify each site's antipode, i.e. for a site with position r, find
/// the site which has position -r. A lattice with a tetrahedron base (q = 3)
/// does not have an antipode for every site.
void QfeLatticeS2::UpdateAntipodes() {
  std::unordered_map<std::string, int> antipode_map;
  for (int s = 0; s < n_sites; s++) {
    // find antipode
    Vec3 anti_r = -r[s];
    std::string key = Vec3ToString(r[s], 6);
    std::string anti_key = Vec3ToString(anti_r, 6);

    if (antipode_map.find(anti_key) != antipode_map.end()) {
      // antipode found in map
      int a = antipode_map[anti_key];
      antipode[s] = a;
      antipode[a] = s;
      antipode_map.erase(anti_key);
    } else {
      // antipode not found yet
      antipode_map[key] = s;
    }
  }

  if (antipode_map.size()) {
    // print error message if there are any unpaired sites
    fprintf(stderr, "no antipode found for %lu/%d sites\n", antipode_map.size(),
            n_sites);
    std::unordered_map<std::string, int>::iterator it;
    for (it = antipode_map.begin(); it != antipode_map.end(); it++) {
      fprintf(stderr, "%04d %s\n", it->second, it->first.c_str());
    }
  }
}

/// @brief Find the circumcenter of face
/// @param f Face id
/// @return Coordinates of face circumcenter
Vec3 QfeLatticeS2::FaceCircumcenter(int f) {
  double sq_edge_1 = EdgeSquared(faces[f].edges[0]);
  double sq_edge_2 = EdgeSquared(faces[f].edges[1]);
  double sq_edge_3 = EdgeSquared(faces[f].edges[2]);

  double w1 = sq_edge_1 * (sq_edge_2 + sq_edge_3 - sq_edge_1);
  double w2 = sq_edge_2 * (sq_edge_3 + sq_edge_1 - sq_edge_2);
  double w3 = sq_edge_3 * (sq_edge_1 + sq_edge_2 - sq_edge_3);

  Vec3 r1 = w1 * r[faces[f].sites[0]];
  Vec3 r2 = w2 * r[faces[f].sites[1]];
  Vec3 r3 = w3 * r[faces[f].sites[2]];

  return (r1 + r2 + r3) / (w1 + w2 + w3);
}

/// @brief Calculate the squared length of link
/// @param l Link id
/// @return Squared length of link
double QfeLatticeS2::EdgeSquared(int l) {
  int s_a = links[l].sites[0];
  int s_b = links[l].sites[1];
  Vec3 dr = r[s_a] - r[s_b];
  return dr.squaredNorm();
}

/// @brief Calculate the length of a link
/// @param l Link id
/// @return Link length
double QfeLatticeS2::EdgeLength(int l) { return sqrt(EdgeSquared(l)); }

/// @brief Calculate the flat area of a triangular face.
/// @param f Face id
/// @return Area of triangular face
double QfeLatticeS2::FlatArea(int f) {
  double a = EdgeLength(faces[f].edges[0]);
  double b = EdgeLength(faces[f].edges[1]);
  double c = EdgeLength(faces[f].edges[2]);
  double area = (a + b + c) * (b + c - a) * (c + a - b) * (a + b - c);
  return 0.25 * sqrt(area);
}

/// @brief Calculate FEM weights based on vertex coordinates.
void QfeLatticeS2::UpdateWeights() {
  // set site weights to zero
  for (int s = 0; s < n_sites; s++) {
    sites[s].wt = 0.0;
  }

  // loop over links to update weights
  for (int l = 0; l < n_links; l++) {
    links[l].wt = 0.0;
    for (int i = 0; i < 2; i++) {
      // find the other two edges of this face
      int f = links[l].faces[i];
      int e = 0;
      while (faces[f].edges[e] != l) e++;
      int e1 = (e + 1) % 3;
      int e2 = (e + 2) % 3;
      int l1 = faces[f].edges[e1];
      int l2 = faces[f].edges[e2];

      // find the area associated with this face
      double sq_edge = EdgeSquared(l);
      double sq_edge_1 = EdgeSquared(l1);
      double sq_edge_2 = EdgeSquared(l2);
      double half_wt = (sq_edge_1 + sq_edge_2 - sq_edge) / (8.0 * FlatArea(f));
      links[l].wt += half_wt;

      // add to the weights of the two sites connected by this link
      sites[links[l].sites[0]].wt += 0.25 * half_wt * sq_edge;
      sites[links[l].sites[1]].wt += 0.25 * half_wt * sq_edge;
    }
  }

  // normalize site weights to 1
  double site_wt_sum = 0.0;
  for (int s = 0; s < n_sites; s++) {
    site_wt_sum += sites[s].wt;
  }

  double site_wt_norm = site_wt_sum / double(n_sites);
  for (int s = 0; s < n_sites; s++) {
    sites[s].wt /= site_wt_norm;
  }

  // set face weights equal to their flat area
  double face_area_sum = 0.0;
  for (int f = 0; f < n_faces; f++) {
    double face_area = FlatArea(f);
    face_area_sum += face_area;
    faces[f].wt = face_area;
  }

  // normalize face areas
  double face_wt_norm = face_area_sum / double(n_faces);
  for (int f = 0; f < n_faces; f++) {
    faces[f].wt /= face_wt_norm;
  }
}

/// @brief Calculate the global effective lattice spacing
/// @return Lattice spacing a/r
double QfeLatticeS2::CalcLatticeSpacing() {
  double area_sum = 0.0;
  for (int f = 0; f < n_faces; f++) {
    area_sum += FlatArea(f);
  }
  double area_mean = area_sum / double(n_sites);
  return sqrt(area_mean);
}

/// @brief Optimize the site weights so that all linear combinations of
/// spherical harmonics invariant under the relevant symmetry group can be
/// integrated exactly up to order @p l_max.
/// @param l_max Maximum spherical harmonic eigenvalue to optimize
void QfeLatticeS2::OptimizeIntegrator(int l_max) {
  // determine which l,m combinations need to be integrated exactly
  std::vector<int> l_relevant;
  std::vector<int> m_relevant;
  int m_spacing = (q == 5) ? 5 : 4;

  for (int l = 0; l <= l_max; l++) {
    // odd l only contributes for tetrahedron
    if ((l % 2) && (q != 3)) continue;

    // number of functions at this l (overcounts for tetrahedron)
    int n_l = (l / 2) + (l / 3) + (l / q) - l + 1;

    for (int i = 0; i < n_l; i++) {
      l_relevant.push_back(l);
      int m = i * m_spacing;

      // odd ell (tetrahedron only)
      if (l % 2) m = ((2 * i + 1) * m_spacing) / 2;

      m_relevant.push_back(m);
    }
  }

  // number of relevant functions
  int n_relevant = l_relevant.size();

  // generate the rectangular matrix
  ComplexMat S = ComplexMat::Zero(n_relevant, n_distinct);
  for (int s = 0; s < n_sites; s++) {
    int id = sites[s].id;
    double theta = acos(r[s].z());
    double phi = atan2(r[s].y(), r[s].x());
    for (int i = 0; i < n_relevant; i++) {
      int l = l_relevant[i];
      int m = m_relevant[i];
      S(i, id) += boost::math::spherical_harmonic(l, m, theta, phi);
    }
  }

  // use the current weights as an initial guess
  ComplexVec x0(n_distinct);
  for (int id = 0; id < n_distinct; id++) {
    int s = distinct_first[id];
    x0(id) = sites[s].wt;
  }

  // the right hand side is the spherical harmonic orthogonality condition
  ComplexVec b = ComplexVec::Zero(n_relevant);
  b(0) = 0.28209479177387814347 * vol;  // n_sites / sqrt(4 pi)

  // compute the solution
  Eigen::LeastSquaresConjugateGradient<ComplexMat> cg;
  cg.compute(S);
  assert(cg.info() == Eigen::Success);
  ComplexVec wt = cg.solveWithGuess(b, x0);

  // apply the improved weights to the sites
  for (int s = 0; s < n_sites; s++) {
    int id = sites[s].id;
    if (s == distinct_first[id]) {
      assert(real(wt(id)) > 0.0);
      // printf("%04d %.12f %.12f\n", id, sites[s].wt, real(wt(id)));
    }
    sites[s].wt = real(wt(id));
  }
}

/// @brief Update spherical harmonic values at each site, up to a maximum l
/// eigenvalue of @p l_max
/// @param l_max Maximum spherical harmonic eigenvalue to calculate
void QfeLatticeS2::UpdateYlm(int l_max) {
  int n_ylm = ((l_max + 1) * (l_max + 2)) / 2;
  using boost::math::spherical_harmonic;

  for (int s = 0; s < n_sites; s++) {
    ylm[s].resize(n_ylm);
    double theta = acos(r[s].z());
    double phi = atan2(r[s].y(), r[s].x());

    for (int i = 0, l = 0, m = 0; i < n_ylm; i++, m++) {
      if (m > l) {
        m = 0;
        l++;
      }
      assert(i < n_ylm);
      ylm[s][i] = spherical_harmonic(l, m, theta, phi);
    }
  }
}

/// @brief Retrieve a pre-calculated spherical harmonic value.
/// @param s Site index
/// @param l Spherical harmonic eigenvalue
/// @param m Spherical harmonic eigenvalue
/// @return Spherical harmonic evaluated at site @p s
Complex QfeLatticeS2::GetYlm(int s, int l, int m) {
  int abs_m = fabs(m);
  assert(abs_m <= l);

  int i = (l * (l + 1)) / 2 + abs_m;
  assert(i < ylm[s].size());
  Complex y = ylm[s][i];

  if (m < 0) {
    y = conj(y);
    if (abs_m & 1) {
      y *= -1;
    }
  }

  return y;
}

/// @brief Calculate a spherical harmonic value (not pre-calculated)
/// @param s Site index
/// @param l Spherical harmonic eigenvalue
/// @param m Spherical harmonic eigenvalue
/// @return Spherical harmonic evaluated at site @p s
Complex QfeLatticeS2::CalcYlm(int s, int l, int m) {
  double theta = acos(r[s].z());
  double phi = atan2(r[s].y(), r[s].x());

  return boost::math::spherical_harmonic(l, m, theta, phi);
}

/// @brief Calculate the cosine of the angle between two sites. This function
/// assumes that the coordinates have been projected onto the unit sphere.
/// @param s1 1st site index
/// @param s2 2nd site index
/// @return cosine of the angle between @p s1 and @p s2
double QfeLatticeS2::CosTheta(int s1, int s2) {
  if (s1 == s2) return 1.0;
  if (antipode[s1] == s2) return -1.0;
  return r[s1].dot(r[s2]);
}

/// @brief Print the cartesian coordinates of the sites. This is helpful for
/// making plots in e.g. Mathematica.
void QfeLatticeS2::PrintCoordinates() {
  printf("{");
  for (int s = 0; s < n_sites; s++) {
    printf("{%.12f,%.12f,%.12f}", r[s].x(), r[s].y(), r[s].z());
    printf("%c\n", s == (n_sites - 1) ? '}' : ',');
  }
}

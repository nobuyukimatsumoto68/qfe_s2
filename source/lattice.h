// lattice.h

#pragma once

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <random>
#include <stack>
#include <string>
#include <unordered_map>
#include <vector>

#include "rng.h"

#define MAX_SITE_NEIGHBORS 20
#define MAX_LINK_FACES 12
#define MAX_FACE_EDGES 4
#define MAX_FACE_CELLS 2
#define MAX_CELL_FACES 4
#define MAX_CELL_SITES 4

struct QfeSite {
  double wt;                          // site weight
  int nn;                             // number of nearest neighbors
  int links[MAX_SITE_NEIGHBORS];      // nearest neighbor links
  int neighbors[MAX_SITE_NEIGHBORS];  // nearest neighbor sites
  int id;
};

struct QfeLink {
  double wt;                  // link weight
  int sites[2];               // sites attached by this link
  int n_faces;                // number of faces that include this link
  int faces[MAX_LINK_FACES];  // faces that include this link
};

struct QfeFace {
  double wt;                       // face weight
  int n_edges;                     // number of links (3 for triangle)
  int n_cells;                     // number of cells that include this face
  int cells[MAX_FACE_CELLS];       // cells around this face
  int edges[MAX_FACE_EDGES];       // links around this face
  int sites[MAX_FACE_EDGES];       // sites around this face
  bool flip_edge[MAX_FACE_EDGES];  // edges which are flipped wrt this face
};

struct QfeCell {
  double wt;
  int n_faces;
  int faces[MAX_CELL_FACES];
  int sites[MAX_CELL_SITES];
};

class QfeLattice {
 public:
  QfeLattice();
  virtual void WriteLattice(FILE* file);
  virtual void WriteSite(FILE* file, int s);
  virtual void WriteLink(FILE* file, int l);
  virtual void WriteFace(FILE* file, int f);
  virtual void WriteCell(FILE* file, int c);
  virtual void ReadLattice(FILE* file);
  virtual void ReadSite(FILE* file, int s);
  virtual void ReadLink(FILE* file, int l);
  virtual void ReadFace(FILE* file, int f);
  virtual void ReadCell(FILE* file, int c);
  void SeedRng(unsigned int seed);
  void InitRect(int Nx, int Ny, double wt_x, double wt_y);
  void InitTriangle(int N, double wt1, double wt2, double wt3);
  void InitTriangle(int Nx, int Ny, double wt1, double wt2, double wt3);
  virtual void ResizeSites(int n_sites);
  virtual void InterpolateSite(int s, int s_a, int s_b, int num, int den);
  virtual void AddDimension(int n_slices);
  int FindLink(int s1, int s2);
  int FindFace(int s1, int s2, int s3);
  int FindCell(int s1, int s2, int s3, int s4);
  int AddLink(int a, int b, double wt);
  int AddFace(int a, int b, int c, double wt = 1.0);
  int AddFace(int a, int b, int c, int d, double wt = 1.0);
  int AddCell(int a, int b, int c, int d, double wt = 1.0);
  void UpdateDistinct();
  void Refine2D(int n_refine);
  void PrintSites();
  void PrintLinks();
  void PrintFaces();
  void PrintCells();
  void CheckConnectivity();
  void CheckConsistency();

  int n_sites;
  int n_links;
  int n_faces;
  int n_cells;
  double vol;

  std::vector<QfeSite> sites;
  std::vector<QfeLink> links;
  std::vector<QfeFace> faces;
  std::vector<QfeCell> cells;

  // symmetrically distinct sites
  std::vector<int> distinct_n_sites;  // number of sites for each distinct id
  std::vector<int> distinct_first;    // representative site for distinct group
  int n_distinct;

  QfeRng rng;
};

QfeLattice::QfeLattice() {
  n_sites = 0;
  n_links = 0;
  n_faces = 0;
  n_cells = 0;
  vol = 0.0;
  n_distinct = 0;
}

void QfeLattice::WriteLattice(FILE* file) {
  fprintf(file, "begin_sites\n");
  fprintf(file, "n_sites %d\n", n_sites);
  for (int s = 0; s < n_sites; s++) {
    WriteSite(file, s);
    fprintf(file, "\n");
  }
  fprintf(file, "end_sites\n");

  fprintf(file, "begin_links\n");
  fprintf(file, "n_links %d\n", n_links);
  for (int l = 0; l < n_links; l++) {
    WriteLink(file, l);
    fprintf(file, "\n");
  }
  fprintf(file, "end_links\n");

  fprintf(file, "begin_faces\n");
  fprintf(file, "n_faces %d\n", n_faces);
  for (int f = 0; f < n_faces; f++) {
    WriteFace(file, f);
    fprintf(file, "\n");
  }
  fprintf(file, "end_faces\n");

  fprintf(file, "begin_cells\n");
  fprintf(file, "n_cells %d\n", n_cells);
  for (int c = 0; c < n_cells; c++) {
    WriteCell(file, c);
    fprintf(file, "\n");
  }
  fprintf(file, "end_cells\n");
}

void QfeLattice::ReadLattice(FILE* file) {
  sites.clear();
  n_sites = 0;
  links.clear();
  n_links = 0;
  faces.clear();
  n_faces = 0;
  vol = 0.0;

  char buffer[200];
  int n;

  while (!feof(file)) {
    fscanf(file, "%s\n", buffer);

    if (strcmp(buffer, "begin_sites") == 0) {
      // read sites
      fscanf(file, "%s %d\n", buffer, &n);
      if (strcmp(buffer, "n_sites")) {
        printf("invalid value \"%s\", expected n_sites\n", buffer);
        return;
      }

      ResizeSites(n);
      vol = 0.0;
      for (int s = 0; s < n_sites; s++) {
        ReadSite(file, s);
        vol += sites[s].wt;
        getc(file);  // read newline character
      }

      fscanf(file, "%s\n", buffer);
      if (strcmp(buffer, "end_sites")) {
        printf("invalid value \"%s\", expected end_sites\n", buffer);
        return;
      }
    } else if (strcmp(buffer, "begin_links") == 0) {
      // read links
      fscanf(file, "%s %d\n", buffer, &n);
      if (strcmp(buffer, "n_links")) {
        printf("invalid value \"%s\", expected n_links\n", buffer);
        return;
      }

      for (int l = 0; l < n; l++) {
        ReadLink(file, l);
        getc(file);  // read newline character
      }

      fscanf(file, "%s\n", buffer);
      if (strcmp(buffer, "end_links")) {
        printf("invalid value \"%s\", expected end_links\n", buffer);
        return;
      }
    } else if (strcmp(buffer, "begin_faces") == 0) {
      // read faces
      fscanf(file, "%s %d\n", buffer, &n);
      if (strcmp(buffer, "n_faces")) {
        printf("invalid value \"%s\", expected n_faces\n", buffer);
        return;
      }

      for (int f = 0; f < n; f++) {
        ReadFace(file, f);
        getc(file);  // read newline character
      }

      fscanf(file, "%s\n", buffer);
      if (strcmp(buffer, "end_faces")) {
        printf("invalid value \"%s\", expected end_faces\n", buffer);
        return;
      }
    } else if (strcmp(buffer, "begin_cells") == 0) {
      // read cells
      fscanf(file, "%s %d\n", buffer, &n);
      if (strcmp(buffer, "n_cells")) {
        printf("invalid value \"%s\", expected n_cells\n", buffer);
        return;
      }

      for (int c = 0; c < n; c++) {
        ReadCell(file, c);
        getc(file);  // read newline character
      }

      fscanf(file, "%s\n", buffer);
      if (strcmp(buffer, "end_cells")) {
        printf("invalid value \"%s\", expected end_cells\n", buffer);
        return;
      }
    }
  }

  UpdateDistinct();
}

void QfeLattice::WriteSite(FILE* file, int s) {
  fprintf(file, "%04d %.16e %04d", s, sites[s].wt, sites[s].id);
}

void QfeLattice::ReadSite(FILE* file, int s) {
  int s_chk;
  fscanf(file, "%d %lf %d", &s_chk, &sites[s].wt, &sites[s].id);
  if (s != s_chk) {
    printf("non-matching site index: %04d %04d\n", s, s_chk);
  }
}

void QfeLattice::WriteLink(FILE* file, int l) {
  fprintf(file, "%04d %.16e %04d %04d", l, links[l].wt, links[l].sites[0],
          links[l].sites[1]);
}

void QfeLattice::ReadLink(FILE* file, int l) {
  int l_chk;
  double wt;
  int s_a, s_b;
  fscanf(file, "%d %lf %d %d", &l_chk, &wt, &s_a, &s_b);
  if (l != l_chk) {
    printf("non-matching link index: %04d %04d\n", l, l_chk);
  }
  AddLink(s_a, s_b, wt);
}

void QfeLattice::WriteFace(FILE* file, int f) {
  fprintf(file, "%04d %.16e %04d %04d %04d", f, faces[f].wt, faces[f].sites[0],
          faces[f].sites[1], faces[f].sites[2]);
}

void QfeLattice::ReadFace(FILE* file, int f) {
  int f_chk;
  double wt;
  int s_a, s_b, s_c;
  fscanf(file, "%d %lf %d %d %d", &f_chk, &wt, &s_a, &s_b, &s_c);
  if (f != f_chk) {
    printf("non-matching face index: %04d %04d\n", f, f_chk);
  }
  AddFace(s_a, s_b, s_c);
}

void QfeLattice::WriteCell(FILE* file, int c) {
  fprintf(file, "%04d %.16e %04d %04d %04d %04d", c, cells[c].wt,
          cells[c].sites[0], cells[c].sites[1], cells[c].sites[2],
          cells[c].sites[3]);
}

void QfeLattice::ReadCell(FILE* file, int c) {
  int c_chk;
  double wt;
  int s_a, s_b, s_c, s_d;
  fscanf(file, "%d %lf %d %d %d %d", &c_chk, &wt, &s_a, &s_b, &s_c, &s_d);
  if (c != c_chk) {
    printf("non-matching cell index: %04d %04d\n", c, c_chk);
  }
  AddCell(s_a, s_b, s_c, s_d);
}

/**
 * @brief Reset and seed the random number generator
 *
 * @param seed Random number generator seed value
 */

void QfeLattice::SeedRng(unsigned int seed) { rng = QfeRng(seed); }

/**
 * @brief Create a flat, rectangular lattice with periodic boundary conditions
 *
 * @param Nx Lattice Size
 * @param Ny Lattice Size
 * @param wtx Link weights
 */

void QfeLattice::InitRect(int Nx, int Ny, double wt1, double wt2) {
  // create sites
  ResizeSites(Nx * Ny);
  vol = double(Nx * Ny);

  // set all site weights to 1.0
  for (int s = 0; s < n_sites; s++) {
    sites[s].wt = 1.0;
    sites[s].nn = 0;
    sites[s].id = 0;
  }

  // create links
  links.clear();

  for (int s = 0; s < n_sites; s++) {
    int x = s % Nx;
    int y = s / Nx;

    // add links in the "forward" direction (2 links per site)
    // each link will end up with 4 neighbors
    int xp1 = (x + 1) % Nx;
    int yp1 = (y + 1) % Ny;
    int f = AddFace(s, xp1 + y * Nx, xp1 + yp1 * Nx, x + yp1 * Nx, 1.0);
    links[faces[f].edges[0]].wt = wt1;
    links[faces[f].edges[1]].wt = wt2;
    links[faces[f].edges[2]].wt = wt1;
    links[faces[f].edges[3]].wt = wt2;
  }
}

/**
 * @brief Creates a flat triangulated lattice with periodic boundary
 * conditions.
 *
 * @param N Lattice size
 * @param wtx The weights given to links in the 3 directions of
 * the triangular lattice
 */

void QfeLattice::InitTriangle(int N, double wt1, double wt2, double wt3) {
  InitTriangle(N, N, wt1, wt2, wt3);
}

/**
 * @brief Creates a flat triangulated lattice with periodic boundary
 * conditions. Lattice dimensions can be asymmetrical.
 *
 * @param Nx Lattice size (x or short direction)
 * @param Ny Lattice size (y or long direction)
 * @param wtx The weights given to links in the 3 directions of
 * the triangular lattice
 */

void QfeLattice::InitTriangle(int Nx, int Ny, double wt1, double wt2,
                              double wt3) {
  // create sites
  ResizeSites(Nx * Ny);
  vol = double(Nx * Ny);

  // set all site weights to 1.0
  for (int s = 0; s < n_sites; s++) {
    sites[s].wt = 1.0;
    sites[s].nn = 0;
    sites[s].id = 0;
  }

  // create links
  links.clear();

  for (int s = 0; s < n_sites; s++) {
    int x = s % Nx;
    int y = s / Nx;

    // add links in the "forward" direction (3 links per site)
    // each link will end up with 6 neighbors
    int xp1 = (x + 1) % Nx;
    int yp1 = (y + 1) % Ny;

    // add "right-handed" faces
    AddFace(s, xp1 + y * Nx, x + yp1 * Nx);
    AddFace(xp1 + y * Nx, xp1 + yp1 * Nx, x + yp1 * Nx);
    int l1 = FindLink(s, xp1 + y * Nx);
    int l2 = FindLink(s, x + yp1 * Nx);
    int l3 = FindLink(xp1 + y * Nx, x + yp1 * Nx);
    links[l1].wt = wt1;
    links[l2].wt = wt2;
    links[l3].wt = wt3;
  }
}

/**
 * @brief Change the number of sites.
 */

void QfeLattice::ResizeSites(int n_sites) {
  this->n_sites = n_sites;
  sites.resize(n_sites);
}

/**
 * @brief Set the position of site @p s by interpolating between sites @p s_a
 * and @p s_b. The parameters @p num and @den define a fraction between 0 and 1
 * that determines how far from site @p s_a to put site @p s, with values of
 * 0 and 1 giving the coordinates of site a and site b, respectively.
 *
 * Subclasses can override this function to set the coordinates of the new
 * point.
 */

void QfeLattice::InterpolateSite(int s, int s_a, int s_b, int num, int den) {
  return;
}

/**
 * @brief Add an extra dimension with @p n_slices slices perpendicular to
 * current lattice.
 */

void QfeLattice::AddDimension(int n_slices) {
  int n_sites_slice = n_sites;  // number of sites per slice

  // create sites for the new slices in the extra dimension
  ResizeSites(n_sites_slice * n_slices);
  vol = 0.0;

  // copy sites from the first slice
  for (int s0 = 0; s0 < n_sites_slice; s0++) {
    QfeSite* site0 = &sites[s0];
    vol += site0->wt * n_slices;

    for (int t = 1; t < n_slices; t++) {
      int s = t * n_sites_slice + s0;
      sites[s].nn = 0;  // add links later
      sites[s].wt = site0->wt;
      sites[s].id = site0->id;
    }
  }

  // duplicate links on the other slices
  int n_links_slice = links.size();
  for (int t = 1; t < n_slices; t++) {
    for (int l = 0; l < n_links_slice; l++) {
      // add links in the other slices
      int s_a = t * n_sites_slice + links[l].sites[0];
      int s_b = t * n_sites_slice + links[l].sites[1];
      AddLink(s_a, s_b, links[l].wt);
    }
  }

  // add links to connect the slices with periodic boundary conditions
  for (int s = 0; s < n_sites; s++) {
    AddLink(s, (n_sites_slice + s) % n_sites, sites[s].wt);
  }
}

/**
 * @brief Finds the link connecting sites @p s1 and @p s2 and returns the link
 * index. Returns -1 if no link exists.
 */

int QfeLattice::FindLink(int s1, int s2) {
  for (int n = 0; n < sites[s1].nn; n++) {
    if (sites[s1].neighbors[n] == s2) return sites[s1].links[n];
  }
  return -1;
}

/**
 * @brief Finds the face connecting sites @p s1, @p s2, and @p s3 and returns
 * the face index. Returns -1 if no face exists.
 */

int QfeLattice::FindFace(int s1, int s2, int s3) {
  int l = FindLink(s1, s2);
  if (l == -1) return -1;

  for (int n = 0; n < links[l].n_faces; n++) {
    int f = links[l].faces[n];
    if (faces[f].n_edges != 3) continue;
    for (int e = 0; e < 3; e++) {
      if (faces[f].sites[e] == s3) return f;
    }
  }
  return -1;
}

int QfeLattice::FindCell(int s1, int s2, int s3, int s4) {
  int f = FindFace(s1, s2, s3);
  if (f == -1) return -1;

  for (int n = 0; n < faces[f].n_cells; n++) {
    int c = faces[f].cells[n];
    if (cells[c].n_faces != 4) continue;
    for (int i = 0; i < 4; i++) {
      if (cells[c].sites[i] == s4) return c;
    }
  }
  return -1;
}

/**
 * @brief Adds a link from site @p a to @p b with weight @p wt. Returns the
 * link index.
 */

int QfeLattice::AddLink(int a, int b, double wt) {
  int l = links.size();  // link index

  QfeLink link;
  link.wt = wt;
  link.sites[0] = a;
  link.sites[1] = b;
  link.n_faces = 0;

  // add site neighbors only if sites are dynamic
  if (a < n_sites) {
    int nn_a = sites[a].nn;
    sites[a].neighbors[nn_a] = b;
    sites[a].links[nn_a] = l;
    sites[a].nn++;
  }

  if (b < n_sites) {
    int nn_b = sites[b].nn;
    sites[b].neighbors[nn_b] = a;
    sites[b].links[nn_b] = l;
    sites[b].nn++;
  }

  links.push_back(link);
  n_links = links.size();
  return l;
}

/**
 * @brief Update the number of sites in each distinct group and find the
 * first site from each distinct group.
 */

void QfeLattice::UpdateDistinct() {
  distinct_n_sites.clear();
  distinct_first.clear();
  n_distinct = 0;

  for (int s = 0; s < n_sites; s++) {
    int id = sites[s].id;
    while (id >= n_distinct) {
      distinct_n_sites.push_back(0);
      distinct_first.push_back(0);
      n_distinct++;
    }

    if (distinct_n_sites[id] == 0) {
      // first site with this id
      distinct_first[id] = s;
    }
    distinct_n_sites[id]++;
    // printf("%04d %d\n", s, id);
  }
}

/**
 * @brief Add a triangular face with corner sites @p a, @p b, and @p c.
 * Returns the face index.
 */

int QfeLattice::AddFace(int a, int b, int c, double wt /*=1.0*/) {
  int f = faces.size();  // face index

  QfeFace face;
  face.wt = wt;
  face.n_edges = 3;
  face.n_cells = 0;
  int l;

  face.sites[0] = a;
  face.sites[1] = b;
  face.sites[2] = c;

  // link a to b
  l = FindLink(a, b);
  face.flip_edge[2] = false;
  if (l == -1) {
    l = AddLink(a, b, 1.0);
  } else if (links[l].sites[0] == b) {
    face.flip_edge[2] = true;
  }
  face.edges[2] = l;
  links[l].faces[links[l].n_faces] = f;
  links[l].n_faces++;

  // link b to c
  l = FindLink(b, c);
  face.flip_edge[0] = false;
  if (l == -1) {
    l = AddLink(b, c, 1.0);
  } else if (links[l].sites[0] == c) {
    face.flip_edge[0] = true;
  }
  face.edges[0] = l;
  links[l].faces[links[l].n_faces] = f;
  links[l].n_faces++;

  // link c to a
  l = FindLink(c, a);
  face.flip_edge[1] = false;
  if (l == -1) {
    l = AddLink(c, a, 1.0);
  } else if (links[l].sites[1] == a) {
    face.flip_edge[1] = true;
  }
  face.edges[1] = l;
  links[l].faces[links[l].n_faces] = f;
  links[l].n_faces++;

  faces.push_back(face);
  n_faces = faces.size();
  return f;
}

/**
 * @brief Add a rectangular face with corner sites @p a, @p b, @p c, and @p d.
 * Returns the face index.
 */

int QfeLattice::AddFace(int a, int b, int c, int d, double wt /*=1.0*/) {
  int f = faces.size();  // face index
  QfeFace face;
  face.wt = wt;
  face.n_edges = 4;
  face.n_cells = 0;
  int l;

  face.sites[0] = a;
  face.sites[1] = b;
  face.sites[2] = c;
  face.sites[3] = d;

  // link a to b
  l = FindLink(a, b);
  face.flip_edge[3] = false;
  if (l == -1) {
    l = AddLink(a, b, 1.0);
  } else if (links[l].sites[1] == b) {
    face.flip_edge[3] = true;
  }
  face.edges[3] = l;
  links[l].faces[links[l].n_faces] = f;
  links[l].n_faces++;

  // link b to c
  l = FindLink(b, c);
  face.flip_edge[0] = false;
  if (l == -1) {
    l = AddLink(b, c, 1.0);
  } else if (links[l].sites[1] == c) {
    face.flip_edge[0] = true;
  }
  face.edges[0] = l;
  links[l].faces[links[l].n_faces] = f;
  links[l].n_faces++;

  // link c to d
  l = FindLink(c, d);
  face.flip_edge[1] = false;
  if (l == -1) {
    l = AddLink(c, d, 1.0);
  } else if (links[l].sites[1] == d) {
    face.flip_edge[1] = true;
  }
  face.edges[1] = l;
  links[l].faces[links[l].n_faces] = f;
  links[l].n_faces++;

  // link d to a
  l = FindLink(d, a);
  face.flip_edge[2] = false;
  if (l == -1) {
    l = AddLink(d, a, 1.0);
  } else if (links[l].sites[1] == a) {
    face.flip_edge[2] = true;
  }
  face.edges[2] = l;
  links[l].faces[links[l].n_faces] = f;
  links[l].n_faces++;

  faces.push_back(face);
  n_faces = faces.size();
  return f;
}

/**
 * @brief Add a tetrahedal cell with corner sites @p a, @p b, @p c, and @p d.
 * Returns the cell index.
 */

int QfeLattice::AddCell(int a, int b, int c, int d, double wt) {
  int cell_id = cells.size();  // cell index
  QfeCell cell;
  cell.wt = wt;
  cell.n_faces = 4;

  // site is opposite face with same index
  cell.sites[0] = a;
  cell.sites[1] = b;
  cell.sites[2] = c;
  cell.sites[3] = d;

  int f;

  // create face bcd
  f = FindFace(b, c, d);
  if (f == -1) {
    f = AddFace(b, c, d, 1.0);
  }
  cell.faces[0] = f;
  faces[f].cells[faces[f].n_cells] = cell_id;
  faces[f].n_cells++;

  // create face acd
  f = FindFace(a, c, d);
  if (f == -1) {
    f = AddFace(a, c, d, 1.0);
  }
  cell.faces[1] = f;
  faces[f].cells[faces[f].n_cells] = cell_id;
  faces[f].n_cells++;

  // create face abd
  f = FindFace(a, b, d);
  if (f == -1) {
    f = AddFace(a, b, d, 1.0);
  }
  cell.faces[2] = f;
  faces[f].cells[faces[f].n_cells] = cell_id;
  faces[f].n_cells++;

  // create face abc
  f = FindFace(a, b, c);
  if (f == -1) {
    f = AddFace(a, b, c, 1.0);
  }
  cell.faces[3] = f;
  faces[f].cells[faces[f].n_cells] = cell_id;
  faces[f].n_cells++;

  cells.push_back(cell);
  n_cells = cells.size();
  return n_cells;
}

/**
 * @brief Refine every face in a 2-dimensional lattice by splitting each link
 * into @p n_refine sublinks and partitioning the face appropriately.
 */

void QfeLattice::Refine2D(int n_refine) {
  if (n_refine < 2) return;

  // copy the old links and faces
  std::vector<QfeLink> old_links = links;
  std::vector<QfeFace> old_faces = faces;

  // remove all links and faces
  links.clear();
  n_links = 0;
  faces.clear();
  n_faces = 0;
  for (int s = 0; s < n_sites; s++) {
    sites[s].nn = 0;
  }

  int n_old_sites = n_sites;  // offset of new sites
  int n_new_sites = 0;

  // create n-1 new sites per edge
  n_new_sites += old_links.size() * (n_refine - 1);

  // create (n-1)(n-2)/2 interior sites per face
  n_new_sites += (old_faces.size() * (n_refine - 1) * (n_refine - 2)) / 2;

  ResizeSites(n_old_sites + n_new_sites);
  vol = double(n_sites);
  for (int s = n_old_sites; s < n_sites; s++) {
    sites[s].wt = 1.0;
    sites[s].nn = 0;
    sites[s].id = 0;
  }

  // map from an ordered pair of old sites to an array of new sites running
  // along the old edge between them
  std::unordered_map<std::string, std::vector<int>> edge_sites;

  int s = n_old_sites;
  for (int l = 0; l < old_links.size(); l++) {
    // corner sites for this edge
    int s_a = old_links[l].sites[0];
    int s_b = old_links[l].sites[1];

    std::vector<int> s_edge;
    s_edge.push_back(s_a);
    for (int n = 1; n < n_refine; n++) {
      InterpolateSite(s, s_a, s_b, n, n_refine);
      sites[s].id = (n_refine - abs(2 * n - n_refine)) / 2;
      s_edge.push_back(s);
      s++;
    }
    s_edge.push_back(s_b);

    // edge sites from a to b
    char key[50];
    sprintf(key, "%d_%d", s_a, s_b);
    edge_sites[key] = s_edge;

    // edge sites from b to a
    sprintf(key, "%d_%d", s_b, s_a);
    std::reverse(s_edge.begin(), s_edge.end());
    edge_sites[key] = s_edge;
  }

  // refine interior of old faces
  for (int f = 0; f < old_faces.size(); f++) {
    int s_corner[3];
    s_corner[0] = old_faces[f].sites[0];
    s_corner[1] = old_faces[f].sites[1];
    s_corner[2] = old_faces[f].sites[2];

    std::vector<int> e_outer[3];
    std::vector<int> e_inner[3];
    char key[50];
    sprintf(key, "%d_%d", s_corner[0], s_corner[1]);
    e_outer[0] = edge_sites[key];
    sprintf(key, "%d_%d", s_corner[1], s_corner[2]);
    e_outer[1] = edge_sites[key];
    sprintf(key, "%d_%d", s_corner[2], s_corner[0]);
    e_outer[2] = edge_sites[key];

    // distinct id of the first site in the first layer
    int layer_id = n_refine / 2 + 1;

    while (e_outer[0].size() >= 3) {
      int outer_size = e_outer[0].size();
      int inner_size = outer_size - 3;

      // total number of sites in inner loop
      int n_inner_loop = (inner_size - 1) * 3;
      if (inner_size == 0) {
        n_inner_loop = 0;
      }
      if (inner_size == 1) {
        // single center point
        n_inner_loop = 1;
      }

      // first site in the inner layer
      int s_inner = s;

      // add faces to connect inner and outer sites
      for (int e = 0; e < 3; e++) {
        int ep1 = (e + 1) % 3;  // next edge
        int em1 = (e + 2) % 3;  // previous edge

        // sites for interpolating
        int s_a = e_outer[em1][outer_size - 2];
        int s_b = e_outer[ep1][1];

        // corner triangle
        AddFace(e_outer[e][0], e_outer[e][1], s_a);

        e_inner[e].clear();

        for (int i = 0; i < inner_size; i++) {
          // connect to the previous inner site
          int sm1 = s;
          if (i == 0) {
            // when i=0, the previous site is on the previous outer edge
            sm1 = s_a;
          } else if (e == 2 && i == (inner_size - 1)) {
            // for the very last point, loop back around to the first point
            s = s_inner;
          } else if (inner_size != 1) {
            // if there is more than one inner site go to the next site
            s++;
          }

          // add two faces connecting this inner site
          AddFace(e_outer[e][i + 1], s, sm1);
          AddFace(e_outer[e][i + 1], e_outer[e][i + 2], s);

          e_inner[e].push_back(s);
          InterpolateSite(s, s_a, s_b, i + 1, inner_size + 1);
          sites[s].id =
              layer_id + (inner_size - 1 - abs(2 * i - (inner_size - 1))) / 2;
        }
      }

      // add center triangle for special cases
      if (inner_size == 0) {
        AddFace(e_outer[0][1], e_outer[1][1], e_outer[2][1]);
      } else if (inner_size == 2) {
        AddFace(e_inner[0][0], e_inner[1][0], e_inner[2][0]);
      }

      // go to the next layer on the interior of the face
      s = s_inner + n_inner_loop;
      e_outer[0] = e_inner[0];
      e_outer[1] = e_inner[1];
      e_outer[2] = e_inner[2];

      layer_id += (inner_size - 1) / 2 + 1;
      // printf("layer_id: %d\n", layer_id);
    }
  }
}

/**
 * @brief Print a list of sites with their weights and neighbors.
 */

void QfeLattice::PrintSites() {
  printf("\n*** sites ***\n");
  for (int s = 0; s < sites.size(); s++) {
    printf("%04d", s);
    printf(" %.12f", sites[s].wt);
    for (int n = 0; n < sites[s].nn; n++) {
      printf(" %04d", sites[s].neighbors[n]);
    }
    printf("\n");
  }
  printf("*************\n");
}

/**
 * @brief Print a list of links with their weights and attached sites.
 */

void QfeLattice::PrintLinks() {
  printf("\n*** links ***\n");
  for (int l = 0; l < links.size(); l++) {
    printf("%04d", l);
    printf(" %.12f", links[l].wt);
    printf(" %04d", links[l].sites[0]);
    printf(" %04d", links[l].sites[1]);
    for (int n = 0; n < links[l].n_faces; n++) {
      printf(" %04d", links[l].faces[n]);
    }
    printf("\n");
  }
  printf("*************\n");
}

void QfeLattice::PrintFaces() {
  printf("\n*** faces ***\n");
  for (int f = 0; f < faces.size(); f++) {
    printf("%04d", f);
    printf(" %.12f", faces[f].wt);
    for (int n = 0; n < faces[f].n_edges; n++) {
      printf(" %04d", faces[f].edges[n]);
    }
    printf("\n");
  }
  printf("*************\n");
}

void QfeLattice::PrintCells() {
  printf("\n*** cells ***\n");
  for (int c = 0; c < cells.size(); c++) {
    printf("%04d", c);
    printf(" %.12f", cells[c].wt);
    for (int n = 0; n < cells[c].n_faces; n++) {
      printf(" %04d", cells[c].faces[n]);
    }
    printf("\n");
  }
  printf("*************\n");
}

/**
 * @brief Check that all lattice sites are connected.
 */

void QfeLattice::CheckConnectivity() {
  printf("\n*** connectivity check ***\n");
  printf("n_sites: %d\n", n_sites);

  // keep track of which sites are connected
  std::vector<bool> is_connected(sites.size());

  // create the stack
  std::stack<int> stack;

  // start with site 0
  stack.push(0);
  is_connected[0] = true;

  int n_connected = 0;

  while (stack.size() != 0) {
    n_connected++;
    int s = stack.top();
    stack.pop();
    QfeSite* site = &sites[s];

    for (int n = 0; n < site->nn; n++) {
      int s = site->neighbors[n];
      if (is_connected[s]) continue;
      is_connected[s] = true;
      stack.push(s);
    }
  }

  printf("connected sites: %d\n", n_connected);
  printf("disconnected sites: %d\n", int(sites.size()) - n_connected);
  printf("**************************\n");
}

/**
 * @brief Check that the neighbor lists match the link sites.
 */

void QfeLattice::CheckConsistency() {
  printf("\n*** consistency check ***\n");

  int n_inconsistent = 0;

  // make sure each site's neighbor table is consistent with its links
  for (int l = 0; l < links.size(); l++) {
    int s_a = links[l].sites[0];
    int s_b = links[l].sites[1];
    int n;

    // check 1st site
    for (n = 0; n < sites[s_a].nn; n++) {
      if (sites[s_a].links[n] == l) break;
    }

    if (n == sites[s_a].nn) {
      // link not found
      printf("link %04d not found in neighbor table for site %04d\n", l, s_a);
      n_inconsistent++;
    } else if (sites[s_a].neighbors[n] != s_b) {
      printf("site %04d neighbor %04d mismatch (link %04d, site %04d)\n", s_a,
             sites[s_a].neighbors[n], l, s_b);
      n_inconsistent++;
    }

    // check 2nd site
    for (n = 0; n < sites[s_b].nn; n++) {
      if (sites[s_b].links[n] == l) break;
    }

    if (n == sites[s_b].nn) {
      // link not found
      printf("link %04d not found in neighbor table for site %04d\n", l, s_b);
      n_inconsistent++;
    } else if (sites[s_b].neighbors[n] != s_a) {
      printf("site %04d neighbor %04d mismatch (link %04d, site %04d)\n", s_b,
             sites[s_b].neighbors[n], l, s_a);
      n_inconsistent++;
    }
  }
  printf("%d inconsistencies found\n", n_inconsistent);
  printf("*************************\n");
}

// ising.h

#pragma once

#include <cmath>
#include <stack>
#include <string>
#include <unordered_map>
#include <vector>

#include "lattice.h"

/// @brief Ising spin model
class QfeIsing {
 public:
  QfeIsing(QfeLattice* lattice, double beta);
  void WriteField(FILE* file);
  void ReadField(FILE* file);
  double Action();
  double MeanSpin();
  void HotStart();
  void ColdStart();
  double Metropolis();
  int WolffUpdate();
  int SWUpdate();
  int FindSWRoot(int s);

  QfeLattice* lattice;
  std::vector<double> spin;  // Z2 field
  double beta;               // bare coupling

  std::vector<bool> is_clustered;  // keeps track of which sites are clustered
  std::vector<int> wolff_cluster;  // array of clustered sites
  std::vector<int> sw_root;        // root for each site
};

/// @brief Ising model constructor
/// @param lattice The Ising spins are defined on the lattice sites
/// @param beta spin coupling coefficient (use link weights for local couplings)
QfeIsing::QfeIsing(QfeLattice* lattice, double beta) {
  this->lattice = lattice;
  this->beta = beta;
  spin.resize(lattice->n_sites);
  is_clustered.resize(lattice->n_sites);
  sw_root.resize(lattice->n_sites);
}

/// @brief Write the Ising spins to a data file. Each spin value is packed into
/// a single bit to reduce file size.
/// @param file Data file
void QfeIsing::WriteField(FILE* file) {
  int buf_size = (spin.size() + 7) / 8;
  std::vector<unsigned char> spin_buf(buf_size, 0);

  char mask = 1;
  for (int i = 0; i < spin.size(); i++) {
    if (spin[i] == -1) {
      spin_buf[i / 8] |= mask;
    }
    mask <<= 1;
    if (!mask) mask = 1;
  }

  fwrite(spin_buf.data(), 1, buf_size, file);
}

/// @brief Read the Ising spins from a data file
/// @param file Data file
void QfeIsing::ReadField(FILE* file) {
  int buf_size = (spin.size() + 7) / 8;
  std::vector<unsigned char> spin_buf(buf_size, 0);
  fread(spin_buf.data(), 1, buf_size, file);

  ColdStart();
  char mask = 1;
  for (int i = 0; i < spin.size(); i++) {
    if (spin_buf[i / 8] & mask) {
      spin[i] = -1;
    }
    mask <<= 1;
    if (!mask) mask = 1;
  }
}

/// @brief Calculate the Ising model action
/// @return The action per spin
double QfeIsing::Action() {
  double action = 0.0;

  // sum over links
  for (int l = 0; l < lattice->n_links; l++) {
    QfeLink* link = &lattice->links[l];
    int a = link->sites[0];
    int b = link->sites[1];
    action -= beta * spin[a] * spin[b] * link->wt;
  }

  return action / lattice->vol;
}

/// @brief Compute the mean spin value, weighted by site weights.
/// @return The mean spin (i.e. magnetization)
double QfeIsing::MeanSpin() {
  double m = 0.0;
  for (int s = 0; s < lattice->n_sites; s++) {
    m += spin[s] * lattice->sites[s].wt;
  }
  return m / lattice->vol;
}

/// @brief Randomize all spins
void QfeIsing::HotStart() {
  for (int s = 0; s < lattice->n_sites; s++) {
    if (lattice->rng.RandBool()) {
      spin[s] = 1.0;
    } else {
      spin[s] = -1.0;
    }
  }
}

/// @brief Set all spins to +1
void QfeIsing::ColdStart() { std::fill(spin.begin(), spin.end(), 1.0); }

/// @brief Metropolis update algorithm (i.e. Rosenbluth-Teller)
/// ref: N. Metropolis, et al., J. Chem. Phys. 21, 1087 (1953).
/// @return The acceptance rate
double QfeIsing::Metropolis() {
  int accept = 0;
  for (int s = 0; s < lattice->n_sites; s++) {
    double delta_S = 0.0;

    // sum over links connected to this site
    QfeSite* site = &lattice->sites[s];
    for (int n = 0; n < site->nn; n++) {
      int l = site->links[n];
      double link_wt = lattice->links[l].wt;
      delta_S += beta * spin[site->neighbors[n]] * link_wt;
    }
    delta_S *= 2.0 * spin[s];

    // metropolis algorithm
    if (delta_S <= 0.0 || lattice->rng.RandReal() < exp(-delta_S)) {
      spin[s] *= -1.0;
      accept++;
    }
  }
  return double(accept) / double(lattice->n_sites);
}

/// @brief Wolff cluster update algorithm
/// ref: U. Wolff, Phys. Rev. Lett. 62, 361 (1989).
/// @return Cluster size
int QfeIsing::WolffUpdate() {
  // remove all sites from the cluster
  std::fill(is_clustered.begin(), is_clustered.end(), false);
  wolff_cluster.clear();

  // create the stack
  std::stack<int> stack;

  // choose a random site and add it to the cluster
  int s = lattice->rng.RandInt(0, lattice->n_sites - 1);
  wolff_cluster.push_back(s);
  is_clustered[s] = true;
  stack.push(s);

  while (stack.size() != 0) {
    s = stack.top();
    stack.pop();

    // flip the spin
    double value = spin[s];
    spin[s] = -value;

    // try to add neighbors
    QfeSite* site = &lattice->sites[s];
    for (int n = 0; n < site->nn; n++) {
      int l = site->links[n];
      double link_wt = lattice->links[l].wt;
      s = site->neighbors[n];

      // skip if the site is already clustered
      if (is_clustered[s]) continue;

      // check if link is clustered
      double rate = -2.0 * beta * value * spin[s] * link_wt;
      if (rate >= 0.0 || lattice->rng.RandReal() < exp(rate)) continue;

      wolff_cluster.push_back(s);
      is_clustered[s] = true;
      stack.push(s);
    }
  }

  return wolff_cluster.size();
}

/// @brief Swendsen-Wang update algorithm
/// ref: R.H. Swendsen and J.S. Wang, Phys. Rev. Lett. 58, 86 (1987)
/// @return The total number of clusters
int QfeIsing::SWUpdate() {
  // each site begins in its own cluster
  std::iota(std::begin(sw_root), std::end(sw_root), 0);

  // loop over all links
  for (int l = 0; l < lattice->n_links; l++) {
    int s1 = lattice->links[l].sites[0];
    int s2 = lattice->links[l].sites[1];
    double link_wt = lattice->links[l].wt;

    // check if link is clustered
    double rate = -2.0 * beta * spin[s1] * spin[s2] * link_wt;
    if (rate >= 0.0 || lattice->rng.RandReal() < exp(rate)) continue;

    // find the root node for each site
    int r1 = FindSWRoot(s1);
    int r2 = FindSWRoot(s2);

    if (r1 == r2) continue;
    int r = std::min(r1, r2);
    sw_root[r1] = r;
    sw_root[r2] = r;
  }

  std::unordered_map<int, int> sw_map;
  std::vector<bool> is_flipped;
  int n_clusters = 0;
  for (int s = 0; s < lattice->n_sites; s++) {
    // find the root node for this site
    int r = FindSWRoot(s);

    if (sw_map.find(r) == sw_map.end()) {
      sw_map[r] = n_clusters;
      is_flipped.push_back(lattice->rng.RandReal() > 0.5);
      n_clusters++;
    }

    int c = sw_map[r];

    // flip half the clusters
    if (is_flipped[c]) spin[s] = -spin[s];
  }

  return n_clusters;
}

/// @brief Find the root node for a Swendsen-Wang cluster. Also update the trail
/// from the starting node to the root node to all be in the same cluster
/// @param s Starting node
/// @return Root node for this cluster
int QfeIsing::FindSWRoot(int s) {
  int root = sw_root[s];

  // find the root
  while (root != sw_root[root]) root = sw_root[root];

  // update the trail from s to the root
  while (s != root) {
    int old_root = sw_root[s];
    sw_root[s] = root;
    s = old_root;
  }

  return root;
}

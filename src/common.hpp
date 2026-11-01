#pragma once
#include <bits/stdc++.h>
#include <boost/unordered/unordered_flat_map.hpp>
#include <boost/unordered/unordered_flat_set.hpp>
using LL = long long;
using VI = std::vector<int>;
using VVI = std::vector<VI>;
using PII = std::pair<int, int>;
using PID = std::pair<int, double>;
using PD_PII = std::pair<double, PII>;
using hclock = std::chrono::steady_clock;
template <class K, class V> using Map = boost::unordered_flat_map<K, V>;
template <class K> using Set = boost::unordered_flat_set<K>;
inline double seconds(hclock::time_point t) {
  return std::chrono::duration<double>(hclock::now() - t).count();
}
static_assert(sizeof(int) == 4 && sizeof(PII) == 8,
              "binary format requires 32-bit endpoints");
struct Vertex {
  VI nbr;
  int degree = 0;
};
struct Graph {
  int vertices_num_ = 0;
  LL m = 0;
  std::vector<Vertex *> vertices_;
  std::vector<PII> edges;
  Graph(const std::string &path) {
    std::ifstream f(path, std::ios::binary);
    int64_t n64, m64;
    f.read((char *)&n64, 8);
    f.read((char *)&m64, 8);
    if (!f || n64 <= 0 || n64 > INT_MAX || m64 < 0 || m64 > INT_MAX)
      throw std::runtime_error("bad graph header");
    vertices_num_ = n64;
    m = m64;
    edges.resize(m);
    f.read((char *)edges.data(), m * sizeof(PII));
    if (!f)
      throw std::runtime_error("truncated graph");
    vertices_.resize(vertices_num_);
    VI deg(vertices_num_);
    if (!std::is_sorted(edges.begin(), edges.end()) ||
        std::adjacent_find(edges.begin(), edges.end()) != edges.end())
      throw std::runtime_error(
          "input edges must be sorted and unique; run preprocess.py first");
    for (auto [u, v] : edges) {
      if (u < 0 || v < 0 || u >= vertices_num_ || v >= vertices_num_ || u >= v)
        throw std::runtime_error("invalid edge");
      deg[u]++;
      deg[v]++;
    }
    for (int u = 0; u < vertices_num_; u++) {
      vertices_[u] = new Vertex;
      vertices_[u]->nbr.reserve(deg[u]);
    }
    for (auto [u, v] : edges) {
      vertices_[u]->nbr.push_back(v);
      vertices_[v]->nbr.push_back(u);
    }
    for (auto x : vertices_)
      x->degree = x->nbr.size();
    order_nbrs();
  }
  ~Graph() {
    for (auto p : vertices_)
      delete p;
  }
  void order_nbrs() {
    for (auto p : vertices_)
      std::sort(p->nbr.begin(), p->nbr.end());
  }
};
struct DisjointSetUnion {
  VI fa;
  explicit DisjointSetUnion(int n) : fa(n + 1) {
    std::iota(fa.begin(), fa.end(), 0);
  }
  int find(int x) { return x == fa[x] ? x : (fa[x] = find(fa[x])); }
  void joint(int u, int v) {
    u = find(u);
    v = find(v);
    if (u != v)
      fa[v] = u;
  }
};
inline int num_common_nbr(const VI &x, const VI &y) {
  int ix = 0, iy = 0, ans = 0;
  while (ix < (int)x.size() && iy < (int)y.size()) {
    if (x[ix] < y[iy])
      ++ix;
    else if (x[ix] == y[iy])
      ++ans, ++ix, ++iy;
    else
      ++iy;
  }
  return ans;
}
inline PII min_pair(int u, int v) { return {std::min(u, v), std::max(u, v)}; }
inline uint64_t key(int a, int b) {
  if (a > b)
    std::swap(a, b);
  return (uint64_t(uint32_t(a)) << 32) | uint32_t(b);
}
struct Encoding {
  LL blocks = 0, se = 0, cp = 0, cm = 0, loops = 0, cost = 0, native = 0;
};
inline Encoding encode(Graph &g, const VI &sid, bool signed_mode,
                       const std::string &out) {
  int n = g.vertices_num_, k = *std::max_element(sid.begin(), sid.end()) + 1;
  VI sz(k);
  VVI members(k);
  for (int u = 0; u < n; u++) {
    sz[sid[u]]++;
    members[sid[u]].push_back(u);
  }
  Map<uint64_t, int> cnt;
  cnt.reserve(g.m);
  for (auto [u, v] : g.edges)
    cnt[key(sid[u], sid[v])]++;
  Encoding s;
  for (int x : sz)
    if (x)
      s.blocks++;
  std::vector<PII> se, cp, cm;
  Set<uint64_t> selected;
  for (auto [kv, e] : cnt) {
    int a = kv >> 32, b = uint32_t(kv);
    LL t = a == b ? LL(sz[a]) * (sz[a] - 1) / 2 : LL(sz[a]) * sz[b];
    bool use = signed_mode ? (1 + t - e < e) : (e == t && e > 1);
    if (use) {
      s.se++;
      s.cm += t - e;
      s.loops += (a == b);
      if (!out.empty()) {
        se.emplace_back(a, b);
        selected.insert(kv);
      }
    } else
      s.cp += e;
    LL unit = a == b ? 1 : 2;
    s.native += signed_mode ? std::min(2LL * e, unit + 2 * (t - e))
                            : ((e == t && t > 0) ? unit : 2LL * e);
  }
  s.cost = s.se + s.cp + s.cm;
  if (!out.empty()) {
    auto parent = std::filesystem::path(out).parent_path();
    if (!parent.empty())
      std::filesystem::create_directories(parent);
    std::ofstream pf(out + ".partition.bin", std::ios::binary);
    pf.write((char *)sid.data(), sid.size() * sizeof(int));
    pf.close();
    if (!pf)
      throw std::runtime_error("cannot write partition output");
    for (auto [u, v] : g.edges)
      if (!selected.contains(key(sid[u], sid[v])))
        cp.emplace_back(u, v);
    if (signed_mode)
      for (auto [a, b] : se)
        for (int u : members[a])
          for (int v : members[b]) {
            if (a == b && v <= u)
              continue;
            if (!std::binary_search(g.vertices_[u]->nbr.begin(),
                                    g.vertices_[u]->nbr.end(), v))
              cm.push_back(min_pair(u, v));
          }
    auto save = [&](std::string ext, std::vector<PII> &a) {
      std::sort(a.begin(), a.end());
      std::ofstream f(out + ext, std::ios::binary);
      f.write((char *)a.data(), a.size() * sizeof(PII));
      f.close();
      if (!f)
        throw std::runtime_error("cannot write summary output");
    };
    save(".se.bin", se);
    save(".cp.bin", cp);
    save(".cm.bin", cm);
    if ((LL)cp.size() != s.cp || (LL)cm.size() != s.cm)
      throw std::runtime_error("encoding count mismatch");
  }
  return s;
}
inline void printenc(const Encoding &e, LL m) {
  std::cout << "\"blocks\":" << e.blocks << ",\"superedges\":" << e.se
            << ",\"positive_corrections\":" << e.cp
            << ",\"negative_corrections\":" << e.cm
            << ",\"superloops\":" << e.loops << ",\"cost\":" << e.cost
            << ",\"relative_size\":" << (m ? double(e.cost) / m : 0.0)
            << ",\"native_opt_cost\":" << e.native
            << ",\"native_relative_size\":"
            << (m ? double(e.native) / (2 * m) : 0.0);
}

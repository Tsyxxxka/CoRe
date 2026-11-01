#include "common.hpp"

namespace core {
constexpr int THREE_FOCUS_INTERVAL = 4;
constexpr int MAX_FOCUS_VERTICES = 3, MAX_REPARTITION_ATOMS = 2 * MAX_FOCUS_VERTICES;
struct ConstructionStats {
  LL merges = 0, merge_evals = 0, merge_cache_hits = 0, windows = 0,
     commits = 0;
  LL skipped = 0, dfs_states = 0, focus_edge_scans = 0, parent_entries = 0;
  LL relocation_rounds = 0, relocation_candidates = 0, relocation_attempts = 0,
     relocation_commits = 0, relocation_gain = 0;
};

struct ScoreRatio {
  uint32_t num = 0, den = 1;
  ScoreRatio() = default;
  ScoreRatio(LL p, LL q) : num(uint32_t(p)), den(uint32_t(q)) {
#ifndef NDEBUG
    if (p < 0 || p > UINT32_MAX || q <= 0 || q > UINT32_MAX)
      throw std::runtime_error("ratio out of range");
#endif
  }
};
inline int compare(ScoreRatio a, ScoreRatio b) {
  uint64_t x = uint64_t(a.num) * b.den, y = uint64_t(b.num) * a.den;
  return (x > y) - (x < y);
}
struct SamplingConfig {
  int direct_samples = 64, low_degree_pivots = 6, pivot_fanout = 32, top_k = 4, count_mode = 0;
};
struct CandidateCountStats {
  LL pool_vertices = 0, direct_count_vertices = 0, bulk_count_vertices = 0, direct_count_pairs = 0,
     binary_count_pairs = 0, ub_pruned = 0, scan_probes = 0, binary_probes = 0,
     wedge_probes = 0;
};
struct CandidateMatrix {
  int width = 0;
  VI values;
  void reset(int n, int k, int fill) {
    width = k;
    values.assign(size_t(n) * k, fill);
  }
  std::span<int> operator[](int u) {
    return {values.data() + size_t(u) * width, size_t(width)};
  }
  std::span<const int> operator[](int u) const {
    return {values.data() + size_t(u) * width, size_t(width)};
  }
};
struct CandidateData {
  CandidateMatrix candidate_sets, common_neighbor_counts;
  VI vertex_order;
  int top_k = 4;
  CandidateCountStats stats;
};
struct RankedCandidate {
  ScoreRatio score;
  int v = 0, common_neighbors = 0;
};
inline bool better(const RankedCandidate &a, const RankedCandidate &b) {
  int z = compare(a.score, b.score);
  return z ? z > 0 : a.v > b.v;
}
inline int binary_depth(unsigned n) { return std::bit_width(n); }
inline bool binary_member(const VI &a, int x, LL &probes) {
  int lo = 0, hi = int(a.size()) - 1;
  while (lo <= hi) {
    int mid = lo + (hi - lo) / 2, y = a[mid];
    probes++;
    if (x == y)
      return true;
    if (y < x)
      lo = mid + 1;
    else
      hi = mid - 1;
  }
  return false;
}
CandidateData CandidateGeneration(const Graph &g, const SamplingConfig &cfg) {
  int n = g.vertices_num_;
  CandidateData out;
  out.top_k = cfg.top_k;
  out.candidate_sets.reset(n, out.top_k, -1);
  out.common_neighbor_counts.reset(n, out.top_k, 0);
  out.vertex_order.resize(n);
  std::iota(out.vertex_order.begin(), out.vertex_order.end(), 0);
  VI seen(n, -1), adjmark(n, -1), codeg(n), degree(n), pool, piv;
  for (int u = 0; u < n; u++)
    degree[u] = g.vertices_[u]->degree;
  std::vector<ScoreRatio> quality(n);
  std::vector<RankedCandidate> heap;
  for (int u = 0; u < n; u++) {
    const VI &nu = g.vertices_[u]->nbr;
    LL du = nu.size();
    pool.clear();
    heap.clear();
    for (int x : nu)
      adjmark[x] = u;
    auto add = [&](int v) {
      if (v != u && seen[v] != u) {
        seen[v] = u;
        pool.push_back(v);
      }
    };
    int take = std::min<int>(cfg.direct_samples, nu.size());
    for (int i = 0; i < take; i++)
      add(nu[size_t(i) * nu.size() / take]);
    piv = nu;
    int np = std::min<int>(cfg.low_degree_pivots, piv.size());
    std::partial_sort(
        piv.begin(), piv.begin() + np, piv.end(), [&](int a, int b) {
          return degree[a] != degree[b] ? degree[a] < degree[b] : a < b;
        });
    for (int j = 0; j < np; j++) {
      const VI &nx = g.vertices_[piv[j]]->nbr;
      int q = std::min<int>(cfg.pivot_fanout, nx.size());
      for (int i = 0; i < q; i++)
        add(nx[size_t(i) * nx.size() / q]);
    }
    out.stats.pool_vertices += pool.size();
    LL W = 0, Q = 0;
    for (int x : nu)
      W += degree[x];
    for (int v : pool)
      Q += std::min<LL>(degree[v], du * binary_depth(unsigned(degree[v])));
    bool bulk =
        cfg.count_mode == 2 || (cfg.count_mode == 0 && W + LL(pool.size()) < Q);
    if (bulk) {
      out.stats.bulk_count_vertices++;
      for (int v : pool)
        codeg[v] = 0;
      for (int x : nu)
        for (int v : g.vertices_[x]->nbr) {
          out.stats.wedge_probes++;
          if (seen[v] == u)
            codeg[v]++;
        }
    } else
      out.stats.direct_count_vertices++;
    for (int v : pool) {
      LL dv = degree[v];
      bool adjacent = adjmark[v] == u;
      if (!bulk && (int)heap.size() == out.top_k) {
        LL cbar = std::min(du, dv) - int(adjacent);
        RankedCandidate upper{{cbar, std::max(1LL, du + dv - cbar)}, v, 0};
        if (!better(upper, heap.front())) {
          out.stats.ub_pruned++;
          continue;
        }
      }
      int c = 0;
      if (bulk)
        c = codeg[v];
      else {
        out.stats.direct_count_pairs++;
        const VI &nv = g.vertices_[v]->nbr;
        if (dv <= du * binary_depth(unsigned(dv))) {
          for (int x : nv) {
            out.stats.scan_probes++;
            c += adjmark[x] == u;
          }
        } else {
          out.stats.binary_count_pairs++;
          for (int x : nu)
            c += binary_member(nv, x, out.stats.binary_probes);
        }
      }
      if (c == 0 && !adjacent)
        continue;
      RankedCandidate rank{{c, du + dv - c}, v, c};
      if ((int)heap.size() < out.top_k) {
        heap.push_back(rank);
        std::push_heap(heap.begin(), heap.end(), better);
      } else if (better(rank, heap.front())) {
        std::pop_heap(heap.begin(), heap.end(), better);
        heap.back() = rank;
        std::push_heap(heap.begin(), heap.end(), better);
      }
    }
    std::sort(heap.begin(), heap.end(), better);
    for (int j = 0; j < (int)heap.size(); j++) {
      out.candidate_sets[u][j] = heap[j].v;
      out.common_neighbor_counts[u][j] = heap[j].common_neighbors;
    }
    quality[u] = heap.empty() ? ScoreRatio{} : heap.front().score;
  }
  std::sort(out.vertex_order.begin(), out.vertex_order.end(), [&](int a, int b) {
    int z = compare(quality[a], quality[b]);
    return z ? z > 0 : a < b;
  });
  return out;
}
struct CachedRejectedMerge {
  uint64_t av = 0, bv = 0;
};
struct ExternalRelationRow {
  int id;
  std::array<int, MAX_REPARTITION_ATOMS> e{};
};
struct RepartitionAtom {
  int parent = -1, vertex = -1, size = 0;
}; 

struct SupernodePartition {
  const Graph &g;
  bool allow_negative_corrections;
  int n;
  VI supernode_of, available_supernode_ids, focus_vertex_index, external_supernode_index, supernode_size, internal_edge_count;
  std::vector<LL> degree_volume;
  VVI supernode_vertices;
  std::vector<Map<int, int>> relation_edge_counts;
  std::vector<uint64_t> partition_revision;
  Map<uint64_t, CachedRejectedMerge> failed;
  bool use_gain_cache = true, fast_gain_evaluation = true;
  VVI complete_neighbors;
  std::vector<uint64_t> complete_neighbors_revision;
  void prepare_complete_neighbors(int a) {
    if (complete_neighbors_revision[a] == partition_revision[a])
      return;
    auto &row = complete_neighbors[a];
    row.clear();
    LL sa = supernode_size[a];
    for (auto [x, e] : relation_edge_counts[a])
      if (x != a && LL(e) == sa * supernode_size[x])
        row.push_back(x);
    complete_neighbors_revision[a] = partition_revision[a];
  }
  LL cost = 0;
  ConstructionStats stats;
  std::vector<LL> complete_relation_savings;
  VI internal_edge_cache;
  bool complete_savings_current = true;
  std::vector<ExternalRelationRow> ext; 
  SupernodePartition(const Graph &gg, VI init, bool signed_mode, bool cache = true)
      : g(gg), allow_negative_corrections(signed_mode), n(gg.vertices_num_), supernode_of(std::move(init)),
        focus_vertex_index(n, -1), external_supernode_index(n, -1), supernode_size(n, 0), internal_edge_count(n, 0),
        degree_volume(n, 0), supernode_vertices(n), relation_edge_counts(n), partition_revision(n, 1), use_gain_cache(cache) {
    bool singleton = true;
    for (int u = 0; u < n; u++) {
      supernode_vertices[supernode_of[u]].push_back(u);
      supernode_size[supernode_of[u]]++;
      degree_volume[supernode_of[u]] += g.vertices_[u]->degree;
      singleton &= (supernode_of[u] == u);
    }
    for (int a = n - 1; a >= 0; a--)
      if ((supernode_size[a] == 0))
        available_supernode_ids.push_back(a);
    if (singleton) {
      
      for (int u = 0; u < n; u++) {
        relation_edge_counts[u].reserve(g.vertices_[u]->nbr.size());
        for (int v : g.vertices_[u]->nbr)
          relation_edge_counts[u].emplace(v, 1);
      }
    } else
      for (auto [u, v] : g.edges) {
        int a = supernode_of[u], b = supernode_of[v];
        relation_edge_counts[a][b]++;
        if (a != b)
          relation_edge_counts[b][a]++;
        else
          internal_edge_count[a]++;
      }
    complete_neighbors.resize(n);
    complete_neighbors_revision.assign(n, 0);
    complete_relation_savings.assign(n, 0);
    internal_edge_cache = internal_edge_count;
    if (!singleton)
      for (int a = 0; a < n; a++)
        if (supernode_size[a])
          for (auto [x, e] : relation_edge_counts[a])
            if (x != a && LL(e) == LL(supernode_size[a]) * supernode_size[x])
              complete_relation_savings[a] += LL(e) - 1;
    if (use_gain_cache)
      failed.reserve(size_t(n) * 4);
    cost = total();
  }
  LL f(LL cap, LL e) const {
#ifndef NDEBUG
    if (e < 0 || e > cap)
      throw std::runtime_error("edge count outside block capacity");
#endif
    return allow_negative_corrections ? std::min(e, 1 + cap - e) : ((e == cap && cap > 0) ? 1 : e);
  }
  int get(int a, int b) const {
    if (a == b)
      return internal_edge_count[a];
    auto it = relation_edge_counts[a].find(b);
    return it == relation_edge_counts[a].end() ? 0 : it->second;
  }
  LL total() const {
    LL result = 0;
    for (int a = 0; a < n; a++)
      if (!(supernode_size[a] == 0))
        for (auto [b, e] : relation_edge_counts[a])
          if (b >= a) {
            LL sa = supernode_size[a], sb = supernode_size[b];
            result += f(a == b ? sa * (sa - 1) / 2 : sa * sb, e);
          }
    return result;
  }
  std::pair<LL, LL> evaluate_merge_gain(int a, int b) {
    if (a > b)
      std::swap(a, b);
    uint64_t k = key(a, b);
    if (use_gain_cache) {
      auto it = failed.find(k);
      if (it != failed.end() && it->second.av == partition_revision[a] &&
          it->second.bv == partition_revision[b]) {
        stats.merge_cache_hits++;
        return {0, 1};
      }
    }
    stats.merge_evals++;
    LL sa = supernode_size[a], sb = supernode_size[b], sw = sa + sb;
    int aa = get(a, a), bb = get(b, b), ab = get(a, b);
    LL before =
        f(sa * (sa - 1) / 2, aa) + f(sb * (sb - 1) / 2, bb) + f(sa * sb, ab);
    LL gain;
    if (!allow_negative_corrections && fast_gain_evaluation) {
      LL wab = (LL(ab) == sa * sb) ? LL(ab) - 1 : 0;
      LL ia = f(sa * (sa - 1) / 2, aa), ib = f(sb * (sb - 1) / 2, bb);
      before = ia + degree_volume[a] - 2LL * aa - complete_relation_savings[a] + ib + degree_volume[b] -
               2LL * bb - complete_relation_savings[b] - f(sa * sb, ab);
      gain = ia + ib + f(sa * sb, ab) - f(sw * (sw - 1) / 2, LL(aa) + bb + ab) -
             complete_relation_savings[a] - complete_relation_savings[b] + 2 * wab;
      int small = a, other = b;
      if (relation_edge_counts[small].size() > relation_edge_counts[other].size())
        std::swap(small, other);
      for (auto [x, e] : relation_edge_counts[small])
        if (x != a && x != b && LL(e) == LL(supernode_size[small]) * supernode_size[x] &&
            LL(get(other, x)) == LL(supernode_size[other]) * supernode_size[x])
          gain += sw * supernode_size[x] - 1;
    } else if (!allow_negative_corrections) {
      gain = before - f(sw * (sw - 1) / 2, LL(aa) + bb + ab);
      for (auto [x, e] : relation_edge_counts[a])
        if (x != a && x != b) {
          LL sx = supernode_size[x], cap = sa * sx;
          before += (e == cap ? 1 : e);
          if (e == cap) {
            gain -= cap - 1;
            if (get(b, x) == sb * sx)
              gain += sw * sx - 1;
          }
        }
      for (auto [x, e] : relation_edge_counts[b])
        if (x != a && x != b) {
          LL cap = sb * LL(supernode_size[x]);
          before += (e == cap ? 1 : e);
          if (e == cap)
            gain -= cap - 1;
        }
    } else {
      LL after = f(sw * (sw - 1) / 2, LL(aa) + bb + ab);
      for (auto [x, e] : relation_edge_counts[a])
        if (x != a && x != b) {
          int eb = get(b, x);
          LL sx = supernode_size[x];
          before += f(sa * sx, e) + f(sb * sx, eb);
          after += f(sw * sx, LL(e) + eb);
        }
      for (auto [x, e] : relation_edge_counts[b])
        if (x != a && x != b && !relation_edge_counts[a].contains(x)) {
          LL sx = supernode_size[x];
          before += f(sb * sx, e);
          after += f(sw * sx, e);
        }
      gain = before - after;
    }
    if (use_gain_cache && gain <= 0)
      failed[k] = {partition_revision[a], partition_revision[b]};
    return {gain, before};
  }
  void merge(int a, int b, LL gain) {
    if (supernode_size[a] < supernode_size[b])
      std::swap(a, b);
    int internal = get(a, a) + get(b, b) + get(a, b);
    
    LL sa = supernode_size[a], sb = supernode_size[b], sw = sa + sb, new_penalty = 0;
    for (auto [x, ea] : relation_edge_counts[a])
      if (x != a && x != b) {
        int eb = get(b, x);
        LL sx = supernode_size[x];
        LL olda = LL(ea) == sa * sx ? LL(ea) - 1 : 0,
           oldb = LL(eb) == sb * sx ? LL(eb) - 1 : 0;
        LL now = LL(ea) + eb == sw * sx ? LL(ea) + eb - 1 : 0;
        complete_relation_savings[x] += now - olda - oldb;
        new_penalty += now;
        partition_revision[x]++;
      }
    for (auto [x, eb] : relation_edge_counts[b])
      if (x != a && x != b && !relation_edge_counts[a].contains(x)) {
        if (LL(eb) == sb * supernode_size[x])
          complete_relation_savings[x] -= LL(eb) - 1;
        partition_revision[x]++;
      }
    complete_relation_savings[a] = new_penalty;
    complete_relation_savings[b] = 0;
    partition_revision[a]++;
    partition_revision[b]++;
    internal_edge_cache[a] = internal;
    internal_edge_cache[b] = 0;
    for (auto [x, e] : relation_edge_counts[b])
      if (x != a && x != b) {
        relation_edge_counts[a][x] += e;
        relation_edge_counts[x].erase(b);
        relation_edge_counts[x][a] = relation_edge_counts[a][x];
      }
    relation_edge_counts[a].erase(b);
    if (internal)
      relation_edge_counts[a][a] = internal;
    else
      relation_edge_counts[a].erase(a);
    relation_edge_counts[b].clear();
    for (int u : supernode_vertices[b]) {
      supernode_of[u] = a;
      supernode_vertices[a].push_back(u);
    }
    internal_edge_count[a] = internal;
    internal_edge_count[b] = 0;
    supernode_size[a] += supernode_size[b];
    supernode_size[b] = 0;
    degree_volume[a] += degree_volume[b];
    degree_volume[b] = 0;
    supernode_vertices[b].clear();
    available_supernode_ids.push_back(b);
    cost -= gain;
    stats.merges++;
  }
  void PositiveCoarsen(const CandidateData &candidates) {
    for (int pass = 0;; pass++) {
      LL start = cost, accepted = 0;
      for (int u : candidates.vertex_order) {
        int a = supernode_of[u], best = -1;
        LL bg = 0;
        ScoreRatio bs{};
        VI tried(candidates.top_k);
        int nt = 0;
        for (int j = 0; j < candidates.top_k; j++) {
          int v = candidates.candidate_sets[u][j];
          if (v < 0)
            break;
          int b = supernode_of[v];
          if (a == b)
            continue;
          bool seen = false;
          for (int h = 0; h < nt; h++)
            seen |= (tried[h] == b);
          if (seen)
            continue;
          tried[nt++] = b;
          auto [gain, before] = evaluate_merge_gain(a, b);
          if (gain > 0) {
            ScoreRatio score{gain, std::max(1LL, before)};
            if (compare(score, bs) > 0) {
              bs = score;
              best = b;
              bg = gain;
            }
          }
        }
        if (best >= 0) {
          merge(a, best, bg);
          accepted++;
        }
      }
      std::cerr << "exact merging pass=" << pass << " before=" << start
                << " after=" << cost << " accepted=" << accepted << '\n';
      if (!accepted)
        break;
    }
  }

  void relaxed_merging(const CandidateData &candidates) {
    bool cache = use_gain_cache;
    use_gain_cache = false;
    cost = total();
    auto started = hclock::now();
    struct Entry {
      ScoreRatio score;
      int u, v;
    };
    auto less = [](const Entry &a, const Entry &b) {
      int z = compare(a.score, b.score);
      return z ? z < 0 : std::tie(a.u, a.v) < std::tie(b.u, b.v);
    };
    std::vector<Entry> heap;
    heap.reserve(size_t(n) * candidates.top_k);
    for (int u = 0; u < n; u++)
      for (int j = 0; j < candidates.top_k; j++) {
        int v = candidates.candidate_sets[u][j];
        if (v < 0)
          break;
        if (u > v &&
            std::find(candidates.candidate_sets[v].begin(), candidates.candidate_sets[v].end(), u) != candidates.candidate_sets[v].end())
          continue;
        int common_neighbor_count = candidates.common_neighbor_counts[u][j],
            adj = std::binary_search(g.vertices_[u]->nbr.begin(),
                                     g.vertices_[u]->nbr.end(), v);
        LL den = LL(g.vertices_[u]->degree) + g.vertices_[v]->degree;
        LL num = common_neighbor_count + adj;
        if (num > 0)
          heap.push_back({{num, den}, std::min(u, v), std::max(u, v)});
      }
    size_t pairs = heap.size();
    LL popped = 0, requeued = 0, accepted = 0;
    std::make_heap(heap.begin(), heap.end(), less);
    while (!heap.empty()) {
      std::pop_heap(heap.begin(), heap.end(), less);
      auto item = heap.back();
      heap.pop_back();
      popped++;
      int a = supernode_of[item.u], b = supernode_of[item.v];
      if (a == b)
        continue;
      auto [gain, before] = evaluate_merge_gain(a, b);
      LL ab = f(LL(supernode_size[a]) * supernode_size[b], get(a, b));
      if (gain + ab <= 0)
        continue;
      ScoreRatio score{gain + ab, std::max(1LL, before + ab)};
      if (compare(score, item.score) < 0) {
        item.score = score;
        heap.push_back(item);
        std::push_heap(heap.begin(), heap.end(), less);
        requeued++;
      } else {
        merge(a, b, gain);
        accepted++;
      }
    }
    std::cerr << "relaxed merging cost=" << cost << " pairs=" << pairs
              << " pops=" << popped << " requeues=" << requeued
              << " accepted=" << accepted << " seconds=" << seconds(started)
              << '\n';
    use_gain_cache = cache;
    failed.clear();
    cost = total();
  }
  int external(int x) {
    int pos = external_supernode_index[x];
    if (pos >= 0)
      return pos;
    pos = ext.size();
    external_supernode_index[x] = pos;
    ext.push_back({x, {}});
    return pos;
  }
  void clear_workspace() {
    for (auto &row : ext)
      external_supernode_index[row.id] = -1;
    ext.clear();
  }

  bool search_local_partition(std::span<const int> input_focus) {
    if (input_focus.empty())
      return false;
    if (input_focus.size() > MAX_FOCUS_VERTICES)
      throw std::runtime_error("at most three focus vertices");
    std::array<int, MAX_FOCUS_VERTICES> fb{}, ob{};
    std::copy(input_focus.begin(), input_focus.end(), fb.begin());
    std::sort(fb.begin(), fb.begin() + input_focus.size());
    int nf =
        std::unique(fb.begin(), fb.begin() + input_focus.size()) - fb.begin();
    std::span<const int> focus(fb.data(), nf);
    for (int i = 0; i < nf; i++)
      ob[i] = supernode_of[focus[i]];
    std::sort(ob.begin(), ob.begin() + nf);
    int np = std::unique(ob.begin(), ob.begin() + nf) - ob.begin();
    std::span<const int> old(ob.data(), np);
    LL nverts = 0, volume = 0;
    for (int a : old) {
      nverts += supernode_size[a];
      
      volume += degree_volume[a];
    }
    
    RepartitionAtom atoms[MAX_REPARTITION_ATOMS];
    LL atom_volume[MAX_REPARTITION_ATOMS] = {};
    int r = 0, oldmask[3] = {}, residual[3] = {-1, -1, -1};
    int fa[3], fp[3];
    for (int p = 0; p < np; p++) {
      int a = old[p], num = 0;
      LL peeled_volume = 0;
      for (int i = 0; i < nf; i++)
        if (supernode_of[focus[i]] == a) {
          fa[i] = r;
          fp[i] = p;
          oldmask[p] |= 1 << r;
          atom_volume[r] = g.vertices_[focus[i]]->degree;
          peeled_volume += atom_volume[r];
          atoms[r++] = {a, focus[i], 1};
          num++;
        }
      if ((int)supernode_size[a] > num) {
        residual[p] = r;
        oldmask[p] |= 1 << r;
        atom_volume[r] = degree_volume[a] - peeled_volume;
        atoms[r++] = {a, -1, int(supernode_size[a]) - num};
      }
    }
    if (r < 2)
      return false;
    if (r > MAX_REPARTITION_ATOMS)
      throw std::runtime_error("too many atoms");
    stats.windows++;
    clear_workspace();
    int eij[MAX_REPARTITION_ATOMS][MAX_REPARTITION_ATOMS] = {}, toparent[3][3] = {}, fedge[3][3] = {};
    auto oldpos = [&](int a) {
      for (int p = 0; p < np; p++)
        if (old[p] == a)
          return p;
      return -1;
    };
    
    for (int p = 0; p < np; p++)
      if (residual[p] >= 0)
        for (auto [x, e] : relation_edge_counts[old[p]]) {
          stats.parent_entries++;
          if (oldpos(x) < 0)
            ext[external(x)].e[residual[p]] = e;
        }
    for (int i = 0; i < nf; i++)
      focus_vertex_index[focus[i]] = i;
    
    for (int i = 0; i < nf; i++) {
      for (int v : g.vertices_[focus[i]]->nbr) {
        stats.focus_edge_scans++;
        int j = focus_vertex_index[v];
        if (j >= 0)
          fedge[i][j] = 1;
        int p = oldpos(supernode_of[v]);
        if (p >= 0)
          toparent[i][p]++;
        else {
          int x = supernode_of[v];
          auto &row = ext[external(x)];
          row.e[fa[i]]++;
          int rest = residual[fp[i]];
          if (rest >= 0)
            row.e[rest]--;
        }
      }
    }
    for (int u : focus)
      focus_vertex_index[u] = -1;
    auto setedge = [&](int a, int b, int e) {
      if (a > b)
        std::swap(a, b);
      if (e < 0)
        throw std::runtime_error("negative peeled count");
      eij[a][b] = e;
    };
    for (int i = 0; i < nf; i++)
      for (int j = i + 1; j < nf; j++)
        setedge(fa[i], fa[j], fedge[i][j]);
    for (int i = 0; i < nf; i++)
      for (int p = 0; p < np; p++)
        if (residual[p] >= 0) {
          int e = toparent[i][p];
          for (int j = 0; j < nf; j++)
            if (fp[j] == p)
              e -= fedge[i][j];
          setedge(fa[i], residual[p], e);
        }
    for (int p = 0; p < np; p++)
      if (residual[p] >= 0) {
        int e = get(old[p], old[p]);
        for (int i = 0; i < nf; i++)
          if (fp[i] == p)
            e -= toparent[i][p];
        for (int i = 0; i < nf; i++)
          if (fp[i] == p)
            for (int j = i + 1; j < nf; j++)
              if (fp[j] == p)
                e += fedge[i][j];
        setedge(residual[p], residual[p], e);
        for (int q = p + 1; q < np; q++)
          if (residual[q] >= 0) {
            int ec = get(old[p], old[q]);
            for (int i = 0; i < nf; i++) {
              if (fp[i] == p)
                ec -= toparent[i][q];
              if (fp[i] == q)
                ec -= toparent[i][p];
            }
            for (int i = 0; i < nf; i++)
              if (fp[i] == p)
                for (int j = 0; j < nf; j++)
                  if (fp[j] == q)
                    ec += fedge[i][j];
            setedge(residual[p], residual[q], ec);
          }
      }
    int N = 1 << r;
    LL sz[64] = {}, inside[64] = {}, unary[64] = {};
    for (int s = 1; s < N; s++) {
      int i = __builtin_ctz((unsigned)s), rest = s & (s - 1);
      sz[s] = sz[rest] + atoms[i].size;
      inside[s] = inside[rest] + eij[i][i];
      for (int bits = rest; bits; bits &= bits - 1) {
        int j = __builtin_ctz((unsigned)bits);
        inside[s] += eij[std::min(i, j)][std::max(i, j)];
      }
      unary[s] = f(sz[s] * (sz[s] - 1) / 2, inside[s]);
    }
    if (!allow_negative_corrections) {

      LL weights[64] = {}, number[64] = {}, edge_mass[6] = {}, mass[64] = {};
      for (const auto &row : ext) {
        LL sx = supernode_size[row.id];
        int mask = 0;
        for (int i = 0; i < r; i++) {
          edge_mass[i] += row.e[i];
          if (row.e[i] == LL(atoms[i].size) * sx)
            mask |= 1 << i;
        }
        weights[mask] += sx;
        number[mask]++;
      }
      for (int i = 0; i < r; i++)
        for (int s = 0; s < N; s++)
          if (!(s & (1 << i))) {
            weights[s] += weights[s | (1 << i)];
            number[s] += number[s | (1 << i)];
          }
      for (int s = 1; s < N; s++) {
        int i = __builtin_ctz((unsigned)s);
        mass[s] = mass[s & (s - 1)] + edge_mass[i];
        unary[s] += mass[s] - sz[s] * weights[s] + number[s];
      }
    } else {
      for (const auto &row : ext) {
        LL es[64] = {}, sx = supernode_size[row.id];
        for (int s = 1; s < N; s++) {
          int i = __builtin_ctz((unsigned)s);
          es[s] = es[s & (s - 1)] + row.e[i];
          unary[s] += f(sz[s] * sx, es[s]);
        }
      }
    }
    auto paircost = [&](int a, int b) {
      return f(sz[a] * sz[b], inside[a | b] - inside[a] - inside[b]);
    };
    LL before = 0;
    for (int i = 0; i < np; i++) {
      before += unary[oldmask[i]];
      for (int j = 0; j < i; j++)
        before += paircost(oldmask[i], oldmask[j]);
    }
    LL best = before;
    int bestmask[6], nbest = np, chosen[6];
    std::copy(oldmask, oldmask + np, bestmask);
    auto dfs = [&](auto &&self, int rem, LL acc, int depth) -> void {
      stats.dfs_states++;
      if (!rem) {
        if (acc < best) {
          best = acc;
          nbest = depth;
          std::copy(chosen, chosen + depth, bestmask);
        }
        return;
      }
      if (acc >= best)
        return;
      int bit = rem & -rem, rest = rem ^ bit;
      for (int sub = rest;; sub = (sub - 1) & rest) {
        int s = sub | bit;
        LL delta = unary[s];
        for (int i = 0; i < depth; i++)
          delta += paircost(s, chosen[i]);
        if (acc + delta < best) {
          chosen[depth] = s;
          self(self, rem ^ s, acc + delta, depth + 1);
        }
        if (!sub)
          break;
      }
    };
    dfs(dfs, N - 1, 0, 0);
    if (best >= before)
      return false;
    
    VVI newnodes(nbest);
    for (int k = 0; k < nbest; k++)
      for (int bits = bestmask[k]; bits; bits &= bits - 1) {
        int i = __builtin_ctz((unsigned)bits);
        if (atoms[i].vertex >= 0)
          newnodes[k].push_back(atoms[i].vertex);
        else
          for (int u : supernode_vertices[atoms[i].parent])
            if (!std::binary_search(focus.begin(), focus.end(), u))
              newnodes[k].push_back(u);
      }
    for (auto &row : ext) {
      int x = row.id;
      partition_revision[x]++;
      for (int a : old) {
        if (complete_savings_current) {
          int e = get(a, x);
          if (e == LL(supernode_size[a]) * supernode_size[x])
            complete_relation_savings[x] -= LL(e) - 1;
        }
        relation_edge_counts[x].erase(a);
      }
    }
    for (int a : old) {
      relation_edge_counts[a].clear();
      supernode_vertices[a].clear();
      supernode_size[a] = 0;
      internal_edge_count[a] = 0;
      degree_volume[a] = 0;
      available_supernode_ids.push_back(a);
      partition_revision[a]++;
    }
    int ids[6];
    for (int k = 0; k < nbest; k++) {
      if (available_supernode_ids.empty())
        throw std::runtime_error("no free block ID");
      int a = available_supernode_ids.back();
      available_supernode_ids.pop_back();
      ids[k] = a;
      partition_revision[a]++;
      supernode_size[a] = int(newnodes[k].size());
      internal_edge_count[a] = int(inside[bestmask[k]]);
      degree_volume[a] = 0;
      for (int bits = bestmask[k]; bits; bits &= bits - 1)
        degree_volume[a] += atom_volume[__builtin_ctz((unsigned)bits)];
      supernode_vertices[a] = std::move(newnodes[k]);
      for (int u : supernode_vertices[a])
        supernode_of[u] = a;
      if (complete_savings_current) {
        complete_relation_savings[a] = 0;
        internal_edge_cache[a] = int(inside[bestmask[k]]);
      }
      if (inside[bestmask[k]])
        relation_edge_counts[a][a] = int(inside[bestmask[k]]);
    }
    
    for (int k = 0; k < nbest; k++) {
      int a = ids[k], mask = bestmask[k];
      for (const auto &row : ext) {
        int e = 0;
        for (int bits = mask; bits; bits &= bits - 1)
          e += row.e[__builtin_ctz((unsigned)bits)];
        if (e > 0) {
          relation_edge_counts[a][row.id] = e;
          relation_edge_counts[row.id][a] = e;
          if (complete_savings_current && LL(e) == LL(supernode_size[a]) * supernode_size[row.id]) {
            complete_relation_savings[a] += LL(e) - 1;
            complete_relation_savings[row.id] += LL(e) - 1;
          }
        }
      }
      for (int j = 0; j < k; j++) {
        int other = bestmask[j];
        LL e = inside[mask | other] - inside[mask] - inside[other];
        if (e > 0) {
          relation_edge_counts[a][ids[j]] = int(e);
          relation_edge_counts[ids[j]][a] = int(e);
          if (complete_savings_current && e == LL(supernode_size[a]) * supernode_size[ids[j]]) {
            complete_relation_savings[a] += e - 1;
            complete_relation_savings[ids[j]] += e - 1;
          }
        }
      }
    }
    cost -= before - best;
    stats.commits++;
    if (n <= 31 && cost != total())
      throw std::runtime_error("incremental repartition cost mismatch");
    return true;
  }

  void relocate_vertices(const CandidateData &candidates, int rounds = 1,
                     bool natural = false) {
    if (allow_negative_corrections)
      throw std::runtime_error(
          "vertex relocation requires NCF mode");
    VVI &gamma = complete_neighbors;
    std::vector<uint64_t> &gv = complete_neighbors_revision;
    VI abtag(n, -1), abedge(n, 0);
    complete_savings_current = true;
    VI hit(n), tag(n, -1), used(n, -1), which;
    std::vector<LL> weight(n);
    VI number(n);
    VI touched;
    touched.reserve(128);
    which.reserve(256);
    int stamp = 0;
    auto saving = [](LL cap, LL e) -> LL {
      return cap > 0 && cap == e ? cap - 1 : 0;
    };
    auto ready = [&](int a) { prepare_complete_neighbors(a); };
    for (int round = 0; round < rounds; round++) {
      LL begin = cost, attempts = 0, applied = 0, candidates_count = 0;
      auto started = hclock::now();
      for (int step = 0; step < n; step++) {
        int pos = round % 2 ? n - 1 - step : step;
        int u = natural ? pos : candidates.vertex_order[pos], a = supernode_of[u], sa = supernode_size[a];
        if (stamp == INT_MAX - 1) {
          std::fill(tag.begin(), tag.end(), -1);
          std::fill(used.begin(), used.end(), -1);
          std::fill(abtag.begin(), abtag.end(), -1);
          stamp = 0;
        }
        ++stamp;
        touched.clear();
        which.clear();
        for (int v : g.vertices_[u]->nbr) {
          int x = supernode_of[v];
          if (tag[x] != stamp) {
            tag[x] = stamp;
            hit[x] = 0;
            touched.push_back(x);
          }
          hit[x]++;
        }
        auto h = [&](int x) { return tag[x] == stamp ? hit[x] : 0; };
        auto add = [&](int b, LL sx, int num) {
          if (b == a || (supernode_size[b] == 0))
            return;
          if (used[b] != stamp) {
            used[b] = stamp;
            weight[b] = 0;
            number[b] = 0;
            which.push_back(b);
          }
          weight[b] += sx;
          number[b] += num;
        };
        const int ha = h(a);
        LL ar = sa - 1;
        LL ra = saving(ar * (ar - 1) / 2, LL(internal_edge_cache[a]) - ha) -
                saving(LL(sa) * (sa - 1) / 2, internal_edge_cache[a]);
        LL pu = ar > 0 && ha == ar ? ar - 1 : 0;
        for (auto [x, e] : relation_edge_counts[a])
          if (x != a) {
            abtag[x] = stamp;
            abedge[x] = e;
            LL sx = supernode_size[x], remain = LL(e) - h(x);
            ra += saving(ar * sx, remain) - saving(LL(sa) * sx, e);
            if (ar > 0 && ha == ar && remain == ar * sx)
              add(x, ar, 1);
          }
        for (int x : touched)
          if (x != a && hit[x] == (int)supernode_size[x]) {
            LL sx = supernode_size[x];
            pu += sx - 1;
            ready(x);
            for (int b : gamma[x])
              add(b, sx, 1);
            add(x, 0, 0); 
          }
        LL best = std::max(0LL, ra + pu);
        int target = -1;
        for (int b : which) {
          candidates_count++;
          LL sb = supernode_size[b];
          int hb = h(b), ab = abtag[b] == stamp ? abedge[b] : 0,
              bb = internal_edge_cache[b];
          LL gain = ra - complete_relation_savings[b] + saving(LL(sa) * sb, ab) -
                    saving(ar * sb, LL(ab) - hb) +
                    saving(sb * (sb + 1) / 2, LL(bb) + hb) -
                    saving(sb * (sb - 1) / 2, bb) + (sb + 1) * weight[b] -
                    number[b];

          if (gain > best ||
              (gain == best && gain > 0 && target >= 0 && b < target)) {
            best = gain;
            target = b;
          }
        }

        if (best > 0) {
          LL old = cost;
          attempts++;
          bool ok = target < 0 ? search_local_partition(std::array<int, 1>{u})
                               : search_local_partition(std::array<int, 2>{u, supernode_vertices[target][0]});
          if (ok) {
            if (old - cost < best)
              throw std::runtime_error(
                  "relocation predicted gain exceeds exact local repartition gain");
            applied++;
          }
        }
      }
      stats.relocation_rounds++;
      stats.relocation_candidates += candidates_count;
      stats.relocation_attempts += attempts;
      stats.relocation_commits += applied;
      stats.relocation_gain += begin - cost;
      std::cerr << "relocation pass=" << round << " before=" << begin
                << " after=" << cost << " candidates=" << candidates_count
                << " attempts=" << attempts << " commits=" << applied
                << " seconds=" << seconds(started) << '\n';
      if (cost != total())
        throw std::runtime_error("relocation cost mismatch");
      for (int a = 0; a < n; a++)
        if (!(supernode_size[a] == 0)) {
          LL p = 0;
          for (auto [x, e] : relation_edge_counts[a])
            if (x != a && e == LL(supernode_size[a]) * supernode_size[x])
              p += LL(e) - 1;
          if (p != complete_relation_savings[a] || get(a, a) != internal_edge_cache[a])
            throw std::runtime_error("relocation cached savings mismatch");
        }
      if (!applied)
        break;
    }
    
  }
  void local_repartition_pass(const CandidateData &candidates, int pass, int triple_period) {
    LL start = cost, comm0 = stats.commits;
    int idx = 0;
    for (int step = 0; step < n; step++) {
      int u = candidates.vertex_order[(pass % 2) ? n - 1 - step : step], first = -1;
      std::array<int, 3> focus{u};
      int nf = 1;
      for (int j = 0; j < candidates.top_k; j++) {
        int v = candidates.candidate_sets[u][(j + pass) % candidates.top_k];
        if (v >= 0 && supernode_of[v] != supernode_of[u]) {
          first = v;
          break;
        }
      }
      if (first < 0)
        first = candidates.candidate_sets[u][pass % candidates.top_k];
      if (first >= 0)
        focus[nf++] = first;
      if (triple_period > 0 && (idx++ % triple_period) == 0) {
        int v = candidates.candidate_sets[u][(pass + 1) % candidates.top_k];
        if (v >= 0)
          focus[nf++] = v;
      }
      search_local_partition(std::span<const int>(focus.data(), nf));
    }
    std::cerr << "repartition pass=" << pass << " before=" << start
              << " after=" << cost << " commits=" << stats.commits - comm0
              << '\n';
    if (cost != total())
      throw std::runtime_error("whole-cost mismatch after local pass");
  }
};
struct MergingTiming {
  double relaxed_seconds = 0;
  double exact_seconds = 0;
};

MergingTiming SupernodeMerge(SupernodePartition &partition,
                             const CandidateData &candidates) {
  MergingTiming time;
  auto started = hclock::now();
  partition.relaxed_merging(candidates);
  time.relaxed_seconds = seconds(started);
  started = hclock::now();
  partition.PositiveCoarsen(candidates);
  time.exact_seconds = seconds(started);
  return time;
}

double SupernodeRepartition(SupernodePartition &partition,
                            const CandidateData &candidates) {
  auto started = hclock::now();
  for (int phase = 0; phase < 2; ++phase) {
    partition.local_repartition_pass(candidates, phase, THREE_FOCUS_INTERVAL);
    partition.PositiveCoarsen(candidates);
  }
  return seconds(started);
}

double VertexRelocation(SupernodePartition &partition,
                         const CandidateData &candidates) {
  auto started = hclock::now();
  partition.relocate_vertices(candidates, 1, true);
  return seconds(started);
}

void print_construction_stats(const ConstructionStats &s) {
  std::cout << "{\"merges\":" << s.merges
            << ",\"merge_evals\":" << s.merge_evals
            << ",\"merge_cache_hits\":" << s.merge_cache_hits
            << ",\"windows\":" << s.windows << ",\"commits\":" << s.commits
            << ",\"skipped_windows\":" << s.skipped
            << ",\"dfs_states\":" << s.dfs_states
            << ",\"focus_edge_scans\":" << s.focus_edge_scans
            << ",\"parent_entries\":" << s.parent_entries
            << ",\"relocation_rounds\":" << s.relocation_rounds
            << ",\"relocation_candidates\":" << s.relocation_candidates
            << ",\"relocation_attempts\":" << s.relocation_attempts
            << ",\"relocation_commits\":" << s.relocation_commits
            << ",\"relocation_gain\":" << s.relocation_gain << '}';
}
}

int main(int argc, char **argv) {
  try {
    if (argc != 3)
      throw std::runtime_error("usage: core graph.bin output_prefix");

    auto total_started = hclock::now();
    auto started = total_started;
    Graph graph(argv[1]);
    int vertex_count = graph.vertices_num_;
    double read_seconds = seconds(started);

    started = hclock::now();
    auto candidate_data = core::CandidateGeneration(graph, core::SamplingConfig{});
    double candidate_seconds = seconds(started);

    VI supernode_assignment(vertex_count);
    std::iota(supernode_assignment.begin(), supernode_assignment.end(), 0);
    started = hclock::now();
    core::SupernodePartition partition(graph, supernode_assignment, false, true);
    partition.fast_gain_evaluation = true;
    double state_init_seconds = seconds(started);

    auto merging_times = core::SupernodeMerge(partition, candidate_data);
    LL merge_only_cost = partition.cost;
    double repartition_seconds = core::SupernodeRepartition(partition, candidate_data);
    LL pre_relocation_cost = partition.cost;
    double relocation_seconds = core::VertexRelocation(partition, candidate_data);

    if (partition.cost != partition.total())
      throw std::runtime_error("bad final cost");
    supernode_assignment = partition.supernode_of;
    VI compact_supernode_id(vertex_count, -1);
    int next_supernode_id = 0;
    for (int vertex = 0; vertex < vertex_count; vertex++) {
      int current = supernode_assignment[vertex];
      if (compact_supernode_id[current] < 0)
        compact_supernode_id[current] = next_supernode_id++;
      supernode_assignment[vertex] = compact_supernode_id[current];
    }

    started = hclock::now();
    auto encoding = encode(graph, supernode_assignment, false, argv[2]);
    double encode_write_seconds = seconds(started);
    if (encoding.cost != partition.cost)
      throw std::runtime_error("independent encoding mismatch");

    const core::SamplingConfig config;
    std::cout << std::setprecision(12)
              << "{\"algorithm\":\"CoRe\",\"mode\":\"positive\",\"n\":"
              << vertex_count << ",\"m\":" << graph.m
              << ",\"focus_interval\":" << core::THREE_FOCUS_INTERVAL
              << ",\"top_k\":" << config.top_k
              << ",\"direct_samples\":" << config.direct_samples
              << ",\"low_degree_pivots\":" << config.low_degree_pivots
              << ",\"pivot_fanout\":" << config.pivot_fanout
              << ",\"count_mode\":" << config.count_mode << ',';
    printenc(encoding, graph.m);
    std::cout << ",\"read_seconds\":" << read_seconds
              << ",\"candidate_seconds\":" << candidate_seconds
              << ",\"state_init_seconds\":" << state_init_seconds
              << ",\"relaxed_merge_seconds\":" << merging_times.relaxed_seconds
              << ",\"exact_merge_seconds\":" << merging_times.exact_seconds
              << ",\"repartition_seconds\":" << repartition_seconds
              << ",\"relocation_seconds\":" << relocation_seconds
              << ",\"pre_relocation_cost\":" << pre_relocation_cost
              << ",\"encode_write_seconds\":" << encode_write_seconds
              << ",\"build_seconds\":"
              << candidate_seconds + state_init_seconds +
                     merging_times.relaxed_seconds + merging_times.exact_seconds +
                     repartition_seconds + relocation_seconds + encode_write_seconds
              << ",\"total_seconds\":" << seconds(total_started)
              << ",\"merge_only_cost\":" << merge_only_cost
              << ",\"stats\":";
    core::print_construction_stats(partition.stats);
    const auto &candidate_stats = candidate_data.stats;
    std::cout << ",\"count_stats\":{\"pool_vertices\":" << candidate_stats.pool_vertices
              << ",\"direct_count_vertices\":" << candidate_stats.direct_count_vertices
              << ",\"bulk_count_vertices\":" << candidate_stats.bulk_count_vertices
              << ",\"direct_count_pairs\":" << candidate_stats.direct_count_pairs
              << ",\"binary_count_pairs\":" << candidate_stats.binary_count_pairs
              << ",\"ub_pruned\":" << candidate_stats.ub_pruned
              << ",\"scan_probes\":" << candidate_stats.scan_probes
              << ",\"binary_probes\":" << candidate_stats.binary_probes
              << ",\"wedge_probes\":" << candidate_stats.wedge_probes << "}}\n";
  } catch (const std::exception &error) {
    std::cerr << "ERROR: " << error.what() << '\n';
    return 1;
  }
}

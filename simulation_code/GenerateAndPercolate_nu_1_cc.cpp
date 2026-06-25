#include <iostream>
#include <vector>
#include <array>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <random>
#include <algorithm>
#include <set>

// RNG utilities
std::mt19937 rng(std::random_device{}());
std::uniform_real_distribution<> uni(0.0, 1.0);
double rand_uniform() { return uni(rng); }

// Hash pair for unordered_set
struct PairHash {
    template <typename T1, typename T2>
    std::size_t operator() (const std::pair<T1,T2>& p) const {
        return std::hash<T1>()(p.first) ^ std::hash<T2>()(p.second);
    }
};
// And for triangles
struct Array3Hash {
    std::size_t operator()(const std::array<int, 3>& arr) const {
        std::size_t h1 = std::hash<int>()(arr[0]);
        std::size_t h2 = std::hash<int>()(arr[1]);
        std::size_t h3 = std::hash<int>()(arr[2]);
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

inline std::pair<int, int> ordered_pair(int i, int j) {
    return (i < j) ? std::make_pair(i, j) : std::make_pair(j, i);
}

// Parameters
const int N = 1000000;
const double lambda_edges = 4.5;
const double lambda_triangles = 1.0;

// Network data
std::vector<std::vector<int>> adjacency(N);
std::vector<std::array<int, 3>> triangles;

// Percolation state
std::unordered_set<std::pair<int, int>, PairHash> occupied_edges;
std::vector<bool> occupied_nodes(N, false);
std::unordered_map<int, std::vector<std::array<int, 3>>> triangles_per_node;

// --- Network Generation ---
void generate_random_graph() {
    std::poisson_distribution<> edge_dist(lambda_edges);
    std::poisson_distribution<> triangle_dist(lambda_triangles);

    std::vector<int> stubs;
    for (int i = 0; i < N; ++i) {
        int degree = edge_dist(rng);
        for (int d = 0; d < degree; ++d)
            stubs.push_back(i);
    }
    std::shuffle(stubs.begin(), stubs.end(), rng);
    for (size_t i = 0; i + 1 < stubs.size(); i += 2) {
        int u = stubs[i];
        int v = stubs[i + 1];
        if (u != v) {
            adjacency[u].push_back(v);
            adjacency[v].push_back(u);
        }
    }

    std::vector<int> t_stubs;
    for (int i = 0; i < N; ++i) {
        int num_tri = triangle_dist(rng);
        for (int d = 0; d < num_tri; ++d)
            t_stubs.push_back(i);
    }
    std::shuffle(t_stubs.begin(), t_stubs.end(), rng);
    for (size_t t = 0; t + 2 < t_stubs.size(); t += 3) {
        int i = t_stubs[t];
        int j = t_stubs[t+1];
        int k = t_stubs[t+2];
        if (i != j && i != k && j != k) {
            triangles.push_back({i, j, k});
            triangles_per_node[i].push_back({i, j, k});
            triangles_per_node[j].push_back({i, j, k});
            triangles_per_node[k].push_back({i, j, k});
        }
    }
}

// --- Dynamic Percolation ---
int dynamic_percolation(int seed_node, double T1, double T2) {
    std::queue<int> q;
    q.push(seed_node);
    occupied_nodes[seed_node] = true;
    int size = 1;

    // We track exposure for every node (id) within every triangle (ordered triplet of nodes)
    std::unordered_map<std::array<int, 3>, std::set<int>, Array3Hash> triangle_exposure;
    std::unordered_map<std::array<int, 3>, std::set<int>, Array3Hash> triangle_transmission;
    auto triangle_key = [](const std::array<int,3>& tri) {
        std::array<int, 3> sorted_tri = tri;
        std::sort(sorted_tri.begin(), sorted_tri.end());
        return std::to_string(sorted_tri[0]) + "-" + std::to_string(sorted_tri[1]) + "-" + std::to_string(sorted_tri[2]);
    };
    auto triangle_hash = [&](const std::array<int, 3>& tri) {
        return std::hash<std::string>()(triangle_key(tri));
    };

    // Until we run out of new nodes
    while (!q.empty()) {

        // Transmit around node u
        int u = q.front(); q.pop();

        // Single edges
        for (int v : adjacency[u]) {
            if (!occupied_nodes[v] && rand_uniform() < T1) {
                occupied_edges.insert(ordered_pair(u, v));
                occupied_nodes[v] = true;
                q.push(v);
                ++size;
            }
        }

        // Triangles
        for (const auto& tri : triangles_per_node[u]) {
            std::vector<std::pair<int, int>> edges = {
                {tri[0], tri[1]},
                {tri[1], tri[2]},
                {tri[2], tri[0]}
            };
            for (auto [x, y] : edges) {
                if (occupied_edges.count(ordered_pair(x, y))) continue; // already occupied
                if (x != u && y != u) continue; // only proceed for edge connected to u
                int other = (x == u) ? y : x; //who is the other node?
                double T = T1;
                if(triangle_exposure[tri].count(other) == 1) T = T2; //was the other node already exposed once?
                if (!occupied_nodes[other] && rand_uniform() < T) {
                    //if(T == T2) std::cout << "Second exposure for " << other << " in " << tri[0] <<","<< tri[1] <<","<< tri[2] << std::endl;
                    //else std::cout << "First exposure for " << other << " in " << tri[0] <<","<< tri[1] <<","<< tri[2] << std::endl;
                    occupied_edges.insert(ordered_pair(x, y));
                    occupied_nodes[other] = true;
                    q.push(other);
                    ++size;
                }
                triangle_exposure[tri].insert(other); //track exposure
            }
        }
    }
    return size;
}

// --- Main ---
int main(int argc, const char *argv[]) {

    double T1_min = atof(argv[1]);
    double T1_max = atof(argv[2]);
    double number_of_t_vals = atof(argv[3]);
    double t1_step = (T1_max - T1_min)/number_of_t_vals;
    const int number_of_runs = atof(argv[4]);

    
    for (double T1 = T1_min; T1 <=T1_max; T1 = T1 + t1_step){
    double T2 = 1.0-pow(1.0-T1,10.0);
    
    int overall_cluster = 0;
    int significantly_large = 0;
    for (size_t i = 0; i < number_of_runs; ++i) {
        if (i % 100 == 0) {
            std::fill(adjacency.begin(), adjacency.end(), std::vector<int>());
            triangles.clear();
            triangles_per_node.clear();
            generate_random_graph();
        }
        
        std::fill(occupied_nodes.begin(), occupied_nodes.end(), false);
        occupied_edges.clear();
        int seed = rng() % N;
        int cluster = dynamic_percolation(seed, T1, T2);
//            if (cluster >= 0.01*N){
//                overall_cluster += cluster;
//            }
//            if (cluster > 0.01*N){
//                ++significantly_large;
//            }
//            if (i % 1000 == 0){
//                std::cout << 100.0*i/number_of_runs << "% complete.\n";
//            }
        std::cout << T1 << ' ' << cluster << std::endl;
    }
    
    // Output
    //std::cout << "T1 = " << T1 << std::endl;
    //std::cout << "Mean Cascade Size = "<<1.0*overall_cluster/(N*number_of_runs) << std::endl;
    //std::cout << "Probability of a large cascade: " << 1.0*significantly_large/number_of_runs << std::endl;
//        if (significantly_large == 0) {
//            significantly_large = 1;
//        }
    //std::cout << T1 << ' ' << 1.0*overall_cluster/(N*significantly_large) << std::endl;
    
}
 //if((1.0*cluster)/(1.0*N) > 0.02) std::cout << 1 << std::endl;
//else std::cout << 0 << std::endl;

    return 0;
}

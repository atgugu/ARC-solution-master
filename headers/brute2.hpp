//#include <ext/pb_ds/assoc_container.hpp>
//using namespace __gnu_pbds;

double now();

#include <random>
#include <chrono>

struct State {
  vImage vimg;
  int depth;
  bool isvec;

  State() = default;
  State(vImage_&& vimg_, bool isvec_, int depth_)
    : vimg(std::move(vimg_)), isvec(isvec_), depth(depth_) {}

  ull hash() const {
    ull r = isvec ? 1 : 0;
    ull multiplier = 123413491;
    size_t i = 0;
    size_t n = vimg.size();

    // Process 4 elements per iteration when possible
    for (; i + 3 < n; i += 4) {
      r = r * multiplier + hashImage(vimg[i]);
      r = r * multiplier + hashImage(vimg[i + 1]);
      r = r * multiplier + hashImage(vimg[i + 2]);
      r = r * multiplier + hashImage(vimg[i + 3]);
      multiplier *= 31;  // Update multiplier to reduce collision risk
    }

    // Process any remaining elements
    for (; i < n; ++i) {
      r = r * multiplier + hashImage(vimg[i]);
    }

    return r;
  }
};


#include "efficient.hpp"



struct Node {
  State state;
  //vImage vimg;
  //bool isvec;
  vector<pair<int,int>> child;
  //int depth;
  int par, pfi;
  bool freed, ispiece;
  Node() {
    freed = ispiece = false;
  }
  Node(vImage_ vimg_, bool isvec_, int depth_, int par_ = -1, int pfi_ = -1) : par(par_), pfi(pfi_) {
    freed = ispiece = false;
    state.vimg = vimg_;
    state.isvec = isvec_;
    state.depth = depth_;
  }
};

struct Functions3 {
  vector<int> listed, cost;
  vector<string> names;
  vector<function<bool(const State&,State&)>> f_list;


  void add(const string& name, int cost, const function<bool(const State&,State&)>&func, int list, double prob = 1.0f);
  void add(string name, int cost, const function<Image(Image_)>&f, int list = 1, double prob = 1.0f);
  void add(string name, int cost, const function<vImage(Image_)>&f, int list = 1, double prob = 1.0f);
  void add(string name, int cost, const function<Image(vImage_)>&f, int list = 1, double prob = 1.0f);
  void add(string name, int cost, const function<vImage(vImage_)>&f, int list = 1, double prob = 1.0f);
  void add(const vector<point>&sizes, string name, int cost, const function<Image(Image_,Image_)>&f, int list = 1, double prob = 1.0f);
  string getName(int fi);
  int findfi(string name);
};



struct DAG {
  Functions3 funcs;
  //vector<Node> node;
  TinyNodeBank tiny_node;
  int givens;
  point target_size;
  //gp_hash_table<ull, int> hashi;
  TinyHashMap hashi;
  vector<int> binary;
  int add(const State&nxt, bool force = false);
  Image getImg(int nodei);
  void build();
  
  void buildBinary();
  void initial(Image_ test_in, const vector<pair<Image,Image>>&train, vector<point> sizes, int ti);
  // void benchmark();
  int applyFunc(int curi, int fi, const State&state);
  int applyFunc(int curi, int fi);
  void applyFunc(string name, bool vec);
  void applyFuncs(vector<pair<string,int>> names, bool vec);
};


struct Pieces;
vector<DAG> brutePieces2(Image_ test_in, const vector<pair<Image,Image>>&train, vector<point> out_sizes);

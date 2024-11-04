#include "precompiled_stl.hpp"
#include <chrono>
using namespace std;
#include "utils.hpp"
#include "read.hpp"
#include "normalize.hpp"
#include "core_functions.hpp"
#include "image_functions.hpp"
#include "image_functions2.hpp"

#include "visu.hpp"

#include "brute2.hpp"
#include "pieces.hpp"

#include "timer.hpp"
#include <unordered_set>

extern int MAXDEPTH, print_nodes;

//double build_f_time = 0, apply_f_time = 0;
//double real_f_time = 0;

double now() {
  ll t = chrono::high_resolution_clock::now().time_since_epoch().count();
  static ll time0 = 0;
  if (time0 == 0) time0 = t;
  return (t-time0)*1e-9;
}
//double now() { return chrono::steady_clock::now().time_since_epoch().count()*1e-9;}

// Timer build_f_time, apply_f_time, real_f_time, add_time, find_child_time, add_child_time, hash_time, map_time, total_time;
// Timer state_time;

void Functions3::add(const string& name, int cost_, const function<bool(const State&,State&)>&func, int list) {
  //if (cost_ != 10) cout << name << endl;
  //assert(cost_ == 10);
  if (list) listed.push_back(names.size());
  names.push_back(name);
  f_list.push_back(func);
  cost.push_back(cost_);
}

void Functions3::add(string name, int cost, const function<Image(Image_)>&f, int list) { //list = 1
  auto func = [f](const State& cur, State& nxt) {

    //if (cur.isvec) return false;

    nxt.vimg.resize(cur.vimg.size());
    nxt.isvec = cur.isvec;

    int area = 0;
    for (int i = 0; i < cur.vimg.size(); ++i) {
      // real_f_time.start();
      nxt.vimg[i] = f(cur.vimg[i]);
      // real_f_time.stop();

      area += nxt.vimg[i].w*nxt.vimg[i].h;
      if (area > MAXPIXELS) return false;
    }
    return true;
  };
  add(name, cost, func, list);
}

void Functions3::add(string name, int cost, const function<vImage(Image_)>&f, int list) { //list = 1
  const int buffer = 5;
  auto func = [f,cost](const State& cur, State& nxt) {
    if (cur.isvec || cur.depth+cost+buffer > MAXDEPTH) return false;
    // real_f_time.start();
    nxt.vimg = f(cur.vimg[0]);
    // real_f_time.stop();
    nxt.isvec = true;
    return true;
  };
  add(name, cost, func, list);
}

void Functions3::add(string name, int cost, const function<Image(vImage_)>&f, int list) { //list = 1
  auto func = [f](const State& cur, State& nxt) {
    if (!cur.isvec) return false;
    nxt.vimg.resize(1);
    // real_f_time.start();
    nxt.vimg[0] = f(cur.vimg);
    // real_f_time.stop();
    nxt.isvec = false;
    return true;
  };
  add(name, cost, func, list);
}

void Functions3::add(string name, int cost, const function<vImage(vImage_)>&f, int list) { //list = 1
  auto func = [f](const State& cur, State& nxt) {
    if (!cur.isvec) return false;
    // real_f_time.start();
    nxt.vimg = f(cur.vimg);
    // real_f_time.stop();
    nxt.isvec = true;
    return true;
  };
  add(name, cost, func, list);
}

void Functions3::add(const vector<point>&sizes, string name, int cost, const function<Image(Image_,Image_)>&f, int list) { //list = 1
  int szi = 0;
  for (point sz : sizes) {
    Image arg2 = core::empty(sz);
    auto func = [f,arg2](const State& cur, State& nxt) {

      if (cur.isvec) return false;

      nxt.vimg.resize(cur.vimg.size());

      int area = 0;
      for (int i = 0; i < cur.vimg.size(); ++i) {
	// real_f_time.start();
	nxt.vimg[i] = f(cur.vimg[i], arg2);
	// real_f_time.stop();

	area += nxt.vimg[i].w*nxt.vimg[i].h;
	if (area > MAXPIXELS) return false;
      }
      nxt.isvec = cur.isvec;
      return true;
    };
    add(name+" "+to_string(szi++), cost, func, list);
  }
}

string Functions3::getName(int fi) {
  assert(fi >= 0 && fi < names.size());
  return names[fi];
}
int Functions3::findfi(string name) {
  int fi = find(names.begin(), names.end(), name)-names.begin();
  if (fi == names.size()) {
    cerr << name << " is not a known function" << endl;
    assert(0);
  }
  return fi;
}


Functions3 initFuncs3(const vector<point>&sizes, const std::unordered_map<int, int> &colorMap) {
  Functions3 funcs;

  // Unary

  //invert is filterCol(img, 0)
  for (int c = 0; c < 10;++c)
    funcs.add("filterCol "+to_string(c), 6, [c](Image_ img) {return filterCol(img, c);});
  for (int c = 1; c < 10;++c)
    funcs.add("eraseCol "+to_string(c), 6,
	      [c](Image_ img) {return eraseCol(img, c);});

  for (int c = 1; c < 10;++c)
    funcs.add("colShape "+to_string(c), 7,
	      [c](Image_ img) {return colShape(img, c);}, 0);

  funcs.add("compress", 8, [](Image_ img) {return compress(img);});
  funcs.add("getPos", 8, getPos);
  funcs.add("getSize0", 8, getSize0);
  funcs.add("getSize", 8, getSize);
  funcs.add("hull0", 10, hull0);
  funcs.add("hull", 8, hull);
  funcs.add("toOrigin", 10, toOrigin);
  funcs.add("Fill", 9, Fill);
  funcs.add("interior", 8, interior);
  funcs.add("interior2", 8, interior2);
  funcs.add("border", 9, border);
  funcs.add("center", 7, center);
  funcs.add("majCol", 8, majCol);

  funcs.add("greedyFillBlack", 10, [](Image_ img) {return greedyFillBlack(img);});
  funcs.add("greedyFillBlack2", 10, [](Image_ img) {return greedyFillBlack2(img);});
  funcs.add("highlightEdges", 10, [colorMap](Image_ img) {return highlightEdges(img, colorMap);});
  funcs.add("maskByColorMap", 10, [colorMap](Image_ img) {return maskByColorMap(img, colorMap);});
  for (const auto &c : colorMap) {
      int backgroundColor = c.first;
      funcs.add("replaceBackground" + std::to_string(backgroundColor), 10, 
                [colorMap, backgroundColor](Image_ img) {
                    return replaceBackground(img, backgroundColor, colorMap);
                });}

  for (int i = 1; i < 11; ++i)
    funcs.add("rigid "+to_string(i), 8,
	      [i](Image_ img) {return rigid(img, i);});
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b){
      int comp = 10;
      if (b==0) comp = 2;
      funcs.add("count "+to_string(a)+" "+to_string(b), 10,
		[a,b](Image_ img) {return count(img, a, b);});
    }
  for (int i = 0; i < 15; ++i)
    funcs.add("smear "+to_string(i), 9,
	      [i](Image_ img) {return smear(img, i);});

  for (int i = 0; i < 4; ++i)
  funcs.add("diagonalSmear "+to_string(i), 9,
    [i](Image_ img) {return diagonalSmear(img, i);});
  // for (int i = 0; i < 8; ++i)
  // for (int j = 0; j < 1; ++j)
  // funcs.add("shiftImage "+to_string(i) + to_string(j), 5,
  //   [i, j](Image_ img) {return shiftImage(img, i, j, false);});

  funcs.add("dilate1", 8, dilate1);
  funcs.add("erode1", 8, erode1);
  funcs.add("dilate2", 8, dilate2);
  funcs.add("erode2", 8, erode2);
  funcs.add("dilate3", 9, dilate3);
  funcs.add("erode3", 9, erode3);

  // funcs.add("morphOpening", 8, morphOpening);
  // funcs.add("morphClosing", 8, morphClosing);
  // for (int i = 2; i < 10; ++i)
  // funcs.add("detectPatterns "+to_string(i), 9,
  //   [i](Image_ img) {return detectPatterns(img, i);});

  funcs.add("removeGrid", 10, removeGrid);
  funcs.add("detectVerticalStripes2", 10, [](Image_ img) { return detectVerticalStripes(img, 2); });
  funcs.add("detectVerticalStripes3", 10, [](Image_ img) { return detectVerticalStripes(img, 3); });
  funcs.add("detectVerticalStripes4", 10, [](Image_ img) { return detectVerticalStripes(img, 4); });
  funcs.add("detectVerticalStripes5", 10, [](Image_ img) { return detectVerticalStripes(img, 5); });
  funcs.add("detectDiagonalPattern2", 10, [](Image_ img) { return detectDiagonalPattern(img, 2); });
  funcs.add("detectDiagonalPattern3", 10, [](Image_ img) { return detectDiagonalPattern(img, 3); });
  funcs.add("detectHorizontalStripes2", 10, [](Image_ img) { return detectHorizontalStripes(img, 2); });
  funcs.add("detectHorizontalStripes3", 10, [](Image_ img) { return detectHorizontalStripes(img, 3); });
  funcs.add("detectCrossPattern3", 10, [](Image_ img) { return detectCrossPattern(img, 3); });
  funcs.add("detectCrossPattern5", 10, [](Image_ img) { return detectCrossPattern(img, 5); });

  funcs.add("detectCheckerboardPattern2", 10, [](Image_ img) { return detectCheckerboardPattern(img, 2); });
  funcs.add("detectCheckerboardPattern3", 10, [](Image_ img) { return detectCheckerboardPattern(img, 3); });
  funcs.add("detectCheckerboardPattern5", 10, [](Image_ img) { return detectCheckerboardPattern(img, 5); });
  // funcs.add("trainAndPredict", 10, [train](Image_ img ) { return trainAndPredict(img, train); });
  // funcs.add("extractAndApplyConstantShapes", 10, [train](Image_ img ) { return extractAndApplyConstantShapes(img, train); });

  funcs.add("detectRepeatingPattern", 10, [](Image_ img) {
      Image pattern;
      int offsetX, offsetY;
      bool hasPattern = detectRepeatingPattern(img, pattern, offsetX, offsetY);
      if (hasPattern) {
          // Return the pattern image
          return pattern;
      } else {
          // Return an empty image
          return Image{{0, 0}, {0, 0}, {}};
      }
  });

  funcs.add("applyColorMapping", 10, [colorMap](Image_ img) {
      return applyColorMapping(img, colorMap);
  });

  funcs.add("trimToContent", 10, trimToContent);
  funcs.add("detectRepeatingPatternWithHole1", 10,
	      [](Image_ img) {return detectRepeatingPatternWithHole(img, true);});
  funcs.add("detectRepeatingPatternWithHole1", 10,
	      [](Image_ img) {return detectRepeatingPatternWithHole(img, false);});

  // funcs.add("detectTranslation1DPattern", 10, detectTranslation1DPattern);
  // funcs.add("detectTranslationPattern", 10, detectTranslationPattern);
// funcs.add("enforceRotationalSymmetry90", 10, enforceRotationalSymmetry90);
// funcs.add("enforceRotationalSymmetry180", 10, enforceRotationalSymmetry180);

  // for (int id = 0; id < 4; ++id)
  //   funcs.add("gridFilter "+to_string(id), 10,
	//       [id](Image_ img) {return gridFilter(img, id, id);});

  funcs.add("mostCommonShape", 10, mostCommonShape);
  funcs.add("largestShape", 10, largestShape);
  funcs.add("smallestShape", 10, smallestShape);
  funcs.add("mostCommonColorShape", 10, mostCommonColorShape);
  funcs.add("enclosedShapes", 10, enclosedShapes);
  funcs.add("symmetricShape", 10, symmetricShape);

  // funcs.add("repairRotationalSymmetry", 10, repairRotationalSymmetry);

  // crashes
  // for (int id = 2; id <3; ++id)
  //   funcs.add("compressPatches "+to_string(id), 10,
	//       [id](Image_ img) {return compressPatches(img,id);});

  funcs.add("compress2", 10, compress2);
  funcs.add("compress3", 10, compress3);
  funcs.add("compressH", 10, [](Image_ img) { return compressHV(img, true); });
  funcs.add("compressV", 10, [](Image_ img) { return compressHV(img, false); });
  funcs.add("outlineShapes", 10, outlineShapes);
  funcs.add("diagonalBridge", 10, diagonalBridge);
  funcs.add("fillIslands", 10, fillIslands);
  funcs.add("reverseSizes", 10, reverseSizes);
  funcs.add("compressSymmetry", 10, compressSymmetry);
  funcs.add("bridgeGaps", 10, bridgeGaps);
  funcs.add("connectNearestShapes", 10, connectNearestShapes);
  funcs.add("connectFarthestShapes", 10, connectFarthestShapes);

  for (int id = 3; id < 5; ++id)
    funcs.add("denseRegionShape "+to_string(id), 10,
	      [id](Image_ img) {return denseRegionShape(img,id);});

  for (int id = 0; id < 5; ++id)
    funcs.add("connect "+to_string(id), 10,
	      [id](Image_ img) {return connect(img,id);});
  for (int id = 1; id < 4; ++id)
    funcs.add("removeNoise "+to_string(id), 10,
	      [id](Image_ img) {return removeNoise(img,id);});

  for (int id = 0; id < 6; ++id)
    funcs.add("rearrangeShapes "+to_string(id), 8,
	      [id](Image_ img) {return rearrangeShapes(img,id);});

  for (int id : {0,1})
    funcs.add("spreadCols "+to_string(id), 10,
	      [id](Image_ img) {return spreadCols(img, id);});

  for (int id = 0; id < 4; ++id)
    funcs.add("half "+to_string(id), 8,
	      [id](Image_ img) {return half(img, id);});

  for (int id = 1; id < 5; ++id)
    funcs.add("zoomIn "+to_string(id), 8,
	      [id](Image_ img) {return zoomIn(img, id);});

  for (int id = 2; id <= 5; ++id)
    funcs.add("upscaleImage "+to_string(id), 8,
	      [id](Image_ img) {return upscaleImage(img, id);});

  for (int id = 2; id <= 3; ++id)
    funcs.add("downscaleImage "+to_string(id), 10,
	      [id](Image_ img) {return downscaleImage(img, id);});

  for (int id = 2; id <= 4; ++id)
    funcs.add("stretchImageH "+to_string(id), 8,
	      [id](Image_ img) {return stretchImage(img, id, 0);});

  for (int id = 2; id <= 4; ++id)
    funcs.add("stretchImageV "+to_string(id), 8,
	      [id](Image_ img) {return stretchImage(img, id, 1);});

  // for (int id = 0; id <= 4; ++id)
  //   funcs.add("shuffleRows "+to_string(id), 8,
	//       [id](Image_ img) {return shuffleRowsOrColumns(img,true, id);});

  // for (int id = 0; id <= 4; ++id)
  //   funcs.add("shuffleColumns "+to_string(id), 8,
	//       [id](Image_ img) {return shuffleRowsOrColumns(img,false, id);});

    funcs.add("makeBorder", 8,
	    [](Image_ img) {return makeBorder(img, 1);});

  for (int id : {0,1})
    funcs.add("makeBorder2 "+to_string(id), 8,
	      [id](Image_ img) {return makeBorder2(img, id);});

  for (int id = -3; id < 3; ++id)
    funcs.add("circularShiftColumn "+to_string(id), 10,
	      [id](Image_ img) {return circularShift(img, false, id);});

    for (int id = -3; id < 3; ++id)
    funcs.add("circularShiftRow "+to_string(id), 10,
	      [id](Image_ img) {return circularShift(img, true, id);});

  // for (int id = 0; id < 4; ++id)
  //   funcs.add("diagonalGravity "+to_string(id), 10,
	//       [id](Image_ img) {return diagonalGravity(img,id);});
  // for (int ammount = 0; ammount < 3; ++ammount)
  // for (int id = 0; id < 8; ++id)
  //   funcs.add("shiftImageCrop "+to_string(id), 5,
	//       [id, ammount](Image_ img) {return shiftImageCrop(img, id, ammount);});

  // // for (int ammount = 0; ammount < 3; ++ammount)
  // // for (int id = 0; id < 8; ++id)
  // //   funcs.add("shiftImageExpand "+to_string(id), 5,
  // //       [id, ammount](Image_ img) {return shiftImageExpand(img, id, ammount);});


  for (int dy = -3; dy <= 3; ++dy) {
    for (int dx = -3; dx <= 3; ++dx) {
      funcs.add("Move "+to_string(dx)+" "+to_string(dy), 8,
		[dx,dy](Image_ img) {return Move(img, Pos(dx,dy));}, 0);
    }
  }

  // Binary
  funcs.add(sizes, "embed", 10, embed);
  funcs.add(sizes, "wrap", 10, wrap);
  funcs.add(sizes, "broadcast", 10, [](Image_ a, Image_ b) {return broadcast(a,b);});
  funcs.add(sizes, "repeat 0",  10, [](Image_ a, Image_ b) {return repeat(a,b);});
  funcs.add(sizes, "repeat 1",  10, [](Image_ a, Image_ b) {return repeat(a,b,1);});
  funcs.add(sizes, "repeat 2",  10, [](Image_ a, Image_ b) {return repeat(a,b,2);});
  funcs.add(sizes, "repeat 3",  10, [](Image_ a, Image_ b) {return repeat(a,b,3);});
  funcs.add(sizes, "repeat 4",  10, [](Image_ a, Image_ b) {return repeat(a,b,4);});
  funcs.add(sizes, "extend2",  10, [](Image_ a, Image_ b) {return extend2(a,b);});

  // funcs.add(sizes, "repeat 5",  10, [](Image_ a, Image_ b) {return repeat(a,b,5);});
  // funcs.add(sizes, "repeat 6",  10, [](Image_ a, Image_ b) {return repeat(a,b,6);});
  // funcs.add(sizes, "repeat 7",  10, [](Image_ a, Image_ b) {return repeat(a,b,7);});
  // funcs.add(sizes, "repeat 8",  10, [](Image_ a, Image_ b) {return repeat(a,b,8);});
  // funcs.add(sizes, "repeat 9",  10, [](Image_ a, Image_ b) {return repeat(a,b,9);});

  funcs.add(sizes, "mirror 0",  10, [](Image_ a, Image_ b) {return mirror(a,b);});
  funcs.add(sizes, "mirror 1",  10, [](Image_ a, Image_ b) {return mirror(a,b,1);});
  funcs.add(sizes, "mirror 2",  10, [](Image_ a, Image_ b) {return mirror(a,b,2);});
  funcs.add(sizes, "mirror 3",  10, [](Image_ a, Image_ b) {return mirror(a,b,3);});
  funcs.add(sizes, "mirror 4",  10, [](Image_ a, Image_ b) {return mirror(a,b,4);});
  // funcs.add(sizes, "mirror 5",  10, [](Image_ a, Image_ b) {return mirror(a,b,5);});
  // funcs.add(sizes, "mirror 6",  10, [](Image_ a, Image_ b) {return mirror(a,b,6);});
  // funcs.add(sizes, "mirror 7",  10, [](Image_ a, Image_ b) {return mirror(a,b,7);});
  // funcs.add(sizes, "mirror 8",  10, [](Image_ a, Image_ b) {return mirror(a,b,8);});
  // funcs.add(sizes, "mirror 9",  10, [](Image_ a, Image_ b) {return mirror(a,b,9);});
  funcs.add(sizes, "ringSmear",  10, [](Image_ a, Image_ b) {return ringSmear(a,b);});
  // crash
  // funcs.add(sizes, "diagonalSmear1",  10, [](Image_ a, Image_ b) {return diagonalSmear(a,b,1);});
  // funcs.add(sizes, "diagonalSmear2",  10, [](Image_ a, Image_ b) {return diagonalSmear(a,b,2);});
  // funcs.add(sizes, "diagonalSmear3",  10, [](Image_ a, Image_ b) {return diagonalSmear(a,b,3);});
  // funcs.add(sizes, "diagonalSmear4",  10, [](Image_ a, Image_ b) {return diagonalSmear(a,b,4);});


  //Split
  funcs.add("cut",       8, [](Image_ img) {return cut(img);});
  funcs.add("splitCols", 7, [](Image_ img) {return splitCols(img);});
  funcs.add("splitAll",     8, splitAll);
  funcs.add("splitColumns", 8, splitColumns);
  funcs.add("splitRows",    8, splitRows);
  funcs.add("insideMarked", 10, insideMarked);
  for (int id = 0; id < 4; ++id)
    funcs.add("gravity "+to_string(id), 10,
	      [id](Image_ img) {return gravity(img,id);});

  //Join
  for (int id = 0; id < 14; ++id)
    funcs.add("pickMax "+to_string(id), 8,
	      [id](vImage_ v) {return pickMax(v,id);});
  for (int id = 0; id < 1; ++id)
    funcs.add("pickUnique "+to_string(id), 8,
	      [id](vImage_ v) {return pickUnique(v,id);});

  funcs.add("composeGrowing", 10, composeGrowing);
  funcs.add("stackLine", 7, stackLine);
  for (int id = 0; id < 3; ++id) //consider going to 4
    funcs.add("myStack "+to_string(id), 10,
	      [id](vImage_ v) {return myStack(v,id);}); //


  //Vector
  for (int id = 0; id < 14; ++id)
    funcs.add("pickMaxes "+to_string(id), 8,
	      [id](vImage_ v) {return pickMaxes(v,id);});
  for (int id = 0; id < 14; ++id)
    funcs.add("pickNotMaxes "+to_string(id), 8,
	      [id](vImage_ v) {return pickNotMaxes(v,id);});

  //funcs.add("smear",    [](Image_ a, Image_ b) {return smear(a,b,6);});

  //funcs.add(insideMarked); //only do once at depth 0

  // Smear diagonals?

  // outerProducts

  //for (int id = 0; id < 4; id++)
  //  funcs.add([id](Image_ img) {return gravity(img, id);});

  // Image makeBorder(Image_ img, int bcol = 1);
  // Image makeBorder2(Image_ img, Image_ bord);
  // Image greedyFillBlack(Image_ img, int N = 3);
  // Image extend2(Image_ img, Image_ room);
  // Image replaceTemplate(Image_ in, Image_ need_, Image_ marked_, int overlapping = 0, int rigids = 0);
  // Image swapTemplate(Image_ in, Image_ a, Image_ b, int rigids = 0);

  // funcs.add("heuristicCut", heuristicCut);

  return funcs;
}

Image DAG::getImg(int ni) {
  return tiny_node.getImg(ni);
  //assert(tiny_node.getImg(ni) == node[ni].vimg[0]);
  //return node[ni].state.vimg[0];
  /*
  //cout << nodei << endl;
  assert(nodei >= 0 && nodei < (int)node.size());
  assert(!node[nodei].isvec);
  if (node[nodei].pfi == embed1fi) {
    assert(funcs.f_list[embed1fi](node[nodei], tmp_node));
    //assert(tmp_node.vimg == node[nodei].vimg);
    return tmp_node.vimg[0];
  } else {
    assert(node[nodei].vimg.size());
    return node[nodei].vimg[0];
    }*/
}


int DAG::add(const State&nxt, bool force) { //force = false
  // hash_time.start();
  ull h = nxt.hash();
  // hash_time.stop();
  const int nodes = tiny_node.size();
  // map_time.start();
  auto [nodei,inserted] = hashi.insert(h,nodes);
  //auto [it,inserted] = hashi.insert({h,nodes}); int nodei = it->second;
  //assert(inserted == tiny_inserted);
  //assert(nodei == tiny_nodei);

  // map_time.stop();

  // add_time.start();
  if (inserted || force) {
    bool ispiece = !nxt.isvec;
    if (!nxt.isvec && target_size != point{-1,-1})
      ispiece &= (nxt.vimg[0].p == point{0,0} && nxt.vimg[0].sz == target_size);
    /*{
      Node n;
      n.state = nxt;
      n.ispiece = ispiece;
      n.freed = false;
      node.push_back(n);
      }*/
    tiny_node.append(nxt, ispiece);
  }
  // add_time.stop();
  return nodei;
}


void DAG::build() {
  // build_f_time.start();

  for (int curi = 0; curi < tiny_node.size(); ++curi) {
    const int depth = tiny_node[curi].depth;
    if (depth+1 > MAXDEPTH) continue;

    //vector<pair<int,int>> child;
    State nxt;
    // state_time.start();
    State cur_state = tiny_node.getState(curi);
    // state_time.stop();
    for (int fi : funcs.listed) {
      // here modify the cost based on depth
      nxt.depth = depth+funcs.cost[fi];
      if (nxt.depth > MAXDEPTH) continue;
      if (funcs.f_list[fi](cur_state, nxt)) {
	int newi = add(nxt);
	//child.emplace_back(fi, newi);
	tiny_node.addChild(curi, fi, newi);
      } else {
	tiny_node.addChild(curi, fi, -1);
	//child.emplace_back(fi, -1);
      }

    }
    //node[curi].child = child;
  }

  // build_f_time.stop();
}

void DAG::initial(Image_ test_in, const vector<pair<Image,Image>>&train, vector<point> sizes, int ti) {
  if (sizes.size() > 1)
    target_size = sizes[1];
  else
    target_size = point{-1,-1};

  Image in = ti < train.size() ? train[ti].first : test_in;

  add(State({in}, false, 0), true);

  //Output sizes
  for (point sz : sizes)
    add(State({core::empty(sz)}, false, 10), true);

  // Outputs of other trains
  for (int tj = 0; tj < train.size(); ++tj)
    add(State({ti != tj ? train[tj].second : core::empty(train[tj].second.sz)}, false, 10), true);

  //add(State({greedyFillBlack2(in)}, false, 10), true);

  //filterCol?

  givens = tiny_node.size();
}

//time each function
// void DAG::benchmark() {
//   vector<pair<double,int>> v;
//   const unsigned int tinynodesize =tiny_node.size();
//   for (int fi : funcs.listed) {
//     double start_time = now();
//     State nxt;
//     for (int i = 0; i < tinynodesize; ++i) {
//       funcs.f_list[fi](tiny_node.getState(i), nxt);
//     }
//     double elapsed = now()-start_time;
//     v.emplace_back(elapsed, fi);
//   }
//   sort(v.begin(), v.end());
//   for (auto [t,fi] : v)
//     printf("%.1f ms - %s\n", t*1e3, funcs.getName(fi).c_str());
// }

int DAG::applyFunc(int curi, int fi, const State&state) {

  // find_child_time.start();
  //auto it = lower_bound(node[curi].child.begin(), node[curi].child.end(), make_pair(fi,-1));
  int it2 = tiny_node.getChild(curi, fi);
  // find_child_time.stop();
  if (it2 != TinyChildren::None) {//it != node[curi].child.end() && it->first == fi) {
    //if (it2 != it->second) cout << it2 << ' ' << it->second << endl;
    //assert(it2 == it->second);
    return it2;//it->second;
  }
  //assert(it2 == -2);

  State nxt;
  nxt.depth = tiny_node[curi].depth+funcs.cost[fi];
  //nxt.par = curi;

  int newi = -1;
  bool ok = false;
  // apply_f_time.start();
  ok = funcs.f_list[fi](state, nxt);
  // apply_f_time.stop();

  if (ok) {
    //nxt.pfi = fi;
    newi = add(nxt);
  }

  // add_child_time.start();
  tiny_node.addChild(curi, fi, newi);
  //node[curi].child.emplace_back(fi, newi);
  //sort(node[curi].child.begin(), node[curi].child.end());
  // add_child_time.stop();
  return newi;
}

int DAG::applyFunc(int curi, int fi) {
  // state_time.start();
  State state = tiny_node.getState(curi);
  // state_time.stop();

  return applyFunc(curi, fi, state);
}

void DAG::applyFunc(string name, bool vec) {
  int fi = funcs.findfi(name);

  const int start_nodes = tiny_node.size();

  for (int curi = 0; curi < start_nodes; ++curi) {
    if (tiny_node[curi].isvec == vec) applyFunc(curi, fi);
  }
}

void DAG::applyFuncs(vector<pair<string,int>> names, bool vec) {
  vector<pair<int,int>> fis;
  for (auto [name,id] : names) {
    fis.emplace_back(funcs.findfi(name), id);
  }

  vector<int> parid(tiny_node.size(),-1);
  const unsigned int tinynodesize = tiny_node.size();
  const unsigned int fissize = fis.size();
  for (int curi = 0; curi < tinynodesize; ++curi) {
    if (tiny_node[curi].isvec != vec) continue;
    // state_time.start();
    State state = tiny_node.getState(curi);
    // state_time.stop();

    for (int i = 0; i < fissize; ++i) {
      auto [fi,id] = fis[i];
      if (id <= parid[curi]) continue;
      if (applyFunc(curi, fi, state) == parid.size()) {
        parid.push_back(id);
        assert(parid.size() == tiny_node.size());
      }
    }
  }
}

void DAG::buildBinary() {
  int fis = *max_element(funcs.listed.begin(), funcs.listed.end())+1;
  binary.assign(fis*fis, -1);
  vector<State> state(fis);
  vector<int> active(fis), memi(fis);
  for (int fi : funcs.listed) {
    int curi = tiny_node.getChild(0, fi);
    if (curi >= 0) {
      active[fi] = 1;
      state[fi] = tiny_node.getState(curi);
      memi[fi] = curi;
    }
  }
  //ham
  //#pragma omp parallel for
  for (int fa : funcs.listed) {
    if (!active[fa]) continue;
    for (int fb : funcs.listed) {
      if (!active[fb]) continue;

      if (state[fa].isvec || state[fb].isvec) continue;
      State nxt;
      nxt.vimg = {align(state[fa].vimg[0], state[fb].vimg[0])};
      nxt.depth = 2; //TODO
      nxt.isvec = false;
      binary[fa*fis+fb] = add(nxt);
      /*if (binary[fa*fis+fb] != memi[fa]) {
	cout << binary[fa*fis+fb] << ' ' << memi[fa] << ' ' << tiny_node.size() << endl;
      }
      assert(binary[fa*fis+fb] == memi[fa]);*/
    }
  }
}

std::vector<int> getAllColors(const Image_ &test_in, const std::vector<std::pair<Image, Image>> &train) {
    std::unordered_set<int> colorSet;

    // Collect colors from the test_in image
    for (int i = 0; i < test_in.h; ++i) {
        for (int j = 0; j < test_in.w; ++j) {
            int color = test_in(i, j);
            if (color >= 0 && color <= 9) {
                colorSet.insert(color);
            }
        }
    }

    // Collect colors from each image in the train vector
    for (const auto &[input, output] : train) {
        for (int i = 0; i < input.h; ++i) {
            for (int j = 0; j < input.w; ++j) {
                int color = input(i, j);
                if (color >= 0 && color <= 9) {
                    colorSet.insert(color);
                }
            }
        }
        for (int i = 0; i < output.h; ++i) {
            for (int j = 0; j < output.w; ++j) {
                int color = output(i, j);
                if (color >= 0 && color <= 9) {
                    colorSet.insert(color);
                }
            }
        }
    }

    // Convert the set of colors to a sorted vector
    std::vector<int> colors(colorSet.begin(), colorSet.end());
    std::sort(colors.begin(), colors.end());
    return colors;
}

std::unordered_map<int, int> getColorMapping(const std::vector<std::pair<Image, Image>> &train) {
    std::unordered_map<int, int> colorMap;

    // Iterate over each input-output pair in the train vector
    for (const auto &[input, output] : train) {
        for (int i = 0; i < input.h; ++i) {
            for (int j = 0; j < input.w; ++j) {
                int inputColor = input(i, j);
                int outputColor = output(i, j);

                // Only map valid colors between 0 and 9
                if (inputColor >= 0 && inputColor <= 9 && outputColor >= 0 && outputColor <= 9) {
                    // If this color already has a mapping, ensure it is consistent
                    if (colorMap.count(inputColor)) {
                        // If an inconsistent mapping is found, this could indicate an error
                        if (colorMap[inputColor] != outputColor) {
                            // Handle inconsistency (this can be logged, replaced, or ignored)
                            // For now, we keep the original mapping
                        }
                    } else {
                        // If the color is not mapped yet, add the mapping
                        colorMap[inputColor] = outputColor;
                    }
                }
            }
        }
    }

    return colorMap;
}
vector<DAG> brutePieces2(Image_ test_in, const vector<pair<Image,Image>>&train, vector<point> out_sizes) {
  const int print = 0;
  const size_t trainsize = train.size();
  vector<DAG> dag(trainsize+1);

  int all_train_out_mask = 0, and_train_out_mask = ~0;
  for (int ti = 0; ti < trainsize; ++ti)
    and_train_out_mask &= core::colMask(train[ti].second);

  for (int ti = 0; ti <= trainsize; ++ti) {
    vector<point> sizes;
    if (ti < trainsize)
      sizes.push_back(train[ti].first.sz);
    else
      sizes.push_back(test_in.sz);

    if (out_sizes.size())
      sizes.push_back(out_sizes[ti]);
    std::unordered_map<int, int> colors = getColorMapping(train);
    dag[ti].funcs = initFuncs3(sizes, colors);

    dag[ti].initial(test_in, train, sizes, ti);

    // total_time.start();

    // double start_time = now();
    dag[ti].build();
    // if (print) cout << now()-start_time << endl;
    //dag[ti].buildBinary();
    //if (print) cout << now()-start_time << endl;
    dag[ti].applyFunc("composeGrowing", 1);
    // if (print) cout << now()-start_time << endl;

    if (sizes.size() > 1) {
      vector<pair<string,int>> toapply;
      toapply.emplace_back("toOrigin",0);
      for (int c = 1; c <= 5;++c)
	if (and_train_out_mask>>c&1)
	  toapply.emplace_back("colShape "+to_string(c),1);
      toapply.emplace_back("embed 1",2);
      dag[ti].applyFuncs(toapply, 0);
      /*
      dag[ti].applyFunc("toOrigin", 0);
      if (print) cout << now()-start_time << endl;

      int mask;
      if (ti < train.size()) {
	mask = core::colMask(train[ti].second);
	all_train_out_mask |= mask;
      } else {
	mask = all_train_out_mask;
      }

      for (int c = 1; c <= 5;++c)
	if (and_train_out_mask>>c&1)
	dag[ti].applyFunc("colShape "+to_string(c), 0);
      if (print) cout << now()-start_time << endl;

      if (ti < train.size())
	deducePositions(dag[ti], train[ti].second);
      if (print) cout << now()-start_time << endl;

      dag[ti].applyFunc("embed 1", 0);
      */
      // if (print) cout << now()-start_time << endl;

      /*if (ti < train.size())
	deducePositions(dag[ti], train[ti].second);

	if (print) cout << now()-start_time << endl;*/

  //     // total_time.stop();
  //     if(print){
  //       // total_time.print("Total time");
  //       // build_f_time.print("Build f time");
  //       // apply_f_time.print("Apply f time");
  //       // real_f_time .print("Real f time ");

  //       // add_time.print("Add time");
  //       // find_child_time.print("Find child");
  //       // add_child_time.print("Add child");
  //       // hash_time.print("Hash");
  //       // map_time.print("Map");

  //       // state_time.print("getState");
  //     }
  //     //exit(0);
  //     /*FILE*fp = fopen("images.txt", "w");
  //     for (Node&n : dag[ti].node) {
	// for (Image_ img : n.vimg) {
	//   fprintf(fp, "%d %d %d %d\n", img.x, img.y, img.w, img.h);
	//   for (char c : img.mask)
	//     fprintf(fp, "%c", '0'+c);
	//   fprintf(fp, "\n");
	// }
  //     }
  //     exit(0);*/
  //     //dag[ti].freeAll();
  } 
    //dag[ti].benchmark();
    //exit(0);

    /*
    for (Node&n : dag[ti].node) {
      if (!n.isvec && n.pfi == dag[ti].embed1fi) {
	n.vimg.clear();
	n.vimg.shrink_to_fit();
      }
    }
    */

    /*if (ti < train.size()) {
      for (Node&n : dag[ti].node) {
	if (n.par > -1 && (n.isvec || n.img[0].w > 30 || n.img[0].h > 30)) {// || n.img[0].p != point{0,0} || n.img[0].sz != given_sizes[ti][1])) {
	  n.img.clear();
	  n.img.shrink_to_fit();
	}
      }
      dag[ti].hashi.clear();
      }*/
  }


  // if (out_sizes.size() && print_nodes) {
  //   cout << "Dag sizes: ";
  //   for (int i = 0; i < dag.size(); ++i)
  //     cout << dag[i].tiny_node.size() << " ";
  //   cout << endl;
  // }

  return dag;
}

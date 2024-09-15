/*#include <vector>
#include <cassert>
#include <iostream>
#include <tuple>
#include <functional>*/
#include "precompiled_stl.hpp"
using namespace std;

#include "utils.hpp"
#include "core_functions.hpp"

namespace core {
  unsigned short colMask(Image_ img) {
    unsigned short mask = 0;
    for (unsigned short i = 0; i < img.h; ++i)
      for (unsigned short j = 0; j < img.w; ++j)
	mask |= 1<<img(i,j);
    return mask;
  }
  unsigned short countCols(Image_ img, int include0) {//include0 = 0
    unsigned short mask = colMask(img);
    if (!include0) mask = mask&~1;
    return __builtin_popcount(mask);
  }
  unsigned short count(Image_ img) {
    unsigned short ans = 0;
    for (unsigned short i = 0; i < img.h; ++i)
      for (unsigned short j = 0; j < img.w; ++j)
	ans += img(i,j) > 0;
    return ans;
  }


  Image full(point p, point sz, int filling) {//filling = 1
    Image ret;
    ret.p = p;
    ret.sz = sz;
    ret.mask.assign(ret.h*ret.w, filling);
    return ret;
  }

  Image empty(point p, point sz) {
    return full(p, sz, 0);
  }

  Image full(point sz, int filling) {//filling = 1
    Image ret;
    ret.p = {0,0};
    ret.sz = sz;
    ret.mask.assign(ret.h*ret.w, filling);
    return ret;
  }

  Image empty(point sz) {
    return full(point{0,0}, sz, 0);
  }


  bool isRectangle(Image_ a) {
    return count(a) == a.w*a.h;
  }

  void countComponents_dfs(Image&img, unsigned short r, unsigned short c) {
    img(r,c) = 0;
    for (unsigned short nr = r-1; nr <= r+1; ++nr)
      for (unsigned short nc = c-1; nc <= c+1; ++nc)
	if (nr >= 0 && nr < img.h && nc >= 0 && nc < img.w && img(nr,nc))
	  countComponents_dfs(img,nr,nc);
  }

  unsigned short countComponents(Image img) {
    unsigned short  ans = 0;
    for (unsigned short i = 0; i < img.h; ++i) {
      for (unsigned short j = 0; j < img.w; ++j) {
	if (img(i,j)) {
	  countComponents_dfs(img,i,j);
	  /*function<void(int,int)> dfs = [&](int r, int c) {
	    if (r < 0 || r >= img.h || c < 0 || c >= img.w || !img(r,c)) return;
	    img(r,c) = 0;
	    for (int nr : {r-1,r,r+1})
	      for (int nc : {c-1,c,c+1})
		dfs(nr,nc);
	  };
	  dfs(i,j);*/
	  ++ans;
	}
      }
    }
    return ans;
  }


  char majorityCol(Image_ img, int include0) { //include0 = 0
    short cnt[10] = {};
    for (short i = 0; i < img.h; ++i)
      for (short j = 0; j < img.w; ++j) {
	char c = img(i,j);
	if (c >= 0 && c < 10)
	  cnt[c]++;
      }
    if (!include0) cnt[0] = 0;
    short ret = 0;
    short ma = cnt[ret];
    for (short c = 1; c < 10; ++c) {
      if (cnt[c] > ma) {
	ma = cnt[c];
	ret = c;
      }
    }
    return ret;
  }
  Image subImage(Image_ img, point p, point sz) {
    assert(p.x >= 0 && p.y >= 0 && p.x+sz.x <= img.w && p.y+sz.y <= img.h && sz.x >= 0 && sz.y >= 0);
    Image ret;
    ret.p = p+img.p;
    ret.sz = sz;
    ret.mask.resize(ret.w*ret.h);
    for (int i = 0; i < ret.h; ++i)
      for (int j = 0; j < ret.w; ++j)
	ret(i,j) = img(i+p.y, j+p.x);
    return ret;
  }

  vector<pair<Image,int>> splitCols(Image_ img, int include0) { //include0 = 0
    vector<pair<Image,int>> ret;
    int mask = colMask(img);
    for (char c = !include0; c < 10; ++c) {
      if (mask>>c&1) {
	Image s = img;
	for (int i = 0; i < s.h; ++i)
	  for (int j = 0; j < s.w; ++j)
	    s(i,j) = s(i,j) == c;
	ret.emplace_back(s, c);
      }
    }
    return ret;
  }
};

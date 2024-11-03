#include "precompiled_stl.hpp"
using namespace std;
#include "utils.hpp"
#include "core_functions.hpp"
#include "image_functions.hpp"

#include "visu.hpp"
#include "brute2.hpp"
#include "pieces.hpp"
#include "compose2.hpp"

#include "deduce_op.hpp"

pair<Image,Image> iOuterProductSI(Image_ img, unsigned int w, unsigned int h) {
    const unsigned int imghh = img.h / h;
    const unsigned int imgww = img.w / w;
    if (img.w * img.h <= 0 || img.w % w || img.h % h) return {badImg, badImg};
    
    Image big = core::full({static_cast<int>(imgww), static_cast<int>(imghh)}, -1);
    Image small = core::full({static_cast<int>(w), static_cast<int>(h)}, -1);

    for (int ii = 0; ii < imghh; ++ii) {
        const unsigned int iih = ii * h;
        for (int jj = 0; jj < imgww; ++jj) {
            const unsigned int jjw = jj * w;
            bool all0 = true;
            for (unsigned int i = 0; i < h && all0; ++i) {
                const unsigned int iihi = iih + i;
                for (unsigned int j = 0; j < w && all0; ++j) {
                    if (img(iihi, jjw + j)) all0 = false;
                }
            }

            big(ii, jj) = !all0;

            if (!all0) {
                for (unsigned int i = 0; i < h; ++i) {
                    const unsigned int iihi = iih + i;
                    for (unsigned int j = 0; j < w; ++j) {
                        char& a = small(i, j);
                        const unsigned char b = img(iihi, jjw + j);
                        if (a != -1 && a != b) return {badImg, badImg};
                        a = b;
                    }
                }
            }
        }
    }
    return {big, small};
}


pair<Image,Image> iOuterProductIS(Image_ img, unsigned int w, unsigned int h) {
  if (img.w*img.h <= 0 || img.w%w || img.h%h) return {badImg,badImg};
  const unsigned int imghh = img.h/h;
  const unsigned int imgww = img.w/w;
  Image big = core::full({static_cast<int>(imgww),static_cast<int>(imghh)},-1);
  Image small = core::full({static_cast<int>(w),static_cast<int>(h)},-1);

  for (unsigned int ii = 0; ii < imghh; ++ii) {
    const unsigned int iih = ii*h;
    for (unsigned int jj = 0; jj < imgww; ++jj) {
      const unsigned int jjw = jj*w;
      unsigned int mask = 0;
      for (unsigned int i = 0; i < h; ++i){
        const unsigned int iihi = iih+i;
	for (unsigned int j = 0; j < w; ++j)
	  mask |= 1<<img(iihi,jjw+j);}

      if (__builtin_popcount(mask&~1) > 1) return {badImg,badImg};
      big(ii,jj) = 31-__builtin_clz(mask);
      if (big(ii,jj)) {
	for (unsigned int i = 0; i < h; ++i) {
    const unsigned int iihi =iih+i;
	  for (unsigned int j = 0; j < w; ++j) {
	    char& a = small(i,j);
	    const unsigned char b = img(iihi,jjw+j) > 0;
	    if (a != -1 && a != b) return {badImg,badImg};
	    a = b;
	  }
	}
      }
    }
  }
  return {big, small};
}




deduceOuterProduct::deduceOuterProduct(const vector<pair<Image,Image>>& train) {

  auto score = [](vector<pair<Image,Image>> pa, int fi) {
    double ans = 0;
    for (int k : {0,1}) {
      vector<Image> imgs;
      int allequal = 1;
      for (int ti = 0; ti < pa.size(); ++ti) {
	Image img = k ? pa[ti].second : pa[ti].first;
	if (img.w*img.h <= 0) return 1e9;

	imgs.push_back(img);
	if (imgs[0] != imgs.back()) allequal = 0;
      }
      if (allequal && imgs.size() > 1) continue;
      //ham
      //#pragma omp parallel for
      for (Image_ img : imgs) {
	const int cols = __builtin_popcount(core::colMask(img)&~1);
	if (cols <= 1 && core::isRectangle(img)) {
	  ans += log(img.w+1)+log(img.h+1);
	} else if (cols <= 1) {
	  ans += log(2)*img.w*img.h;
	} else {
	  ans += log(10)*img.w*img.h;
	}
      }
      //ans *= 1-1e-5*1;//(fi*2-1);
    }
    return ans;
  };


  int minw = 1e9, minh = 1e9;
  for (const auto& [in,out] : train) {
    minw = min(minw, out.w);
    minh = min(minh, out.h);
  }

  /*{
    Image img = train[0].second;
    auto [a,b] = iOuterProductSI(img,img.w,img.h);
    assert(reconstruct(a,b) == img);
    exit(0);
    }*/

  rec_funci = -1;
  double best_score = 1e9;
  for (int h = 1; h <= minh; ++h) {
    for (int w = 1; w <= minw; ++w) {
      for (int l : {0,1}) {
	for (int k : {0,1}) {
	  auto f = k ? iOuterProductSI : iOuterProductIS;
	  vector<pair<Image,Image>> is;
	  int bad = 0;
	  for (const auto& [in,out] : train) {
	    int sw = w, sh = h;
	    if (l) {
	      if (out.w%w || out.h%h) bad = 1;
	      sw = out.w/w;
	      sh = out.h/h;
	    }
	    is.push_back(f(out, sw, sh));
	  }
	  const double entropy = score(is,k);
	  if (entropy < best_score) {
	    best_score = entropy;
	    rec_funci = k;
	    train_targets = is;
	  }
	}
      }
    }
  }
  assert(rec_funci != -1);
  const unsigned int trainsize = train.size();
  for (int k : {0,1}) {
    auto f = k ? iOuterProductSI : iOuterProductIS;
    vector<double> best_at(trainsize, 1e9);
    vector<pair<Image,Image>> best_single(trainsize, {badImg,badImg});
    //ham
    //#pragma omp parallel for
    for (unsigned int ti = 0; ti < trainsize; ++ti) {
      Image target = train[ti].second;
      for (unsigned int h = 1; h <= target.h; ++h) {
	for (unsigned int w = 1; w <= target.w; ++w) {
	  const auto is = f(target, w, h);
	  const double entropy = score({is},k);
	  if (entropy < best_at[ti]) {
	    best_at[ti] = entropy;
	    best_single[ti] = is;
	  }
	}
      }
    }
    const double entropy = score(best_single,k);
    if (entropy < best_score) {
      best_score = entropy;
      rec_funci = k;
      train_targets = best_single;
    }
  }

  assert(rec_funci != -1);
  assert(train_targets.size() == trainsize);
  //ham
  //#pragma omp parallel for
  for (int ti = 0; ti < trainsize; ++ti) {
    Image a, b;
    tie(a,b) = train_targets[ti];
    //print(a);
    //print(b);
    //print(train[ti].second);
    assert(reconstruct(a,b) == train[ti].second);
    //cout << "OK" << endl;
    }
}


Image deduceOuterProduct::reconstruct(Image_ a, Image_ b) {
  auto f = rec_funci ? outerProductSI : outerProductIS;
  return f(a,b);
}

extern int MAXDEPTH;

void addDeduceOuterProduct(Pieces&pieces, const vector<pair<Image,Image>>& train, vector<Candidate>&cands) {
  deduceOuterProduct deduce_op(train);

  int interestings = 0;
  for (const auto& [in,out] : deduce_op.train_targets) {
    if (core::count(in) > 1 && core::count(out) > 1) ++interestings;
  }
  const unsigned trainsize = train.size();
  const unsigned trainsize1 = train.size() + 1;

  if (interestings*2 < trainsize) return;

  vImage a, b;
  for (int k : {0,1}) {
    vImage&cand = k ? b : a;
    int best_match = -1;

    auto add = [&](vImage_ vi) {
      int matches = 0;
      for (int i = 0; i < trainsize; ++i) {
	Image_ target = k ? deduce_op.train_targets[i].second : deduce_op.train_targets[i].first;
	matches += (vi[i] == target);
      }
      if (matches > best_match) {
	best_match = matches;
	cand = vi;
      }
    };

    for (int pi = 0; pi < pieces.mem.size(); pi += pieces.dag.size()) {
      const int*ind = &pieces.mem[pi];
      //TODO: Use hashes to compare instead of full images
      vImage imgs;
      for (int i = 0; i <= train.size(); ++i) {
	if (pieces.dag[i].tiny_node[ind[i]].isvec) continue;
	Image_ img = pieces.dag[i].getImg(ind[i]);
	imgs.push_back(img);
	imgs.back().p = point{0,0};
      }
      if (imgs.size() == trainsize1)
	add(imgs);
    }
    for (auto [x,y] : deduce_op.train_targets) {
      add(vImage(trainsize1, (k ? y : x)));
    }
  }

  {
    assert(a.size() == train.size()+1);
    assert(b.size() == train.size()+1);

    // TODO: Use correct depths
    vImage imgs;
    for (int i = 0; i <= train.size(); ++i)
      imgs.push_back(deduce_op.reconstruct(a[i],b[i]));
    cands.emplace_back(imgs, 2, MAXDEPTH, MAXDEPTH*2);
  }

  /*
  visu.next(to_string(si));
  for (auto [in,out] : train)
    visu.add(in,out);

  visu.next(to_string(si)+"!");
  for (auto [in,out] : deduce_op.train_targets) {
    visu.add(in,out);
    }*/
}

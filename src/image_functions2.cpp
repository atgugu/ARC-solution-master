#include "precompiled_stl.hpp"
using namespace std;
#include "utils.hpp"
#include "core_functions.hpp"
#include "image_functions.hpp"
#include "spec.hpp"
#include "read.hpp"
#include "normalize.hpp"
#include <vector>
#include <queue>
#include <algorithm>
#include <unordered_set>
#include <functional>
#include <tuple>
#include <cmath>


vImage splitAll(Image_ img) {
  vector<Image> ret;
  Image done = core::empty(img.p,img.sz);
  for (int i = 0; i < img.h; ++i) {
    for (int j = 0; j < img.w; ++j) {
      if (!done(i,j)) {
	Image toadd = core::empty(img.p,img.sz);
	function<void(int,int,int)> dfs = [&](int r, int c, int col) {
	  if (r < 0 || r >= img.h || c < 0 || c >= img.w || img(r,c) != col || done(r,c)) return;
	  toadd(r,c) = img(r,c)+1;
	  done(r,c) = 1;
	  for (int d = 0; d < 4; ++d)
	    dfs(r+(d==0)-(d==1),c+(d==2)-(d==3),col);
	};
	dfs(i,j,img(i,j));
	toadd = compress(toadd);
	for (int i = 0; i < toadd.h; ++i) {
	  for (int j = 0; j < toadd.w; ++j) {
	    toadd(i,j) = max(0, toadd(i,j)-1);
	  }
	}
	//Image i = interior(toadd);
	//if (!core::count(compose(img,i,2)))
	if (core::count(toadd))
	  ret.push_back(toadd);
      }
    }
  }
  return ret;
}


Image eraseCol(Image img, int col) {
  for (int i = 0; i < img.h; ++i)
    for (int j = 0; j < img.w; ++j)
      if (img(i,j) == col) img(i,j) = 0;
  return img;
}

Image removeGrid(const Image& img) {
    if (img.w <= 1 || img.h <= 1) {
        return img;  // Return original if too small for a grid
    }

    std::vector<int> rowGridLines, colGridLines;

    // Detect potential row grid lines
    for (int i = 1; i < img.h - 1; ++i) {
        bool isGridLine = true;
        for (int j = 1; j < img.w; ++j) {
            if (img(i, j) != img(i, j - 1)) {
                isGridLine = false;
                break;
            }
        }
        if (isGridLine) {
            rowGridLines.push_back(i);
        }
    }

    // Detect potential column grid lines
    for (int j = 1; j < img.w - 1; ++j) {
        bool isGridLine = true;
        for (int i = 1; i < img.h; ++i) {
            if (img(i, j) != img(i - 1, j)) {
                isGridLine = false;
                break;
            }
        }
        if (isGridLine) {
            colGridLines.push_back(j);
        }
    }

    // Check if grid lines are regularly spaced
    if (rowGridLines.size() > 1) {
        int rowSpacing = rowGridLines[1] - rowGridLines[0];
        for (size_t i = 1; i < rowGridLines.size(); ++i) {
            if (rowGridLines[i] - rowGridLines[i - 1] != rowSpacing) {
                return img;  // Irregular spacing indicates no grid
            }
        }
    } else {
        rowGridLines.clear();
    }

    if (colGridLines.size() > 1) {
        int colSpacing = colGridLines[1] - colGridLines[0];
        for (size_t j = 1; j < colGridLines.size(); ++j) {
            if (colGridLines[j] - colGridLines[j - 1] != colSpacing) {
                return img;  // Irregular spacing indicates no grid
            }
        }
    } else {
        colGridLines.clear();
    }

    // If no regular grid was found, return the original image
    if (rowGridLines.empty() && colGridLines.empty()) {
        return img;
    }

    // Create a new image by excluding detected grid lines
    int newHeight = img.h - rowGridLines.size();
    int newWidth = img.w - colGridLines.size();
    Image result = core::empty({img.p.x, img.p.y}, {newWidth, newHeight});

    int newRow = 0;
    for (int i = 0; i < img.h; ++i) {
        if (std::find(rowGridLines.begin(), rowGridLines.end(), i) != rowGridLines.end()) {
            continue;  // Skip grid row
        }
        int newCol = 0;
        for (int j = 0; j < img.w; ++j) {
            if (std::find(colGridLines.begin(), colGridLines.end(), j) != colGridLines.end()) {
                continue;  // Skip grid column
            }
            result(newRow, newCol) = img(i, j);
            ++newCol;
        }
        ++newRow;
    }

    return result;
}

Image detectRepeatingPatternWithHole(
    Image_ img,
    bool fillHole  // Option to fill the hole in the entire image or return only the missing pattern
) {
    int imgHeight = img.h;
    int imgWidth = img.w;

    // Find the hole in the image (pixels with color 0)
    std::vector<point> holePixels;
    for (int i = 0; i < imgHeight; ++i) {
        for (int j = 0; j < imgWidth; ++j) {
            if (img(i, j) == 0) {
                holePixels.push_back({i, j});
            }
        }
    }

    if (holePixels.empty()) {
        // No hole found, return the original image as-is
        return img;
    }

    // Try different tile sizes starting from 1x1 up to half the image size
    for (int tileHeight = 1; tileHeight <= imgHeight / 2; ++tileHeight) {
        for (int tileWidth = 1; tileWidth <= imgWidth / 2; ++tileWidth) {

            // Extract the tile pattern from known (non-hole) pixels
            Image tile = {{0, 0}, {tileWidth, tileHeight}, std::vector<char>(tileWidth * tileHeight, -1)};
            bool incompleteTile = false;

            for (int i = 0; i < tileHeight; ++i) {
                for (int j = 0; j < tileWidth; ++j) {
                    int imgRow = i;
                    int imgCol = j;

                    if (imgRow < imgHeight && imgCol < imgWidth) {
                        char pixelValue = img(imgRow, imgCol);
                        if (pixelValue != 0) {
                            tile(i, j) = pixelValue;
                        } else {
                            incompleteTile = true;
                        }
                    } else {
                        incompleteTile = true;
                    }
                }
            }

            if (incompleteTile) {
                // Cannot construct a complete tile from this position, skip
                continue;
            }

            // Now, verify if this tile can generate the entire image (ignoring holes)
            bool patternMatches = true;
            for (int i = 0; i < imgHeight; ++i) {
                for (int j = 0; j < imgWidth; ++j) {
                    char expectedPixel = tile(i % tileHeight, j % tileWidth);
                    char actualPixel = img(i, j);
                    if (actualPixel != 0 && actualPixel != expectedPixel) {
                        patternMatches = false;
                        break;
                    }
                }
                if (!patternMatches) break;
            }

            if (patternMatches) {
                // The tile pattern repeats across the image (excluding holes)

                // Determine the bounding box of the hole
                int minRow = holePixels[0].x;
                int maxRow = holePixels[0].x;
                int minCol = holePixels[0].y;
                int maxCol = holePixels[0].y;

                for (const auto& p : holePixels) {
                    if (p.x < minRow) minRow = p.x;
                    if (p.x > maxRow) maxRow = p.x;
                    if (p.y < minCol) minCol = p.y;
                    if (p.y > maxCol) maxCol = p.y;
                }

                int holeHeight = maxRow - minRow + 1;
                int holeWidth = maxCol - minCol + 1;

                if (fillHole) {
                    // Option 1: Return the entire image with the hole filled
                    Image filledImg = img;
                    for (const auto& p : holePixels) {
                        int i = p.x;
                        int j = p.y;
                        char expectedPixel = tile(i % tileHeight, j % tileWidth);
                        filledImg(i, j) = expectedPixel;
                    }
                    return filledImg;

                } else {
                    // Option 2: Return only an image of the size of the hole, filled with the missing pattern
                    Image missingPat = {{0, 0}, {holeWidth, holeHeight}, std::vector<char>(holeWidth * holeHeight)};
                    for (int i = 0; i < holeHeight; ++i) {
                        for (int j = 0; j < holeWidth; ++j) {
                            missingPat(i, j) = tile((minRow + i) % tileHeight, (minCol + j) % tileWidth);
                        }
                    }
                    return missingPat;
                }
            }
        }
    }

    // No repeating pattern found that matches the image with the hole
    return Image{{0, 0}, {0, 0}, {}};  // Return an empty image to indicate no pattern found
}
Image enforceSymmetry(Image_ img, unsigned int symmetryType = 0) {
    Image result = img;
    if (symmetryType ==  0) {
        for (int i = 0; i < img.h / 2; ++i) {
            for (int j = 0; j < img.w; ++j) {
                result(img.h - 1 - i, j) = img(i, j);
            }
        }
    } else if (symmetryType == 1) {
        for (int i = 0; i < img.h; ++i) {
            for (int j = 0; j < img.w / 2; ++j) {
                result(i, img.w - 1 - j) = img(i, j);
            }
        }
    } else if (symmetryType == 2) {
        for (int i = 0; i < img.h; ++i) {
            for (int j = 0; j < img.w; ++j) {
                if (i != j && i < img.h && j < img.w) {
                    result(j, i) = img(i, j);
                }
            }
        }
    }
    return result;
}

// Helper function to upscale for morphological comparison
Image upscaleImageHelper(const Image& img, int factor) {
    int newWidth = img.w * factor;
    int newHeight = img.h * factor;
    Image upscaled = core::empty({img.p.x, img.p.y}, {newWidth, newHeight});

    for (int i = 0; i < img.h; ++i) {
        for (int j = 0; j < img.w; ++j) {
            char pixel = img(i, j);
            for (int y = i * factor; y < (i + 1) * factor; ++y) {
                for (int x = j * factor; x < (j + 1) * factor; ++x) {
                    upscaled(y, x) = pixel;
                }
            }
        }
    }
    return upscaled;
}

// Function to check if two images are morphologically similar
bool morphologicallySimilar(const Image& img1, const Image& img2) {
    if (img1.w != img2.w || img1.h != img2.h) return false;
    int similarCount = 0, totalPixels = img1.w * img1.h;
    for (int i = 0; i < img1.h; ++i) {
        for (int j = 0; j < img1.w; ++j) {
            if (img1(i, j) == img2(i, j)) {
                similarCount++;
            }
        }
    }
    double similarityRatio = static_cast<double>(similarCount) / totalPixels;
    return similarityRatio > 0.95;  // Consider morphologically similar if 95% or more pixels match
}

Image downscaleImage(const Image& img, int factor) {
    // Determine the new dimensions
    int newWidth = img.w / factor;
    int newHeight = img.h / factor;

    // Ensure the dimensions remain at least 1x1
    if (newWidth < 1 || newHeight < 1) {
        return img;  // Downscale too aggressive, return original image
    }

    // Initialize the downscaled image
    Image downscaled = core::empty({img.p.x, img.p.y}, {newWidth, newHeight});

    // Populate the downscaled image by taking the most common color in each block
    for (int i = 0; i < newHeight; ++i) {
        for (int j = 0; j < newWidth; ++j) {
            std::unordered_map<char, int> colorCount;
            int startY = i * factor;
            int startX = j * factor;
            
            // Count colors in the current block
            for (int y = startY; y < std::min(startY + factor, img.h); ++y) {
                for (int x = startX; x < std::min(startX + factor, img.w); ++x) {
                    colorCount[img(y, x)]++;
                }
            }
            
            // Find the most common color in this block
            char dominantColor = 0;
            int maxCount = 0;
            for (const auto& [color, count] : colorCount) {
                if (count > maxCount) {
                    maxCount = count;
                    dominantColor = color;
                }
            }
            downscaled(i, j) = dominantColor;
        }
    }

    // Check if the downscaled image retains morphology by upscaling it back and comparing
    Image reupscaled = upscaleImageHelper(downscaled, factor);
    if (morphologicallySimilar(img, reupscaled)) {
        return downscaled;
    } else {
        return img;  // If morphology changes significantly, return the original image
    }
}



bool detectRepeatingPattern(Image_ img, Image& pattern, int& offsetX, int& offsetY) {
    // Try different tile sizes starting from 1x1 up to half the image size
    for (int tileHeight = 1; tileHeight <= img.h / 2; ++tileHeight) {
        for (int tileWidth = 1; tileWidth <= img.w / 2; ++tileWidth) {
            bool patternFound = true;
            // Extract the pattern from the top-left corner
            Image tile = {{0, 0}, {tileWidth, tileHeight}, vector<char>(tileWidth * tileHeight)};
            for (int i = 0; i < tileHeight; ++i) {
                for (int j = 0; j < tileWidth; ++j) {
                    tile(i, j) = img(i, j);
                }
            }
            // Check if the image can be reconstructed by repeating this tile
            for (int i = 0; i < img.h; ++i) {
                for (int j = 0; j < img.w; ++j) {
                    char expectedPixel = tile(i % tileHeight, j % tileWidth);
                    if (img(i, j) != expectedPixel) {
                        patternFound = false;
                        break;
                    }
                }
                if (!patternFound) break;
            }
            if (patternFound) {
                pattern = tile;
                offsetX = 0;
                offsetY = 0;
                return true;
            }
        }
    }
    return false;
}

Image gridFilter(Image_ img, int cellHeight, int cellWidth) {
    int numRows = img.h / cellHeight;
    int numCols = img.w / cellWidth;
    Image result = {{0, 0}, {numCols, numRows}, vector<char>(numRows * numCols)};

    for (int row = 0; row < numRows; ++row) {
        for (int col = 0; col < numCols; ++col) {
            // Extract cell
            std::unordered_map<char, int> colorCount;
            for (int i = 0; i < cellHeight; ++i) {
                for (int j = 0; j < cellWidth; ++j) {
                    int x = row * cellHeight + i;
                    int y = col * cellWidth + j;
                    if (x < img.h && y < img.w) {
                        char color = img(x, y);
                        colorCount[color]++;
                    }
                }
            }
            // Determine the mode color in the cell
            char modeColor = 0;
            int maxCount = 0;
            for (const auto& [color, count] : colorCount) {
                if (count > maxCount) {
                    maxCount = count;
                    modeColor = color;
                }
            }
            result(row, col) = modeColor;
        }
    }
    return result;
}


// Looks for 4 corners
vImage insideMarked(Image_ in) {
  vector<Image> ret;
  for (int i = 0; i+1 < in.h; ++i) {
    for (int j = 0; j+1 < in.w; ++j) {
      for (int h = 1; i+h+1 < in.h; ++h) {
	for (int w = 1; j+w+1 < in.w; ++w) {
	  char col = in(i,j);
	  if (!col) continue;
	  int ok = 1;
	  for (int k = 0; k < 4; ++k) {
	    int x = j+k%2*w, y = i+k/2*h;
	    for (int d = 0; d < 4; ++d) {
	      if ((d != 3-k) == (in(y+d/2,x+d%2) != col)) {
		ok = 0;
		goto fail;
	      }
	    }
	  }
	fail:
	  if (ok) {
	    Image inside = invert(core::full(point{j+1,i+1}, point{w,h}));
	    ret.push_back(compose(inside,in,3));
	  }
	}
      }
    }
  }
  return ret;
}

Image makeBorder(Image_ img, int bcol = 1) {
  Image ret = hull0(img);
  for (int i = 0; i < ret.h; ++i) {
    for (int j = 0; j < ret.w; ++j) {
      if (img(i,j) == 0) {
	int ok = 1;
	for (int ni : {i-1,i,i+1}) {
	  for (int nj : {j-1,j,j+1}) {
	    if (img.safe(ni,nj) == 0) {
	      ok = 0;
	      break;
	    }
	  }
    if(ok == 0) break;
	}
	if (ok) {
	  ret(i,j) = bcol;
	}
      }
    }
  }
  return ret;
}


Image makeBorder2(Image_ img, int usemaj = 1) {
  int bcol = 1;
  if (usemaj) bcol = core::majorityCol(img);

  point rsz = img.sz+point{2,2};
  if (max(rsz.x,rsz.y) > MAXSIDE || rsz.x*rsz.y > MAXAREA) return badImg;
  Image ret = core::full(img.p-point{1,1}, rsz, bcol);
  for (int i = 0; i < img.h; ++i)
    for (int j = 0; j < img.w; ++j)
      ret(i+1,j+1) = img(i,j);
  return ret;
}

Image makeBorder2(Image_ img, Image_ bord) {
  int bcol = core::majorityCol(bord);
  point rsz = img.sz+bord.sz+bord.sz;
  if (max(rsz.x,rsz.y) > MAXSIDE || rsz.x*rsz.y > MAXAREA) return badImg;
  Image ret = core::full(img.p-bord.sz, rsz, bcol);
  for (int i = 0; i < img.h; ++i)
    for (int j = 0; j < img.w; ++j)
      ret(i+bord.h,j+bord.w) = img(i,j);
  return ret;
}


//Delete black rows / cols
Image compress2(Image_ img) {
  vector<int> row(img.h), col(img.w);
  for (int i = 0; i < img.h; ++i)
    for (int j = 0; j < img.w; ++j)
      if (img(i,j)) row[i] = col[j] = 1;
  vector<int> rows, cols;
  for (int i = 0; i < img.h; ++i) if (row[i]) rows.push_back(i);
  for (int j = 0; j < img.w; ++j) if (col[j]) cols.push_back(j);
  Image ret = core::empty(point{(int)cols.size(), (int)rows.size()});
  for (int i = 0; i < ret.h; ++i)
    for (int j = 0; j < ret.w; ++j)
      ret(i,j) = img(rows[i],cols[j]);
  return ret;
}


//Group single color rectangles
Image compress3(Image_ img) {
    const int height = img.h, width = img.w;
    if (width * height <= 0) return badImg;

    std::vector<int> row(height, 0), col(width, 0);
    row[0] = col[0] = 1;

    // Mark rows and columns that have unique pixels
    for (int i = 1; i < height; ++i) {
        for (int j = 1; j < width; ++j) {
            if (img(i, j) != img(i - 1, j)) row[i] = 1;
            if (img(i, j) != img(i, j - 1)) col[j] = 1;
        }
    }

    // Collect the indices of the marked rows and columns
    std::vector<int> rows, cols;
    rows.reserve(height);
    cols.reserve(width);

    for (int i = 0; i < height; ++i) {
        if (row[i]) rows.push_back(i);
    }
    for (int j = 0; j < width; ++j) {
        if (col[j]) cols.push_back(j);
    }

    // Create the compressed image
    const int compressedHeight = rows.size();
    const int compressedWidth = cols.size();
    Image ret = core::empty(point{compressedWidth, compressedHeight});

    for (int i = 0; i < compressedHeight; ++i) {
        const int rowIdx = rows[i];
        for (int j = 0; j < compressedWidth; ++j) {
            ret(i, j) = img(rowIdx, cols[j]);
        }
    }

    return ret;
}


Image greedyFill(Image &ret, std::vector<std::pair<int, std::vector<int>>> &pieces, Spec &done, int bw, int bh, int &donew) {
    // Sort pieces by size descending for better filling priority
    std::sort(pieces.rbegin(), pieces.rend());

    const int dw = ret.w - bw + 1, dh = ret.h - bh + 1;
    if (dw < 1 || dh < 1)
        return badImg;

    std::vector<int> dones(dw * dh, -1);
    std::priority_queue<std::tuple<unsigned int, int, int>> pq;

    // Function to recalculate priorities for each grid position efficiently
    auto recalc = [&](int i, int j) {
        unsigned int count = 0;
        for (int y = 0; y < bh; ++y)
            for (int x = 0; x < bw; ++x)
                count += done(i + y, j + x);

        if (count != dones[i * dw + j]) {
            dones[i * dw + j] = count;
            pq.emplace(count, j, i);
        }
    };

    // Initialize priority queue by calculating initial counts for each window in the image
    for (int i = 0; i + bh <= ret.h; ++i)
        for (int j = 0; j + bw <= ret.w; ++j)
            recalc(i, j);

    while (!pq.empty()) {
        auto [ds, j, i] = pq.top();
        pq.pop();

        // Ignore outdated priority values
        if (ds != dones[i * dw + j])
            continue;

        bool placedPattern = false;

        // Place the largest valid pattern that fits in the current region
        for (auto &[cnt, mask] : pieces) {
            bool canPlace = true;

            for (int y = 0; y < bh && canPlace; ++y) {
                const int iy = i + y;
                for (int x = 0; x < bw; ++x) {
                    if (done(iy, j + x) && ret(iy, j + x) != mask[y * bw + x]) {
                        canPlace = false;
                        break;
                    }
                }
            }

            // If pattern can be placed, apply it to the output image and update `done`
            if (canPlace) {
                for (int y = 0; y < bh; ++y) {
                    const int iy = i + y;
                    for (int x = 0; x < bw; ++x) {
                        if (!done(iy, j + x)) {
                            done(iy, j + x) = donew;
                            ret(iy, j + x) = mask[y * bw + x];
                        }
                    }
                }

                // Update the relevant region around the newly placed pattern
                const int ibh = min(static_cast<int>(i + bh), dh);
                const int jbw = std::min(static_cast<int>(j + bw), dw);

                for (int y = std::max(static_cast<int>(i - bh + 1), int(0)); y < ibh; ++y)
                    for (int x = std::max(static_cast<int>(j - bw + 1), int(0)); x < jbw; ++x)
                        recalc(y, x);

                placedPattern = true;
                break;
            }
        }

        if (!placedPattern) {
            // No suitable pattern was found for the current position, exit with bad image
            return badImg;
        }
    }

    return ret;
}
Image greedyFillBlack(Image_ img, int N = 3) {
  Image ret = core::empty(img.p, img.sz);
  Spec done;
  done.sz = img.sz;
  done.mask.assign(done.w*done.h,0);

  int donew = 1e6;
  //#pragma omp parallel for
  for (int i = 0; i < ret.h; ++i) {
    for (int j = 0; j < ret.w; ++j) {
      if (img(i,j)) {
	ret(i,j) = img(i,j);
	done(i,j) = donew;
      }
    }
  }

  map<vector<int>,int> piece_cnt;
  vector<int> mask;
  const int bw = N, bh = N;
  mask.reserve(bw*bh);

  for (int r = 0; r < 8; ++r) {
    Image rot = rigid(img,r);
    for (int i = 0; i+bh <= rot.h; ++i) {
      for (int j = 0; j+bw <= rot.w; ++j) {
	mask.resize(0);
	int ok = 1;
	for (int y = 0; y < bh; ++y){
    const int iy = i+y;
	  for (int x = 0; x < bw; ++x) {
	    char c = rot(iy,j+x);
	    mask.push_back(c);
	    if (!c) ok = 0;
	  }}
	if (ok)
	  ++piece_cnt[mask];
      }
    }
  }
  vector<pair<int,vector<int>>> piece;
  for (auto&[p,c] : piece_cnt)
    piece.emplace_back(c,p);

  return greedyFill(ret, piece, done, bw, bh, donew);
}

Image greedyFillBlack2(Image_ img, int N = 3) {
  Image filled = greedyFillBlack(img, N);
  return compose(filled, img, 4);
}

Image extend2(Image_ img, Image_ room) {
  Image ret = core::empty(room.p, room.sz);
  Spec done;
  done.sz = room.sz;
  done.mask.assign(done.w*done.h,0);

  point d = room.p-img.p;
  int donew = 1e6;
  //ham
  //#pragma omp parallel for
  for (int i = 0; i < ret.h; ++i) {
    for (int j = 0; j < ret.w; ++j) {
      const int x = j+d.x;
      const int y = i+d.y;
      if (x >= 0 && y >= 0 && x < img.w && y < img.h) {
	ret(i,j) = img(y,x);
	done(i,j) = donew;
      }
    }
  }

  map<vector<int>,int> piece_cnt;
  vector<int> mask;
  const int bw = 3, bh = 3;
  mask.reserve(bw*bh);
  //ham
  //#pragma omp parallel for
  for (int r = 0; r < 8; ++r) {
    Image rot = rigid(img,r);
    for (int i = 0; i+bh <= rot.h; ++i) {
      for (int j = 0; j+bw <= rot.w; ++j) {
	mask.resize(0);
	for (int y = 0; y < bh; ++y)
	  for (int x = 0; x < bw; ++x)
	    mask.push_back(rot(i+y,j+x));
	++piece_cnt[mask];
      }
    }
  }
  vector<pair<int,vector<int>>> piece;
  for (auto&[p,c] : piece_cnt)
    piece.emplace_back(c,p);

  return greedyFill(ret, piece, done, bw, bh, donew);
}



Image connect(Image_ img, int id) {
  assert(id >= 0 && id < 3);
  Image ret = core::empty(img.p, img.sz);

  //Horizontal
  if (id == 0 || id == 2) {
    for (int i = 0; i < img.h; ++i) {
      int last = -1, lastc = -1;
      for (int j = 0; j < img.w; ++j) {
	if (img(i,j)) {
	  if (img(i,j) == lastc) {
	    for (int k = last+1; k < j; ++k)
	      ret(i,k) = lastc;
	  }
	  lastc = img(i,j);
	  last = j;
	}
      }
    }
  }

  //Vertical

  if (id == 1 || id == 2) {

    for (int j = 0; j < img.w; ++j) {
      int last = -1, lastc = -1;
      for (int i = 0; i < img.h; ++i) {
	if (img(i,j)) {
	  if (img(i,j) == lastc) {
	    for (int k = last+1; k < i; ++k)
	      ret(k,j) = lastc;
	  }
	  lastc = img(i,j);
	  last = i;
	}
      }
    }
  }

  return ret;
}

Image replaceTemplate(Image_ in, Image_ need_, Image_ marked_, int overlapping = 0, int rigids = 0) {
  if (marked_.sz != need_.sz) return badImg;
  if (need_.w*need_.h <= 0) return in;

  const int rots = rigids ? 8 : 1;
  vector<Image> needr(rots), markedr(rots);
  //ham
  //#pragma omp parallel for
  for (int r = 0; r < rots; ++r) {
    needr[r] = rigid(need_,r);
    markedr[r] = rigid(marked_,r);
  }

  Image ret = in;
  //ham
  //#pragma omp parallel for
  for (int r = 0; r < rots; ++r) {
    Image_ need = needr[r];
    Image_ marked = markedr[r];

    for (int i = 0; i+need.h <= ret.h; ++i) {
      for (int j = 0; j+need.w <= ret.w; ++j) {
	int ok = 1;
  //#pragma omp parallel for
	for (int y = 0; y < need.h; ++y)
	  for (int x = 0; x < need.w; ++x)
	    if ((overlapping ? in : ret)(i+y,j+x) != need(y,x)) ok = 0;

	if (overlapping == 2) {
    //#pragma omp parallel for
	  for (int y = -1; y <= need.h; ++y) {
	    for (int x = -1; x <= need.w; ++x) {
	      if (x >= 0 && y >= 0 && x < need.w && y < need.h) continue;

	      char nn = need(clamp(y,0,need.h-1),
			     clamp(x,0,need.w-1));
	      if (nn && nn == in.safe(i+y,j+x)) ok = 0;
	    }
	  }
	}

	if (ok) {
	  for (int y = 0; y < need.h; ++y)
	    for (int x = 0; x < need.w; ++x)
	      ret(i+y,j+x) = marked(y,x);
	}
      }
    }
  }
  return ret;
}



Image swapTemplate(Image_ in, Image_ a, Image_ b, int rigids = 0) {
  if (a.sz != b.sz) return badImg;
  if (a.w*a.h <= 0) return in;

  const int rots = rigids ? 8 : 1;
  vector<Image> ar(rots), br(rots);
  for (int r = 0; r < rots; ++r) {
    ar[r] = rigid(a,r);
    br[r] = rigid(b,r);
  }
  Image done = hull0(in), ret = in;
  //ham
  //#pragma omp parallel for
  for (int k : {0,1}) {
    //#pragma omp parallel for
    for (int r = 0; r < rots; ++r) {
      Image_ need = k ? ar[r] : br[r];
      Image_ to   = k ? br[r] : ar[r];

      for (int i = 0; i+need.h <= ret.h; ++i) {
	for (int j = 0; j+need.w <= ret.w; ++j) {

	  int ok = 1;
	  for (int y = 0; y < need.h; ++y)
	    for (int x = 0; x < need.w; ++x)
	      if (done(i+y,j+x) || ret(i+y,j+x) != need(y,x)) { ok = 0; break; } //ham
	  if (ok) {
	    for (int y = 0; y < need.h; ++y) {
	      for (int x = 0; x < need.w; ++x) {
		ret(i+y,j+x) = to(y,x);
		done(i+y,j+x) = 1;
	      }
	    }
	  }
	}

      }
    }
  }
  return ret;
}


Image spreadCols(Image img, int skipmaj = 0) {
  int skipcol = -1;
  if (skipmaj)
    skipcol = core::majorityCol(img);

  Image done = hull0(img);
  queue<tuple<int,int,int>> q;
  for (int i = 0; i < img.h; ++i) {
    for (int j = 0; j < img.w; ++j) {
      if (img(i,j)) {
	if (img(i,j) != skipcol)
	  q.emplace(j,i,img(i,j));
	done(i,j) = 1;
      }
    }
  }
  while (q.size()) {
    auto [j,i,c] = q.front();
    q.pop();
    for (int d = 0; d < 4; ++d) {
      int ni = i+(d==0)-(d==1);
      int nj = j+(d==2)-(d==3);
      if (ni >= 0 && nj >= 0 && ni < img.h && nj < img.w && !done(ni,nj)) {
	img(ni,nj) = c;
	done(ni,nj) = 1;
	q.emplace(nj,ni,c);
      }
    }
  }
  return img;
}


Image dilation(Image_ img) {
    Image result = img;
    vector<point> directions = {{-1,0},{1,0},{0,-1},{0,1}};
    for (int i = 0; i < img.h; ++i) {
        for (int j = 0; j < img.w; ++j) {
            if (img(i, j) != 0) {
                for (const auto& dir : directions) {
                    int ni = i + dir.x;
                    int nj = j + dir.y;
                    if (ni >= 0 && ni < img.h && nj >= 0 && nj < img.w) {
                        result(ni, nj) = img(i, j);
                    }
                }
            }
        }
    }
    return result;
}

Image erosion(Image_ img) {
    Image result = img;
    vector<point> directions = {{-1,0},{1,0},{0,-1},{0,1}};
    for (int i = 0; i < img.h; ++i) {
        for (int j = 0; j < img.w; ++j) {
            if (img(i, j) != 0) {
                for (const auto& dir : directions) {
                    int ni = i + dir.x;
                    int nj = j + dir.y;
                    if (img.safe(ni, nj) == 0) {
                        result(i, j) = 0;
                        break;
                    }
                }
            }
        }
    }
    return result;
}

Image detectCrossPattern(Image_ img, unsigned int size = 3) {
    if (img.h < size || img.w < size || size < 3) return img;

    Image pattern = {{0, 0}, img.sz, std::vector<char>(img.h * img.w, 0)};
    int mid = size / 2;

    for (int i = mid; i <= img.h - mid - 1; ++i) {
        for (int j = mid; j <= img.w - mid - 1; ++j) {
            int val = img(i, j);
            bool isCross = true;
            for (int k = -mid; k <= mid; ++k) {
                if (img(i + k, j) != val || img(i, j + k) != val) {
                    isCross = false;
                    break;
                }
            }
            if (isCross) {
                for (int k = -mid; k <= mid; ++k) {
                    pattern(i + k, j) = val;
                    pattern(i, j + k) = val;
                }
            }
        }
    }
    return pattern;
}
Image detectHorizontalStripes(Image_ img, unsigned int size = 2) {
    if (img.h < size) return img;

    Image pattern = {{0, 0}, img.sz, std::vector<char>(img.h * img.w, 0)};
    for (int i = 0; i <= img.h - size; ++i) {
        for (int j = 0; j < img.w; ++j) {
            bool isStripe = true;
            int val = img(i, j);
            for (int k = 1; k < size; ++k) {
                if (img(i + k, j) != val) {
                    isStripe = false;
                    break;
                }
            }
            if (isStripe) {
                for (int k = 0; k < size; ++k) {
                    pattern(i + k, j) = val;
                }
            }
        }
    }
    return pattern;
}
Image detectVerticalStripes(Image_ img, unsigned int size = 2) {
    if (img.w < size) return img;

    Image pattern = {{0, 0}, img.sz, std::vector<char>(img.h * img.w, 0)};
    for (int j = 0; j <= img.w - size; ++j) {
        for (int i = 0; i < img.h; ++i) {
            bool isStripe = true;
            int val = img(i, j);
            for (int k = 1; k < size; ++k) {
                if (img(i, j + k) != val) {
                    isStripe = false;
                    break;
                }
            }
            if (isStripe) {
                for (int k = 0; k < size; ++k) {
                    pattern(i, j + k) = val;
                }
            }
        }
    }
    return pattern;
}
struct pair_hash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2>& pair) const noexcept {
        std::size_t h1 = std::hash<T1>{}(pair.first);
        std::size_t h2 = std::hash<T2>{}(pair.second);
        return h1 ^ (h2 << 1);
    }
};

Image detectTranslationPattern(Image_ img) {
    int h = img.h;
    int w = img.w;

    // Find the smallest vertical and horizontal shifts where the pattern repeats
    int bestR = h;
    int bestS = w;
    for (int r = 1; r < h; ++r) {
        bool matches = true;
        for (int i = 0; i < h - r; ++i) {
            for (int j = 0; j < w; ++j) {
                if (img(i, j) != img(i + r, j)) {
                    matches = false;
                    break;
                }
            }
            if (!matches) break;
        }
        if (matches) {
            bestR = r;
            break;
        }
    }

    for (int s = 1; s < w; ++s) {
        bool matches = true;
        for (int i = 0; i < h; ++i) {
            for (int j = 0; j < w - s; ++j) {
                if (img(i, j) != img(i, j + s)) {
                    matches = false;
                    break;
                }
            }
            if (!matches) break;
        }
        if (matches) {
            bestS = s;
            break;
        }
    }

    if (bestR == h && bestS == w) {
        // No repeating pattern found
        return img;
    }

    // Create the pattern image
    Image pattern = {{0, 0}, {bestS, bestR}, std::vector<char>(bestR * bestS)};
    for (int i = 0; i < bestR; ++i) {
        for (int j = 0; j < bestS; ++j) {
            pattern(i, j) = img(i, j);
        }
    }
    return pattern;
}
Image enforceRotationalSymmetry90(Image_ img) {
    int h = img.h;
    int w = img.w;
    Image result = img;

    if (h != w) {
        // Rotation symmetry requires a square image
        return img;
    }

    int n = h; // Since h == w
    for (int i = 0; i < (n + 1) / 2; ++i) {
        for (int j = i; j < n - i - 1; ++j) {
            // Rotate the pixels in groups of four
            char temp = result(i, j);
            result(i, j) = result(n - 1 - j, i);
            result(n - 1 - j, i) = result(n - 1 - i, n - 1 - j);
            result(n - 1 - i, n - 1 - j) = result(j, n - 1 - i);
            result(j, n - 1 - i) = temp;
        }
    }
    return result;
}

Image enforceRotationalSymmetry180(Image_ img) {
    int h = img.h;
    int w = img.w;
    Image result = img;

    for (int i = 0; i < (h + 1) / 2; ++i) {
        for (int j = 0; j < w; ++j) {
            int ni = h - 1 - i;
            int nj = w - 1 - j;
            char temp = result(i, j);
            result(i, j) = result(ni, nj);
            result(ni, nj) = temp;
        }
    }
    return result;
}

Image enforceDiagonalSymmetryNESW(Image_ img) {
    int h = img.h;
    int w = img.w;
    Image result = img;

    int minDim = std::min(h, w);
    for (int i = 0; i < minDim; ++i) {
        for (int j = 0; j < w - i - 1; ++j) {
            // Swap pixels symmetrically across the NE-SW diagonal
            char temp = result(i, j);
            result(i, j) = result(h - j - 1, w - i - 1);
            result(h - j - 1, w - i - 1) = temp;
        }
    }
    return result;
}

Image detectTranslation1DPattern(Image_ img) {
    int h = img.h;
    int w = img.w;
    std::vector<std::pair<int, int>> possibleShifts;

    for (int r = -h + 1; r < h; ++r) {
        for (int s = -w + 1; s < w; ++s) {
            if (r == 0 && s == 0) continue;
            bool possible = true;
            std::unordered_map<std::pair<int, int>, char, pair_hash> equivColors;
            for (int i = 0; i < h; ++i) {
                if (!possible) break;
                for (int j = 0; j < w; ++j) {
                    int u = i * s - j * r;
                    int v = (i * r + j * s + 100 * (r * r + s * s)) % (r * r + s * s);
                    char color = img(i, j);
                    std::pair<int, int> key = {u, v};
                    if (equivColors.find(key) == equivColors.end()) {
                        equivColors[key] = color;
                    } else if (equivColors[key] != color) {
                        possible = false;
                        break;
                    }
                }
            }
            if (possible) {
                possibleShifts.emplace_back(r, s);
            }
        }
    }

    if (possibleShifts.empty()) {
        // No repeating pattern found
        return img;
    }

    // Choose the shift with the minimal absolute sum
    std::sort(possibleShifts.begin(), possibleShifts.end(), [](const auto& a, const auto& b) {
        return std::abs(a.first) + std::abs(a.second) < std::abs(b.first) + std::abs(b.second);
    });

    int bestR = possibleShifts[0].first;
    int bestS = possibleShifts[0].second;


    // Build the equivalence classes
    std::unordered_map<std::pair<int, int>, std::vector<std::pair<int, int>>, pair_hash> equivalenceClasses;
    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; ++j) {
            int u = i * bestS - j * bestR;
            int v = (i * bestR + j * bestS + 100 * (bestR * bestR + bestS * bestS)) % (bestR * bestR + bestS * bestS);
            equivalenceClasses[{u, v}].emplace_back(i, j);
        }
    }

    // Build the pattern image
    Image pattern = img;
    for (const auto& [key, positions] : equivalenceClasses) {
        if (positions.size() > 1) {
            char color = img(positions[0].first, positions[0].second);
            for (const auto& pos : positions) {
                pattern(pos.first, pos.second) = color;
            }
        }
    }
    return pattern;
}



Image detectDiagonalPattern(Image_ img, unsigned int size = 2) {
    if (img.h < size || img.w < size) return img;

    Image pattern = {{0, 0}, img.sz, std::vector<char>(img.h * img.w, 0)};
    for (int i = 0; i <= img.h - size; ++i) {
        for (int j = 0; j <= img.w - size; ++j) {
            bool isDiagonal = true;
            int val = img(i, j);
            for (int k = 1; k < size; ++k) {
                if (img(i + k, j + k) != val) {
                    isDiagonal = false;
                    break;
                }
            }
            if (isDiagonal) {
                for (int k = 0; k < size; ++k) {
                    pattern(i + k, j + k) = val;
                }
            }
        }
    }
    return pattern;
}

Image detectCheckerboardPattern(Image_ img, unsigned int size = 2) {
    if (img.h < size || img.w < size) return img;

    Image pattern = {{0, 0}, img.sz, std::vector<char>(img.h * img.w, 0)};
    for (int i = 0; i <= img.h - size; i += size) {
        for (int j = 0; j <= img.w - size; j += size) {
            int val = img(i, j);
            bool isCheckerboard = true;
            for (int x = 0; x < size; ++x) {
                for (int y = 0; y < size; ++y) {
                    if ((x + y) % 2 == 0) {
                        if (img(i + x, j + y) != val) {
                            isCheckerboard = false;
                            break;
                        }
                    } else {
                        if (img(i + x, j + y) == val) {
                            isCheckerboard = false;
                            break;
                        }
                    }
                }
                if (!isCheckerboard) break;
            }
            if (isCheckerboard) {
                for (int x = 0; x < size; ++x) {
                    for (int y = 0; y < size; ++y) {
                        pattern(i + x, j + y) = img(i + x, j + y);
                    }
                }
            }
        }
    }
    return pattern;
}

Image repairRotationalSymmetry(const Image& img) {
    if (img.h <= 0 || img.w <= 0) {
        return img;
    }

    Image result = img;  // Create a copy of the input image

    // Repair rotational symmetry by averaging
    for (int i = 0; i < img.h; ++i) {
        for (int j = 0; j < img.w; ++j) {
            int opposite_i = img.h - i - 1;
            int opposite_j = img.w - j - 1;

            // Calculate the average value between the pixel and its rotationally symmetric counterpart
            int average_value = (img(i, j) + img(opposite_i, opposite_j)) / 2;

            // Update both the original and the symmetric pixel in the result
            result(i, j) = average_value;
            result(opposite_i, opposite_j) = average_value;
        }
    }

    return result;
}

Image applyColorMapping(const Image_ &img, const std::unordered_map<int, int> &colorMap) {
    // Create a copy of the image to modify and return
    Image result = img;

    for (int i = 0; i < img.h; ++i) {
        for (int j = 0; j < img.w; ++j) {
            int originalColor = img(i, j);

            // Check if the color is in the map
            if (colorMap.find(originalColor) != colorMap.end()) {
                // Change the color based on the mapping
                result(i, j) = colorMap.at(originalColor);
            }
            // If not found in map, keep the original color
        }
    }

    return result;
}
Image detectPatterns(Image_ img, unsigned int size = 2) {
    if(img.h < size || img.w < size) return img;

    Image pattern = {{0,0}, img.sz, vector<char>(img.h * img.w, 0)};
    for (int i = 0; i <= img.h - 2; ++i) {
        for (int j = 0; j <= img.w - 2; ++j) {
            int val = img(i, j);
            if (val == img(i, j + 1) && val == img(i + 1, j) && val == img(i + 1, j + 1)) {
                pattern(i, j) = val;
                pattern(i, j + 1) = val;
                pattern(i + 1, j) = val;
                pattern(i + 1, j + 1) = val;
            }
        }
    }
    return pattern;
}

Image detectPatterns2Image2(Image_ img) {
  return detectPatterns(img, 2);
}

Image detectPatterns2Image3(Image_ img) {
  return detectPatterns(img, 3);
}

Image detectPatterns2Image4(Image_ img) {
  return detectPatterns(img, 4);
}

Image detectPatterns2Image5(Image_ img) {
  return detectPatterns(img, 5);
}

Image detectPatterns2Image6(Image_ img) {
  return detectPatterns(img, 6);
}

Image detectPatterns2Image7(Image_ img) {
  return detectPatterns(img, 7);
}

Image detectPatterns2Image8(Image_ img) {
  return detectPatterns(img, 8);
}

Image detectPatterns2Image9(Image_ img) {
  return detectPatterns(img, 9);
}

Image combineShapes(const vector<Image>& shapes, int rows, int cols) {
    // Initialize the result image with dimensions (rows, cols) and filled with 0s
    Image result = {{0, 0}, {cols, rows}, vector<char>(rows * cols, 0)};

    // Iterate through each shape and place its pixels in the appropriate position in `result`
    for (const auto& shape : shapes) {
        for (int i = 0; i < shape.h; ++i) {
            for (int j = 0; j < shape.w; ++j) {
                // Only place the pixel if it is part of the shape (non-zero value in `shape.mask`)
                if (shape(i, j) != 0) {
                    int x = shape.p.x + i;  // Original x-coordinate of the pixel in `result`
                    int y = shape.p.y + j;  // Original y-coordinate of the pixel in `result`

                    // Ensure coordinates are within bounds before assigning
                    if (x >= 0 && x < rows && y >= 0 && y < cols) {
                        result(x, y) = shape(i, j);
                    }
                }
            }
        }
    }

    return result;
}

#include <vector>
#include <algorithm>
#include <unordered_set>

bool detectRepeatingPatternWithHole(
    Image_ img,
    Image& filledImage,     // Output: the entire image with the hole filled
    Image& missingPattern,  // Output: only the missing pattern that fills the hole
    bool fillHole           // Option to fill the hole in the entire image or return only the missing pattern
) {
    int imgHeight = img.h;
    int imgWidth = img.w;

    // Find the hole in the image (pixels with color 0)
    std::vector<point> holePixels;
    for (int i = 0; i < imgHeight; ++i) {
        for (int j = 0; j < imgWidth; ++j) {
            if (img(i, j) == 0) {
                holePixels.push_back({i, j});
            }
        }
    }

    if (holePixels.empty()) {
        // No hole found, the image is complete
        return false;
    }

    // Try different tile sizes starting from 1x1 up to half the image size
    for (int tileHeight = 1; tileHeight <= imgHeight / 2; ++tileHeight) {
        for (int tileWidth = 1; tileWidth <= imgWidth / 2; ++tileWidth) {

            // Extract the tile pattern from known (non-hole) pixels
            Image tile = {{0, 0}, {tileWidth, tileHeight}, std::vector<char>(tileWidth * tileHeight, -1)};
            bool incompleteTile = false;

            for (int i = 0; i < tileHeight; ++i) {
                for (int j = 0; j < tileWidth; ++j) {
                    // Map tile position to image position
                    int imgRow = i;
                    int imgCol = j;

                    if (imgRow < imgHeight && imgCol < imgWidth) {
                        char pixelValue = img(imgRow, imgCol);
                        if (pixelValue != 0) {
                            tile(i, j) = pixelValue;
                        } else {
                            incompleteTile = true;
                        }
                    } else {
                        incompleteTile = true;
                    }
                }
            }

            if (incompleteTile) {
                // Cannot construct a complete tile from this position, skip
                continue;
            }

            // Now, verify if this tile can generate the entire image (ignoring holes)
            bool patternMatches = true;
            for (int i = 0; i < imgHeight; ++i) {
                for (int j = 0; j < imgWidth; ++j) {
                    char expectedPixel = tile(i % tileHeight, j % tileWidth);
                    char actualPixel = img(i, j);
                    if (actualPixel != 0 && actualPixel != expectedPixel) {
                        patternMatches = false;
                        break;
                    }
                }
                if (!patternMatches) break;
            }

            if (patternMatches) {
                // The tile pattern repeats across the image (excluding holes)

                // Now, fill the hole(s)
                Image filledImg = img;
                for (const auto& p : holePixels) {
                    int i = p.x;
                    int j = p.y;
                    char expectedPixel = tile(i % tileHeight, j % tileWidth);
                    filledImg(i, j) = expectedPixel;
                }

                filledImage = filledImg;

                if (fillHole) {
                    // Return the entire image with the hole filled
                    return true;
                } else {
                    // Return only the missing pattern corresponding to the hole
                    // Determine the bounding box of the hole
                    int minRow = holePixels[0].x;
                    int maxRow = holePixels[0].x;
                    int minCol = holePixels[0].y;
                    int maxCol = holePixels[0].y;

                    for (const auto& p : holePixels) {
                        if (p.x < minRow) minRow = p.x;
                        if (p.x > maxRow) maxRow = p.x;
                        if (p.y < minCol) minCol = p.y;
                        if (p.y > maxCol) maxCol = p.y;
                    }

                    int holeHeight = maxRow - minRow + 1;
                    int holeWidth = maxCol - minCol + 1;

                    // Extract the missing pattern from the filled image
                    Image missingPat = {{0, 0}, {holeWidth, holeHeight}, std::vector<char>(holeWidth * holeHeight)};
                    for (int i = 0; i < holeHeight; ++i) {
                        for (int j = 0; j < holeWidth; ++j) {
                            missingPat(i, j) = filledImg(minRow + i, minCol + j);
                        }
                    }

                    missingPattern = missingPat;
                    return true;
                }
            }
        }
    }

    // No repeating pattern found that matches the image with the hole
    return false;
}


Image trimToContent(Image_ img) {
    int minX = img.h, minY = img.w, maxX = -1, maxY = -1;
    for (int i = 0; i < img.h; ++i) {
        for (int j = 0; j < img.w; ++j) {
            if (img(i, j) != 0) {
                minX = std::min(minX, i);
                minY = std::min(minY, j);
                maxX = std::max(maxX, i);
                maxY = std::max(maxY, j);
            }
        }
    }
    if (minX > maxX || minY > maxY) {
        // No content found
        return Image{{0, 0}, {0, 0}, {}};
    }
    int newHeight = maxX - minX + 1;
    int newWidth = maxY - minY + 1;
    Image result = {{0, 0}, {newWidth, newHeight}, vector<char>(newWidth * newHeight)};
    for (int i = 0; i < newHeight; ++i) {
        for (int j = 0; j < newWidth; ++j) {
            result(i, j) = img(minX + i, minY + j);
        }
    }
    return result;
}


vector<Image> extractConnectedComponents(Image_ img) {
    int rows = img.h;
    int cols = img.w;
    vector<bool> visited(rows * cols, false); // 1D visited array
    vector<Image> components;

    // Directions for 4-neighbor connectivity (up, down, left, right)
    vector<point> directions = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};

    // Helper lambda for 2D to 1D index
    auto index = [&](int x, int y) { return x * cols + y; };

    // Iterate through each pixel in the image
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (!visited[index(i, j)] && img(i, j) != 0) {  // Non-zero and unvisited
                char value = img(i, j);

                // Bounding box of the component
                int minX = i, maxX = i, minY = j, maxY = j;

                // Points belonging to this component
                vector<point> points;

                // BFS queue for flood fill
                queue<point> q;
                q.push({i, j});
                visited[index(i, j)] = true;

                while (!q.empty()) {
                    point p = q.front();
                    q.pop();

                    points.push_back(p);

                    // Update bounding box
                    minX = min(minX, p.x);
                    maxX = max(maxX, p.x);
                    minY = min(minY, p.y);
                    maxY = max(maxY, p.y);

                    // Explore neighbors
                    for (const point& dir : directions) {
                        int nx = p.x + dir.x;
                        int ny = p.y + dir.y;

                        // Check bounds and ensure the pixel is part of the same component
                        if (nx >= 0 && nx < rows && ny >= 0 && ny < cols &&
                            !visited[index(nx, ny)] && img(nx, ny) == value) {
                            visited[index(nx, ny)] = true;
                            q.push({nx, ny});
                        }
                    }
                }

                // Calculate width and height of the component
                int compWidth = maxY - minY + 1;
                int compHeight = maxX - minX + 1;

                // Initialize the mask for this component with zeros
                Image component = {{minX, minY}, {compWidth, compHeight}, vector<char>(compWidth * compHeight, 0)};

                // Place the points in the component image's mask
                for (const point& pt : points) {
                    int localX = pt.x - minX;
                    int localY = pt.y - minY;
                    component(localX, localY) = value;
                }

                // Add the component to the list of components
                components.push_back(std::move(component));
            }
        }
    }

    return components;
}


Image mostCommonShape(Image_ img) {
    auto shapes = extractConnectedComponents(img);
    map<int, int> shapeCounts;
    for (const auto& shape : shapes) {
        int size = shape.mask.size();
        shapeCounts[size]++;
    }
    int maxCount = 0;
    int commonSize = 0;
    for (const auto& [size, count] : shapeCounts) {
        if (count > maxCount) {
            maxCount = count;
            commonSize = size;
        }
    }
    vector<Image> commonShapes;
    for (const auto& shape : shapes) {
        if (shape.mask.size() == commonSize) {
            commonShapes.push_back(shape);
        }
    }
    Image result = combineShapes(commonShapes, img.h, img.w);
    return result;
}


Image dilate(Image_ img, int iterations = 1) {
    Image result = img;
    for (int iter = 0; iter < iterations; ++iter) {
        result = dilation(result);
    }
    return result;
}

Image erode(Image_ img, int iterations = 1) {
    Image result = img;
    for (int iter = 0; iter < iterations; ++iter) {
        result = erosion(result);
    }
    return result;
}

Image dilate1(Image_ img) {
    return dilate(img, 1);
}

Image erode1(Image_ img) {
    return erode(img, 1);
}

Image dilate2(Image_ img) {
    return dilate(img, 2);
}

Image erode2(Image_ img) {
    return erode(img, 2);
}

Image dilate3(Image_ img) {
    return dilate(img, 3);
}

Image erode3(Image_ img) {
    return erode(img, 3);
}

Image ringSmear(Image_ base, Image_ room) {
    Image ret = embed(base, hull(room));
    point center = {base.w / 2, base.h / 2};

    for (int r = 0; r < std::max(base.h, base.w) / 2; ++r) {
        for (int i = -r; i <= r; ++i) {
            for (int j = -r; j <= r; ++j) {
                int x = center.x + i;
                int y = center.y + j;
                if (x >= 0 && x < base.w && y >= 0 && y < base.h) {
                    if (room.safe(y, x) && base.safe(y, x)) {
                        ret(y, x) = base(y, x);
                    }
                }
            }
        }
    }

    return ret;
}

Image checkerboardSmear(Image_ base, Image_ room) {
    Image ret = embed(base, hull(room));

    for (int i = 0; i < ret.h; ++i) {
        for (int j = 0; j < ret.w; ++j) {
            if ((i + j) % 2 == 0) {  // Checkerboard condition
                if (room(i, j) && base.safe(i, j)) {
                    ret(i, j) = base(i, j);
                }
            }
        }
    }

    return ret;
}

Image diagonalSmear(Image_ base, Image_ room, int id) {
    assert(id >= 0 && id < 4);
    const int arr[] = {1, 2, 3, 4};  // 1: Top-left to Bottom-right, 2: Bottom-right to Top-left, etc.
    const int mask = arr[id];

    point d = room.p - base.p;
    Image ret = embed(base, hull(room));

    if (mask == 1) {
        // Top-left to Bottom-right
        for (int i = 0; i < ret.h; ++i) {
            char c = 0;
            for (int j = 0; j < ret.w; ++j) {
                if (!room(i, j)) c = 0;
                else if (base.safe(i + d.y, j + d.x)) c = base(i + d.y, j + d.x);
                if (c) ret(i, j) = c;
            }
        }
    } else if (mask == 2) {
        // Bottom-right to Top-left
        for (int i = ret.h - 1; i >= 0; --i) {
            char c = 0;
            for (int j = ret.w - 1; j >= 0; --j) {
                if (!room(i, j)) c = 0;
                else if (base.safe(i + d.y, j + d.x)) c = base(i + d.y, j + d.x);
                if (c) ret(i, j) = c;
            }
        }
    } else if (mask == 3) {
        // Top-right to Bottom-left
        for (int i = 0; i < ret.h; ++i) {
            char c = 0;
            for (int j = ret.w - 1; j >= 0; --j) {
                if (!room(i, j)) c = 0;
                else if (base.safe(i + d.y, j + d.x)) c = base(i + d.y, j + d.x);
                if (c) ret(i, j) = c;
            }
        }
    } else if (mask == 4) {
        // Bottom-left to Top-right
        for (int i = ret.h - 1; i >= 0; --i) {
            char c = 0;
            for (int j = 0; j < ret.w; ++j) {
                if (!room(i, j)) c = 0;
                else if (base.safe(i + d.y, j + d.x)) c = base(i + d.y, j + d.x);
                if (c) ret(i, j) = c;
            }
        }
    }

    return ret;
}



vImage splitColumns(Image_ img) {
  if (img.w*img.h <= 0) return {};
  vector<Image> ret(img.w);
  for (int j = 0; j < img.w; ++j) {
    ret[j].p = {j,0};
    ret[j].sz = {1,img.h};
    ret[j].mask.resize(img.h);
    for (int i = 0; i < img.h; ++i)
      ret[j].mask[i] = img(i,j);
  }
  return ret;
}
vImage splitRows(Image_ img) {
  if (img.w*img.h <= 0) return {};
  vector<Image> ret(img.h);
  for (int i = 0; i < img.h; ++i) {
    ret[i].p = {0,i};
    ret[i].sz = {img.w,1};
    ret[i].mask.resize(img.w);
    for (int j = 0; j < img.w; ++j)
      ret[i].mask[j] = img(i,j);
  }
  return ret;
}


Image half(Image_ img, int id) {
  assert(id >= 0 && id < 4);
  if (id == 0) {
    return core::subImage(img, {0,0}, {img.w/2, img.h});
  } else if (id == 1) {
    return core::subImage(img, {img.w-img.w/2,0}, {img.w/2, img.h});
  } else if (id == 2) {
    return core::subImage(img, {0,0}, {img.w, img.h/2});
  } else if (id == 3) {
    return core::subImage(img, {0,img.h-img.h/2}, {img.w, img.h/2});
  }
  return badImg;
}
Image diagonalSmear(Image_ img, int id) {
    Image ret = img;

    if (id == 0) { // Top-left to bottom-right
        for (int i = 0; i < ret.h; ++i) {
            char c = 0;
            for (int j = 0; j < ret.w; ++j) {
                if (!img(i, j)) c = 0;
                else c = img(i, j);
                if (c) ret(i, j) = c;
            }
        }
    } else if (id == 1) { // Bottom-right to top-left
        for (int i = ret.h - 1; i >= 0; --i) {
            char c = 0;
            for (int j = ret.w - 1; j >= 0; --j) {
                if (!img(i, j)) c = 0;
                else c = img(i, j);
                if (c) ret(i, j) = c;
            }
        }
    } else if (id == 2) { // Top-right to bottom-left
        for (int i = 0; i < ret.h; ++i) {
            char c = 0;
            for (int j = ret.w - 1; j >= 0; --j) {
                if (!img(i, j)) c = 0;
                else c = img(i, j);
                if (c) ret(i, j) = c;
            }
        }
    } else if (id == 3) { // Bottom-left to top-right
        for (int i = ret.h - 1; i >= 0; --i) {
            char c = 0;
            for (int j = 0; j < ret.w; ++j) {
                if (!img(i, j)) c = 0;
                else c = img(i, j);
                if (c) ret(i, j) = c;
            }
        }
    }

    return ret;
}


Image smear(Image_ img, int id) {
  assert(id >= 0 && id < 15);
  const pair<int,int> R = {1,0}, L = {-1,0}, D = {0,1}, U = {0,-1};
  const pair<int,int> X = {1,1}, Y = {-1,-1}, Z = {1,-1}, W = {-1,1};
  const vector<pair<int,int>> d[15] = {{R},{L},{D},{U},
				       {R,L},{D,U},
				       {R,L,D,U},
				       {X},{Y},{Z},{W},
				       {X,Y},{Z,W},
				       {X,Y,Z,W},
				       {R,L,D,U,X,Y,Z,W}};
  const int w = img.w;
  Image ret = img;


  for (auto [dx,dy] : d[id]) {
    int di = dy*w+dx;

    for (int i = 0; i < ret.h; ++i) {
      int step = i == 0 || i == ret.h-1 ? 1 : max(ret.w-1,1);
      for (int j = 0; j < ret.w; j += step) {
	if (i-dy < 0 || j-dx < 0 || i-dy >= img.h || j-dx >= img.w) {
	  int steps = MAXSIDE;
	  if (dx ==-1) steps = min(steps, j+1);
	  if (dx == 1) steps = min(steps, img.w-j);
	  if (dy ==-1) steps = min(steps, i+1);
	  if (dy == 1) steps = min(steps, img.h-i);

	  int ind = i*w+j;
	  int end_ind = ind+steps*di;
	  int c = 0;
	  for (; ind != end_ind; ind += di) {
	    if (img.mask[ind]) c = img.mask[ind];
	    if (c) ret.mask[ind] = c;
	  }
	}
      }
    }
  }

  return ret;
}
Image diagonalGravity(Image_ img, int corner) {
    Image result = img;
    bool moved = true;

    while (moved) {
        moved = false;

        if (corner == 0) { // Top-left
            for (int i = 1; i < img.h; ++i) {
                for (int j = 1; j < img.w; ++j) {
                    if (result(i, j) != 0 && result(i - 1, j - 1) == 0) {
                        std::swap(result(i, j), result(i - 1, j - 1));
                        moved = true;
                    }
                }
            }
        } else if (corner == 1) { // Bottom-right
            for (int i = img.h - 2; i >= 0; --i) {
                for (int j = img.w - 2; j >= 0; --j) {
                    if (result(i, j) != 0 && result(i + 1, j + 1) == 0) {
                        std::swap(result(i, j), result(i + 1, j + 1));
                        moved = true;
                    }
                }
            }
        } else if (corner == 2) { // Top-right
            for (int i = 1; i < img.h; ++i) {
                for (int j = img.w - 2; j >= 0; --j) {
                    if (result(i, j) != 0 && result(i - 1, j + 1) == 0) {
                        std::swap(result(i, j), result(i - 1, j + 1));
                        moved = true;
                    }
                }
            }
        } else if (corner == 3) { // Bottom-left
            for (int i = img.h - 2; i >= 0; --i) {
                for (int j = 1; j < img.w; ++j) {
                    if (result(i, j) != 0 && result(i + 1, j - 1) == 0) {
                        std::swap(result(i, j), result(i + 1, j - 1));
                        moved = true;
                    }
                }
            }
        }
    }

    return result;
}


Image upscaleImage(const Image& img, int scaleFactor) {
    // Calculate new dimensions
    int newWidth = img.w * scaleFactor;
    int newHeight = img.h * scaleFactor;

    // Ensure the upscaled image fits within 30x30
    if (newWidth > 30 || newHeight > 30) {
        // Return original image if it cannot be upscaled within the 30x30 limit
        return img;
    }

    // Create a new image with the upscaled dimensions
    Image upscaled = core::empty({img.p.x, img.p.y}, {newWidth, newHeight});

    // Duplicate each pixel according to the scale factor
    for (int i = 0; i < img.h; ++i) {
        for (int j = 0; j < img.w; ++j) {
            char pixelValue = img(i, j);
            for (int dy = 0; dy < scaleFactor; ++dy) {
                for (int dx = 0; dx < scaleFactor; ++dx) {
                    upscaled(i * scaleFactor + dy, j * scaleFactor + dx) = pixelValue;
                }
            }
        }
    }

    return upscaled;
}

Image stretchImage(const Image& img, int stretchFactor, int direction) {
    int newWidth = img.w;
    int newHeight = img.h;

    // Determine new dimensions based on the direction of the stretch
    if (direction == 0) {
        newWidth = img.w * stretchFactor;
    } else if (direction == 1) {
        newHeight = img.h * stretchFactor;
    }

    // Ensure the stretched image fits within 30x30
    if (newWidth > 30 || newHeight > 30) {
        // Return the original image if the stretched version exceeds the limit
        return img;
    }

    // Create a new image with the stretched dimensions
    Image stretched = core::empty({img.p.x, img.p.y}, {newWidth, newHeight});

    // Populate the stretched image
    for (int i = 0; i < img.h; ++i) {
        for (int j = 0; j < img.w; ++j) {
            char pixelValue = img(i, j);

            // Fill the appropriate region in the stretched image
            if (direction == 0) {
                // Horizontally duplicate pixels
                for (int dx = 0; dx < stretchFactor; ++dx) {
                    stretched(i, j * stretchFactor + dx) = pixelValue;
                }
            } else if (direction == 1) {
                // Vertically duplicate pixels
                for (int dy = 0; dy < stretchFactor; ++dy) {
                    stretched(i * stretchFactor + dy, j) = pixelValue;
                }
            }
        }
    }

    return stretched;
}


Image zoomIn(const Image& img, int id) {
    int newWidth, newHeight;
    int startX = 0, startY = 0;
    Image result;

    // Determine the region based on the id
    switch (id) {
        case 1: // First quarter (top-left)
            newWidth = img.w / 2;
            newHeight = img.h / 2;
            startX = 0;
            startY = 0;
            break;
        case 2: // Second quarter (top-right)
            newWidth = img.w / 2;
            newHeight = img.h / 2;
            startX = img.w / 2;
            startY = 0;
            break;
        case 3: // Third quarter (bottom-left)
            newWidth = img.w / 2;
            newHeight = img.h / 2;
            startX = 0;
            startY = img.h / 2;
            break;
        case 4: // Fourth quarter (bottom-right)
            newWidth = img.w / 2;
            newHeight = img.h / 2;
            startX = img.w / 2;
            startY = img.h / 2;
            break;
        case 5: // Top half
            newWidth = img.w;
            newHeight = img.h / 2;
            startX = 0;
            startY = 0;
            break;
        case 6: // Bottom half
            newWidth = img.w;
            newHeight = img.h / 2;
            startX = 0;
            startY = img.h / 2;
            break;
        case 7: // Left half
            newWidth = img.w / 2;
            newHeight = img.h;
            startX = 0;
            startY = 0;
            break;
        case 8: // Right half
            newWidth = img.w / 2;
            newHeight = img.h;
            startX = img.w / 2;
            startY = 0;
            break;
        default:
            // If an unsupported ID is provided, return an empty image
            return core::empty(img.p, {0, 0});
    }

    // Create the result image and copy the specified region from the original image
    result = core::empty({startX, startY}, {newWidth, newHeight});
    for (int i = 0; i < newHeight; ++i) {
        for (int j = 0; j < newWidth; ++j) {
            result(i, j) = img(startY + i, startX + j);
        }
    }

    return result;
}



Image circularShift(Image_ img, bool isRowShift, int shiftAmount) {
    Image result = img;

    if (isRowShift) { // Shift rows circularly
        for (int i = 0; i < img.h; ++i) {
            for (int j = 0; j < img.w; ++j) {
                int newCol = (j + shiftAmount) % img.w;
                if (newCol < 0) newCol += img.w;
                result(i, newCol) = img(i, j);
            }
        }
    } else { // Shift columns circularly
        for (int j = 0; j < img.w; ++j) {
            for (int i = 0; i < img.h; ++i) {
                int newRow = (i + shiftAmount) % img.h;
                if (newRow < 0) newRow += img.h;
                result(newRow, j) = img(i, j);
            }
        }
    }

    return result;
}

Image mirror2(Image_ a, Image_ line) {
  Image ret;
  if (line.w > line.h) {
    ret = rigid(a,5);
    ret.x = a.x;
    ret.y = line.y*2+line.h-a.y-a.h;
  } else {
    ret = rigid(a,4);
    ret.y = a.y;
    ret.x = line.x*2+line.w-a.x-a.w;
  }
  return ret;
}

vImage gravity(Image_ in, int d) {
  vImage pieces = splitAll(in);
  Image room = hull0(in);
  const int dx = (d==0)-(d==1);
  const int dy = (d==2)-(d==3);

  vImage ret;
  Image out = room;
  sort(pieces.begin(), pieces.end(), [dx,dy](Image_ a, Image_ b) {
      return (a.x-b.x)*dx+(a.y-b.y)*dy > 0;});
  for (Image p : pieces) {
    while (1) {
      p.x += dx;
      p.y += dy;
      const int pxoutx = p.x-out.x;
      const int pyouty = p.y-out.y;
      for (int i = 0; i < p.h; ++i) {
          const int y = i+p.y-out.y;
	for (int j = 0; j < p.w; ++j) {
	  if (p(i,j) == 0) continue;
	  const int x = j+pxoutx;
	  if (x < 0 || y < 0 || x >= out.w || y >= out.h || out(y,x)) {
	    p.x -= dx;
	    p.y -= dy;
	    goto done;
	  }
	}
      }
    }
  done:
    ret.push_back(p);
    out = compose(out, p, 3);
  }
  return ret;
}



Image myStack(vImage_ lens, int id) {
  const int n = lens.size();
  if (!n) return badImg;
  vector<pair<int,int>> order(n);
  for (int i = 0; i < n; ++i) {
    order[i] = {lens[i].w*lens[i].h,i};
  }
  sort(order.begin(), order.end());

  Image out = lens[order[0].second];
  for (int i = 1; i < n; ++i)
    out = myStack(out,lens[order[i].second],id);
  return out;
}

Image stackLine(vImage_ shapes) {
  const int n = shapes.size();
  if (!n) return badImg;
  else if (n == 1) return shapes[0];
  vector<int> xs(n), ys(n);
  for (int i = 0; i < n; ++i) {
    xs[i] = shapes[i].x;
    ys[i] = shapes[i].y;
  }
  sort(xs.begin(), xs.end());
  sort(ys.begin(), ys.end());
  int xmin = 1e9, ymin = 1e9;
  for (int i = 1; i < n; ++i) {
    xmin = min(xmin, xs[i]-xs[i-1]);
    ymin = min(ymin, ys[i]-ys[i-1]);
  }
  int dx = 1, dy = 0;
  if (xmin < ymin) dx = 0, dy = 1;

  vector<pair<int,int>> order(n);
  for (int i = 0; i < shapes.size(); ++i) {
    order[i] = {shapes[i].x*dx+shapes[i].y*dy,i};
  }
  sort(order.begin(), order.end());

  Image out = shapes[order[0].second];
  
  for (int i = 1; i < n; ++i)
    out = myStack(out,shapes[order[i].second],dy);
  return out;
}


Image composeGrowingSlow(vImage_ imgs) {
  const int n = imgs.size();
  if (!n) return badImg;

  vector<pair<int,int>> order(n);
  for (int i = 0; i < n; ++i) {
    order[i] = {core::count(imgs[i]), i};
  }
  sort(order.rbegin(), order.rend());

  Image ret = imgs[order[0].second];
  for (int i = 1; i < n; ++i)
    ret = compose(ret, imgs[order[i].second], 0);
  return ret;
}

// Optimized radix sort for integer keys
void radixSort(std::vector<std::pair<int, int>>& arr) {
    int maxVal = 0;
    for (const auto& p : arr) {
        maxVal = std::max(maxVal, p.first);
    }

    int exp = 1;
    while (maxVal / exp > 0) {
        std::vector<std::vector<std::pair<int, int>>> buckets(10);
        for (const auto& p : arr) {
            buckets[(p.first / exp) % 10].push_back(p);
        }

        int index = 0;
        for (const auto& bucket : buckets) {
            for (const auto& p : bucket) {
                arr[index++] = p;
            }
        }

        exp *= 10;
    }
}

Image composeGrowing(vImage_ imgs) {
    const int n = imgs.size();
    if (!n) return badImg;
    if (n == 1) return imgs[0];

    int minx = 1e9, miny = 1e9, maxx = -1e9, maxy = -1e9;

    for (const Image_& img : imgs) {
        minx = std::min(minx, img.x);
        miny = std::min(miny, img.y);
        maxx = std::max(maxx, img.x + img.w);
        maxy = std::max(maxy, img.y + img.h);
    }

    point rsz = {maxx - minx, maxy - miny};
    if (std::max(rsz.x, rsz.y) > MAXSIDE || rsz.x * rsz.y > MAXAREA || rsz.x <= 0 || rsz.y <= 0)
        return badImg;

    std::vector<std::pair<int, int>> order(n);
    for (int i = 0; i < n; ++i) {
        order[i] = {core::count(imgs[i]), i};
    }

    // Radix sort for integer keys
    int maxVal = 0;
    for (const auto& p : order) {
        maxVal = std::max(maxVal, p.first);
    }

    int exp = 1;
    while (maxVal / exp > 0) {
        std::vector<std::vector<std::pair<int, int>>> buckets(10);
        for (const auto& p : order) {
            buckets[(p.first / exp) % 10].push_back(p);
        }

        int index = 0;
        for (const auto& bucket : buckets) {
            for (const auto& p : bucket) {
                order[index++] = p;
            }
        }

        exp *= 10;
    }

    Image ret = core::empty(point{minx, miny}, rsz);
    for (const auto& [cnt, imgi] : order) {
        const Image_& img = imgs[imgi];
        const int dx = img.x - ret.x;
        const int dy = img.y - ret.y;

        // Direct pixel access
        for (int i = 0; i < img.h; ++i) {
          const int idy = i + dy;
            for (int j = 0; j < img.w; ++j) {
                if (img(i, j)) {
                    ret(idy, j + dx) = img(i, j);
                }
            }
        }
    }

    return ret;
}


Image pickUnique(vImage_ imgs, int id) {
  assert(id == 0);

  const int n = imgs.size();
  if (!n) return badImg;

  //Pick the one with the unique color
  vector<int> mask(n);
  vector<int> cnt(10);
  for (int i = 0; i < n; ++i) {
    mask[i] = core::colMask(imgs[i]);
    for (int c = 0; c < 10; ++c) {
      if (mask[i]>>c&1) cnt[c]++;
    }
  }
  int reti = -1;
  for (int i = 0; i < n; ++i) {
    for (int c = 0; c < 10; ++c) {
      if (mask[i]>>c&1) {
	if (cnt[c] == 1) {
	  if (reti == -1) reti = i;
	  else return badImg;
	}
      }
    }
  }
  if (reti == -1) return badImg;
  return imgs[reti];
}

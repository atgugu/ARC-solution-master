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
#include <array>
#include <vector>
#include <algorithm>
#include <cassert>
#include <utility>
#include <random>
#include <omp.h>

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


Image extractAndApplyConstantShapes(Image_ img, const vector<pair<Image, Image>>& train) {
    const int max_size = 30;
    
    // Initialize the output image as a blank image with all zeros
    Image output = core::empty({img.p.x, img.p.y}, {max_size, max_size});

    // Temporary set to store constant shape positions
    unordered_set<int> constant_shape_positions;

    // Iterate over each training pair
    for (const auto& pair : train) {
        const Image& input = pair.first;
        const Image& target = pair.second;

        // Ensure both images are the same size
        if (input.w != target.w || input.h != target.h) {
            output = img;
            break; // Skip if sizes do not match
        }

        // Identify constant shapes by comparing pixels
        for (int i = 0; i < input.h; ++i) {
            for (int j = 0; j < input.w; ++j) {
                int pos = i * input.w + j;
                if (input(i, j) == target(i, j) && input(i, j) != 0) {
                    constant_shape_positions.insert(pos);
                }
            }
        }
    }

    // Apply constant shapes to the output image based on the input image
    for (int pos : constant_shape_positions) {
        int i = pos / img.w;
        int j = pos % img.w;

        // Only apply if within bounds
        if (i < img.h && j < img.w && i < output.h && j < output.w) {
            output(i, j) = img(i, j); // Transfer the constant shape pixel from input to output
        }
    }

    return output;
}

Image eraseCol(Image img, int col) {
  for (int i = 0; i < img.h; ++i)
    for (int j = 0; j < img.w; ++j)
      if (img(i,j) == col) img(i,j) = 0;
  return img;
}

Image removeGrid(const Image& img) {
    if (img.w <= 1 || img.h <= 1) {
        return badImg; 
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
                return badImg;  // Irregular spacing indicates no grid
            }
        }
    } else {
        rowGridLines.clear();
    }

    if (colGridLines.size() > 1) {
        int colSpacing = colGridLines[1] - colGridLines[0];
        for (size_t j = 1; j < colGridLines.size(); ++j) {
            if (colGridLines[j] - colGridLines[j - 1] != colSpacing) {
                return badImg;  // Irregular spacing indicates no grid
            }
        }
    } else {
        colGridLines.clear();
    }

    // If no regular grid was found, return the original image
    if (rowGridLines.empty() && colGridLines.empty()) {
        return badImg;
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
        return badImg;
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
    return similarityRatio > 0.8;  // Consider morphologically similar if 95% or more pixels match
}

Image downscaleImage(const Image& img, int factor) {
    // Determine the new dimensions
    int newWidth = img.w / factor;
    int newHeight = img.h / factor;

    // Ensure the dimensions remain at least 1x1
    if (newWidth < 1 || newHeight < 1) {
        return badImg;  // Downscale too aggressive, return original image
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
        return badImg;  // If morphology changes significantly, return the original image
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
    // Calculate the number of full cells in each direction
    int numRows = (img.h + cellHeight - 1) / cellHeight;  // Round up to cover entire height
    int numCols = (img.w + cellWidth - 1) / cellWidth;    // Round up to cover entire width

    // Initialize the result image with each cell representing the mode color of that grid
    Image result = {{0, 0}, {numCols, numRows}, std::vector<char>(numRows * numCols, 0)};

    // Process each cell to determine the mode color
    for (int row = 0; row < numRows; ++row) {
        for (int col = 0; col < numCols; ++col) {
            // Count colors in the current cell
            std::unordered_map<char, int> colorCount;
            for (int i = 0; i < cellHeight; ++i) {
                for (int j = 0; j < cellWidth; ++j) {
                    int x = row * cellHeight + i;
                    int y = col * cellWidth + j;
                    
                    // Check if we're within the image bounds
                    if (x < img.h && y < img.w) {
                        char color = img(x, y);
                        colorCount[color]++;
                    }
                }
            }

            // Determine the mode color in this cell
            char modeColor = 0;
            int maxCount = 0;
            for (const auto& [color, count] : colorCount) {
                if (count > maxCount) {
                    maxCount = count;
                    modeColor = color;
                }
            }
            
            // Assign the mode color to the corresponding cell in the result image
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

Image compressPatches(Image_ img, int patchSize = 2) {
    const int height = img.h, width = img.w;
    if (width * height <= 0) return badImg;

    std::vector<std::pair<int, int>> patches;

    // Identify unique patches by comparing pixels within each patch
    for (int i = 0; i < height; i += patchSize) {
        for (int j = 0; j < width; j += patchSize) {
            bool unique = false;
            char referenceColor = img(i, j);

            // Check if all pixels in this patch have the same color
            for (int di = 0; di < patchSize && i + di < height; ++di) {
                for (int dj = 0; dj < patchSize && j + dj < width; ++dj) {
                    if (img(i + di, j + dj) != referenceColor) {
                        unique = true;
                        break;
                    }
                }
                if (unique) break;
            }

            if (unique) {
                patches.emplace_back(i, j);  // Save the starting point of the unique patch
            }
        }
    }

    // Check if any unique patches were identified
    if (patches.empty()) return badImg;

    // Calculate the dimensions of the compressed image
    const int compressedHeight = patches.size() * patchSize;
    const int compressedWidth = patchSize;

    // Ensure compressed dimensions are valid
    if (compressedHeight <= 0 || compressedWidth <= 0) return badImg;

    // Initialize the compressed image
    Image ret = core::empty(point{compressedWidth, compressedHeight});

    // Populate the compressed image with unique patches
    for (size_t i = 0; i < patches.size(); ++i) {
        int startRow = patches[i].first;
        int startCol = patches[i].second;

        for (int di = 0; di < patchSize && (i * patchSize + di) < compressedHeight; ++di) {
            for (int dj = 0; dj < patchSize && dj < compressedWidth; ++dj) {
                if (startRow + di < height && startCol + dj < width) {
                    ret(i * patchSize + di, dj) = img(startRow + di, startCol + dj);
                }
            }
        }
    }

    return ret;
}


Image compressSymmetry(Image_ img) {
    const int height = img.h, width = img.w;
    if (width * height <= 0) return badImg;

    int centerX = width / 2, centerY = height / 2;
    Image ret = core::empty(point{centerX, centerY});

    for (int i = 0; i < centerY; ++i) {
        for (int j = 0; j < centerX; ++j) {
            ret(i, j) = (img(i, j) == img(height - i - 1, j) &&
                         img(i, j) == img(i, width - j - 1)) ? img(i, j) : 0;
        }
    }
    return ret;
}


Image compressHV(Image_ img, bool compressHorizontal) {
    const int height = img.h, width = img.w;
    if (width * height <= 0) return badImg;

    // Mark rows or columns with unique pixels
    std::vector<int> markers(compressHorizontal ? height : width, 0);
    markers[0] = 1;

    for (int i = 1; i < (compressHorizontal ? height : width); ++i) {
        for (int j = 1; j < (compressHorizontal ? width : height); ++j) {
            if ((compressHorizontal && img(i, j) != img(i - 1, j)) ||
                (!compressHorizontal && img(j, i) != img(j, i - 1))) {
                markers[i] = 1;
            }
        }
    }

    // Collect indices to be included in the compressed image
    std::vector<int> indices;
    indices.reserve(compressHorizontal ? height : width);
    for (int i = 0; i < (compressHorizontal ? height : width); ++i) {
        if (markers[i]) indices.push_back(i);
    }

    const int newSize = indices.size();
    Image ret = core::empty(compressHorizontal ? point{width, newSize} : point{newSize, height});

    // Fill the compressed image based on the selected rows or columns
    if (compressHorizontal) {
        for (int i = 0; i < newSize; ++i) {
            int rowIdx = indices[i];
            for (int j = 0; j < width; ++j) {
                ret(i, j) = img(rowIdx, j);
            }
        }
    } else {
        for (int i = 0; i < newSize; ++i) {
            int colIdx = indices[i];
            for (int j = 0; j < height; ++j) {
                ret(j, i) = img(j, colIdx);
            }
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
    // #pragma omp parallel for
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


Image connectFarthestShapes(Image_ img) {
    Image ret = img;  // Start with the original image
    std::unordered_map<char, std::vector<point>> colorPositions;

    // Collect all non-zero pixels by color
    for (int i = 0; i < img.h; ++i) {
        for (int j = 0; j < img.w; ++j) {
            if (img(i, j) != 0) {
                colorPositions[img(i, j)].push_back({i, j});
            }
        }
    }

    // For each color, find and connect the farthest pair of pixels
    for (const auto& [color, positions] : colorPositions) {
        if (positions.size() < 2) continue;  // Skip if there are fewer than two points for this color

        // Initialize variables to store the farthest pair of points
        int maxDistance = 0;
        point p1, p2;

        // Find the farthest pair of points for this color
        for (size_t a = 0; a < positions.size(); ++a) {
            for (size_t b = a + 1; b < positions.size(); ++b) {
                int distance = abs(positions[a].x - positions[b].x) + abs(positions[a].y - positions[b].y);
                if (distance > maxDistance) {
                    maxDistance = distance;
                    p1 = positions[a];
                    p2 = positions[b];
                }
            }
        }

        // Connect the farthest points (p1 and p2) using a Manhattan path
        int x = p1.x, y = p1.y;
        while (x != p2.x || y != p2.y) {
            ret(x, y) = color;
            if (x < p2.x) ++x;
            else if (x > p2.x) --x;
            else if (y < p2.y) ++y;
            else --y;
        }
        ret(x, y) = color;  // Set the endpoint
    }

    return ret;
}


Image connectNearestShapes(Image_ img) {
    Image ret = img;  // Start with the original image
    std::unordered_map<char, std::vector<point>> colorPositions;

    // Collect all non-zero pixels by color
    for (int i = 0; i < img.h; ++i) {
        for (int j = 0; j < img.w; ++j) {
            if (img(i, j) != 0) {
                colorPositions[img(i, j)].push_back({i, j});
            }
        }
    }

    // For each color, connect each pair of pixels with the shortest path
    for (const auto& [color, positions] : colorPositions) {
        for (size_t a = 0; a < positions.size(); ++a) {
            for (size_t b = a + 1; b < positions.size(); ++b) {
                point p1 = positions[a];
                point p2 = positions[b];

                // Connect p1 to p2 using a basic Manhattan path
                int x = p1.x, y = p1.y;
                while (x != p2.x || y != p2.y) {
                    ret(x, y) = color;
                    if (x < p2.x) ++x;
                    else if (x > p2.x) --x;
                    else if (y < p2.y) ++y;
                    else --y;
                }
                ret(x, y) = color;  // Set the endpoint
            }
        }
    }

    return ret;
}

Image fillIslands(Image_ img) {
    Image ret = img;
    std::vector<std::vector<bool>> visited(img.h, std::vector<bool>(img.w, false));

    auto floodFill = [&](int x, int y, char color) {
        std::queue<point> q;
        q.push({x, y});
        visited[x][y] = true;

        while (!q.empty()) {
            point p = q.front();
            q.pop();

            for (const auto& [dx, dy] : std::vector<point>{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}) {
                int nx = p.x + dx, ny = p.y + dy;
                if (nx >= 0 && nx < img.h && ny >= 0 && ny < img.w && !visited[nx][ny] && img(nx, ny) == color) {
                    ret(nx, ny) = color;
                    visited[nx][ny] = true;
                    q.push({nx, ny});
                }
            }
        }
    };

    // Apply flood fill to each unvisited region
    for (int i = 0; i < img.h; ++i) {
        for (int j = 0; j < img.w; ++j) {
            if (!visited[i][j] && img(i, j) != 0) {
                floodFill(i, j, img(i, j));
            }
        }
    }

    return ret;
}


Image diagonalBridge(Image_ img) {
    Image ret = img;

    // Diagonal bridging (top-left to bottom-right)
    for (int i = 0; i < img.h - 1; ++i) {
        for (int j = 0; j < img.w - 1; ++j) {
            if (img(i, j) != 0 && img(i + 1, j + 1) == img(i, j)) {
                ret(i + 1, j) = img(i, j);
                ret(i, j + 1) = img(i, j);
            }
        }
    }

    // Diagonal bridging (top-right to bottom-left)
    for (int i = 0; i < img.h - 1; ++i) {
        for (int j = 1; j < img.w; ++j) {
            if (img(i, j) != 0 && img(i + 1, j - 1) == img(i, j)) {
                ret(i + 1, j) = img(i, j);
                ret(i, j - 1) = img(i, j);
            }
        }
    }

    return ret;
}


Image outlineShapes(Image_ img) {
    Image ret = core::empty(img.p, img.sz);

    // Directions for 8-neighbor connectivity
    std::vector<point> directions = {{-1, -1}, {-1, 0}, {-1, 1},
                                     {0, -1},         {0, 1},
                                     {1, -1}, {1, 0}, {1, 1}};

    // For each pixel, check if it's on the edge of a shape
    for (int i = 0; i < img.h; ++i) {
        for (int j = 0; j < img.w; ++j) {
            if (img(i, j) != 0) {
                bool isEdge = false;
                for (const auto& dir : directions) {
                    int ni = i + dir.x, nj = j + dir.y;
                    if (ni >= 0 && ni < img.h && nj >= 0 && nj < img.w && img(ni, nj) == 0) {
                        isEdge = true;
                        break;
                    }
                }
                if (isEdge) {
                    ret(i, j) = img(i, j);  // Set edge pixel in output image
                }
            }
        }
    }

    return ret;
}

Image bridgeGaps(Image_ img) {
    Image ret = core::empty(img.p, img.sz);

    // Horizontal bridging
    for (int i = 0; i < img.h; ++i) {
        int last = -1, lastc = -1;
        for (int j = 0; j < img.w; ++j) {
            if (img(i, j)) {
                if (img(i, j) == lastc) {
                    for (int k = last + 1; k < j; ++k) {
                        ret(i, k) = lastc;
                    }
                }
                lastc = img(i, j);
                last = j;
            }
        }
    }

    // Vertical bridging
    for (int j = 0; j < img.w; ++j) {
        int last = -1, lastc = -1;
        for (int i = 0; i < img.h; ++i) {
            if (img(i, j)) {
                if (img(i, j) == lastc) {
                    for (int k = last + 1; k < i; ++k) {
                        ret(k, j) = lastc;
                    }
                }
                lastc = img(i, j);
                last = i;
            }
        }
    }

    return ret;
}


Image connect(Image_ img, int id) {
  assert(id >= 0 && id < 5); // Allow diagonal IDs
  Image ret = core::empty(img.p, img.sz);

  // Horizontal
  if (id == 0 || id == 2) {
    for (int i = 0; i < img.h; ++i) {
      int last = -1, lastc = -1;
      for (int j = 0; j < img.w; ++j) {
        if (img(i, j)) {
          if (img(i, j) == lastc) {
            for (int k = last + 1; k < j; ++k)
              ret(i, k) = lastc;
          }
          lastc = img(i, j);
          last = j;
        }
      }
    }
  }

  // Vertical
  if (id == 1 || id == 2) {
    for (int j = 0; j < img.w; ++j) {
      int last = -1, lastc = -1;
      for (int i = 0; i < img.h; ++i) {
        if (img(i, j)) {
          if (img(i, j) == lastc) {
            for (int k = last + 1; k < i; ++k)
              ret(k, j) = lastc;
          }
          lastc = img(i, j);
          last = i;
        }
      }
    }
  }

  // Diagonal (top-left to bottom-right)
  if (id == 3) {
    for (int i = 0; i < img.h; ++i) {
      int lastRow = -1, lastCol = -1, lastc = -1;
      for (int j = 0; j < img.w; ++j) {
        if (i + j < img.h && img(i + j, j)) {
          if (img(i + j, j) == lastc) {
            for (int k = 1; k < i + j - lastRow; ++k)
              ret(lastRow + k, lastCol + k) = lastc;
          }
          lastc = img(i + j, j);
          lastRow = i + j;
          lastCol = j;
        }
      }
    }
  }

  // Diagonal (top-right to bottom-left)
  if (id == 4) {
    for (int i = 0; i < img.h; ++i) {
      int lastRow = -1, lastCol = img.w, lastc = -1;
      for (int j = img.w - 1; j >= 0; --j) {
        if (i + (img.w - 1 - j) < img.h && img(i + (img.w - 1 - j), j)) {
          if (img(i + (img.w - 1 - j), j) == lastc) {
            for (int k = 1; k < i + (img.w - 1 - j) - lastRow; ++k)
              ret(lastRow + k, lastCol - k) = lastc;
          }
          lastc = img(i + (img.w - 1 - j), j);
          lastRow = i + (img.w - 1 - j);
          lastCol = j;
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
  //hamsie
//   #pragma omp parallel for
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
  //hamsie
  #pragma omp parallel for
  for (int k : {0,1}) {
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
    if (img.h < size || img.w < size || size < 3) return badImg;

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
    if (img.h < size) return badImg;

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
    if (img.w < size) return badImg;

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
        return badImg;
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
        return badImg;
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
        return badImg;
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
    if (img.h < size || img.w < size) return badImg;

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
    if (img.h < size || img.w < size) return badImg;

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
        return badImg;
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
    if(img.h < size || img.w < size) return badImg;

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

    return move(components);
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

Image largestShape(Image_ img) {
    auto shapes = extractConnectedComponents(img);
    if (shapes.empty()) {
        return badImg;  // Return an empty image if no shapes are found
    }

    // Initialize a placeholder for the largest shape
    Image largestShape = shapes[0];
    int maxSize = largestShape.mask.size();

    // Iterate through each shape to find the largest by mask size
    for (const auto& shape : shapes) {
        int size = shape.mask.size();
        if (size > maxSize) {
            maxSize = size;
            largestShape = shape;  // Assign the current shape as the largest
        }
    }

    return largestShape;
}

Image mostCommonColorShape(Image_ img) {
    // Count occurrences of each color in the entire image
    std::unordered_map<char, int> colorCount;
    for (int i = 0; i < img.h; ++i) {
        for (int j = 0; j < img.w; ++j) {
            if (img(i, j) != 0) {  // Exclude zero as it's usually background
                colorCount[img(i, j)]++;
            }
        }
    }

    // Check if there are any non-background colors
    if (colorCount.empty()) {
        return badImg;
    }

    // Determine the most common color
    char mostCommonColor = 0;
    int maxCount = 0;
    for (const auto& [color, count] : colorCount) {
        if (count > maxCount) {
            maxCount = count;
            mostCommonColor = color;
        }
    }

    // Extract connected components (shapes)
    auto shapes = extractConnectedComponents(img);

    // Find a shape with the most common color
    for (const auto& shape : shapes) {
        // Calculate the dominant color in this shape
        std::unordered_map<char, int> shapeColorCount;
        for (int i = 0; i < shape.h; ++i) {
            for (int j = 0; j < shape.w; ++j) {
                char pixel = shape(i, j);
                if (pixel != 0) {  // Only count non-zero pixels
                    shapeColorCount[pixel]++;
                }
            }
        }

        // Check if the dominant color in this shape matches the most common color
        char dominantShapeColor = 0;
        int maxShapeColorCount = 0;
        for (const auto& [color, count] : shapeColorCount) {
            if (count > maxShapeColorCount) {
                maxShapeColorCount = count;
                dominantShapeColor = color;
            }
        }

        if (dominantShapeColor == mostCommonColor) {
            return shape;  // Return this shape if its dominant color matches
        }
    }

    return badImg;
}


Image smallestShape(Image_ img) {
    auto shapes = extractConnectedComponents(img);
    if (shapes.empty()) {
        return badImg;  // Return an empty image if no shapes are found
    }

    // Initialize with the first non-zero size shape
    Image smallestShape = shapes[0];
    int minSize = std::numeric_limits<int>::max();
    bool foundNonZeroSize = false;

    for (const auto& shape : shapes) {
        int size = shape.mask.size();
        if (size > 0 && size < minSize) {
            minSize = size;
            smallestShape = shape;
            foundNonZeroSize = true;
        }
    }

    // Return smallest shape only if a valid non-zero size shape was found
    return foundNonZeroSize ? smallestShape : img;
}

Image extractShape(const Image& img, int startX, int startY, int regionSize) {
    // Calculate the width and height for the extracted shape, ensuring it doesn't exceed image bounds
    int actualWidth = std::min(regionSize, img.w - startY);
    int actualHeight = std::min(regionSize, img.h - startX);

    // Initialize the new shape image with the calculated dimensions
    Image shape = {{startX, startY}, {actualWidth, actualHeight}, std::vector<char>(actualWidth * actualHeight, 0)};

    // Copy pixels from the original image to the new shape
    for (int i = 0; i < actualHeight; ++i) {
        for (int j = 0; j < actualWidth; ++j) {
            shape(i, j) = img(startX + i, startY + j);
        }
    }

    return shape;
}

Image denseRegionShape(Image_ img, int regionSize = 3) {
    int maxDensity = 0;
    Image densestShape;
    bool foundDenseRegion = false;

    for (int i = 0; i <= img.h - regionSize; ++i) {
        for (int j = 0; j <= img.w - regionSize; ++j) {
            int density = 0;
            for (int di = 0; di < regionSize; ++di) {
                for (int dj = 0; dj < regionSize; ++dj) {
                    if (img(i + di, j + dj) != 0) density++;
                }
            }

            if (density > maxDensity) {
                maxDensity = density;
                densestShape = extractShape(img, i, j, regionSize);
                foundDenseRegion = true;
            }
        }
    }

    return foundDenseRegion ? densestShape : badImg;
}

bool isEnclosed(const Image& shape, const std::vector<Image>& shapes) {
    // Define the bounding box of the shape
    int minX = shape.p.x;
    int minY = shape.p.y;
    int maxX = minX + shape.h - 1;
    int maxY = minY + shape.w - 1;

    // Check each pixel on the boundary of the shape's bounding box
    for (int i = minX; i <= maxX; ++i) {
        for (int j = minY; j <= maxY; ++j) {
            // Skip non-boundary pixels
            if ((i == minX || i == maxX || j == minY || j == maxY) && shape(i - minX, j - minY) != 0) {
                bool surrounded = false;
                // Check if any other shape surrounds this boundary pixel
                for (const auto& otherShape : shapes) {
                    if (&otherShape != &shape && otherShape.safe(i, j) != 0) {  // Check if another shape overlaps this point
                        surrounded = true;
                        break;
                    }
                }
                if (!surrounded) {
                    return false;  // If any boundary pixel is not surrounded, the shape is not enclosed
                }
            }
        }
    }

    return true;  // All boundary pixels are surrounded by other shapes
}

Image enclosedShapes(Image_ img) {
    auto shapes = extractConnectedComponents(img);
    if (shapes.empty()) {
        return badImg;
    }

    std::vector<Image> enclosedShapes;
    for (const auto& shape : shapes) {
        if (isEnclosed(shape, shapes)) {  // Check if `shape` is fully enclosed by another shape
            enclosedShapes.push_back(shape);
        }
    }

    if (enclosedShapes.empty()) {
        return badImg;
    }

    return combineShapes(enclosedShapes, img.h, img.w);
}
bool isSymmetric(const Image& shape) {
    int height = shape.h;
    int width = shape.w;

    // Check for horizontal symmetry (top and bottom halves mirror each other)
    for (int i = 0; i < height / 2; ++i) {
        for (int j = 0; j < width; ++j) {
            if (shape(i, j) != shape(height - 1 - i, j)) {
                return false;  // Not symmetric horizontally
            }
        }
    }

    // Check for vertical symmetry (left and right halves mirror each other)
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width / 2; ++j) {
            if (shape(i, j) != shape(i, width - 1 - j)) {
                return false;  // Not symmetric vertically
            }
        }
    }

    return true;  // Shape is symmetric both horizontally and vertically
}


Image symmetricShape(Image_ img) {
    auto shapes = extractConnectedComponents(img);
    if (shapes.empty()) {
        std::cout << "No shapes found in the image." << std::endl;
        return badImg;
    }

    Image largestSymmetricShape;
    int maxSize = 0;

    for (const auto& shape : shapes) {
        if (isSymmetric(shape)) {  // Assume `isSymmetric` checks symmetry
            int size = shape.mask.size();
            if (size > maxSize) {
                maxSize = size;
                largestSymmetricShape = shape;
            }
        }
    }

    return maxSize > 0 ? largestSymmetricShape : badImg;
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

Image morphOpening(Image_ img) {
    return dilate(erode(img, 1), 1);
}

Image morphClosing(Image_ img) {
    return erode(dilate(img, 1), 1);
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
  return move(ret);
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
  return move(ret);
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


Image shiftImage(Image_ img, int direction, int amount, bool expandGrid=false) {
    int new_w = img.w, new_h = img.h;
    int dx = 0, dy = 0;

    // Define direction vectors
    switch (direction) {
        case 0: dy = -amount; break; // Up
        case 1: dy = amount; break;  // Down
        case 2: dx = -amount; break; // Left
        case 3: dx = amount; break;  // Right
        case 4: dx = amount; dy = -amount; break; // Diagonal Up-Right
        case 5: dx = amount; dy = amount; break;  // Diagonal Down-Right
        case 6: dx = -amount; dy = -amount; break; // Anti-Diagonal Up-Left
        case 7: dx = -amount; dy = amount; break;  // Anti-Diagonal Down-Left
        default: return badImg; // Invalid direction
    }

    // If expanding the grid, calculate new dimensions
    if (expandGrid) {
        new_w = std::min(MAXSIDE, img.w + abs(dx));
        new_h = std::min(MAXSIDE, img.h + abs(dy));
    }

    // Create a new image with updated dimensions and initialize with 0
    Image ret = core::empty({img.p.x, img.p.y}, {new_w, new_h});

    // Iterate over the original image and shift pixels
    for (int i = 0; i < img.h; ++i) {
        for (int j = 0; j < img.w; ++j) {
            int ni = i + dy;
            int nj = j + dx;

            // Check if the shifted position is within the bounds of the new image
            if (ni >= 0 && ni < ret.h && nj >= 0 && nj < ret.w) {
                ret(ni, nj) = img(i, j);
            }
        }
    }

    return ret;
}






Image padImage(Image_ img, int new_w, int new_h) {
    Image padded_img;
    padded_img.w = new_w;
    padded_img.h = new_h;
    padded_img.mask = vector<char>(new_w * new_h, 0);

    int offset_x = (new_w - img.w) / 2;
    int offset_y = (new_h - img.h) / 2;

    for (int i = 0; i < img.h; ++i) {
        for (int j = 0; j < img.w; ++j) {
            int ni = i + offset_y;
            int nj = j + offset_x;
            if (ni >= 0 && ni < new_h && nj >= 0 && nj < new_w) {
                padded_img(ni, nj) = img(i, j);
            }
        }
    }

    return padded_img;
}

Image rotateImage90(Image_ img) {
    Image rotated_img;
    rotated_img.w = img.h;
    rotated_img.h = img.w;
    rotated_img.mask = vector<char>(rotated_img.w * rotated_img.h, 0);

    for (int i = 0; i < img.h; ++i) {
        for (int j = 0; j < img.w; ++j) {
            rotated_img(j, rotated_img.w - 1 - i) = img(i, j);
        }
    }

    return rotated_img;
}

Image flipImageHorizontal(Image_ img) {
    Image flipped_img = img;

    for (int i = 0; i < img.h; ++i) {
        for (int j = 0; j < img.w / 2; ++j) {
            swap(flipped_img(i, j), flipped_img(i, img.w - 1 - j));
        }
    }

    return flipped_img;
}

Image flipImageVertical(Image_ img) {
    Image flipped_img = img;

    for (int i = 0; i < img.h / 2; ++i) {
        for (int j = 0; j < img.w; ++j) {
            swap(flipped_img(i, j), flipped_img(img.h - 1 - i, j));
        }
    }

    return flipped_img;
}

Image flipImageDiagonal(Image_ img) {
    Image flipped_img;
    flipped_img.w = img.h;
    flipped_img.h = img.w;
    flipped_img.mask = vector<char>(flipped_img.w * flipped_img.h, 0);

    for (int i = 0; i < img.h; ++i) {
        for (int j = 0; j < img.w; ++j) {
            flipped_img(j, i) = img(i, j);
        }
    }

    return flipped_img;
}

Image flipImageAntiDiagonal(Image_ img) {
    Image flipped_img;
    flipped_img.w = img.h;
    flipped_img.h = img.w;
    flipped_img.mask = vector<char>(flipped_img.w * flipped_img.h, 0);

    for (int i = 0; i < img.h; ++i) {
        for (int j = 0; j < img.w; ++j) {
            flipped_img(flipped_img.h - 1 - j, flipped_img.w - 1 - i) = img(i, j);
        }
    }

    return flipped_img;
}

int computeHammingDistance(Image_ img1, Image_ img2) {
    int distance = 0;
    for (int i = 0; i < img1.h * img1.w; ++i) {
        if (img1.mask[i] != img2.mask[i]) {
            ++distance;
        }
    }
    return distance;
}

Image trainAndPredict(Image_ img, const vector<pair<Image, Image>>& train) {
    const int max_size = 30;

    // Prepare augmented training data
    vector<pair<Image, Image>> augmented_train;

    for (const auto& pair : train) {
        Image input = padImage(pair.first, max_size, max_size);
        Image target = padImage(pair.second, max_size, max_size);

        // Add original pair
        augmented_train.push_back({input, target});

        // Rotations (90, 180, 270 degrees)
        for (int r = 0; r < 3; ++r) {
            input = rotateImage90(input);
            target = rotateImage90(target);
            augmented_train.push_back({input, target});
        }

        // Flips: horizontal, vertical, diagonal, anti-diagonal
        vector<Image> flips_input = {
            flipImageHorizontal(input),
            flipImageVertical(input),
            flipImageDiagonal(input),
            flipImageAntiDiagonal(input)
        };
        vector<Image> flips_target = {
            flipImageHorizontal(target),
            flipImageVertical(target),
            flipImageDiagonal(target),
            flipImageAntiDiagonal(target)
        };

        for (int f = 0; f < flips_input.size(); ++f) {
            Image flip_input = flips_input[f];
            Image flip_target = flips_target[f];

            augmented_train.push_back({flip_input, flip_target});

            // Rotations of flips
            for (int r = 0; r < 3; ++r) {
                flip_input = rotateImage90(flip_input);
                flip_target = rotateImage90(flip_target);
                augmented_train.push_back({flip_input, flip_target});
            }
        }
    }

    // Initialize model weights (simulate with a mapping)
    vector<char> model_weights(max_size * max_size, 0);

    int epochs = 3;
    bool perfect_match = false;

    for (int epoch = 0; epoch < epochs && !perfect_match; ++epoch) {
        perfect_match = true;

        for (const auto& pair : augmented_train) {
            const Image& input = pair.first;
            const Image& target = pair.second;

            // Apply model (for simplicity, we use current model_weights)
            Image prediction;
            prediction.w = max_size;
            prediction.h = max_size;
            prediction.mask = vector<char>(max_size * max_size, 0);

            for (int i = 0; i < max_size * max_size; ++i) {
                prediction.mask[i] = model_weights[i];
            }

            // Compute loss (Hamming distance)
            int loss = 0;
            for (int i = 0; i < max_size * max_size; ++i) {
                if (prediction.mask[i] != target.mask[i]) {
                    ++loss;
                    perfect_match = false;
                    // Update model_weights (simple approach)
                    model_weights[i] = target.mask[i];
                }
            }
        }
    }

    // Apply the model to the input image
    Image padded_input = padImage(img, max_size, max_size);

    Image output = core::empty({img.p.x, img.p.y}, {max_size, max_size});

    output.w = max_size;
    output.h = max_size;
    output.mask = vector<char>(max_size * max_size, 0);

    for (int i = 0; i < max_size * max_size; ++i) {
        output.mask[i] = model_weights[i];
    }

    return output;
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


Image rearrangeShapes(Image_ img, int arrangementType) {
    auto shapes = extractConnectedComponents(img);  // Extract all shapes

    if (shapes.empty()) {
        std::cout << "No shapes found in the image." << std::endl;
        return badImg;
    }

    // Sort and arrange shapes based on arrangementType
    switch (arrangementType) {
        case 0:  // Arrange shapes by size (smallest to largest) in a row
            std::sort(shapes.begin(), shapes.end(), [](const Image& a, const Image& b) {
                return a.mask.size() < b.mask.size();
            });
            break;

        case 1:  // Arrange shapes by color in columns
            std::sort(shapes.begin(), shapes.end(), [](const Image& a, const Image& b) {
                char colorA = 0, colorB = 0;
                for (char c : a.mask) if (c != 0) { colorA = c; break; }
                for (char c : b.mask) if (c != 0) { colorB = c; break; }
                return colorA < colorB;
            });
            break;

        case 2:  // Arrange shapes in a grid based on appearance order
            // Already in original order, so no sorting required
            break;

        case 3:  // Arrange shapes in a circular pattern
        case 4:  // Arrange shapes in a square pattern
        case 5:  // Arrange shapes in a triangle pattern
            // No sorting required for circular, square, and triangle patterns
            break;

        default:
            std::cout << "Invalid arrangement type. Please use 0 to 5." << std::endl;
            return badImg;
    }

    int maxWidth = img.w, maxHeight = img.h;
    Image result = core::empty(point{0, 0}, point{maxWidth, maxHeight});

    int x = 0, y = 0;
    int centerX = maxWidth / 2;
    int centerY = maxHeight / 2;
    int radius = std::min(maxWidth, maxHeight) / 4;  // Radius for circle arrangement

    for (size_t index = 0; index < shapes.size(); ++index) {
        const auto& shape = shapes[index];

        // Determine position based on arrangement type
        if (arrangementType == 0) {  // Row arrangement
            if (x + shape.w > maxWidth) {
                x = 0;
                y += shape.h + 1;
            }
        } else if (arrangementType == 1) {  // Column arrangement
            if (y + shape.h > maxHeight) {
                y = 0;
                x += shape.w + 1;
            }
        } else if (arrangementType == 2) {  // Grid arrangement
            int gridSize = std::ceil(std::sqrt(shapes.size()));
            x = (shape.w + 1) * (index % gridSize);
            y = (shape.h + 1) * (index / gridSize);
        } else if (arrangementType == 3) {  // Circle arrangement
            double angle = 2 * M_PI * index / shapes.size();  // Equal angle spacing
            x = centerX + static_cast<int>(radius * cos(angle)) - shape.w / 2;
            y = centerY + static_cast<int>(radius * sin(angle)) - shape.h / 2;
        } else if (arrangementType == 4) {  // Square arrangement
            int side = std::ceil(std::sqrt(shapes.size()));  // Side length for the square layout
            x = centerX + (index % side - side / 2) * (shape.w + 1);
            y = centerY + (index / side - side / 2) * (shape.h + 1);
        } else if (arrangementType == 5) {  // Triangle arrangement
            int row = static_cast<int>(std::floor(std::sqrt(2 * index + 0.25) - 0.5));
            int posInRow = index - (row * (row + 1)) / 2;
            x = centerX + (posInRow - row / 2.0) * (shape.w + 1);
            y = centerY + row * (shape.h + 1);
        }

        // Ensure position is within bounds
        x = std::max(0, std::min(x, maxWidth - shape.w));
        y = std::max(0, std::min(y, maxHeight - shape.h));

        // Place the shape in the result image
        for (int i = 0; i < shape.h; ++i) {
            for (int j = 0; j < shape.w; ++j) {
                if (shape(i, j) != 0) {  // Place only non-zero pixels
                    result(y + i, x + j) = shape(i, j);
                }
            }
        }

        // Update x and y positions based on the arrangement type
        if (arrangementType == 0) {
            x += shape.w + 1;
        } else if (arrangementType == 1) {
            y += shape.h + 1;
        }
    }

    return result;
}

Image resizeShape(const Image& shape, int targetSize) {
    int originalSize = shape.mask.size();
    int newHeight = std::sqrt(targetSize);
    int newWidth = (targetSize + newHeight - 1) / newHeight;  // Ensure enough width to cover the target size

    // Create a resized mask by copying pixels proportionally
    Image resizedShape = core::empty(shape.p, {newWidth, newHeight});
    int index = 0;

    for (int i = 0; i < newHeight; ++i) {
        for (int j = 0; j < newWidth; ++j) {
            if (index < targetSize && index < originalSize && shape.mask[index] != 0) {
                resizedShape(i, j) = shape.mask[index];
            } else {
                resizedShape(i, j) = 0;  // Fill extra space with background (0)
            }
            index++;
        }
    }

    return resizedShape;
}

void rotateSquare(Image& img, int startX, int startY, int squareSize, int times) {
    // Ensure the square region is within the bounds of the image
    if (startX + squareSize > img.h || startY + squareSize > img.w) {
        std::cerr << "Square region goes out of image bounds." << std::endl;
        return;
    }

    for (int time = 0; time < times; ++time)
    for (int i = 0; i < squareSize / 2; ++i) {
        for (int j = i; j < squareSize - i - 1; ++j) {
            // Save the top element
            int temp = img(startX + i, startY + j);

            // Move left to top
            img(startX + i, startY + j) = img(startX + squareSize - 1 - j, startY + i);

            // Move bottom to left
            img(startX + squareSize - 1 - j, startY + i) = img(startX + squareSize - 1 - i, startY + squareSize - 1 - j);

            // Move right to bottom
            img(startX + squareSize - 1 - i, startY + squareSize - 1 - j) = img(startX + j, startY + squareSize - 1 - i);

            // Move top (temp) to right
            img(startX + j, startY + squareSize - 1 - i) = temp;
        }
    }
}


Image rotateSquareFromCenter(const Image& img, int squareSize, int times) {
    // Ensure the square size is valid
    if (squareSize <= 1 || squareSize > img.h || squareSize > img.w) {
        std::cerr << "Invalid square size. It must be greater than 1 and fit within the image dimensions." << std::endl;
        return badImg;  // Return the original image if the square size is invalid
    }

    // Calculate the starting coordinates to center the square
    int centerX = img.h / 2;
    int centerY = img.w / 2;
    int startX = centerX - squareSize / 2;
    int startY = centerY - squareSize / 2;

    // Ensure the square region does not exceed image bounds
    if (startX < 0 || startY < 0 || startX + squareSize > img.h || startY + squareSize > img.w) {
        return badImg;
    }

    // Create a copy of the original image to modify
    Image result = img;

    // Rotate the centered square using the rotateSquare utility function
    rotateSquare(result, startX, startY, squareSize, times);

    return result;
}


Image rotateSquareCorners(const Image& img, int squareSize, int times) {
    // Ensure the square size is valid
    if (squareSize <= 1 || squareSize > img.h || squareSize > img.w) {
        return badImg;  // Return the original image if the square size is invalid
    }

    // Create a copy of the original image to modify
    Image result = img;

    // Rotate each corner of the image using rotateSquare
    // 1. Top-left corner
    rotateSquare(result, 0, 0, squareSize, times);

    // 2. Top-right corner
    rotateSquare(result, 0, result.w - squareSize, squareSize, times);

    // 3. Bottom-left corner
    rotateSquare(result, result.h - squareSize, 0, squareSize, times);

    // 4. Bottom-right corner
    rotateSquare(result, result.h - squareSize, result.w - squareSize, squareSize, times);

    return result;
}


Image reverseSizes(Image_ img) {
    // Extract all shapes from the image
    auto shapes = extractConnectedComponents(img);
    if (shapes.empty()) {
        return badImg;
    }

    // Sort shapes by size (smallest to largest)
    std::vector<Image> sortedShapes = shapes;
    std::sort(sortedShapes.begin(), sortedShapes.end(), [](const Image& a, const Image& b) {
        return a.mask.size() < b.mask.size();
    });

    // Map each shape to its new size by reversing the sorted order
    std::vector<int> targetSizes;
    for (const auto& shape : sortedShapes) {
        targetSizes.push_back(shape.mask.size());
    }
    std::reverse(targetSizes.begin(), targetSizes.end());  // Reverse the order of sizes

    // Create the output image
    Image result = core::empty(img.p, img.sz);

    // Resize each shape to its new target size and place it in its original position
    for (size_t i = 0; i < shapes.size(); ++i) {
        const auto& originalShape = shapes[i];
        int newSize = targetSizes[i];
        Image resizedShape = resizeShape(originalShape, newSize);

        // Place the resized shape in its original location in the result image
        int startX = originalShape.p.x;
        int startY = originalShape.p.y;
        for (int x = 0; x < resizedShape.h; ++x) {
            for (int y = 0; y < resizedShape.w; ++y) {
                if (resizedShape(x, y) != 0) {  // Only place non-zero pixels
                    int targetX = startX + x;
                    int targetY = startY + y;
                    if (targetX < result.h && targetY < result.w) {
                        result(targetX, targetY) = resizedShape(x, y);
                    }
                }
            }
        }
    }

    return result;
}


Image shuffleRowsOrColumns(Image_ img, bool shuffleRows, int seed) {
    // Create a copy of the original image for the result
    Image result = core::empty(img.p, img.sz);

    // Initialize the random number generator with the seed for deterministic behavior
    std::mt19937 rng(seed);

    // Shuffle rows or columns based on `shuffleRows` flag
    if (shuffleRows) {
        // Generate a list of row indices to shuffle
        std::vector<int> rowIndices(img.h);
        std::iota(rowIndices.begin(), rowIndices.end(), 0);  // Fill with {0, 1, ..., img.h - 1}
        std::shuffle(rowIndices.begin(), rowIndices.end(), rng);  // Shuffle the row indices

        // Copy rows in shuffled order
        for (int i = 0; i < img.h; ++i) {
            int srcRow = rowIndices[i];
            for (int j = 0; j < img.w; ++j) {
                result(i, j) = img(srcRow, j);
            }
        }
    } else {
        // Generate a list of column indices to shuffle
        std::vector<int> colIndices(img.w);
        std::iota(colIndices.begin(), colIndices.end(), 0);  // Fill with {0, 1, ..., img.w - 1}
        std::shuffle(colIndices.begin(), colIndices.end(), rng);  // Shuffle the column indices

        // Copy columns in shuffled order
        for (int j = 0; j < img.w; ++j) {
            int srcCol = colIndices[j];
            for (int i = 0; i < img.h; ++i) {
                result(i, j) = img(i, srcCol);
            }
        }
    }

    return result;
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
    if (newWidth > MAXSIDE / 2 || newHeight > MAXSIDE / 2) {
        // Return original image if it cannot be upscaled within the 30x30 limit
        return badImg;
    }

    // Create a new image with the upscaled dimensions
    Image upscaled = core::empty({img.p.x, img.p.y}, {newWidth, newHeight});

    // Duplicate each pixel according to the scale factor
    for (int i = 0; i < img.h; ++i) {
        for (int j = 0; j < img.w; ++j) {
            const char pixelValue = img(i, j);
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
    if (newWidth > MAXSIDE || newHeight > MAXSIDE) {
        // Return the original image if the stretched version exceeds the limit
        return badImg;
    }

    // Create a new image with the stretched dimensions
    Image stretched = core::empty({img.p.x, img.p.y}, {newWidth, newHeight});

    // Populate the stretched image
    for (int i = 0; i < img.h; ++i) {
        for (int j = 0; j < img.w; ++j) {
            const char pixelValue = img(i, j);

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
            return badImg;
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

Image maskByColorMap(const Image& img, const std::unordered_map<int, int>& colorMap) {
    Image mask = core::empty(img.p, img.sz);
    for (int i = 0; i < img.h; ++i) {
        for (int j = 0; j < img.w; ++j) {
            int color = img(i, j);
            if (colorMap.count(color)) {
                mask(i, j) = 1;  // Mark as part of the mask
            } else {
                mask(i, j) = 0;
            }
        }
    }
    return mask;
}



Image replaceBackground(const Image& img, const int backgroundColor, const std::unordered_map<int, int>& colorMap) {
    Image result = img;
    if (colorMap.count(backgroundColor)) {
        int newColor = colorMap.at(backgroundColor);
        for (int i = 0; i < img.h; ++i) {
            for (int j = 0; j < img.w; ++j) {
                if (img(i, j) == backgroundColor) {
                    result(i, j) = newColor;  // Replace background color
                }
            }
        }
    }
    return result;
}

Image highlightEdges(const Image& img, const std::unordered_map<int, int>& colorMap) {
    // Create an empty output image with the same size as the input
    Image result = img;  // Start with a copy of the original image to preserve shape interiors

    // Extract all connected components (shapes) from the image
    auto shapes = extractConnectedComponents(img);

    // Directions for 8-neighbor connectivity to identify edges
    std::vector<point> directions = {{-1, -1}, {-1, 0}, {-1, 1},
                                     {0, -1},          {0, 1},
                                     {1, -1},  {1, 0}, {1, 1}};

    // Process each shape to highlight its edges
    for (const auto& shape : shapes) {
        // Determine the original color of the shape by examining a non-zero pixel in its mask
        int originalColor = 0;
        for (const char& c : shape.mask) {
            if (c != 0) {
                originalColor = c;
                break;
            }
        }

        // Find the mapped edge color for this shape from colorMap
        int edgeColor = colorMap.count(originalColor) ? colorMap.at(originalColor) : originalColor;

        // Iterate through each pixel in the shape to identify and highlight edges
        for (int i = 0; i < shape.h; ++i) {
            for (int j = 0; j < shape.w; ++j) {
                if (shape(i, j) == originalColor) {  // Only process shape pixels
                    bool isEdge = false;

                    // Check if this pixel is on the edge by looking at its 8 neighbors
                    for (const auto& dir : directions) {
                        int ni = shape.p.x + i + dir.x;
                        int nj = shape.p.y + j + dir.y;
                        
                        // If the neighbor pixel is outside the shape or is zero, this is an edge pixel
                        if (ni >= 0 && ni < img.h && nj >= 0 && nj < img.w && img(ni, nj) != originalColor) {
                            isEdge = true;
                            break;
                        }
                    }

                    // Highlight edge in the result image with the mapped edge color
                    if (isEdge) {
                        result(shape.p.x + i, shape.p.y + j) = edgeColor;
                    }
                }
            }
        }
    }

    return result;
}


Image removeNoise(const Image& img, int minSize) {
    // Extract all connected components (shapes) from the image
    auto shapes = extractConnectedComponents(img);
    Image result = core::empty(img.p, img.sz);  // Start with an empty image

    // Iterate through each shape
    for (const auto& shape : shapes) {
        // If the shape size is larger than or equal to minSize, add it to the result image
        if (shape.mask.size() >= minSize) {
            int startX = shape.p.x;
            int startY = shape.p.y;

            // Copy the shape pixels into the result image
            for (int i = 0; i < shape.h; ++i) {
                for (int j = 0; j < shape.w; ++j) {
                    if (shape(i, j) != 0) {  // Only copy non-zero pixels
                        result(startX + i, startY + j) = shape(i, j);
                    }
                }
            }
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



void radixSort(std::vector<std::pair<int, int>>& arr) {
    if (arr.empty()) return;

    // Find the maximum value to determine the number of digits
    int maxVal = 0;
    for (const auto& p : arr) {
        maxVal = std::max(maxVal, p.first);
    }

    int exp = 1;
    std::vector<std::pair<int, int>> output(arr.size());  // Output array for each pass
    const int base = 10;

    while (maxVal / exp > 0) {
        std::array<int, base> count = {0};  // Bucket counters

        // Count occurrences of each digit in the current place
        for (const auto& p : arr) {
            count[(p.first / exp) % base]++;
        }

        // Update count array to positions
        for (int i = 1; i < base; ++i) {
            count[i] += count[i - 1];
        }

        // Build the output array in stable order
        for (int i = arr.size() - 1; i >= 0; --i) {
            int digit = (arr[i].first / exp) % base;
            output[--count[digit]] = arr[i];
        }

        // Copy output back to arr
        std::copy(output.begin(), output.end(), arr.begin());

        // Move to the next digit place
        exp *= base;
    }
}


Image myStack(vImage_ lens, int id) {
  const int n = lens.size();
  if (!n) return badImg;
  vector<pair<int,int>> order(n);
  for (int i = 0; i < n; ++i) {
    order[i] = {lens[i].w*lens[i].h,i};
  }
   sort(order.begin(), order.end());
    // radixSort(order);
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
    // radixSort(order);
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
        std::array<std::vector<std::pair<int, int>>, 10> buckets;
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

        // Direct pixel access with boundary check
        for (int i = 0; i < img.h; ++i) {
            const int idy = i + dy;
            if (idy < 0 || idy >= ret.h) continue; // Check row bounds

            for (int j = 0; j < img.w; ++j) {
                const int idx = j + dx;
                if (idx < 0 || idx >= ret.w) continue; // Check column bounds

                if (img(i, j)) {
                    ret(idy, idx) = img(i, j);
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

Image cutUnion(Image_ img) {
    Image cutRegion = heuristicCut(img);
    vImage components = cut(img, cutRegion);

    Image unionImage = core::empty(img.p, img.sz);

    for (const auto& component : components) {
        for (int i = 0; i < component.h; ++i) {
            for (int j = 0; j < component.w; ++j) {
                if (component(i, j) != 0 && i < img.h && j < img.w) {
                    unionImage(i, j) = component(i, j);
                }
            }
        }
    }

    return unionImage;
}

Image cutIntersection(Image_ img) {
    // Apply heuristic cut to get components
    Image cutRegion = heuristicCut(img);
    vImage components = cut(img, cutRegion);

    // If there are no components, return an empty image
    if (components.empty()) {
        return badImg;
    }

    // Initialize the intersection image with the first component's mask
    Image intersectionImage = core::empty(img.p, img.sz);
    for (int i = 0; i < components[0].h; ++i) {
        for (int j = 0; j < components[0].w; ++j) {
            intersectionImage(i, j) = components[0](i, j);
        }
    }

    // Iterate over each component to compute the intersection
    for (size_t k = 1; k < components.size(); ++k) {
        for (int i = 0; i < img.h; ++i) {
            for (int j = 0; j < img.w; ++j) {
                // Only retain pixels that are present in all components
                if (i < components[k].h && j < components[k].w && components[k](i, j) == 0) {
                    intersectionImage(i, j) = 0;
                }
            }
        }
    }

    return intersectionImage;
}

Image cutDifference(Image_ img) {
    Image cutRegion = heuristicCut(img);
    vImage components = cut(img, cutRegion);

    Image differenceImage = img;

    for (const auto& component : components) {
        for (int i = 0; i < component.h; ++i) {
            for (int j = 0; j < component.w; ++j) {
                if (component(i, j) != 0 && i < img.h && j < img.w) {
                    differenceImage(i, j) = 0;
                }
            }
        }
    }

    return differenceImage;
}

Image cutComposeMultiple(Image_ img, int id) {
    Image cutRegion = heuristicCut(img);
    vImage components = cut(img, cutRegion);

    // Align each component to origin and compose into a single image
    for (auto& component : components) {
        component = toOrigin(component);
    }

    return !components.empty() ? compose(components, id) : img;
}

Image cutFilterByColor(Image_ img, int targetColor) {
    Image cutRegion = heuristicCut(img);
    vImage components = cut(img, cutRegion);

    // Filter components based on the target color
    vImage colorFilteredComponents;
    for (const auto& component : components) {
        if (core::majorityCol(component) == targetColor) {
            colorFilteredComponents.push_back(component);
        }
    }

    // Compose the filtered components back into a single image
    return !colorFilteredComponents.empty() ? compose(colorFilteredComponents, 0) : badImg;
}


namespace core {
  unsigned int colMask(Image_ img);
  unsigned int countCols(Image_ img, int include0 = 0);
  unsigned int count(Image_ img);
  Image full(point p, point sz, int filling = 0);
  Image empty(point p, point sz);
  Image full(point sz, int filling = 0);
  Image empty(point sz);
  bool isRectangle(Image_ a);
  unsigned int countComponents(Image img);
  char majorityCol(Image_ img, int include0 = 0);
  Image subImage(Image_ img, point p, point sz);
  vector<pair<Image,int>> splitCols(Image_ img, int include0 = 0);
};

struct Loader {
  long n;
  string text;
  long long counter;
  long prev;
  short keep_title;
  Loader(long long n_, string text_ = "") {
    n = n_;
    text = text_;
    counter = 0;
    prev = -1;
    keep_title = 0;
  }
  void operator()();
  ~Loader();
};

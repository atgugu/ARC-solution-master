
struct deduceOuterProduct {
  vector<pair<Image,Image>> train_targets;
  int rec_funci;
  deduceOuterProduct(vector<pair<Image,Image>> train);
  Image reconstruct(Image_ a, Image_ b);
};
void addDeduceOutDeduceOuterProduct(Pieces&pieces, vector<pair<Image,Image>> train, vector<Candidate>&cands);

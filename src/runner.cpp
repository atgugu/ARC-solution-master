#include "precompiled_stl.hpp"

using namespace std;

#include "utils.hpp"
#include "core_functions.hpp"
#include "image_functions.hpp"
#include "visu.hpp"
#include "read.hpp"
#include "normalize.hpp"
#include "tasks.hpp"
#include "evals.hpp"

#include "brute2.hpp"

#include "score.hpp"
#include "load.hpp"

#include "deduce_op.hpp"
#include "pieces.hpp"
#include "compose2.hpp"

#include "brute_size.hpp"

#include <thread>

string green(string s) {
  return ("\033[1;32m"+s+"\033[0m");
}
string blue(string s) {
  return ("\033[1;34m"+s+"\033[0m");
}
string yellow(string s) {
  return ("\033[1;33m"+s+"\033[0m");
}
string red(string s) {
  return ("\033[1;31m"+s+"\033[0m");
}


void writeVerdict(int si, string sid, int verdict) {
  
  printf("Task #%2d (%s): ", si, sid.c_str());
  switch (verdict) {
  case 3: cout << green("Correct") << endl; break;
  case 2: cout << yellow("Candidate") << endl; break;
  case 1: cout << blue("Dimensions") << endl; break;
  case 0: cout << red("Nothing") << endl; break;
  default: assert(0);
  }
  
}

int MAXDEPTH = -1; //Argument
int NUMFUNCS = 200;
int MAXSIDE = 100, MAXAREA = 40*40, MAXPIXELS = 40*40*5; //Just default values
unsigned int SEED = 4;
int print_times = 0, print_mem = 0, print_nodes = 0;

void run(int only_sid = -1, int arg = -1) {
    // Default argument processing and feature flag setup
    int no_norm   = (arg >= 10 && arg < 20);
    int add_flips = (arg >= 20 && arg < 40);
    int add_flip_id = (arg >= 30 && arg < 40 ? 7 : 6);
    
    if (arg == -1) arg = 2;
    MAXDEPTH = (arg % 10) * 10;
    NUMFUNCS -= MAXDEPTH;
    int eval = 0;
    
    // Directory and sample selection
    string sample_dir = eval ? "test" : "training";
    int samples = eval ? -1 : -1;

    // Read all samples
    vector<Sample> sample = readAll(sample_dir, samples);
    vector<int> verdict(sample.size());
    Visu visu;

    // Prepare loaders and counters
    int dones = 0;
    Loader load(sample.size());
    
    assert(only_sid < sample.size());

    int scores[4] = {};
    int skips = 0;

    unsigned train_size_max = 0;
    for (const auto& sam: sample) {
        if(sam.train.size() > train_size_max) train_size_max = sam.train.size();
    }

    vector<pair<Image, Image>> train;
    train.reserve(train_size_max);
    vector<point> out_sizes;
    vector<Candidate> cands;
    vector<Candidate> answers;
    for (int si = 0; si < sample.size(); ++si) {
        if (only_sid != -1 && si != only_sid) continue;

        // Progress display
        if (eval) load();
        else if (++dones % 10 == 0) cout << dones << " / " << sample.size() << endl;

        const Sample& s = sample[si];

        // Normalize sample with or without normalization
        Simplifier sim = no_norm ? normalizeDummy(s.train) : normalizeCols(s.train);
        train.clear();

        // Process training data
        for (auto& [in, out] : s.train) {
            train.push_back(sim(in, out));
        }

        // Apply flips if required
        if (add_flips) {
            size_t initial_size = train.size();
            train.reserve(initial_size * 2);
            for (size_t i = 0; i < initial_size;++i) {
                auto [rin, rout] = sim(train[i].first, train[i].second);
                train.push_back({rigid(rin, add_flip_id), rigid(rout, add_flip_id)});
            }
        }

        auto [test_in, test_out] = sim(s.test_in, s.test_out);

        // Feature calculation and dataset constraints
        {
            int insumsz = 0, outsumsz = 0, macols = 0;
            int maxside = 0, maxarea = 0;
            for (auto& [in, out] : s.train) {
                maxside = max({maxside, in.w, in.h, out.w, out.h});
                maxarea = max({maxarea, in.w * in.h, out.w * out.h});
                insumsz += in.w * in.h;
                outsumsz += out.w * out.h;
                macols = max(macols, __builtin_popcount(core::colMask(in)));
            }
            int sumsz = max(insumsz, outsumsz);

            double w[4] = {1.2772523019346949, 0.00655104, 0.70820414, 0.00194519};
            double expect_time3 = w[0] + w[1] * sumsz + w[2] * macols + w[1] * w[2] * sumsz * macols;
            MAXSIDE = maxside * 2;
            MAXAREA = maxarea * 2;
            MAXPIXELS = MAXAREA * 5;
        }

        // Size calculation for brute force methods
        out_sizes.clear();
        out_sizes = bruteSize(test_in, train);

        // Generate candidate pieces
        Pieces pieces;
        {
            double start_time = now();
            vector<DAG> dags = brutePieces2(test_in, train, out_sizes);
            pieces = makePieces2(dags, train, out_sizes);
        }

        // Clear memory to avoid memory bloat
        #pragma omp parallel for
        for (DAG& d : pieces.dag) {
            d.hashi.clear();
            for (TinyNode& n : d.tiny_node.node) {
                n.child.clear();
            }
        }

        // Score and evaluate
        int s1 = (out_sizes.back() == test_out.sz) ? 1 : 0;
        cands.clear();
        cands = composePieces2(pieces, train, out_sizes);
        addDeduceOuterProduct(pieces, train, cands);
        cands = evaluateCands(cands, train);
        int s2 = scoreCands(cands, test_in, test_out);

        // Pick best candidates
        answers.clear();
        answers = cands;
        {
            sort(cands.begin(), cands.end());
            set<ul> seen;
            answers.clear();
            answers.reserve(3 + skips * 3);

            for (const Candidate& cand : cands) {
                ull h = hashImage(cand.imgs.back());
                if (seen.insert(h).second) {
                    answers.push_back(cand);
                    if (answers.size() == 3 + skips * 3) break;
                }
            }
            if (skips * 3 < answers.size()) answers.erase(answers.begin(), answers.begin() + skips * 3);
        }

        // Reconstruct answers
        const unsigned ansSize = answers.size();
        vector<Image> rec_answers;
        vector<double> answer_scores;
        rec_answers.reserve(ansSize);
        answer_scores.reserve(ansSize);

        for (Candidate& cand : answers) {
            rec_answers.push_back(sim.rec(s.test_in, cand.imgs.back()));
            answer_scores.push_back(add_flips ? cand.score / (2 - 1e-5) : cand.score);
        }

        int s3 = scoreAnswers(rec_answers, s.test_in, s.test_out);

        // Visualization and logging
        if (!eval) {
            visu.next(to_string(si) + " - test");
            for (auto& [in, out] : train) visu.add(in, out);
            visu.add(test_in, test_out);
            visu.next(to_string(si) + " - cands");
            for (int i = 0; i < min((int)answers.size(), 5); ++i) {
                visu.add(test_in, answers[i].imgs.back());
            }
        }

        // Verdict calculation
        if (!eval) {
            verdict[si] = (s3) ? 3 : (s2) ? 2 : (s1) ? 1 : 0;
            scores[verdict[si]]++;
            writeVerdict(si, s.id, verdict[si]);
        }

        // Write final answers to CSV
        string fn = "output/answer_" + to_string(only_sid) + "_" + to_string(arg) + ".csv";
        writeAnswersWithScores(s, fn, rec_answers, answer_scores);
    }

    // Final score display
    if (!eval && only_sid == -1) {
        for (int si = 0; si < sample.size(); ++si) {
            Sample& s = sample[si];
            writeVerdict(si, s.id, verdict[si]);
        }

        for (int i = 3; i > 0; --i) scores[i - 1] += scores[i];
        printf("\nTotal: %4d\nPieces: %4d\nCands: %4d\nCorrect: %3d\n", scores[0], scores[1], scores[2], scores[3]);
    }
}
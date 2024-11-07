DEBUG = -Wshadow -Wall -fsanitize=address -fsanitize=undefined -D_GLIBCXX_DEBUG -g -Wno-sign-compare -Wno-shadow -Wno-char-subscripts -Wno-unused-variable
#-fprofile-generate -fprofile-update -mcpu=sandybridge
FLAGS = -std=c++17 -floop-interchange -mavx2 -floop-interchange -floop-unroll-and-jam -mfma -ftree-vectorize -fno-strict-aliasing -foptimize-sibling-calls -fprefetch-loop-arrays -floop-nest-optimize -fipa-pta -fno-strict-overflow -fgcse -fopenmp -flto=auto  -march=native -mtune=native -funroll-loops -finline-functions -fomit-frame-pointer -ffast-math -fno-stack-protector -fno-math-errno -funsafe-math-optimizations -fno-trapping-math -fno-signed-zeros -fassociative-math -freciprocal-math -ffinite-math-only -fexcess-precision=fast -fno-rounding-math -fno-signaling-nans -fcx-limited-range -fno-math-errno

LIBS = -lstdc++ -lstdc++fs

.DEFAULT_GOAL := run

PREC = headers/precompiled_stl.hpp
HEAD = utils.hpp read.hpp visu.hpp image_functions.hpp normalize.hpp core_functions.hpp brute2.hpp timer.hpp efficient.hpp
HEAD_PATH = $(addprefix headers/,$(HEAD)) $(PREC).gch | obj

OBJ = read.o core_functions.o image_functions.o image_functions2.o visu.o normalize.o tasks.o runner.o score.o load.o evals.o brute2.o deduce_op.o pieces.o compose2.o brute_size.o efficient.o

OBJ_PATH = $(addprefix obj/,$(OBJ))

obj:
	mkdir -p obj
	mkdir -p output

$(PREC).gch: $(PREC)
	g++ -c $< $(FLAGS)

obj/%.o: src/%.cpp $(HEAD_PATH)
	g++ -c $< $(FLAGS) -o $@ -I headers

all: $(OBJ_PATH)

run: src/main.cpp $(OBJ_PATH) $(HEAD_PATH) headers/tasks.hpp
	g++ src/main.cpp $(OBJ_PATH) $(FLAGS) $(LIBS) -o run -I headers

count_tasks: src/count_tasks.cpp obj/read.o headers/utils.hpp
	g++ src/count_tasks.cpp obj/read.o $(FLAGS) $(LIBS) -o count_tasks -I headers

# Clean target
clean:
	rm -rf obj/*.o run  store/*
clean_profiles:
	rm -rf profilerdata/*

gen_profile: FLAGS += -fprofile-generate -fprofile-dir=./profilerdata -Ofast
gen_profile: clean_profiles all run count_tasks

use_profile: FLAGS += -fprofile-use -fprofile-dir=./profilerdata -Ofast
use_profile: all run count_tasks

test: FLAGS += -Ofast
test: clean all run count_tasks

save_space: FLAGS += -fprofile-use -fprofile-dir=./profilerdata -Os
save_space: clean all run count_tasks
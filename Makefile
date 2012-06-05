CC = g++
CFLAGS = -Wall -Werror -ansi -pedantic -I/home/ugoc/include -Iinclude

MACHINE = $(shell uname -m)
LDPATH = -L/home/ugoc/lib

vpath %.cpp src
vpath %.o   obj
vpath %.h   include

.PHONY: all clean mk_machine_dir

BAREOBJ = query_vite_runner.o dtw_util.o query_hmm.o
OBJ = $(addprefix obj/, $(BAREOBJ))

TARGET = bin/std_viterbi_dtw


all: CFLAGS += -O2 -DLOGGING

debug: CFLAGS += -g -DLOGGING

all: mk_machine_dir $(TARGET)

debug: $(TARGET)

mk_machine_dir:
	@mkdir -p obj/
	@mkdir -p bin/

bin/std_viterbi_dtw: LDFLAGS += -pthread \
	-l hmmlite \
	-l atlas_wrapper \
	-l atlas_wrapper_type \
	-l lapack \
	-l lapack_atlas \
	-l anguso_arg_parser \
	-l parm \
	-l feature \
	-l segtree \
	-l ugoc_utility

obj/%.o: src/%.cpp include/%.h
	${CC} ${CFLAGS} $< -c -o $@ ${LDFLAGS} $(LDPATH)

bin/std_viterbi_dtw: std_viterbi_dtw.cpp ${OBJ}
	mkdir -p bin
	${CC} ${CFLAGS} $^ -o $@ ${LDFLAGS} $(LDPATH)

clean:
	rm -f bin/std_viterbi_dtw

allclean:
	rm -f bin/std_viterbi_dtw ${OBJ}

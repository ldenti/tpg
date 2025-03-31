CXX=g++
CXXFLAGS=-g -O3 -Wall
LIBS=
LDFLAGS=-lz

CWD:=$(shell pwd)

LIB_DIR:=lib
INC_DIR:=include
EXT_DIR:=ext
LOG_DIR:=logs

LIBHANDLEGRAPH_DIR:=$(EXT_DIR)/libhandlegraph
SDSL_DIR:=$(EXT_DIR)/sdsl-lite
GBWT_DIR:=$(EXT_DIR)/gbwt
GBWTGRAPH_DIR=$(EXT_DIR)/gbwtgraph

LIB_DEPS =
LIB_DEPS += $(LIB_DIR)/libhandlegraph.a
LIB_DEPS += $(LIB_DIR)/libsdsl.a
LIB_DEPS += $(LIB_DIR)/libgbwt.a
LIB_DEPS += $(LIB_DIR)/libgbwtgraph.a

.PHONY: all deps clean clean-all

all: tpg

deps: $(LIB_DIR)/libgbwtgraph.a

logs:
	mkdir -p $(CWD)/$(LOG_DIR)

$(LIB_DIR)/libhandlegraph.a: logs $(LIBHANDLEGRAPH_DIR)/src/include/handlegraph/*.hpp $(LIBHANDLEGRAPH_DIR)/src/*.cpp
	echo "* Building libhandlegraph"
	cd $(LIBHANDLEGRAPH_DIR) && rm -Rf build CMakeCache.txt CMakeFiles && mkdir build && cd build && cmake -DCMAKE_INSTALL_PREFIX=$(CWD) -DCMAKE_INSTALL_LIBDIR=lib .. &> $(CWD)/$(LOG_DIR)/libhandlegraph.log && make &>> $(CWD)/$(LOG_DIR)/libhandlegraph.log && make install &>> $(CWD)/$(LOG_DIR)/libhandlegraph.log

$(LIB_DIR)/libsdsl.a $(LIB_DIR)/libdivsufsort.a $(LIB_DIR)/libdivsufsort64.a &: logs $(SDSL_DIR)/lib/*.cpp $(SDSL_DIR)/include/sdsl/*.hpp
	echo "* Building sdsl-lite"
	cd $(SDSL_DIR) && BUILD_PORTABLE=1 ./install.sh $(CWD) &> $(CWD)/$(LOG_DIR)/sdsl.log

$(LIB_DIR)/libgbwt.a: logs $(LIB_DIR)/libsdsl.a $(LIB_DIR)/libdivsufsort.a $(LIB_DIR)/libdivsufsort64.a $(wildcard $(GBWT_DIR)/src/*.cpp) $(wildcard $(GBWT_DIR)/include/gbwt/*.h)
	echo "* Building gbwt"
	cp -r $(GBWT_DIR)/include/gbwt $(CWD)/$(INC_DIR)/
	cd $(GBWT_DIR) && make clean &> $(CWD)/$(LOG_DIR)/gbwt.log && make &>> $(CWD)/$(LOG_DIR)/gbwt.log && mv lib/libgbwt.a $(CWD)/$(LIB_DIR)

$(LIB_DIR)/libgbwtgraph.a: logs $(LIB_DIR)/libgbwt.a $(LIB_DIR)/libsdsl.a $(LIB_DIR)/libdivsufsort.a $(LIB_DIR)/libdivsufsort64.a $(LIB_DIR)/libhandlegraph.a $(wildcard $(GBWTGRAPH_DIR)/src/*.cpp) $(wildcard $(GBWTGRAPH_DIR)/include/gbwtgraph/*.h)
	echo "* Building gbwtgraph"
	cp -r $(GBWTGRAPH_DIR)/include/gbwtgraph $(CWD)/$(INC_DIR)/
	cd $(GBWTGRAPH_DIR) && make clean &> $(CWD)/$(LOG_DIR)/gbwtgraph.log && make &>> $(CWD)/$(LOG_DIR)/gbwtgraph.log && mv lib/libgbwtgraph.a $(CWD)/$(LIB_DIR) 

%.o: %.cpp $(LIB_DEPS)
	@echo '* Compiling $<'
	$(CXX) $(CFLAGS) -I$(INC_DIR) -o $@ -c $<

tpg: main.o graph.o path.o segments.o misc.o $(LIB_DEPS)
	@echo "* Linking $<"
	$(CXX) $(CFLAGS) $(LIBS) -o $@ $^ $(LDFLAGS)

clean:
	rm -f tpg main.o graph.o segments.o path.o misc.o

clean-all:
	rm -rf tpg main.o graph.o segments.o path.o misc.o $(LIB_DIR) $(INC_DIR)

CPPFLAGS=-O3 -g -Wall -pedantic
CC=gcc

SRC=extract_reassortments.cc test_tree_code.cc mcmc_split_info.cc tree.cc splits.cc util.cc gamma-prob.c build_incompat_graph.cc catalog.cc

giraf: giraf.o extract_reassortments.o mcmc_split_info.o tree.o splits.o util.o dist.o gamma-prob.o build_incompat_graph.o catalog.o
	$(CXX) -o $@ $^

all: giraf

advanced: extract_reassortments mcmc_split_info build_incompat_graph

extract_reassortments: main_extract.o extract_reassortments.o
	$(CXX) -o $@ $^

mcmc_split_info: main_split.o mcmc_split_info.o splits.o tree.o util.o
	$(CXX) -o $@ $^

build_incompat_graph: main_graph.o build_incompat_graph.o dist.o gamma-prob.o splits.o tree.o util.o
	$(CXX) -o $@ $^

test_tree_code: test_tree_code.o splits.o tree.o util.o
	$(CXX) -o $@ $^

depend:
	makedepend -- $(CFLAGS) -- $(SRC)

clean:
	rm -f giraf
	rm -f extract_reassortments mcmc_split_info build_incompat_graph 
	rm -f binomial_invcdf normal_invcdf
	rm -f test_tree_code
	rm -f *.o

# DO NOT DELETE

extract_reassortments.o: label_types.h bigraph.h scored_set.h mica.h timer.h
extract_reassortments.o: options.h
build_incompat_graph.o: tree.h util.h splits.h dist.h options.h
test_tree_code.o: tree.h util.h splits.h
mcmc_split_info.o: tree.h util.h splits.h options.h
tree.o: tree.h util.h
splits.o: splits.h tree.h util.h
util.o: util.h
catalog.o: catalog.h util.h
dist.o: util.h tree.h dist.h
giraf.o: util.h catalog.h timer.h
main_graph.o: timer.h

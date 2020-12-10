#!/bin/bash
g++ -O3 main.cpp Vertex.cpp Graph.cpp NN.cpp -o max-e-k-sat -std=c++11 -larmadillo -lmlpack -lboost_serialization

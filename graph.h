// graph.h
// Jay Patel
//
// Basic graph class using adjacency matrix representation.  Currently
// limited to a graph with at most 100 vertices.
//

#pragma once

#include <iostream>
#include <stdexcept>
#include <vector>
#include <map>
#include <set>

using namespace std;

template<typename VertexT, typename WeightT>
class graph {
 private:
  map<VertexT, vector<pair<VertexT, WeightT>>> adjacency_list;
	vector<VertexT> Vertices;
	int edges;

 public:
  //
  // constructor:
  //
  // Constructs an empty graph.
  //
  graph() {
		edges = 0;
  }

	//
	// remove all vertices and edges
	//
	void clear() {
		adjacency_list.clear();
		Vertices.clear();
		edges = 0;
	}

  //
  // NumVertices
  //
  // Returns the # of vertices currently in the graph.
  //
  int NumVertices() const {
    return adjacency_list.size();
  }

  //
  // NumEdges
  //
  // Returns the # of edges currently in the graph.
  //
  int NumEdges() const {
    return this->edges;
  }

  //
  // addVertex
  //
  // add a new vertex to the graph (with no edges),
  // return false if vertex already exists.
  //
  bool addVertex(VertexT v) {
		if (adjacency_list.find(v) != adjacency_list.end()) {
			return false;
		}

		adjacency_list[v];
		Vertices.push_back(v);

    return true;
  }

  //
  // addEdge
  //
  // add a weighted edge from vertex src to vertex to
  // (overwrite weight value if edge already exists),
  // return false if either vertex doesn’t exist
  //
  bool addEdge(VertexT from, VertexT to, WeightT weight) {
		if (adjacency_list.find(from) == adjacency_list.end()) {
			return false;
		} else if (adjacency_list.find(to) == adjacency_list.end()) {
			return false;
		}

		pair<VertexT, WeightT> p(to, weight);
		
		for (auto &e : adjacency_list[from]) {
			if (e.first == to) {
				e.second = weight;
				return true;
			}
		}

		adjacency_list[from].push_back(p);
		this->edges++;
		return true;
  }

  //
  // getWeight
  //
  // return the weight of the edge from vertex src to vertex to,
  // return false if either vertex doesn’t exist
  //
  bool getWeight(VertexT from, VertexT to, WeightT& weight) const {
    if (adjacency_list.find(from) == adjacency_list.end()) {
			return false;
		} else if (adjacency_list.find(to) == adjacency_list.end()) {
			return false;
		}

		for (auto e : adjacency_list.at(from)) {
			if (e.first == to) {
				weight = e.second;
				return true;
			}
		}

    return false;
  }

  //
  // neighbors
  //
  // given a start vertex v, return a set of its
  // neighbors (all vertices v has an edge to)
  //
  set<VertexT> neighbors(VertexT v) const {
    set<VertexT> S;

		if (adjacency_list.find(v) == adjacency_list.end()) {
			return S;
		}

		for (auto edge : adjacency_list.at(v)) {
			S.insert(edge.first);
		}

    return S;
  }

  //
  // getVertices
  //
  // Returns a vector containing all the vertices currently in
  // the graph.
  //
  vector<VertexT> getVertices() const {
    return this->Vertices;
  }

  //
  // dump
  //
  // Dumps the internal state of the graph for debugging purposes.
  //
  // Example:
  //    graph<string,int>  G(26);
  //    ...
  //    G.dump(cout);  // dump to console
  //
  void dump(ostream& output) const {
    output << "***************************************************" << endl;
    output << "********************* GRAPH ***********************" << endl;

    output << "**Num vertices: " << this->NumVertices() << endl;
    output << "**Num edges: " << this->NumEdges() << endl;

    output << endl;
    output << "**Vertices:" << endl;
    for (int i = 0; i < this->NumVertices(); ++i) {
      output << " " << i << ". " << this->Vertices[i] << endl;
    }

    output << endl;
    output << "**Edges:" << endl;
    for (int i = 0; i < this->NumVertices(); i++) {
      output << this->Vertices[i] << ": ";

      for (auto edges : adjacency_list.at(this->Vertices[i])) {
          output << "(" << Vertices[i] << ", " << edges.first << ", " << edges.second << ") ";
      }
			output << endl;
		}
    output << "**************************************************" << endl;
  }
};

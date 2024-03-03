// tests.cpp - to test graph.h using GoogleTest framework
//             (this is encouraged, but not graded)

#include <gtest/gtest.h>
#include "graph.h"

TEST(graph, test1) {
	struct weight {
		int num;
		string str;
	};

	graph<int, weight> g;
	
	for (int i = 0; i < 5; i++) {
		g.addVertex(i);
	}

	int k = 0;
	for (int i = 5; i >= 0; i--) {
		weight st;

		st.num = k * i;
		st.str = "jay";

		g.addEdge(k, i, st);

		k++;
	}
	
	weight w;
	g.getWeight(0, 5, w);

	ASSERT_EQ(w.str, "jay");
}
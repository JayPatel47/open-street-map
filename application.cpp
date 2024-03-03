// CREATIVE COMPONENT: My creative component of this program has
// enhanced the project. It takes multiple user locations until
// the user types "#". It then stores finds if that user input
// is valid building on the map, if not, then prompts the user
// to input another string. Once all inputs have been taken, it
// stores all the building in vector of type BuildingInfo. And
// here, it follows the same functionality as the application(),
// but with vector of buildings and nodes instead of just two
// buildings. I've made a function called centerBetweenPoints()
// which is different than the one in dist.cpp. It finds the
// center point of all the buildings in the vector.
//
// To run the creative component, just insert different building
// strings and hit the "#" once done. The creative component will
// find the shortest path for all the people to a center point.
//
// application.cpp
// Jay Patel
//

#include <iostream>
#include <iomanip>  /*setprecision*/
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <limits>
#include <queue>
#include <stack>
#include <cmath>

#include "tinyxml2.h"
#include "dist.h"
#include "osm.h"
#include "graph.h"

using namespace std;
using namespace tinyxml2;

const double INF = numeric_limits<double>::max();

//
// centerBetweenPoints() - for creative component
//
// Its a new function I made, inspired by the one in dist.cpp.
// It finds the centroid of multiple coordinates. It takes
// vector(containes multiple buildings) as parameter. It
// returns the coordinates of the centroid.
//
Coordinates centerBetweenPoints(vector<BuildingInfo> usersBuildings) {
  double PI = 3.14159265;
	double x = 0, y = 0, z = 0;

	for (auto e : usersBuildings) {
		// degree to radians
		double latRad = e.Coords.Lat * PI / 180.0;
		double lonRad = e.Coords.Lon * PI / 180.0;
		x += cos(latRad) * cos(lonRad);
		y += cos(latRad) * sin(lonRad);
		z += sin(latRad);
	}

	double centroidX = x / usersBuildings.size();
	double centroidY = y / usersBuildings.size();
	double centroidZ = z / usersBuildings.size();
	
	double centroidLon = atan2(centroidY, centroidX);
	double centroidLat = atan2(centroidZ, sqrt(centroidX * centroidX +
																						 centroidY * centroidY));
	//radians to degree
	centroidLat = centroidLat * 180.0 / PI;
	centroidLon = centroidLon * 180.0 / PI;

	return Coordinates(-1, centroidLat, centroidLon);
}

//
// centralBuilding() - for creative component
//
// It finds the building nearest to the centroid. Takes
// vector of usersBuildings and Buildings as parameter.
// returns the building nearest to the centroid.
//
BuildingInfo centralBuilding(vector<BuildingInfo> usersBuildings,
						     vector<BuildingInfo> Buildings, set<string> S) {
	Coordinates midpoint = centerBetweenPoints(usersBuildings);

	double min = INF;
	BuildingInfo center;

	for (auto b : Buildings) {
		double distance = distBetween2Points(midpoint.Lat, midpoint.Lon,
																				 b.Coords.Lat, b.Coords.Lon);

		if (distance < min && S.count(b.Fullname) == 0) {
			min = distance;
			center = b;
		}
	}

	return center;
}

//
// nearestNode()
//
// This function finds the nearest node from the given building.
// Takes the building, map of nodes and vector of footways as
// parameter. Returns the nearest node.
//
long long nearestNode(BuildingInfo b, map<long long, Coordinates> Nodes,
											vector<FootwayInfo> Footways) {
	long long node;
	double min = INF;

	for (auto e : Footways) {
		for (unsigned int i = 0; i < e.Nodes.size(); i++) {
			Coordinates c = Nodes[e.Nodes[i]];
			double distance = distBetween2Points(c.Lat,
												c.Lon, b.Coords.Lat, b.Coords.Lon);

			if (distance < min) {
				min = distance;
				node = e.Nodes[i];
			}
		}
	}

	return node;
}

//
// findNodes() - for creative component
//
// This function calls the nearestNode() function for all
// building in a vector and stores all the nearest nodes
// in a vector. It takes a vector of buildings, map of
// nodes, and vector of footways as parameter. It returns
// the vector that contains all the nearest nodes.
//
vector<long long> findNodes(vector<BuildingInfo> usersBuildings,
							map<long long, Coordinates> Nodes,
							vector<FootwayInfo> Footways) {
	vector<long long> nearestNodes;
	for (auto n : usersBuildings) {
		long long result = nearestNode(n, Nodes, Footways);
		nearestNodes.push_back(result);
	}
	return nearestNodes;
}

//
// printData() - for creative component
//
// this function prints everthing stored in the
// vector of buildings and nearestNodes.
//
void printData(vector<BuildingInfo> usersBuildings,
	vector<long long> nearestNodes, long long dest,
	BuildingInfo center, map<long long, Coordinates> Nodes, bool flag) {
	if (flag == true) {
		cout << "New destination building:" << endl;
		cout << " " << center.Fullname << endl;
		cout << " (" << center.Coords.Lat << ", " << center.Coords.Lon << ")" << endl;
		cout << "Nearest destination node:" << endl;
		cout << " " << dest << endl;
		cout << " (" << Nodes[dest].Lat << ", " << Nodes[dest].Lon << ")" << endl;
		cout << endl;
	} else {
		int i = 1;
		for (auto b : usersBuildings) {
			cout << "Person "<< i <<"'s point:" << endl;
			cout << " " << b.Fullname << endl;
			cout << " (" << b.Coords.Lat << ", " << b.Coords.Lon << ")" << endl;
			i++;
		}

		cout << "Destination Building:" << endl;
		cout << " " << center.Fullname << endl;
		cout << " (" << center.Coords.Lat << ", " << center.Coords.Lon << ")" << endl;
		cout << endl;

		i = 1;
		for (auto n : nearestNodes) {
			cout << "Nearest P" << i << " node:" << endl;
			cout << " " << n << endl;
			cout << " (" << Nodes[n].Lat << ", " << Nodes[n].Lon << ")" << endl;
			i++;
		}

		cout << "Nearest destination node:" << endl;
		cout << " " << dest << endl;
		cout << " (" << Nodes[dest].Lat << ", " << Nodes[dest].Lon << ")" << endl;
		cout << endl;
	}
}

//
// printPath()
//
// prints the vector given in the parameter.
//
void printPath(vector<long long> path) {
	cout << "Path: " << path[0];
	for (unsigned int i = 1; i < path.size(); i++) {
		cout << "->" << path[i];
	}
	cout << endl;
}

//
// getPath()
//
// This function uses the predecessors map to find the path
// to the destination. It takes map predecessors, map distances,
// and start vertex and end vertex.
//
void getPath(map<long long, long long> predecessors,
						 map<long long, double> distances,
						 long long startV, long long endV) {
	if (distances[endV] >= INF) {
		cout << "Sorry, destination unreachable." << endl;
		return;
	}
	stack<long long> st;
	vector<long long> path;
	long long currentV = endV;
	while (currentV != 0) {
		st.push(currentV);
		currentV = predecessors[currentV];
	}

	while (!st.empty()) {
		currentV = st.top();
		st.pop();
		path.push_back(currentV);
	}

	printPath(path);
}

//
// class prioritize
//
// This class is used by the priority_queue in
// dijkstra's algorithm.
//
class prioritize {
	public:
		bool operator()(const pair<long long, double>& p1,
										const pair<long long, double>& p2) const {
			return p1.second > p2.second;
		}
};

//
// dijkstraAlgo()
//
// This is the algorithm that finds the shortest path from start vertex to
// end vertex. It takes a graph, map of distances, map of predecessors and
// start and end vertex.
//
void dijkstraAlgo(graph<long long, double> G, map<long long, double>& distances,
	map<long long, long long>& predecessors, long long startV, long long endV) {
	priority_queue<
  pair<long long, double>,
  vector<pair<long long, double>>,
  prioritize> unvisitedQueue;
	set<long long> visited;

	for (auto currentV : G.getVertices()) {
		distances[currentV] = INF;
		predecessors[currentV] = 0;
		unvisitedQueue.push(make_pair(currentV, INF));
	}

	distances[startV] = 0;
	unvisitedQueue.push(make_pair(startV, 0));

	while (!unvisitedQueue.empty()) {
		pair<long long, double> currentV = unvisitedQueue.top();
		unvisitedQueue.pop();

		if (distances[currentV.first] == INF) {
			break;
		} else if (visited.count(currentV.first) == 1) {
			continue;
		}

		visited.insert(currentV.first);
		for (auto adjV : G.neighbors(currentV.first)) {
			double weight = 0;
			G.getWeight(currentV.first, adjV, weight);
			double alternativePathDist = distances[currentV.first] + weight;

			if (alternativePathDist < distances[adjV]) {
				distances[adjV] = alternativePathDist;
				predecessors[adjV] = currentV.first;
				unvisitedQueue.push(make_pair(adjV, alternativePathDist));
			}
		}
	}
}

//
// shortestPath2() - for creative component
//
// This function just makes some checks for edge cases when
// running the dijkstraAlgo. It calls dijkstraAlgo(), getPath.
// It returns if path finding was successful or not.
//
bool shortestPath2(graph<long long, double> G,
	vector<long long> userNodes, long long endV) {
	map<long long, double> distances;
	map<long long, long long> predecessors;
	set<string> unreachableBuildings;

	for (unsigned int i = 0; i < userNodes.size(); i++) {
		for (unsigned int j = 0; j < userNodes.size(); j++) {
			dijkstraAlgo(G, distances, predecessors, userNodes[i], userNodes[j]);

			if (distances[userNodes[j]] >= INF) {
				cout << "Sorry, destination unreachable." << endl;
				return true;
			}

			distances.clear();
			predecessors.clear();
		}
	}

	for (auto startV : userNodes) {
		dijkstraAlgo(G, distances, predecessors, startV, endV);
		if (distances[endV] >= INF) {
			cout << "At least one person was unable to reach the destination building.";
			cout << " Finding next closest building..." << endl;
			cout << endl;
			return false;
	  }
		distances.clear();
		predecessors.clear();
	}

	int i = 1;
	for (auto startV : userNodes) {
		dijkstraAlgo(G, distances, predecessors, startV, endV);
		cout << "Person " << i << "'s distance to dest: " << distances[endV]
		<< " miles" << endl;
		getPath(predecessors, distances, startV, endV);
		cout << endl;
		i++;
	}

	return true;
}

//
// searchBuilding()
//
// This function finds if the given string is a building on
// map or not and if so then it returns that building as a reference.
// It returns true if the building was found.
//
bool searchBuilding(BuildingInfo& building,
								vector<BuildingInfo> Buildings, string query) {
	for (auto b : Buildings) {
		if (query == b.Abbrev) {
			building = b;
			return true;
		}
	}

	for (auto b : Buildings) {
		if (b.Fullname.find(query) != string::npos) {
			building = b;
			return true;
		}
	}

	return false;
}

//
// creative()
//
// This function finds the shortest path to a centroid from multiple people at
// different locations.
//
void creative(map<long long, Coordinates>& Nodes, vector<FootwayInfo>& Footways,
    vector<BuildingInfo>& Buildings, graph<long long, double> G) {
	string userBuilding = "";
	vector<BuildingInfo> usersBuildings; // stores al buildings of multiple users
	int c = 2;
	BuildingInfo b;

	cout << endl;
  cout << "Enter person 1's building (partial name or abbreviation), or #> ";
  getline(cin, userBuilding);

	while (userBuilding != "#") {
		if (!searchBuilding(b, Buildings, userBuilding)) {
			cout << "Person 1's building not found" << endl;
			c--;
		} else {
			usersBuildings.push_back(b);
		}

		cout << "Enter person " << c << "'s building (partial name or abbreviation), or #> ";
		getline(cin, userBuilding);

		c++;
	}

	set<string> S;  // set of unreachable Buildings
	bool flag = false;  // to indicate that shortestPath was no successful
	while (!0) {
		BuildingInfo centerB = centralBuilding(usersBuildings, Buildings, S);
		vector<long long> userNodes = findNodes(usersBuildings, Nodes, Footways);
		long long destination = nearestNode(centerB, Nodes, Footways);
		printData(usersBuildings, userNodes, destination, centerB, Nodes, flag);

		if (!shortestPath2(G, userNodes, destination)) {
			S.insert(centerB.Fullname);
			flag = true;
			continue;
		}
		break;
	}
}

//
// buildGraph()
//
// builds the given graph. Takes nodes map and
// graph as parameter. Returns the graph.
//
void buildGraph(graph<long long, double>& G,
							map<long long, Coordinates> Nodes,
							vector<FootwayInfo> Footways) {
	for (auto e : Nodes) {
		G.addVertex(e.first);
	}

	for (auto e : Footways) {
		for (unsigned int i = 0; i < e.Nodes.size() - 1; i++) {
			Coordinates c1 = Nodes[e.Nodes[i]];
			Coordinates c2 = Nodes[e.Nodes[i+1]];

			double distance = distBetween2Points(c1.Lat, c1.Lon, c2.Lat, c2.Lon);

			G.addEdge(e.Nodes[i], e.Nodes[i+1], distance);
			G.addEdge(e.Nodes[i+1], e.Nodes[i], distance);
		}
	}
}

//
// nearestBuilding()
//
// This function finds the nearest Building from the given
// coordinates of a midpoint. It returns that building.
// It takes the vector of Buildings and coordinates of
// midpoint as paameters.
//
BuildingInfo nearestBuilding(vector<BuildingInfo> Buildings,
														 Coordinates midpoint) {
	double min = INF;
	BuildingInfo center;

	for (auto b : Buildings) {
		double distance = distBetween2Points(midpoint.Lat, midpoint.Lon,
													b.Coords.Lat, b.Coords.Lon);

		if (distance < min) {
			min = distance;
			center = b;
		}
	}

	return center;
}

//
// shortestPath()
//
// This function just makes some checks for edge cases when
// running the dijkstraAlgo. It calls dijkstraAlgo(), getPath.
// It returns if path finding was successful or not.
//
bool shortestPath(graph<long long, double> G, long long startV1,
									long long startV2, long long endV) {
	map<long long, double> distances1;
	map<long long, double> distances2;
	map<long long, long long> predecessors1;
	map<long long, long long> predecessors2;
	set<string> unreachableBuildings;

	dijkstraAlgo(G, distances1, predecessors1, startV1, startV2);
	if (distances1[startV2] >= INF) {
		cout << "Sorry, destination unreachable." << endl;
		return true;
	}
	distances1.clear();
	predecessors1.clear();
	dijkstraAlgo(G, distances1, predecessors1, startV1, endV);
	dijkstraAlgo(G, distances2, predecessors2, startV2, endV);

	if (distances1[endV] >= INF || distances2[endV] >= INF) {
		cout << "At least one person was unable to reach the destination building.";
		cout << " Finding next closest building..." << endl;
		cout << endl;
		return false;
	}
	cout << "Person 1's distance to dest: " << distances1[endV]
		<< " miles" << endl;
	getPath(predecessors1, distances1, startV1, endV);
	cout << endl;
	cout << "Person 2's distance to dest: " << distances2[endV]
		<< " miles" << endl;
	getPath(predecessors2, distances2, startV2, endV);
	return true;
}

//
// buildingCenter()
//
// this function finds the nearest building from the coordinates
// of midpoint between two buildings. It takes two BildingInfos
// and a vector of buildings as parameter. It returns the nearest
// building from midpoint.
//
BuildingInfo buildingCenter(BuildingInfo b1, BuildingInfo b2,
							vector<BuildingInfo> Buildings, set<string> S) {
	Coordinates midpoint = centerBetween2Points(b1.Coords.Lat, b1.Coords.Lon,
																							 b2.Coords.Lat, b2.Coords.Lon);

	double min = INF;
	BuildingInfo center;

	for (auto b : Buildings) {
		double distance = distBetween2Points(midpoint.Lat, midpoint.Lon,
																				 b.Coords.Lat, b.Coords.Lon);

		if (distance < min && S.count(b.Fullname) == 0) {
			min = distance;
			center = b;
		}
	}

	return center;
}

//
// printData()
//
// This function all data to the console. It takes two buildingInfos,
// and center buildingInfo, three nodes and a map of nodes as
// parameter.
//
void printData(BuildingInfo b1, BuildingInfo b2, BuildingInfo center,
	long long n1, long long n2, long long dest,
	map<long long, Coordinates> Nodes, bool flag) {
	if (flag == true) {
		cout << "New destination building:" << endl;
		cout << " " << center.Fullname << endl;
		cout << " (" << center.Coords.Lat << ", " << center.Coords.Lon << ")" << endl;
		cout << "Nearest destination node:" << endl;
		cout << " " << dest << endl;
		cout << " (" << Nodes[dest].Lat << ", " << Nodes[dest].Lon << ")" << endl;
		cout << endl;
	} else {
		cout << "Person 1's point:" << endl;
		cout << " " << b1.Fullname << endl;
		cout << " (" << b1.Coords.Lat << ", " << b1.Coords.Lon << ")" << endl;

		cout << "Person 2's point:" << endl;
		cout << " " << b2.Fullname << endl;
		cout << " (" << b2.Coords.Lat << ", " << b2.Coords.Lon << ")" << endl;

		cout << "Destination Building:" << endl;
		cout << " " << center.Fullname << endl;
		cout << " (" << center.Coords.Lat << ", " << center.Coords.Lon << ")" << endl;
		cout << endl;

		cout << "Nearest P1 node:" << endl;
		cout << " " << n1 << endl;
		cout << " (" << Nodes[n1].Lat << ", " << Nodes[n1].Lon << ")" << endl;

		cout << "Nearest P2 node:" << endl;
		cout << " " << n2 << endl;
		cout << " (" << Nodes[n2].Lat << ", " << Nodes[n2].Lon << ")" << endl;

		cout << "Nearest destination node:" << endl;
		cout << " " << dest << endl;
		cout << " (" << Nodes[dest].Lat << ", " << Nodes[dest].Lon << ")" << endl;
		cout << endl;
	}
}

//
// application()
//
// This function is the pivotal function that calls other functions
// to help two people find the shortest path to a midpoint between them.
// It takes a map of Nodes, vector of all footways, vector of all buildings,
// and graph as parameter.
//
void application(
    map<long long, Coordinates>& Nodes, vector<FootwayInfo>& Footways,
    vector<BuildingInfo>& Buildings, graph<long long, double> G) {
  string person1Building, person2Building;

  cout << endl;
  cout << "Enter person 1's building (partial name or abbreviation), or #> ";
  getline(cin, person1Building);

  while (person1Building != "#") {
    cout << "Enter person 2's building (partial name or abbreviation)> ";
    getline(cin, person2Building);
		cout << endl;

		BuildingInfo b1, b2;

		if (!searchBuilding(b1, Buildings, person1Building)) {
			cout << "Person 1's building not found" << endl;
			cout << "Enter person 1's building (partial name or abbreviation), or #> ";
			getline(cin, person1Building);
			continue;
		} else if (!searchBuilding(b2, Buildings, person2Building)) {
			cout << "Person 2's building not found" << endl;
			cout << "Enter person 1's building (partial name or abbreviation), or #> ";
			getline(cin, person1Building);
			continue;
		}
		set<string> S;  // set of unreachable Buildings
		bool flag = false;  // to indicate that shortestPath was no successful
		while (!0) {
			BuildingInfo centerB = buildingCenter(b1, b2, Buildings, S);
			long long p1Node = nearestNode(b1, Nodes, Footways);
			long long p2Node = nearestNode(b2, Nodes, Footways);
			long long destination = nearestNode(centerB, Nodes, Footways);
			printData(b1, b2, centerB, p1Node, p2Node, destination, Nodes, flag);
			if (!shortestPath(G, p1Node, p2Node, destination)) {
				S.insert(centerB.Fullname);
				flag = true;
				continue;
			}
			break;
		}

    cout << endl;
    cout << "Enter person 1's building (partial name or abbreviation), or #> ";
    getline(cin, person1Building);
  }
}

int main() {
  // maps a Node ID to it's coordinates (lat, lon)
  map<long long, Coordinates>  Nodes;
  // info about each footway, in no particular order
  vector<FootwayInfo>          Footways;
  // info about each building, in no particular order
  vector<BuildingInfo>         Buildings;
  XMLDocument                  xmldoc;

  cout << "** Navigating UIC open street map **" << endl;
  cout << endl;
  cout << std::setprecision(8);

  string def_filename = "map.osm";
  string filename;

  cout << "Enter map filename> ";
  getline(cin, filename);

  if (filename == "") {
    filename = def_filename;
  }

  //
  // Load XML-based map file
  //
  if (!LoadOpenStreetMap(filename, xmldoc)) {
    cout << "**Error: unable to load open street map." << endl;
    cout << endl;
    return 0;
  }

  //
  // Read the nodes, which are the various known positions on the map:
  //
  int nodeCount = ReadMapNodes(xmldoc, Nodes);

  //
  // Read the footways, which are the walking paths:
  //
  int footwayCount = ReadFootways(xmldoc, Footways);

  //
  // Read the university buildings:
  //
  int buildingCount = ReadUniversityBuildings(xmldoc, Nodes, Buildings);

  //
  // Stats
  //
  assert(nodeCount == (int)Nodes.size());
  assert(footwayCount == (int)Footways.size());
  assert(buildingCount == (int)Buildings.size());

  cout << endl;
  cout << "# of nodes: " << Nodes.size() << endl;
  cout << "# of footways: " << Footways.size() << endl;
  cout << "# of buildings: " << Buildings.size() << endl;


  graph<long long, double> G;
	buildGraph(G, Nodes, Footways);


  cout << "# of vertices: " << G.NumVertices() << endl;
  cout << "# of edges: " << G.NumEdges() << endl;
  cout << endl;

  //
  // Menu
  //
  string userInput;
  cout << "Enter \"a\" for the standard application or "
        << "\"c\" for the creative component application> ";
  getline(cin, userInput);
  if (userInput == "a") {
    // TO DO: add argument for the graph you make.
    application(Nodes, Footways, Buildings, G);
  } else if (userInput == "c") {
    // TO DO: add arguments
    creative(Nodes, Footways, Buildings, G);
  }
  //
  // done:
  //
  cout << "** Done **" << endl;
  return 0;
}

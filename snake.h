#pragma once

// THIS FILE CONTAINS DECLARATIONS OF SOME GLOBAL VARIABLES AND FUNCTIONS DEFINED EITHER IN "snake.cpp" OR IN "main.cpp"
// THIS FILE MUST BE INCLUDED IN "main.cpp" and "snake.cpp" AS BOTH USE THESE GLOBAL VARIABLES AND FUNCTIONS 
using namespace std;

extern Table2D<RGB> image;  // loaded by GUI in main.cpp
extern bool closedContour;    // a flag indicating if contour was closed - set in snake.cpp
extern vector<Point> contour;   // a list of control points for a "snake" - set in snake.cpp
extern int dP;   // spacing between control points (in pixels) - set in snake.cpp

void addToContour(Point click);
void addToContourLast(Point click);
void addRepulseNudge(Point click);
void addAttractNudge(Point click);
void reset_segm();
double DP_move();
void DP_converge();
double computeElasticEnergyDiscrete(vector<Point> cont);
void calculateExtEnergyTables();

extern bool view; // defined in main.cpp (boolean flag set by a check box)
void draw(Point mouse = Point(-1,-1)); // defined in main.cpp, but it is also called in snake.cpp for visualisation 

const double INFTY=1.e20; 

class MyPoint : Point {
	vector<Point> neighbours;
	double energyValue;
	Point bestPoint; //this will be the point with the best energy. 
	Point pathLocation; //location in a table of where the next MyPoint will be for the lowest energy, after the viterbi algorithm is run.
public:
	MyPoint(){};
	/*
	* Default constructor which takes a starting point and will compute all of its neighbours. 
	*/
	MyPoint(Point start) {
		neighbours.clear();
		neighbours.push_back(start);
		neighbours.push_back(Point(start.x - 1, start.y - 1));
		neighbours.push_back(Point(start.x - 1, start.y));
		neighbours.push_back(Point(start.x - 1, start.y + 1));
		neighbours.push_back(Point(start.x, start.y + 1));
		neighbours.push_back(Point(start.x + 1, start.y + 1));
		neighbours.push_back(Point(start.x + 1, start.y));
		neighbours.push_back(Point(start.x + 1, start.y - 1));
		neighbours.push_back(Point(start.x, start.y - 1));

		pathLocation = Point(0, 0);
		bestPoint = start;
		energyValue = INFTY;
	};
	/**
	* normal constructor that takes a fixed point number to determine which state should be fixed
	*/
	MyPoint(int fixed, Point start) {
		neighbours.clear();
		for(int i = 0; i < 9; i ++) {
			if(fixed == 0) {
				neighbours.push_back(start);
			} else if(fixed == 1) {
				neighbours.push_back(Point(start.x - 1, start.y - 1));
			} else if(fixed == 2) {
				neighbours.push_back(Point(start.x - 1, start.y));
			} else if(fixed == 3) {
				neighbours.push_back(Point(start.x - 1, start.y + 1));
			} else if(fixed == 4) {
				neighbours.push_back(Point(start.x, start.y + 1));
			} else if(fixed == 5) {
				neighbours.push_back(Point(start.x + 1, start.y + 1));
			} else if(fixed == 6) {
				neighbours.push_back(Point(start.x + 1, start.y));
			} else if(fixed == 7) {
				neighbours.push_back(Point(start.x + 1, start.y - 1));
			} else if(fixed == 8) {
				neighbours.push_back(Point(start.x, start.y - 1));
			}
		}

		pathLocation = Point(0, 0);
		bestPoint = start;
		energyValue = INFTY;
	};
	/*
	* for initiating a start point and an energy but a path hasn't been decided yet. This is run at the start of the algorithm
	*/
	MyPoint(Point start, double energy) {
		neighbours.clear();
		neighbours.push_back(start);
		neighbours.push_back(Point(start.x - 1, start.y - 1));
		neighbours.push_back(Point(start.x - 1, start.y));
		neighbours.push_back(Point(start.x - 1, start.y + 1));
		neighbours.push_back(Point(start.x, start.y + 1));
		neighbours.push_back(Point(start.x + 1, start.y + 1));
		neighbours.push_back(Point(start.x + 1, start.y));
		neighbours.push_back(Point(start.x + 1, start.y - 1));
		neighbours.push_back(Point(start.x, start.y - 1));

		pathLocation = Point(0, 0);
		bestPoint = start;
		energyValue = energy;
	};
	/*
	* used when a state is needing to point to another state, that when combined, share the shortest amount of energy 
	*/
	MyPoint(Point start, Point next, double energy){
		neighbours.clear();
		neighbours.push_back(start);
		neighbours.push_back(Point(start.x - 1, start.y - 1));
		neighbours.push_back(Point(start.x - 1, start.y));
		neighbours.push_back(Point(start.x - 1, start.y + 1));
		neighbours.push_back(Point(start.x, start.y + 1));
		neighbours.push_back(Point(start.x + 1, start.y + 1));
		neighbours.push_back(Point(start.x + 1, start.y));
		neighbours.push_back(Point(start.x + 1, start.y - 1));
		neighbours.push_back(Point(start.x, start.y - 1));
		
		pathLocation = next;
		bestPoint = start;
		energyValue = energy;
	};
	vector<Point> getNeighbours() {return neighbours;};
	Point getNextLocation() {return pathLocation;}
	double getEnergyValue() {return energyValue;}
	Point getBestPoint() {return bestPoint;}
};


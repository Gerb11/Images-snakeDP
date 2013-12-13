#include <iostream>     // for cout, rand
#include <vector>       // STL data structures
#include "Basics2D.h"
#include "Table2D.h"
#include "Math2D.h"
#include "Image2D.h"
#include "snake.h"
#include <cmath>


// GLOBAL PARAMETERS AND SPECIALISED DATA TYPES
int dP=4;    // default value for spacing between control points (in pixels)
double localMin = INFTY;
bool push = false;
bool pull = false;
bool closedContour=false; // a flag indicating if contour was closed 

// declarations of global variables 
Table2D<RGB> image; // image is "loaded" from a BMP file by function "image_load" in "main.cpp" 
Table2D<double> imageGradients;
Table2D<double> distanceTransform;

Point nudgePoint;
vector<Point> contour; // list of control points of a "snake"

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////Variables you can tune/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
bool stopShrink = false; //flag indicating to use the keep length internal energy to stop it from shrinking to a point
bool useDT = false;
float alpha = 0.1f; // tuning parameter for the snakes internal energy
float beta = 5.0f; //tuning parameter for the snakes external energy
double nudgeTuner = 1000.0f; //tuning parameter for nudge points (r^2)
double dTPenalty = 15; //tuning parameter for distance transform

double fn(const double t) {double w=0.1; return 1.1/(1.0+w*abs(t));}  // function used for setting "pixel penalties" ("Sensitivity" w is a tuning parameter)

// GUI calls this function when button "Clear" is pressed, or when new image is loaded
// THIS FUNCTION IS FULLY IMPLEMENTED, YOU DO NOT NEED TO CHANGE THE CODE IN IT
void reset_segm()
{
	cout << "resetting 'snake'" << endl;

	// remove all points from the "contour"
	while (!contour.empty()) contour.pop_back();
	closedContour=false;
	localMin = INFTY;

}

// GUI calls this function when a user left-clicks on an image pixel while in "Contour" mode
void addToContour(Point p) 
{
	if (closedContour || dP<=0) return;

	// if contour is empty, append control point p to the "contour", 
	// else (if contour is not empty) append to the "contour" 
	// control points on a straignt line segment from 
	// the end of the current contour to point 'p'
	// The interval between control points is conttroled by parameter 'dP'
	if (contour.empty()) contour.push_back(p);
	else if (!(p==contour.back()))
	{ 
		Point last = contour.back();
		Point v = (p - last); 
		int n = 1 + (int) (v.norm()/dP);
		for (int i=1; i<=n; i++) contour.push_back(last + v*(i/((double)n)));
	}
	cout << "contour size = " << contour.size() << endl;
}

// GUI calls this function when a user right-clicks on an image pixel while in "Contour" mode
void addToContourLast(Point p)
{
	if (closedContour || contour.empty() || dP<=0) return;

	addToContour(p);
	// adding the "closing" interval connecting to the first control point
	addToContour(contour[0]);
	contour.pop_back(); // removes unnecessary copy of point contour[0] at the back!!!!

	closedContour = true;
	draw();
	
}

void addAttractNudge(Point click)
{
	cout << "                     ...attract to pixel p=(" << click.x << "," << click.y << ")" << endl;
	// the code below needs to be replaced, it is for fun only :)
	unsigned i, n = (unsigned) contour.size();
	nudgePoint = click;
	pull = true;
	DP_move();
	pull = false;
}

void addRepulseNudge(Point click)
{
	cout << "                     ...repulse from pixel p=(" << click.x << "," << click.y << ")" << endl;
	unsigned i, n = (unsigned) contour.size();
	nudgePoint = click;
	push = true;
	DP_move();
	push = false;
}

///////////////////////////////////////////////////////////////
// DP_move() is a function that computes one optimal move for 
// the snake's control points using DP (Viterbi algorithm)
double DP_move()
{
	if(closedContour) {
		contour.push_back(contour[0]); //for fixed point calculations
	}
	int n = ((int) contour.size());
	
	vector<vector<MyPoint>> energyTotalPath(n); //table representation
	for(int i = 0; i < n; i++) {
		energyTotalPath[i] = vector<MyPoint>(9);
		for(int j = 0; j < 9; j++) {
			energyTotalPath[i][j] = MyPoint(Point(0,0), INFTY);
		}
	}
	
	vector<Point> firstOrder;

	for(int foo = 0; foo < 9; foo++) {//run algorithm 9 times for each of states of the first point being fixed. 
		for (int i = 0; i < n - 1; i++) {
			MyPoint temp;
			//the following 2 if statements are to set all of the neighbours to the same fixed state for closed snake calculations
			if(i == 0 && closedContour) {
				temp = MyPoint(foo, contour[i]);
			} else {
				temp = MyPoint(contour[i]);
			}
			MyPoint temp1;
			if(i + 1 == n - 1 && closedContour) {
				temp1 = MyPoint(foo, contour[i + 1]);
			} else {
				temp1 = MyPoint(contour[i + 1]);
			}

			vector<Point> neighboursTemp = temp.getNeighbours();
			vector<Point> neighboursTemp1 = temp1.getNeighbours();

			for(int j = 0; j < neighboursTemp1.size(); j++) {
				for(int k = 0; k < neighboursTemp.size(); k ++) {
					Point tempState = neighboursTemp[k];
					Point tempState1 = neighboursTemp1[j];

					if(i == 0 && foo == 0) {//fill first row with k states, as they need to have no energy. Only do this on first pass through
						energyTotalPath [i][k] = MyPoint(tempState, 0);
					}

					firstOrder.clear();
					firstOrder.push_back(tempState);
					firstOrder.push_back(tempState1);

					double energyBetweenPoints = computeElasticEnergyDiscrete(firstOrder);

					if(energyTotalPath[i + 1][j].getEnergyValue() > energyBetweenPoints + energyTotalPath[i][k].getEnergyValue()) { //energy is less between these 2 states than anything that came before
						//creates point with tempstate1 being best energy so far, and Point(i,k) being where the next point is in the table,
						//setting that best energy point to the energy between the paths and the energy at the last point.
						energyTotalPath[i + 1][j] = MyPoint(tempState1, Point(k, i), energyBetweenPoints + energyTotalPath[i][k].getEnergyValue()); 
						if(i == 0) {//for the pass through's of the algorithm that are not the first, and have different fixed states for the first point
							energyTotalPath[i][k] = MyPoint(tempState, 0);
						}
					}
				}
			}
		}
	}

	if(closedContour) {
		contour.pop_back();
	}

	if(contour.empty()) {
		return localMin;
	}

	double cheapest = INFTY;
	int cheapestLocation = 0;
	for(int i = 0; i < 9; i++) {//time to find the cheapest path
		if(energyTotalPath[n - 1][i].getEnergyValue() < cheapest) {
			cheapest = energyTotalPath[n - 1][i].getEnergyValue(); 
			cheapestLocation = i;
		}
	}

	if(cheapest < localMin || push || pull) {
		for(int i = n - 1; i >= 0; i--) {//moving the points
			MyPoint temp = energyTotalPath[i][cheapestLocation];
			contour[i] = temp.getBestPoint();
			cheapestLocation = temp.getNextLocation().x; 
		}
		localMin = cheapest;
	}

	return localMin;
}

///////////////////////////////////////////////////////////////
// DP_converge() is a function that runs DP moves for a snake
// until convergence to a local optima position
void DP_converge()
{
	cout << "Converging until there is a local minimum\n";
	while (true) {
		double eng = localMin;
		if(eng == DP_move()) {
			break;
		}
	}
	cout << "Convergence done\n";
}

/*
* computes, from the list of points given, the internal energy using the discrete method 
* as well as only the elastic term for the internal energy. 
* Uses alpha global var for tuning purposes. Will return a double of the energy for the contour
*/
double computeElasticEnergyDiscrete(vector<Point> cont) {
	int n = (int) cont.size();
	double intEnergy = 0;
	double extEnergy = 0;
	double engChange = 0;
	double tempR = pow(nudgeTuner, 2);

	for(int i = 0; i < n; i ++) {
		Point temp = cont[i];
		Point temp1; 
		if(i != n-1) {
			temp1 = cont[i+1];
			if(stopShrink) {
				intEnergy += (pow((sqrt(pow((temp1.x - temp.x), 2) + pow((temp1.y - temp.y), 2)) - dP), 2));
			} else {
				intEnergy += (pow((temp1.x - temp.x), 2) + pow((temp1.y - temp.y), 2));
			}
		}

		if(useDT) {
			extEnergy -= distanceTransform[temp.x][temp.y];
		} else {
			extEnergy += imageGradients[temp.x][temp.y];
			extEnergy *= beta; //tuning parameter
		}

		engChange += (tempR) / ((pow(temp.x - nudgePoint.x, 2)) + (pow(temp.y - nudgePoint.y, 2)));

	}
	intEnergy *= alpha; //tuning parameter

	if(push) {
		return intEnergy - extEnergy + engChange;
	} else if(pull){
		return intEnergy - extEnergy - engChange;
	} else {
		return intEnergy - extEnergy;
	}
}

void calculateExtEnergyTables() {
	imageGradients.reset(image.getWidth(), image.getHeight(), INFTY);
	distanceTransform.reset(image.getWidth(), image.getHeight(), INFTY);
	for(int i = 1; i < image.getWidth() - 1; i ++) {
		for(int j = 1; j < image.getHeight() - 1; j++) {
			double iY = I(image[i][j + 1]) - I(image[i][j - 1]);
			double iX = I(image[i + 1][j]) - I(image[i - 1][j]);
			iY = pow(iY, 2);
			iX = pow(iX, 2);
			imageGradients[i][j] =  iY + iX;
			distanceTransform[i][j] = fn(imageGradients[i][j]);
		}
	}

	for(int i = 1; i < distanceTransform.getWidth(); i++) {
		for(int j = 1; j < distanceTransform.getHeight(); j ++) {
			distanceTransform[i][j] = min(min((distanceTransform[i - 1][j] + 1) * dTPenalty, (distanceTransform[i][ j - 1] + 1) * dTPenalty), distanceTransform[i][j] * dTPenalty); 
		}
	}

	for(int i = distanceTransform.getWidth() - 2; i >= 0; i--) {
		for(int j = distanceTransform.getHeight() - 2; j >= 0; j--) {
			distanceTransform[i][j] = min(min((distanceTransform[i + 1][j] + 1) * dTPenalty, (distanceTransform[i][j + 1] + 1) * dTPenalty), distanceTransform[i][j] * dTPenalty);
		}
	}
}

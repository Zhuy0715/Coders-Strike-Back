#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;
#define PI 3.14159265

/**
 * Auto-generated code below aims at helping you parse
 * the standard input according to the problem statement.
 **/
class Point {
public:
	int x, y;

  Point(){
    x = 0;
    y = 0;
  }

  Point( int x, int y ){
    x = x;
    y = y;
  }

  Point( const Point& p ){
     x = p.x;
     y = p.y;
   }

	Point& operator=(const Point& p) {
		x = p.x;
		y = p.y;
		return *this;
	}

  bool operator==( const Point& p ){
     return ( x == p.x ) && ( y == p.y );
   }

  bool operator!=( const Point& p ){
     return ( x != p.x ) || ( y != p.y );
   }


	Point operator+(Point& p) {
		return Point(x + p.x, y + p.y);
	}
	Point operator-(Point& p) {
		return Point(x - p.x, y - p.y);
	}

	Point& operator+=(Point& v) {
		x += v.x;
		y += v.y;
		return *this;
	}
	Point& operator-=(Point& v) {
		x -= v.x;
		y -= v.y;
		return *this;
	}

	Point operator+(double s) {
		return Point(x + s, y + s);
	}
	Point operator-(double s) {
		return Point(x - s, y - s);
	}
	Point operator*(double s) {
		return Point(x * s, y * s);
	}
	Point operator/(double s) {
		return Point(x / s, y / s);
	}


	Point& operator+=(double s) {
		x += s;
		y += s;
		return *this;
	}
	Point& operator-=(double s) {
		x -= s;
		y -= s;
		return *this;
	}
	Point& operator*=(double s) {
		x *= s;
		y *= s;
		return *this;
	}
	Point& operator/=(double s) {
		x /= s;
		y /= s;
		return *this;
	}

	void rotate(double angle) {
		double radian = angle * PI /180;
		double cosAngle = cos(radian);
		double sinAngle = sin(radian);
		float tx = x * cosAngle - y * sinAngle;
		float ty = x * sinAngle + y * cosAngle;
		x = (int) tx;
		y = (int) ty;
	}

	float dist2(Point p){
		return pow(p.x - x,2)+ pow(p.y - y,2);
	}
	float dist(Point p){
		return sqrt(dist2(p));
	}
  float length(){
    return sqrt(x*x + y*y);
  }
  Point& normalize() {
    if (length() == 0) return *this;
		*this *= (1.0 / length());
		return *this;
  }

};

bool checkpointSave(vector<Point>& checkpoints, int nextCheckpointX, int nextCheckpointY){
    Point newPoint = Point(nextCheckpointX, nextCheckpointY);
    if(checkpoints.empty()){
      checkpoints.push_back(newPoint);
    }else{
      if(checkpoints.back()!=newPoint){
        if(checkpoints.front()==newPoint && checkpoints.size()>1){
          return true;
        }
        checkpoints.push_back(newPoint);
      }
    }
    return false;
}

int getIndice(vector<Point> checkpoints, int x, int y){
    for (auto it = checkpoints.begin() ; it != checkpoints.end(); ++it){
        if((*it).x==x && (*it).y==y){
            return it - checkpoints.begin();
        }
    }
    return -1;
}

double getMaxDistancePoint(vector<Point> checkpoints){
    double distance = (checkpoints.end()-1)->dist2(*checkpoints.begin());
    double maxDistance = distance;
    for (auto it = checkpoints.begin() ; it != checkpoints.end()-1; ++it){
        distance = (*it).dist2(*(it+1));
        if(distance>maxDistance){
            maxDistance = distance;
        }
    }
    return maxDistance;
}
int main()
{

    // game loop
    bool oneLapFinished = false;
    vector<Point> checkpoints;
    double maxDistance = 0;
    bool boost = false;
    int AngleMin = 1;
    int SlowAngle = 90;
    int SlowRadius = 4*600;
    while (1) {
        int x;
        int y;
        int nextCheckpointX; // x position of the next check point
        int nextCheckpointY; // y position of the next check point
        int nextCheckpointDist; // distance to the next checkpoint
        int nextCheckpointAngle; // angle between your pod orientation and the direction of the next checkpoint
        cin >> x >> y >> nextCheckpointX >> nextCheckpointY >> nextCheckpointDist >> nextCheckpointAngle; cin.ignore();
        int opponentX;
        int opponentY;
        cin >> opponentX >> opponentY; cin.ignore();
        //Save the checkpoints and get the maxDistance for boost
        if(!oneLapFinished){
           oneLapFinished = checkpointSave(checkpoints, nextCheckpointX, nextCheckpointY);
        } else{
            if(maxDistance==0){
                maxDistance = getMaxDistancePoint(checkpoints);
            }
        }
        // Write an action using cout. DON'T FORGET THE "<< endl"
        // To debug: cerr << "Debug messages..." << endl;

        int thrust = 100;
        // You have to output the target position
        // followed by the power (0 <= thrust <= 100)
        // i.e.: "x y thrust"

        if( abs(nextCheckpointAngle) >= AngleMin )
        {
            //search for the steer direction
            Point desiredDirection( nextCheckpointX - x, nextCheckpointY - y );
            desiredDirection = desiredDirection.normalize();

            Point currentDirection = desiredDirection;
            currentDirection.rotate( -nextCheckpointAngle );
            currentDirection = currentDirection.normalize();

            Point steeringDirection = ( desiredDirection - currentDirection );
            steeringDirection = steeringDirection.normalize();
            steeringDirection *= 100;

            nextCheckpointX += steeringDirection.x;
            nextCheckpointY += steeringDirection.y;

            // slow down when angle too big
            if( abs(nextCheckpointAngle) >= SlowAngle)
            {
                thrust = 0;
            }
            else if( nextCheckpointDist < SlowRadius )
            {
                thrust *= 1 - abs(nextCheckpointAngle)/(float)SlowAngle;
            }
        }
        else
        {
            if( !boost && oneLapFinished && nextCheckpointDist + 2000 > maxDistance )
            {
                boost = true;
            }
            else if( nextCheckpointDist < SlowRadius )
            {
                //slow down when pod close to checkpoint
                thrust *= nextCheckpointDist / SlowRadius;
            }
        }
        cout << nextCheckpointX << " " << nextCheckpointY << " ";
          if(boost){
              cout << "BOOST" << endl;
          }else{
              cout << thrust << endl;
          }


    }
}

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
	float x, y;

  Point(){
    this->x = 0;
    this->y = 0;
  }

  Point( float x, float y ){
    this->x = x;
    this->y = y;
  }

  Point( const Point& p ){
     this->x = p.x;
     this->y = p.y;
   }

	Point& operator=(const Point& p) {
		this->x = p.x;
		this->y = p.y;
		return *this;
	}

  bool operator==( const Point& p ){
     return ( this->x == p.x ) && ( this->y == p.y );
   }

  bool operator!=( const Point& p ){
     return ( this->x != p.x ) || ( this->y != p.y );
   }


	Point operator+(Point& p) {
		return Point(this->x + p.x, this->y + p.y);
	}
	Point operator-(Point& p) {
		return Point(this->x - p.x, this->y - p.y);
	}

	Point& operator+=(Point& v) {
		this->x += v.x;
		this->y += v.y;
		return *this;
	}
	Point& operator-=(Point& v) {
		this->x -= v.x;
		this->y -= v.y;
		return *this;
	}

	Point operator+(double s) {
		return Point(this->x + s, this->y + s);
	}
	Point operator-(double s) {
		return Point(this->x - s, this->y - s);
	}
	Point operator*(double s) {
		return Point(this->x * s,this->y * s);
	}
	Point operator/(double s) {
		return Point(this->x / s, this->y / s);
	}


	Point& operator+=(double s) {
		this->x += s;
		this->y += s;
		return *this;
	}
	Point& operator-=(double s) {
		this->x -= s;
		this->y -= s;
		return *this;
	}
	Point& operator*=(double s) {
		this->x *= s;
		this->y *= s;
		return *this;
	}
	Point& operator/=(double s) {
		this->x /= s;
		this->y /= s;
		return *this;
	}

	void rotate(double angle) {
		double radian = angle * PI /180;
		double cosAngle = cos(radian);
		double sinAngle = sin(radian);
		float tx = this->x * cosAngle - this->y * sinAngle;
		float ty = this->x * sinAngle + this->y * cosAngle;
		this->x =  tx;
		this->y =  ty;
	}

	float dist2(Point p){
		return pow(p.x - this->x,2)+ pow(p.y - this->y,2);
	}
	float dist(Point p){
		return sqrt(dist2(p));
	}
  float length(){
    return sqrt(this->x*this->x + this->y*this->y);
  }
  Point normalize() {
      if(length()==0){
          return Point(this->x, this->y);
      }
		return Point(this->x/length(), this->y/length());
  }

};

bool checkpointSave(vector<Point>& checkpoints, int nextCheckpointX, int nextCheckpointY){
    Point newPoint (nextCheckpointX, nextCheckpointY);
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
					 if(maxDistance<nextCheckpointDist){
	             maxDistance = nextCheckpointDist;
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
            Point target (nextCheckpointX - x, nextCheckpointY - y );
            target = target.normalize();

            Point position = target;
            position.rotate( -nextCheckpointAngle );
            position = position.normalize();

            Point velocity = target - position;
            velocity = velocity.normalize();
            velocity *= 100;

            nextCheckpointX += velocity.x;
            nextCheckpointY += velocity.y;

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
            if( !boost && oneLapFinished && (nextCheckpointDist + 2000) > maxDistance )
            {
                boost = true;
            }
            else if( nextCheckpointDist < SlowRadius )
            {
                //slow down when pod close to checkpoint
                thrust *= nextCheckpointDist / (float)SlowRadius;
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

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
 struct Pod{
   Point position;
   float angle;
   Point speed;
   int nextCheckPointId;
 };
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
    vector<Point> checkpoints;
    int laps;
    cin >> laps; cin.ignore();
    int checkpointCount;
    cin >> checkpointCount; cin.ignore();
    cerr << "checkpointCount: " << checkpointCount << endl;
    for (int i = 0; i < checkpointCount; i++) {
        int checkpointX;
        int checkpointY;
        cin >> checkpointX >> checkpointY; cin.ignore();
        Point newPoint(checkpointX,checkpointY);
        checkpoints.push_back(newPoint);
    }
    double maxDistance = getMaxDistancePoint(checkpoints);
    bool boost = false;
    int AngleMin = 1;
    int SlowAngle = 90;
    int SlowRadius = 4*600;

    // game loop
    while (1) {
      Pod pods[2];
        for (int i = 0; i < 2; i++) {
            int x; // x position of your pod
            int y; // y position of your pod
            int vx; // x speed of your pod
            int vy; // y speed of your pod
            int angle; // angle of your pod
            int nextCheckPointId; // next check point id of your pod
            cin >> x >> y >> vx >> vy >> angle >> nextCheckPointId; cin.ignore();
            Point position(x,y);
            Point speed(vx, vy);
            pods[i].position = position;
            pods[i].speed = speed;
            pods[i].angle = angle;
            pods[i].nextCheckPointId = nextCheckPointId;
        }
        for (int i = 0; i < 2; i++) {
            int x2; // x position of the opponent's pod
            int y2; // y position of the opponent's pod
            int vx2; // x speed of the opponent's pod
            int vy2; // y speed of the opponent's pod
            int angle2; // angle of the opponent's pod
            int nextCheckPointId2; // next check point id of the opponent's pod
            cin >> x2 >> y2 >> vx2 >> vy2 >> angle2 >> nextCheckPointId2; cin.ignore();
        }

        // Write an action using cout. DON'T FORGET THE "<< endl"
        // To debug: cerr << "Debug messages..." << endl;


        // You have to output the target position
        // followed by the power (0 <= thrust <= 100)
        // i.e.: "x y thrust"
        int thrust = 100;

        for(int i=0;i<2;i++){
            int nextCheckpointX = checkpoints[pods[i].nextCheckPointId].x;
            int nextCheckpointY = checkpoints[pods[i].nextCheckPointId].y;
            Point target(nextCheckpointX,nextCheckpointY);
          Point destination = target - pods[i].position;
          float rotation = atan2(destination.y, destination.x) *180/PI;
          if(rotation < 0){
              rotation += 360;
          }
          int nextCheckpointAngle = rotation - pods[i].angle;
          cerr << "angle: " << nextCheckpointAngle << endl;
          int nextCheckpointDist = target.dist(pods[i].position);
          cerr << "dist: " << nextCheckpointDist << endl;
          if( abs(nextCheckpointAngle) >= AngleMin )
          {
              //search for the steer direction
              Point target (nextCheckpointX - pods[i].position.x, nextCheckpointY - pods[i].position.y );
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
              if( !boost && (nextCheckpointDist + 2000) > maxDistance )
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
}

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <float.h>

using namespace std;
#define PI 3.14159265
#define podRadius 400.0
#define AngleMin 1
#define SlowAngle 90
#define SlowRadius 4*600

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
     return sqrt(dot());
   }
   Point normalize() {
       if(length()==0){
           return Point(this->x, this->y);
       }
 		return Point(this->x/length(), this->y/length());
   }
   float dot(){
     return this->x*this->x + this->y*this->y;
   }
   float dot(Point p){
     return this->x*p.x + this->y*p.y;
   }

 };
 class Pod{
 public:
   Point position;
   float angle;
   Point speed;
   int nextCheckPointId;
   bool shield = false;
   bool shieldActive = true;
   bool boost = false;
   bool boostActive = true;
   int thrust = 100;
   float radius = 400.0;
   float nextCheckpointAngle = 0;
   int checkpointsPassed = 0;
   float score;
   int mass = 1;

   Pod(){
     this->position = Point(0,0);
     this->angle = 0;
     this->speed = Point(0,0);
     this->nextCheckPointId = 0;
   }
   Pod(Point p, float angle,Point speed, int nextCheckPointId){
     this->position = p;
     this-> angle = angle;
     this->speed = speed;
     this->nextCheckPointId = nextCheckPointId;
   }
   void UpdatePod(float x, float y, float vx, float vy, float angle, int nextCheckPointId){
     this->position = Point(x,y);
     this->angle = angle;
     this->speed = Point(vx, vy);
     this->nextCheckPointId = nextCheckPointId;
   }
   Point InitValue(vector<Point> checkpoints, float maxDistance){
     int nextCheckpointX = checkpoints[this->nextCheckPointId].x;
     int nextCheckpointY = checkpoints[this->nextCheckPointId].y;
     int nextCheckpointDist = target.dist(this->position);
     cerr << "dist: " << nextCheckpointDist << endl;

     Point target(nextCheckpointX,nextCheckpointY);
     Point destination = target - this->position;
     float rotation = atan2(destination.y, destination.x) *180/PI;
     if(rotation < 0){
         rotation += 360;
     }
     int nextCheckpointAngle = rotation - this->angle;
     this->nextCheckpointAngle = nextCheckpointAngle;
     cerr << "angle: " << nextCheckpointAngle << endl;

     if( abs(nextCheckpointAngle) >= AngleMin )
     {
         //search for the steer direction
         Point target (nextCheckpointX - this->position.x, nextCheckpointY - this->position.y );
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
             this->thrust = 0;
         }
         else if( nextCheckpointDist < SlowRadius )
         {
             this->thrust *= 1 - abs(nextCheckpointAngle)/(float)SlowAngle;
         }
     }
     else
     {
         if( !boost && (nextCheckpointDist + 2000) > maxDistance )
         {
             this->boost = true;
         }
         else if( nextCheckpointDist < SlowRadius )
         {
             //slow down when pod close to checkpoint
             this->thrust *= nextCheckpointDist / (float)SlowRadius;
         }
     }
     return Point(nextCheckpointX, nextCheckpointY);
   }

   void Rotate(float angle){
     this->angle = (this->angle + angle)%360;
   }

   void Accelerate(Point target){
     if(!this->shield || !this->shieldActive){
       Point direction(cos(this->angle*PI/180), sin(this->angle*PI/180));
       if(this->boost && this->boostActive){
         this->speed += direction*500;
       }else{
         this->speed += this->thrust*direction;
       }
     }
   }
   void Move(vector<Pod>& pods, vector<Point> checkpoints, int checkpointCount){
     float currentTime = 0.0;
     float maxTime = 1.0;
     while(currentTime<maxTime){
       int indexi = 5;
       int indexj = 5;
       float minTime = CollisionTime(pods, currentTime, maxTime, indexi, indexj);
       this->position += minTime * this->speed;
       if(this->position.dist2(checkpoints[this->nextCheckPointId])<600*600){
         this->nextCheckPointId = (this->nextCheckPointId +1)%checkpointCount;
         this->checkpointsPassed+=1;
       }
       if(indexi!=5 && indexj!=5){
         Rebound(pods[indexi], pods[indexj]);
       }
       currentTime += minTime;
     }
   }
   void Friction(){
     this->speed *= 0.85;
   }
   void Round(){
     this->speed = Point((int) this->speed.x, (int) this->speed.y);
     this->position = Point(round(this->position.x), round(this->position.y));
   }
 };
 // class Checkpoints{
 // public:
 //   Point position;
 //   float radius;
 //   Checkpoints(Point p, float radius = 600){
 //     this-> position = p;
 //     this-> radius = radius;
 //   }
 //
 //   bool collision(Point p1,  Point p2){
 //     Point target = p2 - p1;
 //     float distMinToCheckPoint = abs((p2.x-p1.x)*(p2.y-p1.y) - (p1.x-this->position.x)*(p1.y-this->position.y))/p1.dist(p2);
 //     float distMin = p1.dist2(this->position) - distMinToCheckPoint*distMinToCheckPoint;
 //     if(p1.dist2(p2)<distMin){
 //       return false;
 //     }else{
 //       return true;
 //     }
 //   }
 // };
 float CollisionTime(vector<Pod>& pods, float currentTime, float maxTime, int& indexi, int& indexj){
   float minTime = maxTime - currentTime;
   for(int i =0; i<4; i++){
     for(int j=i+1;j<4;j++){
       Point distance = pod.position - this->position;
       Point differenceSpeed = pod.speed - this->speed;
       float a = differenceSpeed.dot();
       if(a>=0){
         float b = -2*differenceSpeed.dot(distance);
         float c = distance.dot() - 4*this->radius*this->radius;
         float delta = b*b - 4*a*c;
         if(delta >= 0){
           float t = (b-sqrt(delta))/(2*a);
           if(t>=0 && (currentTime+t)<maxTime && t<minTime){
             minTime = t;
             indexi = i;
             indexj = j;
           }
         }
       }
     }
   }
   return minTime;
 }

 void Rebound(Pod& pod1, Pod& pod2){
   if(pod1.shield){
     pod1.mass = 10;
   }
   if(pod2.shield){
     pod2.mass = 10;
   }
   Point vp = pod2.position - pod1.position;
   vp = vp.normalize();

   Point vs = pod2.speed - pod1.speed;

   float m = (pod1.mass*pod2.mass)/(pod1.mass + pod2.mass);
   float k = vs.dot(vp);

   float impulse = clamp(-2*m*k, -120, 120);

   pod1.speed += (-impulse/pod1.mass)* vp;
   pod2.speed += (impulse/pod2.mass)* vp;
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

 void PlayOneTurn(vector<Pod>& pods, vector<Point>& checkpoints, int checkpointCount){
   for(int i=0; i<2;i++){
     //rotate
     if(pods[i].nextCheckpointAngle>18){
       pods[i].rotate(18);
     }else if(pods[i].nextCheckpointAngle<-18){
       pods[i].rotate(-18);
     }else{
       pods[i].rotate(pods[i].nextCheckpointAngle);
     }
     //accurate
     pods[i].Accelerate(checkpoints[pods[i].nextCheckPointId]);
     //move
     pods[i].Move(pods, checkpoints, checkpointCount);
     //Friction
     pods[i].Friction();
     //end
     pods[i].Round();
   }
 }
 int EvaluateScore(vector<Pod> pods, vector<Point> checkpoints, int allCheckpoints){
   for(Pod& p : pods){
     p.score = distanceScore(p, checkpoints);
   }

   int AvancePodIndex = (pods[0].score > pods[1].score) ? 0 : 1;
   Pod& AvancePod = pods[firstPodIndex];
   Pod& BotherPod = pods[1-firstPodIndex];

   Pod& opponent = (pods[2].score > pods[3].score) ? pods[2] : pods[3];

   if(AvancePod.checkpointsPassed > allCheckpoints){
     return DBL_MAX;
   }else if(opponent.checkpointsPassed>allCheckpoints){
     return DBL_MIN;
   }

   int scoreDifference = AvancePod.score - opponent.score;

   return scoreDifference;
 }
 auto distanceScore(Pod& p, vector<Point> checkpoints)
 {
     int coefficient = 30000;
     int distance = dist(p.position, checkpoints[p.nextCheckPointId]);
     return coefficient*p.checkpointsPassed - distance;
 };
void Simulation(vector<Pod>& pods, vector<Point>& checkpoints, int allCheckpoints){
  vector<int> scores;
  scores.reserve(simulationNumber);
  vector<Pod> podsCopy = pods;
  for(int i=0;i<simulationNumber;i++){
    PlayOneTurn(podsCopy, checkpoints);
    int score = EvaluateScore(podsCopy, checkpoints, allCheckpoints);
    scores[i] = score;
  }


}
int main()
{
    vector<Point> checkpoints;
    int laps;
    cin >> laps; cin.ignore();
    int checkpointCount;
    int allCheckpoints = laps*checkpointCount;
    cin >> checkpointCount; cin.ignore();
    for (int i = 0; i < checkpointCount; i++) {
        int checkpointX;
        int checkpointY;
        cin >> checkpointX >> checkpointY; cin.ignore();
        Point newPoint(checkpointX,checkpointY);
        checkpoints.push_back(newPoint);
    }
    double maxDistance = getMaxDistancePoint(checkpoints);
    vector<Pod> pods;
    pods.reserve(4);

    // game loop
    while (1) {

        for (int i = 0; i < 2; i++) {
            int x; // x position of your pod
            int y; // y position of your pod
            int vx; // x speed of your pod
            int vy; // y speed of your pod
            int angle; // angle of your pod
            int nextCheckPointId; // next check point id of your pod
            cin >> x >> y >> vx >> vy >> angle >> nextCheckPointId; cin.ignore();
            pods[i].position = Point(x,y);
            pods[i].speed = Point(vx,vy);
            pods[i].angle = angle;
            if(pods[i].nextCheckPointId != nextCheckPointId){
              pods[i].checkpointsPassed += 1;
            }
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
            pods[i+2].position = Point(x2,y2);
            pods[i+2].speed = Point(vx2,vy2);
            pods[i+2].angle = angle2;
            if(pods[i+2].nextCheckPointId != nextCheckPointId2){
              pods[i+2].checkpointsPassed += 1;
            }
            pods[i+2].nextCheckPointId = nextCheckPointId2;
            cerr << "x :" << x2 << "y: " << y2 << endl;
            cerr << "vx :" << vx2 << "vy: " << vy2 << endl;
            cerr << "next check point: " << nextCheckPointId2 << endl;
        }

        // Write an action using cout. DON'T FORGET THE "<< endl"
        // To debug: cerr << "Debug messages..." << endl;


        // You have to output the target position
        // followed by the power (0 <= thrust <= 100)
        // i.e.: "x y thrust"
        int thrust = 100;

        for(int i=0;i<2;i++){
          pods[i].InitValue(checkpoints, maxDistance);

          cout << nextCheckpointX << " " << nextCheckpointY << " ";
            if(boost){
                cout << "BOOST" << endl;
            }else{
                cout << thrust << endl;
            }
        }

    }
}

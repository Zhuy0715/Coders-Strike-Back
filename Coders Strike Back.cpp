#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <float.h>
#include <stdlib.h>
#include <cstdint>
#include <iomanip>
#include <chrono>
#include <string>
#include <vector>

using namespace std;
#define PI 3.14159265
#define podRadius 400.0
#define AngleMin 1
#define SlowAngle 90
#define SlowRadius 4*600
#define maxRotation 18
#define simulationSteps 4;

/**
 * Auto-generated code below aims at helping you parse
 * the standard input according to the problem statement.
 **/
  inline int fastrand()
 {
     static unsigned int g_seed = 100;
     g_seed = (214013*g_seed+2531011);
     return (g_seed>>16)&0x7FFF;
 }
 inline int rand(int a, int b)
 {
     return (fastrand() % (b-a)) + a;
 }
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

 	Point& operator+=(const Point& v) {
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
   void InitValue(vector<Point> checkpoints, float maxDistance){
     int nextCheckpointX = checkpoints[this->nextCheckPointId].x;
     int nextCheckpointY = checkpoints[this->nextCheckPointId].y;
     Point target(nextCheckpointX,nextCheckpointY);
     int nextCheckpointDist = target.dist(this->position);
     cerr << "dist: " << nextCheckpointDist << endl;


     Point destination = target - this->position;
     float rotation = atan2(destination.y, destination.x) *180/PI;
     if(rotation < 0){
         rotation += 360;
     }
     this-> angle = rotation;
     // int nextCheckpointAngle = rotation - this->angle;
     // this->nextCheckpointAngle = nextCheckpointAngle;
     // cerr << "angle: " << nextCheckpointAngle << endl;
     //
     // if( abs(nextCheckpointAngle) >= AngleMin )
     // {
     //     //search for the steer direction
     //     Point target (nextCheckpointX - this->position.x, nextCheckpointY - this->position.y );
     //     target = target.normalize();
     //
     //     Point position = target;
     //     position.rotate( -nextCheckpointAngle );
     //     position = position.normalize();
     //
     //     Point velocity = target - position;
     //     velocity = velocity.normalize();
     //     velocity *= 100;
     //
     //     nextCheckpointX += velocity.x;
     //     nextCheckpointY += velocity.y;
     //
     //     // slow down when angle too big
     //     if( abs(nextCheckpointAngle) >= SlowAngle)
     //     {
     //         this->thrust = 0;
     //     }
     //     else if( nextCheckpointDist < SlowRadius )
     //     {
     //         this->thrust *= 1 - abs(nextCheckpointAngle)/(float)SlowAngle;
     //     }
     // }
     // else
     // {
     //     if( !boost && (nextCheckpointDist + 2000) > maxDistance )
     //     {
     //         this->boost = true;
     //     }
     //     else if( nextCheckpointDist < SlowRadius )
     //     {
     //         //slow down when pod close to checkpoint
     //         this->thrust *= nextCheckpointDist / (float)SlowRadius;
     //     }
     // }
     // return Point(nextCheckpointX, nextCheckpointY);
   }

 };
struct Movement{
  float rotation;
  float thrust;
  bool shield;
  bool boost;
};
class Step
{
    private:
        vector<Movement> moves = vector<Movement>(2);
    public:
        Movement& operator[](size_t m){ return moves[m]; }
};
void Randomize(Movement& m){
    cerr << "begin random" << endl;
    int i = rand(0,12);
    cerr << i << endl;
    if(i<5){
      m.rotation = -maxRotation + rand(-maxRotation, maxRotation);
    }else if(i<10){
      m.thrust = rand(0,100);
    }else if(i<11){
      m.shield = !m.shield;
    }else{
      m.boost = !m.boost;
    }
  }

class Solution{
public:
  std::vector<Step> podsMovement = vector<Step>(4);
  int score = 0;
  Step& operator[](size_t t){ return podsMovement[t]; }
};
float CollisionTime(vector<Pod>& pods, float currentTime, float maxTime, int& indexi, int& indexj){
   float minTime = maxTime - currentTime;
   for(int i =0; i<4; i++){
     for(int j=i+1;j<4;j++){
       Point distance = pods[i].position - pods[j].position;
       Point differenceSpeed = pods[i].speed - pods[j].speed;
       float a = differenceSpeed.dot();
       if(a>=0){
         float b = -2*differenceSpeed.dot(distance);
         float c = distance.dot() - 4*pods[j].radius*pods[j].radius;
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

   float impulse = clamp(-2*m*k, -120.0f, 120.0f);

   pod1.speed +=  vp * (-impulse/pod1.mass);
   pod2.speed +=  vp * (impulse/pod2.mass);
 }
 void Rotate(std::vector<Pod> pods, Step& m){
   for (int i=0; i<2; i++)
   {
     pods[i].angle = ((int)(pods[i].angle + m[i].rotation))%360;
   }
 }

 void Accelerate(std::vector<Pod> pods, Step& m){
   for (int i=0; i<2; i++)
   {
     if(!m[i].shield || !pods[i].shieldActive){
       Point direction(cos(pods[i].angle*PI/180), sin(pods[i].angle*PI/180));
       if(m[i].boost && pods[i].boostActive){
         pods[i].speed += direction*500;
         pods[i].boost = true;
         pods[i].boostActive = false;
       }else{
         pods[i].speed += direction * m[i].thrust;
       }
     }else{
       pods[i].shield = true;
       pods[i].shieldActive = false;
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
     for(int i=0;i<4;i++){
       pods[i].position += (pods[i].speed * minTime) ;
       if(pods[i].position.dist2(checkpoints[pods[i].nextCheckPointId])<600*600){
         pods[i].nextCheckPointId = (pods[i].nextCheckPointId +1)%checkpointCount;
         pods[i].checkpointsPassed+=1;
       }
     }
     if(indexi!=5 && indexj!=5){
       Rebound(pods[indexi], pods[indexj]);
     }
     currentTime += minTime;
   }
 }
 void Friction(std::vector<Pod> pods){
   for(int i =0; i<4;i++){
     pods[i].speed *= 0.85;
   }
 }
 void Round(std::vector<Pod> pods){
   for(int i=0;i<4;i++){
     pods[i].speed = Point((int) pods[i].speed.x, (int) pods[i].speed.y);
     pods[i].position = Point(round(pods[i].position.x), round(pods[i].position.y));
   }

 }

void PlayOneStep(vector<Pod>& pods, vector<Point>& checkpoints, int checkpointCount, Step& Movement){
    //rotate
    Rotate(pods, Movement);
    //accurate
    Accelerate(pods, Movement);
    //move
    Move(pods, checkpoints, checkpointCount);
    //Friction
    Friction(pods);
    //end
    Round(pods);
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
 auto distanceScore(Pod& p, vector<Point> checkpoints)
 {
     int coefficient = 30000;
     int distance = p.position.dist(checkpoints[p.nextCheckPointId]);
     return coefficient*p.checkpointsPassed - distance;
 };
 double EvaluateScore(vector<Pod> pods, vector<Point> checkpoints, int allCheckpoints){
   for(Pod& p : pods){
     p.score = distanceScore(p, checkpoints);
   }

   int AvancePodIndex = (pods[0].score > pods[1].score) ? 0 : 1;
   Pod& AvancePod = pods[AvancePodIndex];
   Pod& BotherPod = pods[1-AvancePodIndex];

   Pod& opponent = (pods[2].score > pods[3].score) ? pods[2] : pods[3];

   if(AvancePod.checkpointsPassed > allCheckpoints){
     return DBL_MAX;
   }else if(opponent.checkpointsPassed>allCheckpoints){
     return DBL_MIN;
   }

   int scoreDifference = AvancePod.score - opponent.score;

   return scoreDifference;
 }

 void UpdateInput(std::vector<Pod>& pods){
   for (int i = 0; i < 4; i++) {
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
 }
 void ConvertSolutionToOutput(Solution& solution, std::vector<Pod> pods){
   Step m = solution[0];
   for (int i=0; i<2; i++)
   {
       Pod& p = pods[i];
       Movement& mov = m[i];

       // generate coordinates from the rotation
       float angle = ((int) (p.angle + mov.rotation)) % 360;
       float angleRad = angle * PI / 180.f;

       float targetDistance = 10000.f; // compute the target arbitrarily far enough, to avoid rounding errors
       Point direction{ targetDistance*cos(angleRad),targetDistance*sin(angleRad) };
       Point target = p.position + direction;

       cout << round(target.x) << " " << round(target.y) << " ";

       if (mov.shield)
           cout << "SHIELD";
           //p.shieldActive  = false;
       else if (mov.boost)
           cout << "BOOST";
           //p.boostActive = false;
       else
           cout << mov.thrust;

       cout << endl;


   }
 }
class Simulation{
public:
  std::vector<Solution> solutions;
  int solutionNumber = 6;

  Simulation(){
    solutions.resize(solutionNumber);
    InitSimulation();
  }
  Simulation(int solutionNumber){
    this->solutionNumber = solutionNumber;
    solutions.resize(solutionNumber);
    InitSimulation();
  }

  void InitSimulation(){
      cerr << "begin init simulation" << endl;
    for (int s=0; s<solutionNumber; s++){
      for (int m=0; m<simulationSteps; m++){
        for (int i=0; i<2; i++){
            cerr << s << m << i << endl;
            Randomize(solutions[s][m][i]);

        }
      }
    }
  }

  void UpdateNextStep(int s){
      for(int i=1;i<simulationSteps;i++){
        for(int j=0;j<2;j++){
          this->solutions[s][i-1][j] = this->solutions[s][i][j];
          if(i==simulationSteps-1){
            Randomize(solutions[s][i][j]);
          }
        }
      }
  }

  Solution& simulate(std::vector<Pod>& pods, std::vector<Point> checkpoints, int allCheckpoints){
    for(int i=0; i<solutionNumber;i++){
      Solution& s = this->solutions[i];
      UpdateNextStep(i);
      vector<Pod> podsCopy = pods;
      for(int j=0;j<simulationSteps;j++){
        PlayOneStep(podsCopy, checkpoints, allCheckpoints, s.podsMovement[j]);
      }
      int score = EvaluateScore(podsCopy, checkpoints, allCheckpoints);
      s.score = score;
    }
    std::sort( this->solutions.begin(), this->solutions.end(),
               [](const Solution& a, const Solution& b) { return a.score > b.score; }
             );
    return this->solutions[0];
  }
};

int main()
{
    cerr << "begin" << endl;
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
    cerr << "init checkpoints" << endl;
    vector<Pod> pods(4);
    Simulation sim = Simulation();
    int step = 0;
    cerr << "init finish" << endl;
    // game loop
    while (1) {
      // init value ->
      // search for the solution randomly ->
      // play the solution ->
      // evaluate solution ->
      // find the best solution ->
      // convert solution to output

        UpdateInput(pods);
        if(step==0){
          for(int i =0; i<4;i++){
            pods[i].InitValue(checkpoints, maxDistance);
          }
        }
        Solution& sol = sim.simulate(pods, checkpoints, allCheckpoints);
        ConvertSolutionToOutput(sol, pods);
        step ++;


    }
}

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <float.h>


using namespace std;
using namespace std::chrono;

#define SIMULATIONSTEPS 4
#define SOLUTIONCOUNT 6

#define maxThrust 100
#define boostThrust 650
#define maxRotation 18
#define checkpointRadius 600.0
#define podRadius 400.0
#define minImpulse 120.0

#define maxTimeFirstStep 500
#define maxTimeEachStep 75

#define PI 3.14159265
#define eps 0.000001

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

  Point(const Point& p ){
     this->x = p.x;
     this->y = p.y;
   }

 Point& operator=(const Point& p) {
   this->x = p.x;
   this->y = p.y;
   return *this;
 }

  bool operator==( Point& p ){
     return ( this->x == p.x ) && ( this->y == p.y );
   }

  bool operator!=( Point& p ){
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
// -------------------------------------------------------------------

class Pod{
public:
  Point position;
  int angle;
  Point speed;
  int nextCheckPointId;
  int shieldCount = 0;
  bool boostActive = true;
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

};
// -------------------------------------------------------------------

void CountShield(bool shiled, Pod& p)
{
  if (shiled)
    p.shieldCount = 4;
  else if (p.shieldCount > 0)
    --p.shieldCount;
}

int getMass(const Pod& p)
{
    if (p.shieldCount == 4)
        return 10;
    return 1;
}

// -------------------------------------------------------------------

 struct Movement{
   int rotation;
   int thrust;
   bool shield;
   bool boost;
 };
 class Step
 {
     public:
         vector<Movement> moves = vector<Movement>(2);
         Movement& operator[](size_t m){ return moves[m]; }
 };
class Solution
{
    public:
        vector<Step> steps = vector<Step>(SIMULATIONSTEPS);
        Step& operator[](size_t t)             { return steps[t]; }
        int score = 0;
};
// -------------------------------------------------------------------

// https://stackoverflow.com/questions/1640258/need-a-fast-random-generator-for-c
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



class Simulation
{
    private:
        vector<Solution> solutions;
        vector<Point> checkpoints;
        int checkpointCount;
        int allCheckpoints;

    public:
        void InitSimulation()
        {
            InitSolutions();
            BoostAtBegin();
        }
        const vector<Point>& Checkpoints() const { return checkpoints; }

        void InitCheckpoints()
        {
            int laps;
            cin >> laps; cin.ignore();
            cin >> checkpointCount; cin.ignore();

            checkpoints.reserve(checkpointCount);
            for (int i=0; i<checkpointCount; i++)
            {
                int checkpointX, checkpointY;
                cin >> checkpointX >> checkpointY; cin.ignore();
                checkpoints[i] = Point(checkpointX, checkpointY);
            }

            allCheckpoints = laps*checkpointCount;
        }

        Solution& simulate(const vector<Pod>& pods, int time)
        {
          auto start = high_resolution_clock::now();
          auto now = high_resolution_clock::now();
          auto duration = duration_cast<milliseconds>(now -start);
          bool timeout = (duration.count() >= time);
            for (int i=0; i<SOLUTIONCOUNT; ++i)
            {
                Solution& s = solutions[i];
                UpdateNextStep(s);
                vector<Pod> podsCopy = pods;
                for (int i=0; i<SIMULATIONSTEPS; i++)
                {
                    PlayOneStep(podsCopy, s[i]);
                }
                s.score = EvaluateScore(podsCopy);
            }

          while(!timeout){
            for (int i=0; i<SOLUTIONCOUNT; ++i)
            {
                Solution& s = solutions[i+SOLUTIONCOUNT];
                s = solutions[i];
                UpdateRandomStep(s);
                vector<Pod> podsCopy = pods;
                for (int i=0; i<SIMULATIONSTEPS; i++)
                {
                    PlayOneStep(podsCopy, s[i]);
                }
                s.score = EvaluateScore(podsCopy);
            }

            std::sort( solutions.begin(), solutions.end(),
                       [](Solution& a, Solution& b) { return a.score > b.score; }
                     );

             now = high_resolution_clock::now();
             duration = duration_cast<milliseconds>(now -start);
             timeout = (duration.count() >= time);
          }

            return solutions[0];
        }

    private:
      void PlayOneStep(vector<Pod>& pods, Step& step) const
      {
          // reference in expert rules
          Rotate(pods, step);
          Accelerate(pods, step);
          Move(pods);
          Friction(pods);
          Round(pods);
      }

      void Rotate(std::vector<Pod>& pods, Step& m) const{
          for (int i=0; i<2; i++)
          {
              Pod& p = pods[i];
              Movement& mov = m[i];
              p.angle = (p.angle + mov.rotation)%360;
          }
      }
       void Accelerate(std::vector<Pod>& pods, Step& step) const{
          for (int i=0; i<2; i++)
          {
              Pod& p = pods[i];
              Movement& m = step[i];
              CountShield(m.shield, p);
              if(p.shieldCount>0){
                continue;
              }
              const float angleRad = p.angle * PI / 180.0;
              Point direction{cos(angleRad), sin(angleRad)};
              if(m.boost && p.boostActive){
                p.speed += direction *500.0f ;
                p.boostActive = false;
              }else{
                p.speed +=  direction * m.thrust ;
              }
          }
      }
      void Move(vector<Pod>& pods) const
      {
          float currentTime = 0.0;
          float maxTime = 1.0;
          while (currentTime < maxTime)
          {
              Pod* p1 = nullptr;
              Pod* p2 = nullptr;
              float minTime = maxTime - currentTime;

              for (int i=0; i<4; i++)
              {
                  for (int j=i+1; j<4; j++)
                  {
                      const float collisionTime = CollisionTime(pods[i], pods[j]);
                      if ( (currentTime+collisionTime < maxTime) && (collisionTime < minTime) )
                      {
                          minTime = collisionTime;
                          p1 = &pods[i];
                          p2 = &pods[j];
                      }
                  }
              }

              for (Pod& p : pods)
              {
                  p.position += p.speed * minTime ;

                  if (p.position.dist2(checkpoints[p.nextCheckPointId]) < checkpointRadius*checkpointRadius)
                  {
                      p.nextCheckPointId = (p.nextCheckPointId + 1) % checkpointCount;
                      ++p.checkpointsPassed;
                  }
              }

              if (p1 != nullptr && p2 != nullptr)
              {
                  Rebound(*p1, *p2);
              }

              currentTime+=minTime;
          }
      }
      void Friction(vector<Pod>& pods) const
      {
          for (Pod& p : pods)
          {
              p.speed *= 0.85;
          }
      }
      void Round(vector<Pod>& pods) const
      {
          for (Pod& p : pods)
          {
              p.speed = Point{ (int) p.speed.x, (int) p.speed.y };
              p.position = Point{round(p.position.x), round(p.position.y)};
          }
      }
        void InitSolutions()
        {
            solutions.resize(2*SOLUTIONCOUNT);

            for (int s=0; s<SOLUTIONCOUNT; s++)
                for (int t=0; t<SIMULATIONSTEPS; t++)
                    for (int i=0; i<2; i++)
                        Randomize(solutions[s][t][i]);
        }

        void BoostAtBegin()
        {
            const float d = checkpoints[0].dist2(checkpoints[1]); // size>1
            float boostDistanceThreshold = 3000.0 * 3000.0;
            if (d < boostDistanceThreshold)
                return;

            for (int i=0; i<2; i++)
            {
                for (int s=0; s<SOLUTIONCOUNT; ++s)
                    solutions[s][0][i].boost = true;
            }
        }

        void Randomize(Movement& m, bool random = false) const{
          if(!random){
            int r = rand(-2*maxRotation, 3*maxRotation);
            if (r > 2*maxRotation)
                m.rotation = 0;
            else
                m.rotation = clamp(r, -maxRotation, maxRotation);
            r = rand(-0.5f * maxThrust, 2*maxThrust);
            m.thrust = clamp(r, 0, maxThrust);
            if((rand(0,10)>6)){
              m.shield = !m.shield;
            }
            if((rand(0,10)>6)){
              m.boost = !m.boost;
            }
          }else{
            int i = rand(0,12);
            //cerr << i << endl;
            if(i<5){
              int r = rand(-2*maxRotation, 3*maxRotation);
              if (r > 2*maxRotation)
                  m.rotation = 0;
              else
                  m.rotation = clamp(r, -maxRotation, maxRotation);
            }else if(i<10){
              int r = rand(-0.5f * maxThrust, 2*maxThrust);
              m.thrust = clamp(r, 0, maxThrust);
            }else if(i<11){
              if((rand(0,10)>6)){
              m.shield = !m.shield;
            }
            }else{
              if((rand(0,10)>6)){
              m.boost = !m.boost;
            }
            }
          }
        }

        void UpdateNextStep(Solution& s) const
        {
            for (int t=1; t<SIMULATIONSTEPS; t++){
              for (int i=0; i<2; i++){
                s[t-1][i] = s[t][i];
                if(t==SIMULATIONSTEPS-1){
                  Movement& m = s[t][i];
                  Randomize(m);
                }
              }
            }
        }

        void UpdateRandomStep(Solution& s) const
        {
            int k = rand(0, SIMULATIONSTEPS);
            Movement& m = s[k][k%2];
            Randomize(m,true);
        }

        int ComputeScore(Solution& sol, vector<Pod>& pods) const
        {
            vector<Pod> podsCopy = pods;
            for (int i=0; i<SIMULATIONSTEPS; i++)
            {
                PlayOneStep(pods, sol[i]);
            }
            sol.score = EvaluateScore(podsCopy);
            return sol.score;
        }
        int distanceScore(Pod& p) const
        {
          int coefficient = 20000;
          int distance = p.position.dist(checkpoints[p.nextCheckPointId]);
          return coefficient*p.checkpointsPassed - distance;
        }
        int EvaluateScore(vector<Pod> pods) const{
            for (Pod& p : pods)
                p.score = distanceScore(p);

            int AvancePodIndex = (pods[0].score > pods[1].score) ? 0 : 1;
            Pod& AvancePod = pods[AvancePodIndex];
            Pod& BotherPod = pods[1-AvancePodIndex];

            Pod& opponent = (pods[2].score > pods[3].score) ? pods[2] : pods[3];

            if (AvancePod.checkpointsPassed > allCheckpoints){
              return INT8_MAX; // I win
            }else if(opponent.checkpointsPassed>allCheckpoints){
              return -INT8_MAX; // opponent wins
            }

            // check the difference between us
            int scoreDifference = AvancePod.score - opponent.score;

            // try to bother my oppenent
            Point opponentCheckpoint = checkpoints[opponent.nextCheckPointId];
            int botherScore = - BotherPod.position.dist(opponentCheckpoint);

            return 2*scoreDifference + botherScore; // advance is better than bother
        }

        float CollisionTime(Pod& p1, Pod& p2) const
        {
            Point distance = p2.position - p1.position;
            Point differenceSpeed = (p2.speed - p1.speed);

            const float a = differenceSpeed.dot();
            if (a < eps)
                return INT8_MAX;

            const float b = -2.0*distance.dot(differenceSpeed);
            const float c = distance.dot() - 4.0*podRadius*podRadius;

            const float delta = b*b - 4.0*a*c;
            if (delta < 0)
                return INT8_MAX;

            const float t = (b - sqrt(delta)) / (2.0 * a);
            if (t <= eps)
                return INT8_MAX;

            return t;
        }

        void Rebound(Pod& a, Pod& b) const
        {
            // https://en.wikipedia.org/wiki/Elastic_collision#Two-dimensional_collision_with_two_moving_objects
            const float mA = getMass(a);
            const float mB = getMass(b);

            Point vp = b.position - a.position;
            vp = vp.normalize();

            Point differenceSpeed = b.speed - a.speed;

            const float m = (mA * mB) / (mA + mB);
            const float k = differenceSpeed.dot(vp);

            const float impulse = clamp((-2.0 * m * k), -minImpulse, minImpulse);

            a.speed += vp * (-impulse/mA);
            b.speed += vp * (impulse/mB);
        }
};

void UpdateInput(Pod& p){
      int x; // x position of your pod
      int y; // y position of your pod
      int vx; // x speed of your pod
      int vy; // y speed of your pod
      int angle; // angle of your pod
      int nextCheckPointId; // next check point id of your pod
      cin >> x >> y >> vx >> vy >> angle >> nextCheckPointId; cin.ignore();
      p.position = Point(x,y);
      p.speed = Point(vx,vy);
      p.angle = angle;
      if(p.nextCheckPointId != nextCheckPointId){
        p.checkpointsPassed += 1;
      }
      p.nextCheckPointId = nextCheckPointId;
}

void ConvertSolutionToOutput(Solution& solution, vector<Pod>& pods)
{
    Step m = solution[0];
   for (int i=0; i<2; i++)
   {
       Pod& p = pods[i];
       Movement& mov = m[i];

       float angle = ((int) (p.angle + mov.rotation)) % 360;
       float angleRad = angle * PI / 180.0;

       float targetDistance = 10000.0;
       Point direction{ targetDistance*cos(angleRad),targetDistance*sin(angleRad) };
       Point target = p.position + direction;

       cout << round(target.x) << " " << round(target.y) << " ";

       if (mov.shield)
           cout << "SHIELD";
       else if (mov.boost)
           cout << "BOOST";
       else
           cout << mov.thrust;

       cout << endl;


   }
}

 void UpdatePodsStatus(Solution& solution, vector<Pod>& pods)
 {
    Step m = solution[0];
     for (int i=0; i<2; i++)
     {
       Pod& p = pods[i];
       Movement& mov = m[i];
       CountShield(mov.shield, p);

       if (p.shieldCount == 0 && mov.boost)
           p.boostActive = false;
     }
 }

void InitAngle(Pod& p, Point target)
{
    // int nextCheckpointX = checkpoints[this->nextCheckPointId].x;
    // int nextCheckpointY = checkpoints[this->nextCheckPointId].y;
    // Point target(nextCheckpointX,nextCheckpointY);

    Point destination = target - p.position;
    float rotation = atan2(destination.y, destination.x) *180/PI;
    if(rotation < 0){
        rotation += 360;
    }
    p.angle = rotation;
}

int main()
{
    Simulation simulation;
    cerr << "init" << endl;
    simulation.InitCheckpoints();
    simulation.InitSimulation();
    vector<Pod> pods(4);

    int step=0;
    while (1)
    {
      // init value ->
      // search for the solution randomly ->
      // play the solution ->
      // evaluate solution ->
      // find the best solution ->
      // convert solution to output
        for (int i=0; i<4; i++)
        {
            UpdateInput(pods[i]);
            if (step == 0)
                InitAngle(pods[i],simulation.Checkpoints()[0]);
        }

        int availableTime = (step==0) ? maxTimeFirstStep : maxTimeEachStep;
        float timeCoefficient = 0.95;

        Solution& s = simulation.simulate(pods, timeCoefficient*availableTime);
        ConvertSolutionToOutput(s, pods);
        UpdatePodsStatus(s, pods);

        ++step;
    }
}

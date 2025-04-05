#include <bits/stdc++.h>
#include <cmath>
using namespace std;

const double launchAngle = M_PI/4; //Radians
const double initPressure = 413685; //Pascals
const double barrelRadius = 0.038; //Meters
const double barrelLen = 0.4; // Meters
const double shirtMass = 0.15; //Kilograms
const double tankVol = 0.0005; //Meters^3
const double initPos = 0.15; //Meters
const double atmPressure = 101325; //Pascals
const double airSpecificHeatRatio = 1.4;
const double valveFlowCoeff = 10;
const double initTemp = 300; //Kelvin
const double dt = 0.0001; //Seconds
double weight = 9.8*sin(launchAngle)*shirtMass;

//Prototype functions
bool inside(double t);
bool critical(double t);
double acceleration(double t);
double velocity(double t);
double displacement(double t);
double barrelVol(double t);
double molesAirBarrel(double t);
double molesAirTank(double t);
double tankPressure(double t);
double barrelPressure(double t);
double barrelTemp(double t);
double barrelDensity(double t);
double tankDensity(double t);
double volumetricFlowRate(double t);
double airForce(double t);

//Inside barrel?
bool inside(double t){ 
    if (0<=displacement(t) && displacement(t) < barrelLen) return 1;
    else return 0;
}
//Critical flow?
bool critical(double t){
    if (tankPressure(t)-2*barrelPressure(t)>=0) return 1;
    else return 0;
}
double acceleration(double t){ //Meters per second^2
    return (inside(t)*airForce(t)-weight)/shirtMass;
}
double velocity(double t){ //Meters per second
    if(t <= 0) return 0;
    else{
        return velocity(t-dt)+acceleration(t-dt)*dt;
    }
}
double displacement(double t){ //Meters
    if(t<=0) return initPos;
    else{
        return displacement(t-dt)+velocity(t-dt)*dt;
    }
}
double barrelVol(double t){ //Meters^3
    return inside(t)*max(M_PI*barrelRadius*barrelRadius*displacement(t),0.0);
}
double molesAirBarrel(double t){ //Moles
    if(t<=0){
        return barrelVol(0)/(8310*initTemp);
    }
    else{
        return molesAirBarrel(t-dt)+0.000271*barrelDensity(t-dt)*volumetricFlowRate(t-dt)*dt;
    }
}
double molesAirTank (double t){ //Moles
    if(t<=0){
        return initPressure*tankVol/(8310*initTemp*atmPressure);
    }
    else{
        return molesAirTank(t-dt)-0.000271*tankDensity(t-dt)*volumetricFlowRate(t-dt)*dt;
    }
}
double tankPressure(double t){ //Pascals
    return 8310*molesAirTank(t)*initTemp/(atmPressure*tankVol);
}
double barrelPressure(double t){ //Pascals
    return 8310*molesAirBarrel(t)*barrelTemp(t)/(atmPressure*barrelVol(t));
}
double barrelTemp(double t){ //Kelvin
    return initTemp*pow(barrelVol(0)/barrelVol(t), airSpecificHeatRatio-1);
}
double barrelDensity(double t){ //Kilograms per meter^3
    return barrelPressure(t)/(287*barrelTemp(t));
}
double tankDensity(double t){ //Kilograms per meter^3
    return tankPressure(t)/(287*initTemp);
}
double volumetricFlowRate(double t){ //Standard Cubic Feet per Hour
    if (critical(t)){
        return 816*valveFlowCoeff*tankPressure(t)/sqrt(1.8*barrelTemp(t));
    }
    else{
        return 962*valveFlowCoeff*sqrt((tankPressure(t)*tankPressure(t)+barrelPressure(t)*barrelPressure(t))/(1.8*barrelTemp(t)));
    }
}
double airForce(double t){ //Newtons
    return inside(t)*barrelPressure(t)-atmPressure*M_PI*barrelRadius*barrelRadius;
}
int main(){
    double t = 0;
    while(t <= 3){
        cout << t << ' ' << velocity(t) << endl;
        t+=0.01;
    }
}
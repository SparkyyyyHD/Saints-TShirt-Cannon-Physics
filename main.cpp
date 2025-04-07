#include <bits/stdc++.h>
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
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
const double dt = 0.01; //Seconds
double weight = 9.8*sin(launchAngle)*shirtMass;

class memoizer{
    public:
        void cache(double key, double value){
            int scaledKey = static_cast<int>(round(key * 100000));
            m[scaledKey] = value;
        }
        bool inMap(double key){
            int scaledKey = static_cast<int>(round(key * 100000));
            if(m.count(scaledKey)) return true;
            return false;
        }
        double getVal(double key){
            int scaledKey = static_cast<int>(round(key * 100000));
            return m[scaledKey];
        }
    private:
        map<int,double> m;
};

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

//Memoizer
memoizer memAcceleration;
memoizer memVelocity;
memoizer memDisplacement;
memoizer memBarrelVol;
memoizer memMolesAirBarrel;
memoizer memMolesAirTank;
memoizer memTankPressure;
memoizer memBarrelPressure;
memoizer memBarrelTemp;
memoizer memBarrelDensity;
memoizer memTankDensity;
memoizer memVolumetricFlowRate;
memoizer memAirForce;



//Inside barrel?
bool inside(double t){ 
    if (0<=displacement(t) && displacement(t) < barrelLen) return 1;
    else return 0;
}
//Critical flow?
bool critical(double t){
    if (0.528*tankPressure(t) > barrelPressure(t)) return 1;
    else return 0;
}
double acceleration(double t){ //Meters per second^2
    if (memAcceleration.inMap(t)) return memAcceleration.getVal(t);
    memAcceleration.cache(t,(inside(t)*airForce(t)-weight)/shirtMass);
    return memAcceleration.getVal(t);
}
double velocity(double t){ //Meters per second
    if (memVelocity.inMap(t)) return memVelocity.getVal(t);
    if(t <= 0) return 0;
    else{
        memVelocity.cache(t,velocity(t-dt)+acceleration(t-dt)*dt);
        return memVelocity.getVal(t);
    }
}
double displacement(double t){ //Meters
    if (memDisplacement.inMap(t)) return memDisplacement.getVal(t);
    if(t<=0) return initPos;
    else{
        memDisplacement.cache(t,displacement(t-dt)+velocity(t-dt)*dt);
        return memDisplacement.getVal(t);
    }
}
double barrelVol(double t){ //Meters^3
    if (memBarrelVol.inMap(t)) return memBarrelVol.getVal(t);
    memBarrelVol.cache(t, M_PI*barrelRadius*barrelRadius*displacement(t));
    return memBarrelVol.getVal(t);
}
double molesAirBarrel(double t){ //Moles n = PV/RT 8.314J/(mol*K)
    if (memMolesAirBarrel.inMap(t)) return memMolesAirBarrel.getVal(t);
    if(t<=0){
        return atmPressure*barrelVol(0)/(8.314*initTemp);
    }
    else{
        memMolesAirBarrel.cache(t,molesAirBarrel(t-dt)+34.6*barrelDensity(t-dt)*volumetricFlowRate(t-dt)*dt); //rho*Q = m' 34.6mol/kg
        return memMolesAirBarrel.getVal(t);
    }
}
double molesAirTank (double t){ //Moles
    if (memMolesAirTank.inMap(t)) return memMolesAirTank.getVal(t);
    if(t<=0){
        return initPressure*tankVol/(8.314*initTemp);
    }
    else{
        memMolesAirTank.cache(t, molesAirTank(t-dt)-34.6*tankDensity(t-dt)*volumetricFlowRate(t-dt)*dt);
        return memMolesAirTank.getVal(t);
    }
}
double tankPressure(double t){ //Pascals P = nRT/V
    if (memTankPressure.inMap(t)) return memTankPressure.getVal(t);
    if (t<=0) return initPressure;
    memTankPressure.cache(t, 8.314*molesAirTank(t)*initTemp/tankVol);
    return memTankPressure.getVal(t);
}
double barrelPressure(double t){ //Pascals
    if (memBarrelPressure.inMap(t)) return memBarrelPressure.getVal(t);
    if (t<=0) return atmPressure;
    memBarrelPressure.cache(t, 8.314*molesAirBarrel(t)*barrelTemp(t)/barrelVol(t));
    return memBarrelPressure.getVal(t);
}
double barrelTemp(double t){ //Kelvin Adiabatic Expansion
    if (memBarrelTemp.inMap(t)) return memBarrelTemp.getVal(t);
    memBarrelTemp.cache(t, initTemp*pow(barrelVol(0)/barrelVol(t), airSpecificHeatRatio-1));
    return memBarrelTemp.getVal(t);
}
double barrelDensity(double t){ //Kilograms per meter^3 p/RT 287.05J/(kg*K)
    if (memBarrelDensity.inMap(t)) return memBarrelDensity.getVal(t);
    memBarrelDensity.cache(t,barrelPressure(t-dt)/(287.05*barrelTemp(t-dt)));
    return memBarrelDensity.getVal(t);
}
double tankDensity(double t){ //Kilograms per meter^3
    if (memTankDensity.inMap(t)) return memTankDensity.getVal(t);
    memTankDensity.cache(t,tankPressure(t-dt)/(287.05*initTemp));
    return memTankDensity.getVal(t);
}
double volumetricFlowRate(double t){ //meters^3 per second 0.0000078658 SCFH/(m^3/s) 1.8 R/K 0.000145038 PSI/Pa 0.0000000210100000296 PSI^2/Pa^2
    if (memVolumetricFlowRate.inMap(t)) return memVolumetricFlowRate.getVal(t);
    if (critical(t-dt)){
        memVolumetricFlowRate.cache(t,0.0000078658*816*valveFlowCoeff*0.000145038*tankPressure(t-dt)/sqrt(1.8*barrelTemp(t-dt)));
    }
    else{
        memVolumetricFlowRate.cache(t,0.0000078658*962*valveFlowCoeff*sqrt((0.0000000210100000296*(tankPressure(t-dt)*tankPressure(t-dt)+barrelPressure(t-dt)*barrelPressure(t-dt)))/(1.8*barrelTemp(t-dt))));
    }
    return memVolumetricFlowRate.getVal(t);
}
double airForce(double t){ //Newtons a=F/m
    if (memAirForce.inMap(t)) return memAirForce.getVal(t);
    memAirForce.cache(t, inside(t-dt)*(barrelPressure(t-dt)-atmPressure)*M_PI*barrelRadius*barrelRadius);
    return memAirForce.getVal(t);
}

int main(){
    double t = 0;
    while(t <= 3){
        cout << t << ' ' << velocity(t) << endl;
        t+=0.01;
    }
}
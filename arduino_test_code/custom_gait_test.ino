#include <Servo.h>
#include "Kinematics.h"

const int kneePin = 9;
const int hipPin = 10;
const float frequency = 0.003;

Servo knee;
Servo hip;

const float kneeOffset = 90.0;
const float hipOffset = 45;
const float kneeRange = 60.0;
const float hipRange = 60.0; 

float timer;

void setup() {
  knee.attach(kneePin);
  hip.attach(hipPin);

  knee.write(kneeOffset);
  hip.write(hipOffset);

  timer = millis();
}

void loop() {
  // put your main code here, to run repeatedly:
  delay(50);

  float t = millis() - timer;
  
  hip.write(hipOffset + hipRange/2*cos(frequency*t));
  knee.write(kneeOffset + kneeRange/2*sin(frequency*t));
}

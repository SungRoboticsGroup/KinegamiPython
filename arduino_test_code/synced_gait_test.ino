#include <Servo.h>

const float frequency = 0.0025;

Servo servos[14];

int servoPins[14] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

float servoOffsets[14] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

const float kneeRange = 60.0;
const float hipRange = 60.0; 
const float kneeOffset = 90.0;
const float hipOffset = 0;

float timer;

void moveServos(float t) {
  for (int i = 0; i < 6; i++) {
    int hipIndex = i*2;
    int kneeIndex = i*2 + 1;
    servos[hipIndex].write(servoOffsets[hipIndex] + hipOffset + hipRange/2*cos(frequency*t + i*PI/6));
    servos[kneeIndex].write(servoOffsets[kneeIndex] + kneeOffset + kneeRange/2*sin(frequency*t + i*PI/6));
  }
}
void setup() {
  moveServos(0.0);

  timer = millis();
}

void loop() {
  delay(50);
  
  moveServos(millis() - timer);
}

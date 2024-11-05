#include <Servo.h>

const float frequency = 0.0025;

Servo servos[14];

int servoPins[14] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

float servoOffsets[14] = {0, 90, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
float servoSigns[14] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

float L1 = 8.5; //hipLength; //TODO: measure
float L2 = 8.5;// kneeLength; // TODO: measure
float radius = 4; //TODO: decide

float timer;

float theta1(float x, float y) {
  float beta = atan(y/x);
  float alpha = acos((L1*L1 - L2*L2 + x*x + y*y)/(2*L1*sqrt(x*x + y*y)));
  return (beta - alpha) * 180 / PI;
}

float theta2(float x, float y) {
  float gamma = acos((L1*L1 + L2*L2 - x*x - y*y)/(2*L1*L2));
  return (PI - gamma) * 180 / PI;
}

void moveServos(float t) {
  for (int i = 0; i < 6; i++) {
    int hipIndex = i*2;
    int kneeIndex = i*2 + 1;
    float x = radius*cos(frequency*t + i*PI/6);
    float y = radius*sin(frequency*t + i*PI/6);
    if (y < 0) {
      y = 0;
    }

    Serial.println(x, y);
    y = L1 + L2 - y - 1; //subtract one for tolerance
    
    // x = 0.41289672736166194;
    // y = 11.017077534966651;

    Serial.println(theta1(x,y));
    Serial.println(theta2(x,y));

    float angle1 = theta1(x,y);
    if (x < 0) {
      angle1 += 180;
    }

    servos[hipIndex].write(servoSigns[hipIndex] * (angle1 + servoOffsets[hipIndex]));
    servos[kneeIndex].write(servoSigns[kneeIndex] * (theta2(x,y) + servoOffsets[kneeIndex]));
  }
}
void setup() {
  Serial.begin(9600);
  Serial.println("HEY");

  for (int i = 0; i < 14; i++) {
    servos[i].attach(servoPins[i]);
  }
  
  moveServos(0.0);

  timer = millis();

  delay(2000);

}

void loop() {
  delay(50);
  
  moveServos(millis() - timer);
}

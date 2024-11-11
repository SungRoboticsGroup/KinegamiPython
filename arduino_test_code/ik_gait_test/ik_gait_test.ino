#include <Servo.h>

const float frequency = 0.0025f;

Servo servos[14];

int servoPins[14] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

float servoOffsets[14] = {0, 90, 0, 90, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
float servoSigns[14] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

bool highFreq[6] = {true,true,true,true,true,true};
int counts[6] = {0,0,0,0,0,0};

float L1 = 8.5; //cm
float L2 = 8.5;
float radius = 4; 

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
    float x,y;

    if (highFreq[i]) {
      //semicircle part of motion
      x = radius*cos(frequency*t + i*PI/6.f);
      y = radius*sin(frequency*t + i*PI/6.f);
    } else {
      //dragging along the floor (five times as long)
      x = radius*cos(frequency*t / 5.f + i*PI/30.f - 6*PI/5.f * counts[i]);
      y = radius*sin(frequency*t / 5.f + i*PI/30.f - 6*PI/5.f * counts[i]);
    }

    //switch between phases of motion
    if (highFreq[i] && (y < 0)) {
      highFreq[i] = false;
      counts[i] = (counts[i] + 1) % 6;
      x = radius*cos(frequency*t / 5.f + i*PI/30.f - 6*PI/5.f * counts[i]);
      y = radius*sin(frequency*t / 5.f + i*PI/30.f - 6*PI/5.f * counts[i]);
    } else if (!highFreq[i] && (y > 0)) {
      highFreq[i] = true;
      x = radius*cos(frequency*t + i*PI/6);
      y = radius*sin(frequency*t + i*PI/6);
    }

    if (y < 0) {
      y = 0;
    }

    // Serial.println(x, y);
    y = L1 + L2 - y - 1; //subtract one for tolerance, so it can drag on the floor
    

    // Serial.println(theta1(x,y));
    // Serial.println(theta2(x,y));

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

// ===================================================================
// Auto-generated from: beach_pepe.session
// Servos: 5   Configurations: 4
//
// 270-degree servos: physical range [0, 270]
//   Servo library signal range [0, 180] maps linearly to [0, 270].
//   Joint states are expressed in [-135, 135] (0 = center).
//   Conversion:  signal = (SIGN * state + 135) * 180.0 / 270.0
// ===================================================================

#include <Servo.h>

#define NUM_SERVOS  5
#define NUM_CONFIGS 4
#define NUM_SEGMENTS (NUM_CONFIGS - 1)

// -- Pin assignments (edit these to match your wiring) ----------------
const int SERVO_PINS[NUM_SERVOS] = {
    3,  // servo 0: joint 0 (CoaxialRDS3225)
    4,  // servo 1: joint 4 (TransverseRDS3225)
    5,  // servo 2: joint 5 (TransverseRDS3225)
    6,  // servo 3: joint 7 (TransverseRDS3225)
    7  // servo 4: joint 9 (TransverseRDS3225)
};

// -- Direction signs (+1 normal, -1 reversed) -----------------------
const int SERVO_SIGN[NUM_SERVOS] = {
   1,  // servo 0
   1,  // servo 1
   1,  // servo 2
   1,  // servo 3
   1  // servo 4
};

// -- Configuration sequence (degrees, range [-135, 135]) -------------
// Each row is one configuration; columns are servos 0 .. N-1.
// Edit these values to adjust the dance.
const float CONFIGS[NUM_CONFIGS][NUM_SERVOS] = {
  {    -0.00,     0.00,     0.00,     0.00,     0.00 },  // config 0
  {   -30.00,    35.00,    35.00,     0.00,     0.00 },  // config 1
  {    30.00,     0.00,     0.00,    30.00,    30.00 },  // config 2
  {     0.00,     0.00,     0.00,     0.00,     0.00 }  // config 3
};

// -- Segment durations (milliseconds) --------------------------------
// Duration between config[i] and config[i+1].
const unsigned long SEGMENT_MS[NUM_SEGMENTS] = {
    1000, // segment 0 -> 1
    1000, // segment 1 -> 2
    1000 // segment 2 -> 3
};

// -- Servo objects ----------------------------------------------------
Servo servos[NUM_SERVOS];

// -- Convert joint state [-135, 135] -> servo signal [0, 180] --------
int stateToSignal(int servoIdx, float stateDeg) {
  // Apply direction sign
  float effective = SERVO_SIGN[servoIdx] * stateDeg;
  // Shift from [-135,135] to [0,270], then scale to [0,180]
  float signal = (effective + 135.0) * 180.0 / 270.0;
  // Clamp to valid servo range
  if (signal < 0.0)   signal = 0.0;
  if (signal > 180.0) signal = 180.0;
  return (int)(signal + 0.5);  // round to nearest int
}

void setup() {
  Serial.begin(9600);
  for (int i = 0; i < NUM_SERVOS; i++) {
    servos[i].attach(SERVO_PINS[i]);
  }
  // Move to initial configuration
  for (int i = 0; i < NUM_SERVOS; i++) {
    servos[i].write(stateToSignal(i, CONFIGS[0][i]));
  }
  delay(500);  // let servos reach start position
}

void loop() {
  for (int seg = 0; seg < NUM_SEGMENTS; seg++) {
    unsigned long segStart = millis();
    unsigned long duration = SEGMENT_MS[seg];

    while (true) {
      unsigned long elapsed = millis() - segStart;
      if (elapsed >= duration) break;

      float t = (float)elapsed / (float)duration;  // 0.0 -> 1.0

      for (int i = 0; i < NUM_SERVOS; i++) {
        float val = CONFIGS[seg][i] * (1.0 - t) + CONFIGS[seg + 1][i] * t;
        servos[i].write(stateToSignal(i, val));
      }

      delay(20);  // ~50 Hz update rate
    }

    // Snap to exact end-of-segment position
    for (int i = 0; i < NUM_SERVOS; i++) {
      servos[i].write(stateToSignal(i, CONFIGS[seg + 1][i]));
    }
  }
  // Loop back to start (sequence repeats)
}

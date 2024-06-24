/*!
 * @file  colorview.ino
 * @brief Gets the ambient light color
 * @copyright   Copyright (c) 2010 DFRobot Co.Ltd (http://www.dfrobot.com)
 * @license     The MIT License (MIT)
 * @author      PengKaixing(kaixing.peng@dfrobot.com)
 * @version     V1.0.0
 * @date        2022-03-16
 * @url         https://github.com/DFRobot/DFRobot_TCS34725
 */
#include "DFRobot_TCS34725.h"

//DFRobot_TCS34725 tcs = DFRobot_TCS34725(&Wire, TCS34725_ADDRESS,TCS34725_INTEGRATIONTIME_50MS, TCS34725_GAIN_1X);
DFRobot_TCS34725 tcs = DFRobot_TCS34725(&Wire, TCS34725_ADDRESS,TCS34725_INTEGRATIONTIME_154MS, TCS34725_GAIN_1X);

float R,G,B,C;  // running average of RGBW values
float const f = 0.2;
uint16_t clear, red, green, blue;
int averages = 7;

void setup() 
{
  Serial.begin(115200);
  delay(2000);

  while(!tcs.begin())  // fails if ID not TCS34725, but TCS34723/TCS34727 otherwise works
  {
    Serial.println("No TCS34725 found ... check your connections");
    delay(1000);
  }

  tcs.getRGBC(&red, &green, &blue, &clear); // get initial values
  R=red;
  G=green;
  B=blue;
  C=clear;
}

void loop() {
  R=0; G=0; B=0; C=0;
  for (int i=0;i<averages;i++) {
    tcs.getRGBC(&red, &green, &blue, &clear);
    R += red;
    G += green;
    B += blue;
    C += clear;
  }
  R = R / (float) averages;
  G = G / (float) averages;
  B = B / (float) averages;
  C = C / (float) averages;

  red = (uint) R;
  green = (uint) G;
  blue = (uint) B;
  clear = (uint) C;

  uint16_t cct = tcs.calculateColortemperature(red, green,blue);
  //uint16_t lux = tcs.calculateLux(red,green,blue);
  
  Serial.print(R,3); Serial.print(", "); 
  Serial.print(G,3); Serial.print(", ");
  Serial.print(B,3); Serial.print(", ");
  Serial.print(C,3); Serial.print(", "); 
  Serial.print(cct);
  // Serial.print(", "); 
  //Serial.print(lux);

  Serial.println();
  // delay(1000);
}


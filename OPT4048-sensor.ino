/*
Example 3 - Color Warmth

This example shows how to retrieve the Corelated Color Temperature (CCT) from the OPT4048.
If you need more information on CCT, check out this Wikipedia article:
https://en.wikipedia.org/wiki/Correlated_color_temperature

Written by Elias Santistevan @ SparkFun Electronics, July 2023

Products:
    Qwiic 1x1: https://www.sparkfun.com/products/22638
    Qwiic Mini: https://www.sparkfun.com/products/22639

Repository:
    https://github.com/sparkfun/SparkFun_OPT4048_Arduino_Library

SparkFun code, firmware, and software is released under the MIT 
License	(http://opensource.org/licenses/MIT).
*/

#include "SparkFun_OPT4048.h"
#include <Wire.h>

SparkFun_OPT4048 myColor;

void setup()
{
    Serial.begin(115200);
    Serial.println("OPT4048 Example 3 Basic Color Warmth.");

    Wire.begin();
    delay(2000);

    if (!myColor.begin()) {
        Serial.println("OPT4048 not detected- check wiring or that your I2C address is correct!");
        while (1) ;
    }

    myColor.setBasicSetup();
    myColor.setRange(RANGE_AUTO);  // or RANGE_144LUX

    // Serial.println("Ready to go!");
}


void loop()
{

    // bool bright = myColor.getTooBrightFlag();
    sfe_color_t c = myColor.getAllADC();
    Serial.print(c.red);
    Serial.print(",");
    Serial.print(c.green);
    Serial.print(",");
    Serial.print(c.blue);
    Serial.print(",");
    Serial.print(c.white);
    Serial.print(",");

    Serial.print(myColor.getLux(),1);
    Serial.print(",");
    Serial.print(myColor.getCCT(),1);
    Serial.print(",");
    Serial.print(myColor.getCIEx(),3);
    Serial.print(",");
    Serial.println(myColor.getCIEy(),3);


    // Delay time is set to the conversion time * number of channels
    // You need three channels for color sensing @ 200ms conversion time = 600ms.
    delay(600);
}

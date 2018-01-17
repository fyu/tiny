/* Source file for the R3 light class */



/* Include files */

#include "R3Graphics.h"



/* Class type definitions */

RN_CLASS_TYPE_DEFINITIONS(R3Light);



/* Public variables */

RNScalar R3ambient_light_intensity = 0.2;
RNRgb R3ambient_light_color(1.0, 1.0, 1.0);



/* Public functions */

int 
R3InitLight()
{
    /* Return success */
    return TRUE;
}



void 
R3StopLight()
{
}



R3Light::
R3Light(void)
  : id(0)
{
}



R3Light::
R3Light(const R3Light& light)
    : active(light.active),
      intensity(light.intensity),
      color(light.color),
      id(-1)
{
}



R3Light::
R3Light(const RNRgb& color,
	RNScalar intensity,
	RNBoolean active)
    : active(active),
      intensity(intensity),
      color(color),
      id(-1)
{
}



R3Light::
~R3Light(void)
{
}



void R3Light::
SetActive(RNBoolean active)
{
    // Set active
    this->active = active;
}



void R3Light::
SetIntensity(RNScalar intensity)
{
    // Set intensity
    this->intensity = intensity;
}



void R3Light::
SetColor(const RNRgb& color)
{
    // Set color
    this->color = color;
}




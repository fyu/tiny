/* Include file for the R3 light class */



/* Initialization functions */

int R3InitLight();
void R3StopLight();



/* Class definition */

class R3Light {
    public:
        // Constructor functions
	R3Light(void);
	R3Light(const R3Light& light);
        R3Light(const RNRgb& color, RNScalar intensity = 1.0, RNBoolean active = TRUE);
        virtual ~R3Light(void);

	// Property functions/operations
	const RNBoolean IsActive(void) const;
  	const RNScalar Intensity(void) const;
  	const RNRgb& Color(void) const;

	// Manipulation functions/operations
	virtual void SetActive(RNBoolean active);
  	virtual void SetIntensity(RNScalar intensity);
  	virtual void SetColor(const RNRgb& color);

	// Reflection evaluation functions
	virtual RNRgb Reflection(const R3Brdf& brdf, const R3Point& eye, 
	    const R3Point& point, const R3Vector& normal) const = 0;
	virtual RNRgb DiffuseReflection(const R3Brdf& brdf, 
	    const R3Point& point, const R3Vector& normal) const = 0;
	virtual RNRgb SpecularReflection(const R3Brdf& brdf, const R3Point& eye, 
	    const R3Point& point, const R3Vector& normal) const = 0;

	// Draw functions/operations
        virtual void Draw(int i) const = 0;

	// Class type definitions
	RN_CLASS_TYPE_DECLARATIONS(R3Light);

    private:
	RNBoolean active;
	RNScalar intensity;
	RNRgb color;
        int id;
};



/* Public variables */

extern RNScalar R3ambient_light_intensity;
extern RNRgb R3ambient_light_color;



/* Inline functions */

inline const RNBoolean R3Light::
IsActive(void) const
{
    // Return status
    return active;
}



inline const RNScalar R3Light::
Intensity(void) const
{
    // Return intensity 
    return intensity;
}



inline const RNRgb& R3Light::
Color(void) const
{
    // Return color 
    return color;
}




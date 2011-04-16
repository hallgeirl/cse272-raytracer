#ifndef CSE272_EMISSIVE_H_INCLUDED
#define CSE272_EMISSIVE_H_INCLUDED

#include "Miro.h"
#include "Vector3.h"

class Emissive : public Material
{
public:

	Emissive()
		:m_power(0)
	{}

    virtual Vector3 shade(const Ray& ray, const HitInfo& hit, const Scene& scene) const
	{
		return m_power;
	}

	virtual void setPower(const float power) {m_power = power;}

	bool         isReflective() const { return false; }
	bool         isRefractive() const { return false; }
	virtual bool isDiffuse() const    { return false; }

	Vector3         getReflection() const       {return Vector3(0);}
	Vector3         getRefraction() const       {return Vector3(0);}
	virtual Vector3 getDiffuse() const          {return Vector3(1.f);}

private:
	float m_power;
};

#endif //CSE272_EMISSIVE_H_INCLUDED

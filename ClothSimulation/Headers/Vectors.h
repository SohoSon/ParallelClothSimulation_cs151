#pragma once

#include <math.h>
#include <iostream>

struct Vec2
{
    double x;
    double y;
    
    Vec2(void)
    {
        x = 0.0;
        y = 0.0;
    }
    Vec2(double x0, double y0)
    {
        x = x0;
        y = y0;
    }
    ~Vec2() {}
    
    Vec2 operator+(Vec2 v)
    {
        return Vec2(x+v.x, y+v.y);
    }
    Vec2 operator-(Vec2 v)
    {
        return Vec2(x-v.x, y-v.y);
    }

    void operator+=(Vec2 v)
    {
        x += v.x;
        y += v.y;
    }
    void operator-=(Vec2 v)
    {
        x -= v.x;
        y -= v.y;
    }
    
};

struct Vec3
{
    double x;
    double y;
    double z;
    
    Vec3(void)
    {
        x = 0.0;
        y = 0.0;
        z = 0.0;
    }
    Vec3(double x0, double y0, double z0)
    {
		x = x0;
		y = y0;
		z = z0;
	}
	~Vec3(){}
    
    static Vec3 cross(Vec3 v1, Vec3 v2)
    {
        return Vec3(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x);
    }
    static double dot(Vec3 v1, Vec3 v2)
    {
        return (v1.x*v2.x + v1.y*v2.y + v1.z*v2.z);
    }
    static double dist(Vec3 v1, Vec3 v2)
    {
        return sqrt(pow(v1.x - v2.x, 2) + pow(v1.y - v2.y, 2) + pow(v1.z - v2.z, 2));
    }
    
    Vec3 minus()
    {
        return Vec3(-x, -y, -z);
    }
    
	Vec3 operator+(Vec3 v)
	{
		return Vec3(x+v.x, y+v.y, z+v.z);
	}
	Vec3 operator-(Vec3 v)
	{
		return Vec3(x-v.x, y-v.y, z-v.z);
	}

    // TODO
    Vec3 operator--() const
    {
        return Vec3(-x, -y, -z);
    }

	void operator+=(Vec3 v)
	{
		x += v.x;
		y += v.y;
		z += v.z;
	}
	void operator-=(Vec3 v)
	{
		x -= v.x;
		y -= v.y;
		z -= v.z;
	}

    // new normalized, because diff is a temporary object, cannot call normalize
    Vec3 normalized() const {
        float len = length();
        if (len > 1e-8) {
            return Vec3(x / len, y / len, z / len);
        }
        return Vec3(0, 0, 0);
    }
    
	Vec3 operator*(double n)
	{
		return Vec3(x*n, y*n, z*n);
	}
	Vec3 operator/(double n)
	{
		return Vec3(x/n, y/n, z/n);
	}
    
    bool operator ==(Vec3 &v)
    {
        return x == v.x && y == v.y&&z == v.z;
    }
    bool operator!=(Vec3 &v)
    {
        return x != v.x || y != v.y || z != v.z;
    }

    double length() const
    {
        return sqrt(x*x + y*y + z*z);
    }
    
    void normalize()
    {
        double w = length();
        if (w < 0.00001) return;
        
        x /= w;
        y /= w;
        z /= w;
    }
	
	void setZeroVec(){
        x = 0.0;
        y = 0.0;
        z = 0.0;
	}

    void setInitVec(){
        x = 0.0;
        y = -9.8 / 100;
        z = 0.0;
	}
    //overloading << operator
    friend std::ostream& operator<<(std::ostream& os, const Vec3& vec) {
        os << "(" << vec.x << ", " << vec.y << ", " << vec.z << ")";
        return os;
    }
};



    // void normalize()
    // {
    //     double w = length();
    //     if (w < 0.00001) return;
        
    //     x /= w;
    //     y /= w;
    //     z /= w;
    // }
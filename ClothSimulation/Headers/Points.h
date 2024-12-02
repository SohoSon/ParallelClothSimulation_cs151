#pragma once

#include "Vectors.h"

struct Vertex
{
public:
    Vec3 position;
    Vec3 normal;
    
    Vertex() {}
    Vertex(Vec3 pos)
    {
        position = pos;
    }
    ~Vertex() {}
};

class alignas(64) Node
{
public:
    double  mass;           // In this project it will always be 1
    bool    isFixed;        // Use to pin the cloth
    Vec2    texCoord;       // Texture coord
    Vec3    normal;         // For smoothly shading
	Vec3	position;
    Vec3    velocity;
    Vec3    force;
	Vec3	acceleration;
    int     SpringCount;
    // std::atomic<int>     TempSpringCount;
    int     TempSpringCount;


public:
    Node(void) {
        mass = 1.0;
        isFixed = false;
        velocity.setZeroVec();
        force.setInitVec();
        acceleration.setZeroVec();
        SpringCount = 0;
        TempSpringCount = 0;
    }
	Node(Vec3 pos)
    {
        mass = 1.0;
        isFixed = false;
        position = pos;
        velocity.setZeroVec();
        force.setInitVec();
        acceleration.setZeroVec();
        SpringCount = 0;
        TempSpringCount = 0;
    }
	~Node(void) {}

	void addForce(Vec3 f)
	{
        force += f;
	}

	void integrate(double timeStep) // Only non-fixed nodes take integration
	{
		if (!isFixed) // Verlet integration
		{
            acceleration = force/mass;
            velocity += acceleration*timeStep;
            position += velocity*timeStep;
        }
        force.setInitVec();
	}
};

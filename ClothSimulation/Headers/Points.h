#pragma once

#include "Vectors.h"
//#include "Spring.h"
#include <vector>

class Spring;

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
    //for node-based parallel
    std::vector<Spring*> connectedSprings; // 新增：存储与节点相连的弹簧

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
        #pragma omp atomic
        force.x += f.x;
        #pragma omp atomic
        force.y += f.y;
        #pragma omp atomic
        force.z += f.z;
	}

	void integrate(double timeStep) // Only non-fixed nodes take integration
	{
		if (!isFixed) // Verlet integration
		{
            acceleration = force/mass; // 计算加速度
            velocity += acceleration*timeStep; // 更新速度
            position += velocity*timeStep; // 更新位置
        }
        force.setInitVec();
	}
    //for node-based parallel
    void addSpring(Spring* spring) {
        //std::cout << "Adding spring to node: " << position << std::endl;
        connectedSprings.push_back(spring);
    }
};

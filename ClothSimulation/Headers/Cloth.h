#pragma once

#include <vector>
#include <omp.h>
#include <chrono>
#include "Spring.h"
#include "Rigid.h"
#include <iomanip> 
class Cloth
{
public:
    const int nodesDensity = 20;
    const int iterationFreq = 100;
    const double structuralCoef = 4000.0;
    const double shearCoef = 50.0;
    const double bendingCoef = 400.0;
    
    enum DrawModeEnum{
        DRAW_NODES,
        DRAW_LINES,
        DRAW_FACES
    };
    DrawModeEnum drawMode = DRAW_FACES;
    
    Vec3 clothPos;
    
    int width, height;
    int nodesPerRow, nodesPerCol;
    
    std::vector<Node*> nodes;
	std::vector<Spring*> springs;
	std::vector<Node*> faces;
    
    Vec2 pin1;
    Vec2 pin2;
    
    // add parallel control variables
    bool enable_parallel = true;  // control whether to enable parallel
    int num_threads = 4;          // control thread number  

	Cloth(Vec3 pos, Vec2 size)
	{
        clothPos = pos;
        width = size.x;
        height = size.y;
        init();

       
	}
	~Cloth()
	{ 
		for (int i = 0; i < nodes.size(); i++) { delete nodes[i]; }
		for (int i = 0; i < springs.size(); i++) { delete springs[i]; }
		nodes.clear();
		springs.clear();
		faces.clear();
        
	}
 
public:
    Node* getNode(int x, int y) { return nodes[y*nodesPerRow+x]; }
    Vec3 computeFaceNormal(Node* n1, Node* n2, Node* n3)
    {
        return Vec3::cross(n2->position - n1->position, n3->position - n1->position);
    }
    
    void pin(Vec2 index, Vec3 offset) // Pin cloth's (x, y) node with offset
    {
        if (!(index.x < 0 || index.x >= nodesPerRow || index.y < 0 || index.y >= nodesPerCol)) {
            getNode(index.x, index.y)->position += offset;
            getNode(index.x, index.y)->isFixed = true;
        }
    }
    void unPin(Vec2 index) // Unpin cloth's (x, y) node
    {
        if (!(index.x < 0 || index.x >= nodesPerRow || index.y < 0 || index.y >= nodesPerCol)) {
            getNode(index.x, index.y)->isFixed = false;
        }
    }
    
	void init()
	{
        nodesPerRow = width * nodesDensity;
        nodesPerCol = height * nodesDensity;
        
        pin1 = Vec2(0, 0);
        pin2 = Vec2(nodesPerRow-1, 0);
        
        /** Add nodes **/
        printf("Init cloth with %d nodes\n", nodesPerRow*nodesPerCol);
        for (int i = 0; i < nodesPerRow; i ++) {
            for (int j = 0; j < nodesPerCol; j ++) {
                /** Create node by position **/
                Node* node = new Node(Vec3((double)j/nodesDensity, -((double)i/nodesDensity), 0));
                /** Set texture coordinates **/
                node->texCoord.x = (double)j/(nodesPerRow-1);
                node->texCoord.y = (double)i/(1-nodesPerCol);
                /** Add node to cloth **/
                nodes.push_back(node);
                
                printf("\t[%d, %d] (%f, %f, %f) - (%f, %f)\n", i, j, node->position.x, node->position.y, node->position.z, node->texCoord.x, node->texCoord.y);
            }
            std::cout << std::endl;
        }
        
        /** Add springs **/
        for (int i = 0; i < nodesPerRow; i ++) {
            for (int j = 0; j < nodesPerCol; j ++) {
                /** Structural **/
                if (i < nodesPerRow-1) springs.push_back(new Spring(getNode(i, j), getNode(i+1, j), structuralCoef));
                if (j < nodesPerCol-1) springs.push_back(new Spring(getNode(i, j), getNode(i, j+1), structuralCoef));
                /** Shear **/
                if (i < nodesPerRow-1 && j < nodesPerCol-1) {
                    springs.push_back(new Spring(getNode(i, j), getNode(i+1, j+1), shearCoef));
                    springs.push_back(new Spring(getNode(i+1, j), getNode(i, j+1), shearCoef));
                }
                /** Bending **/
                if (i < nodesPerRow-2) springs.push_back(new Spring(getNode(i, j), getNode(i+2, j), bendingCoef));
                if (j < nodesPerCol-2) springs.push_back(new Spring(getNode(i, j), getNode(i, j+2), bendingCoef));
            }
        }
        
        pin(pin1, Vec3(1.0, 0.0, 0.0));
        pin(pin2, Vec3(-1.0, 0.0, 0.0));
        
		/** Triangle faces **/
        for (int i = 0; i < nodesPerRow-1; i ++) {
            for (int j = 0; j < nodesPerCol-1; j ++) {
                // Left upper triangle
                faces.push_back(getNode(i+1, j));
                faces.push_back(getNode(i, j));
                faces.push_back(getNode(i, j+1));
                // Right bottom triangle
                faces.push_back(getNode(i+1, j+1));
                faces.push_back(getNode(i+1, j));
                faces.push_back(getNode(i, j+1));
            }
        }
	}
	
	void computeNormal()
	{
        /** Reset nodes' normal **/
        Vec3 normal(0.0, 0.0, 0.0);
        for (int i = 0; i < nodes.size(); i ++) {
            nodes[i]->normal = normal;
        }
        /** Compute normal of each face **/
        for (int i = 0; i < faces.size()/3; i ++) { // 3 nodes in each face
            Node* n1 = faces[3*i+0];
            Node* n2 = faces[3*i+1];
            Node* n3 = faces[3*i+2];
            
            // Face normal
            normal = computeFaceNormal(n1, n2, n3);
            // Add all face normal
            n1->normal += normal;
            n2->normal += normal;
            n3->normal += normal;
        }
        
        for (int i = 0; i < nodes.size(); i ++) {
            nodes[i]->normal.normalize();
        }
	}
	
	void addForce(Vec3 f)
	{		 
		for (int i = 0; i < nodes.size(); i++)
		{
			nodes[i]->addForce(f);
		}
	}

    void simulation(double timeStep, Vec3 gravity, Ground* ground, Ball* ball) {
        // init
        double total_start = omp_get_wtime();
        double simulation_time = 0.0;
        // for(int i = 0; i < nodes.size(); i++) {
        //     nodes[i]->TempSpringCount = nodes[i]->SpringCount;
        // }
        // spring force calculation
        {
            double start = omp_get_wtime();
            if (enable_parallel) {
                #pragma omp parallel for num_threads(num_threads) schedule(static, 256)
                for (int i = 0; i < springs.size(); i++) {
                    springs[i]->applyInternalForce(timeStep);
                 
                    // omp_set_lock(&node_locks[springs[i]->node1]);
                    // omp_set_lock(&node_locks[springs[i]->node2]);
                    springs[i]->node1->TempSpringCount--;
                    springs[i]->node2->TempSpringCount--;
                //     omp_unset_lock(&node_locks[springs[i]->node1]);
                // omp_unset_lock(&node_locks[springs[i]->node2]);
                    if (springs[i]->node1->TempSpringCount <= 0) {
                        springs[i]->node1->TempSpringCount = springs[i]->node1->SpringCount;
                        springs[i]->node1->integrate(timeStep);
                        /** Ground collision **/
                        if (getWorldPos(springs[i]->node1).y < ground->position.y) {
                            springs[i]->node1->position.y = ground->position.y - clothPos.y + 0.01;
                            springs[i]->node1->velocity = springs[i]->node1->velocity * ground->friction;
                        }
                        
                        /** Ball collision **/
                        Vec3 distVec = getWorldPos(springs[i]->node1) - ball->center;
                        double distLen = distVec.length();
                        double safeDist = ball->radius*1.05;
                        if (distLen < safeDist) {
                            distVec.normalize();
                            setWorldPos(springs[i]->node1, distVec*safeDist+ball->center);
                            springs[i]->node1->velocity = springs[i]->node1->velocity * ball->friction;
                        }
                    }
                    if (springs[i]->node2->TempSpringCount <= 0) {
                        springs[i]->node2->integrate(timeStep);
                        springs[i]->node2->TempSpringCount = springs[i]->node2->SpringCount;
                        /** Ground collision **/
                        if (getWorldPos(springs[i]->node2).y < ground->position.y) {
                            springs[i]->node2->position.y = ground->position.y - clothPos.y + 0.01;
                            springs[i]->node2->velocity = springs[i]->node2->velocity * ground->friction;
                        }
                        
                        /** Ball collision **/
                        Vec3 distVec = getWorldPos(springs[i]->node2) - ball->center;
                        double distLen = distVec.length();
                        double safeDist = ball->radius*1.05;
                        if (distLen < safeDist) {
                            distVec.normalize();
                            setWorldPos(springs[i]->node2, distVec*safeDist+ball->center);
                            springs[i]->node2->velocity = springs[i]->node2->velocity * ball->friction;
                        }
                    }
                }
            } else {
                for (int i = 0; i < springs.size(); i++) {
                    springs[i]->applyInternalForce(timeStep);
                    
                    springs[i]->node1->TempSpringCount--;
                    springs[i]->node2->TempSpringCount--;
                 
                    if (springs[i]->node1->TempSpringCount == 0) {
                        springs[i]->node1->integrate(timeStep);
                        springs[i]->node1->TempSpringCount = springs[i]->node1->SpringCount;
                        /** Ground collision **/
                        if (getWorldPos(springs[i]->node1).y < ground->position.y) {
                            springs[i]->node1->position.y = ground->position.y - clothPos.y + 0.01;
                            springs[i]->node1->velocity = springs[i]->node1->velocity * ground->friction;
                        }
                        
                        /** Ball collision **/
                        Vec3 distVec = getWorldPos(springs[i]->node1) - ball->center;
                        double distLen = distVec.length();
                        double safeDist = ball->radius*1.05;
                        if (distLen < safeDist) {
                            distVec.normalize();
                            setWorldPos(springs[i]->node1, distVec*safeDist+ball->center);
                            springs[i]->node1->velocity = springs[i]->node1->velocity * ball->friction;
                        }
                    }
                    if (springs[i]->node2->TempSpringCount == 0) {
                        springs[i]->node2->integrate(timeStep);
                        springs[i]->node2->TempSpringCount = springs[i]->node2->SpringCount;
                        /** Ground collision **/
                        if (getWorldPos(springs[i]->node2).y < ground->position.y) {
                            springs[i]->node2->position.y = ground->position.y - clothPos.y + 0.01;
                            springs[i]->node2->velocity = springs[i]->node2->velocity * ground->friction;
                        }
                        
                        /** Ball collision **/
                        Vec3 distVec = getWorldPos(springs[i]->node2) - ball->center;
                        double distLen = distVec.length();
                        double safeDist = ball->radius*1.05;
                        if (distLen < safeDist) {
                            distVec.normalize();
                            setWorldPos(springs[i]->node2, distVec*safeDist+ball->center);
                            springs[i]->node2->velocity = springs[i]->node2->velocity * ball->friction;
                        }
                    }
                }
            }
            simulation_time = omp_get_wtime() - start;
        }

        double total_time = omp_get_wtime() - total_start;

        // stats
        #pragma omp critical
        {
            static int frame_count = 0;
            static double last_time = 0.0;
            static double fps = 0.0;
            
            frame_count++;
            double current_time = omp_get_wtime();
            
            // update FPS per second
            if (current_time - last_time >= 1.0) {
                fps = frame_count / (current_time - last_time);
                frame_count = 0;
                last_time = current_time;
                
                // only output performance data when FPS is updated
                std::cout << "performance data:\n"
                        << "FPS: " << std::fixed << std::setprecision(2) << fps << "\n"
                        << "total time: " << total_time * 1000 << " ms\n"
                        << "spring calc time: " << simulation_time * 1000 << " ms\n"
                        //<< "Current threads: " << omp_get_max_threads() << "\n"
                        << "Current threads: " << num_threads << "\n"
                        << "total nodes: " << nodes.size() << "\n"
                        << "total springs: " << springs.size() << "\n"
                        << "------------------------\n";
            }
        }
    }

	// void integrate(double airFriction, double timeStep)
	// {
    //     for (int i = 0; i < nodes.size(); i++) {
    //         nodes[i]->integrate(timeStep);
    //     }
        
	// }
	
    Vec3 getWorldPos(Node* n) { return clothPos + n->position; }
    void setWorldPos(Node* n, Vec3 pos) { n->position = pos - clothPos; }
    
	void collisionResponse(Ground* ground, Ball* ball)
	{
        for (int i = 0; i < nodes.size(); i++)
        {
            /** Ground collision **/
            if (getWorldPos(nodes[i]).y < ground->position.y) {
                nodes[i]->position.y = ground->position.y - clothPos.y + 0.01;
                nodes[i]->velocity = nodes[i]->velocity * ground->friction;
            }
            
            /** Ball collision **/
            Vec3 distVec = getWorldPos(nodes[i]) - ball->center;
            double distLen = distVec.length();
            double safeDist = ball->radius*1.05;
            if (distLen < safeDist) {
                distVec.normalize();
                setWorldPos(nodes[i], distVec*safeDist+ball->center);
                nodes[i]->velocity = nodes[i]->velocity*ball->friction;
            }
        }
	}
};

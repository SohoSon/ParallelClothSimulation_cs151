#pragma once

#include <vector>
#include <omp.h>
#include <chrono>
#include "Spring.h"
#include "Rigid.h"
#include <iomanip> 
#include <unordered_map>
#include <cmath>

class Cloth
{
private:

    int hashFunction(int x, int y, int z) const {
        return x * 73856093 ^ y * 19349663 ^ z * 83492791;
    }

public:
    const int nodesDensity = 10;
    const int iterationFreq = 100;
    const double structuralCoef = 4000.0;
    const double shearCoef = 50.0;
    const double bendingCoef = 400.0;
    const double collisionDistance = 0.02; 
    bool handleCollisions = false;
    bool groundRender_visible = false;
    bool ballRender_visible = false;
    enum DrawModeEnum{
        DRAW_NODES,
        DRAW_LINES,
        DRAW_FACES
    };
    DrawModeEnum drawMode = DRAW_LINES;
    
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
    std::vector<Node*>& getNodes() { return nodes; }
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
        std::cout << "Initializing cloth..." << std::endl;
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
                if (i < nodesPerRow-1) {
                    Spring* spring = new Spring(getNode(i, j), getNode(i+1, j), structuralCoef);
                    springs.push_back(spring);
                }
                if (j < nodesPerCol-1) {
                    Spring* spring = new Spring(getNode(i, j), getNode(i, j+1), structuralCoef);
                    springs.push_back(spring);
                }
                /** Shear **/
                if (i < nodesPerRow-1 && j < nodesPerCol-1) {
                    Spring* spring1 = new Spring(getNode(i, j), getNode(i+1, j+1), shearCoef);
                    Spring* spring2 = new Spring(getNode(i+1, j), getNode(i, j+1), shearCoef);
                    springs.push_back(spring1);
                    springs.push_back(spring2);
                }
                /** Bending **/
                if (i < nodesPerRow-2) {
                    Spring* spring = new Spring(getNode(i, j), getNode(i+2, j), bendingCoef);
                    springs.push_back(spring);
                }
                if (j < nodesPerCol-2) {
                    Spring* spring = new Spring(getNode(i, j), getNode(i, j+2), bendingCoef);
                    springs.push_back(spring);
                }
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
        std::cout << "Cloth initialization complete." << std::endl;
	}
	void resetCloth() {
        std::cout << "Resetting cloth..." << std::endl;
        for (int i = 0; i < nodes.size(); i++) { delete nodes[i]; }
        for (int i = 0; i < springs.size(); i++) { delete springs[i]; }
        nodes.clear();
        springs.clear();
        faces.clear();
        init();
        Vec3 initForce(10.0, 40.0, 20.0);
        addForce(initForce);
        std::cout << "Cloth reset complete." << std::endl;
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

void node_based_simulation(double timeStep, Vec3 gravity, Ground* ground, Ball* ball, float cellSize) {
    double total_start = omp_get_wtime();
    double simulation_time = 0.0;
    //std::unordered_map<int, std::vector<Node*>> spatialHash;
    
        double thickness = 0.01;
        double thickness2 = thickness * thickness;
        double maxDist = 1;
        double maxDist2 = maxDist * maxDist;
    //float cellSize = 0.1f; 
    std::vector<std::vector<int>> adjacencyList(nodes.size());
    // 
    if (handleCollisions) {
        std::unordered_map<int, std::vector<int>> spatialHash;
        spatialHash.reserve(nodes.size());
        spatialHash.max_load_factor(0.7);

        for (size_t i = 0; i < nodes.size(); i++) {
            Node* n = nodes[i];
            if (n->isFixed) continue;
            
            int gridX = (int)std::floor(n->position.x / cellSize);
            int gridY = (int)std::floor(n->position.y / cellSize);
            int gridZ = (int)std::floor(n->position.z / cellSize);
            int hashKey = hashFunction(gridX, gridY, gridZ);
            spatialHash[hashKey].push_back((int)i);
        }


    


    // define a function to get the nodes in the bucket corresponding to the given coordinates
    auto getBucketNodes = [&](int x, int y, int z) -> const std::vector<int>& {
        int key = hashFunction(x, y, z);
        auto it = spatialHash.find(key);
        if (it != spatialHash.end()) {
            return it->second;
        }
        static const std::vector<int> empty;
        return empty;
    };

    // for each node, perform the "query" operation
    #pragma omp parallel for num_threads(num_threads)
    for (size_t i = 0; i < nodes.size(); i++) {
        Node* n0 = nodes[i];
        if (n0->isFixed) continue;

        int gridX = (int)std::floor(n0->position.x / cellSize);
        int gridY = (int)std::floor(n0->position.y / cellSize);
        int gridZ = (int)std::floor(n0->position.z / cellSize);

        // search range, including the current grid and adjacent grids
        // you can increase or decrease the search range according to your needs
        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy++) {
                for (int dz = -1; dz <= 1; dz++) {
                    const std::vector<int>& bucketNodes = getBucketNodes(gridX+dx, gridY+dy, gridZ+dz);
                    for (int idx : bucketNodes) {
                        if ((size_t)idx == i) continue; // 
                        Node* n1 = nodes[idx];
                        if (n1->isFixed) continue;
                        
                        double dx = n1->position.x - n0->position.x;
                        double dy = n1->position.y - n0->position.y;
                        double dz = n1->position.z - n0->position.z;
                        double dist2 = dx*dx + dy*dy + dz*dz;

                        if (dist2 < maxDist2) {
                            // n1 is a neighbor of n0
                            adjacencyList[i].push_back(idx);
                        }
                    }
                }
            }
        }
    }
 // output adjacencyList information
    // for (size_t i = 0; i < adjacencyList.size(); i++) {
    //     std::cout << "Node " << i << " has neighbors: ";
    //     for (int neighbor : adjacencyList[i]) {
    //         std::cout << neighbor << " ";
    //     }
    //     std::cout << std::endl;
    // }
    }

    
    if (enable_parallel) {
        #pragma omp parallel for num_threads(num_threads)
        for (int i = 0; i < nodes.size(); i++) {
            Node* node = nodes[i];
            
            // Reset forces
            //node->force = Vec3(0.0, -9.8/100, 0.0);

            // Calculate spring forces
            //std::cout << "connectedSprings: " << node->connectedSprings.size() << std::endl;
            for (Spring* spring : node->connectedSprings) {
                Vec3 force = spring->computeForce(node);
                //std::cout << "Spring " << i << ": Force = " << force << std::endl;
                node->addForce(force);
            }

        }
        #pragma omp parallel for num_threads(num_threads)
        for (int i = 0; i < nodes.size(); i++) {
            Node* node = nodes[i];
            node->integrate(timeStep);

            // Ground collision
            if (groundRender_visible) {
                if (getWorldPos(node).y < ground->position.y) {
                    node->position.y = ground->position.y - clothPos.y + 0.01;
                    node->velocity = node->velocity * ground->friction;
                }
            }

            // Ball collision
            if (ballRender_visible) {
                Vec3 distVec = getWorldPos(node) - ball->center;
                double distLen = distVec.length();
                double safeDist = ball->radius * 1.1;
                if (distLen < safeDist) {
                    distVec.normalize();
                    setWorldPos(node, distVec * safeDist + ball->center);
                    node->velocity = node->velocity * ball->friction;
                }
            }

        if (handleCollisions) {
            Node* n0 = node;
            if (n0->isFixed) continue;
            for (int j : adjacencyList[i]) {
                Node* n1 = nodes[j];
                if (n1->isFixed) continue;
                
                // check if the distance between two nodes is greater than thickness
                Vec3 diff = n1->position - n0->position;
                double dist = diff.length(); // calculate the distance between two nodes
                double dist2 = dist*dist; // calculate the square of the distance between two nodes
                if (dist2 > thickness2 || dist2 < 1e-12) continue;
                
                // position correction

                // // if you need to compare with restPositions (choose logic)sd
                // double restDist2 = (restPositions[j] - restPositions[i]).lengthSquared();
                // if (dist2 > restDist2) {
                //     // original logic: if the current distance is greater than the initial distance, do not process
                //     continue;
                // }
                // // if restDist2 is smaller than thicknessSquared, then minDist = sqrt(restDist2)
                double minDist = thickness; // minimum distance
                // if (restDist2 < thicknessSquared) {
                //     minDist = std::sqrt(restDist2);
                // }
                diff = diff * (1.0/dist);  // normalize
                double penetration = (minDist - dist); // calculate the penetration depth
                Vec3 correction = diff * (0.5 * penetration);
                n0->position -= correction;
                n1->position += correction;
            }
        }
    }



    } else {
        //#pragma omp parallel for num_threads(num_threads)
        for (int i = 0; i < nodes.size(); i++) {
            Node* node = nodes[i];
            
            // Reset forces
            //node->force = Vec3(0.0, -9.8/100, 0.0);

            // Calculate spring forces
            //std::cout << "connectedSprings: " << node->connectedSprings.size() << std::endl;
            for (Spring* spring : node->connectedSprings) {
                Vec3 force = spring->computeForce(node);
                //std::cout << "Spring " << i << ": Force = " << force << std::endl;
                node->addForce(force);
            }

        }
        //#pragma omp parallel for num_threads(num_threads)
        for (int i = 0; i < nodes.size(); i++) {
            Node* node = nodes[i];
            node->integrate(timeStep);

            // Ground collision
            if (groundRender_visible) {
                if (getWorldPos(node).y < ground->position.y) {
                    node->position.y = ground->position.y - clothPos.y + 0.01;
                    node->velocity = node->velocity * ground->friction;
                }
            }

            // Ball collision
            if (ballRender_visible) {
                Vec3 distVec = getWorldPos(node) - ball->center;
                double distLen = distVec.length();
                double safeDist = ball->radius * 1.1;
                if (distLen < safeDist) {
                    distVec.normalize();
                    setWorldPos(node, distVec * safeDist + ball->center);
                    node->velocity = node->velocity * ball->friction;
                }
            }

        //     // Self collision
        //     for (auto& cell : spatialHash) {
        //         std::vector<Node*>& cellNodes = cell.second;
        //         for (size_t i = 0; i < cellNodes.size(); ++i) {
        //             for (size_t j = i + 1; j < cellNodes.size(); ++j) {
        //                 Node* a = cellNodes[i];
        //                 Node* b = cellNodes[j];

        //                 Vec3 diff = b->position - a->position;
        //                 float dist = diff.length();
        //                 const int maxIterations = 10;
        //                 int iteration = 0;

        //                 while (dist < collisionDistance && iteration < maxIterations) {
        //                     float tolerance = 1e-4f;
        //                     Vec3 correction = diff.normalized() * (collisionDistance - dist + tolerance);

        //                     float springStiffness = 500.0f;
        //                     Vec3 springForce = correction * springStiffness;

        //                     if (!a->isFixed) {
        //                         a->addForce(springForce);
        //                     }
        //                     if (!b->isFixed) {
        //                         b->addForce(--springForce);
        //                     }

        //                     diff = b->position - a->position;
        //                     dist = diff.length();
        //                     iteration++;
        //                 }
        //             }
        //         }
        //     } //
        }
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
                      << "FPS: " << std::fixed << std::setprecision(2) << fps / iterationFreq << "\n"
                      << "total time: " << total_time * 1000 << " ms\n"
                      << "simulation time: " << simulation_time * 1000 << " ms\n"
                      << "Current threads: " << num_threads << "\n"
                      << "total nodes: " << nodes.size() << "\n"
                      << "total springs: " << springs.size() << "\n"
                      << "------------------------\n";
        }
    }
}


//simulation
void spring_based_simulation(double timeStep, Vec3 gravity, Ground* ground, Ball* ball) {
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
                    #pragma omp atomic
                    springs[i]->node1->TempSpringCount--;
                    #pragma omp atomic
                    springs[i]->node2->TempSpringCount--;
                //     omp_unset_lock(&node_locks[springs[i]->node1]);
                // omp_unset_lock(&node_locks[springs[i]->node2]);
                    if (springs[i]->node1->TempSpringCount <= 0) {
                        springs[i]->node1->TempSpringCount = springs[i]->node1->SpringCount;
                        springs[i]->node1->integrate(timeStep);
                        //std::cout << "Node " << i << ": Position = " << springs[i]->node1->position << ", Velocity = " << springs[i]->node1->velocity << std::endl;
                        /** Ground collision **/
                        if (groundRender_visible) { 
                            if (getWorldPos(springs[i]->node2).y < ground->position.y) {
                                springs[i]->node2->position.y = ground->position.y - clothPos.y + 0.01;
                                springs[i]->node2->velocity = springs[i]->node2->velocity * ground->friction;
                            }
                        }
                        
                        /** Ball collision **/
                        if (ballRender_visible) {
                            Vec3 distVec = getWorldPos(springs[i]->node2) - ball->center;
                            double distLen = distVec.length();
                            double safeDist = ball->radius*1.1;
                            if (distLen < safeDist) {
                                distVec.normalize();
                                setWorldPos(springs[i]->node2, distVec*safeDist+ball->center);
                                springs[i]->node2->velocity = springs[i]->node2->velocity * ball->friction;
                            }
                        }
                    }
                    if (springs[i]->node2->TempSpringCount <= 0) {
                        springs[i]->node2->integrate(timeStep);
                        springs[i]->node2->TempSpringCount = springs[i]->node2->SpringCount;
                        /** Ground collision **/
                        if (groundRender_visible) { 
                            if (getWorldPos(springs[i]->node2).y < ground->position.y) {
                                springs[i]->node2->position.y = ground->position.y - clothPos.y + 0.01;
                                springs[i]->node2->velocity = springs[i]->node2->velocity * ground->friction;
                            }
                        }
                        
                        /** Ball collision **/
                        if (ballRender_visible) {
                            Vec3 distVec = getWorldPos(springs[i]->node2) - ball->center;
                            double distLen = distVec.length();
                            double safeDist = ball->radius*1.1;
                            if (distLen < safeDist) {
                                distVec.normalize();
                                setWorldPos(springs[i]->node2, distVec*safeDist+ball->center);
                                springs[i]->node2->velocity = springs[i]->node2->velocity * ball->friction;
                            }
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
                        if (groundRender_visible) { 
                            if (getWorldPos(springs[i]->node2).y < ground->position.y) {
                                springs[i]->node2->position.y = ground->position.y - clothPos.y + 0.01;
                                springs[i]->node2->velocity = springs[i]->node2->velocity * ground->friction;
                            }
                        }
                        
                        /** Ball collision **/
                        if (ballRender_visible) {
                            Vec3 distVec = getWorldPos(springs[i]->node2) - ball->center;
                            double distLen = distVec.length();
                            double safeDist = ball->radius*1.1;
                            if (distLen < safeDist) {
                                distVec.normalize();
                                setWorldPos(springs[i]->node2, distVec*safeDist+ball->center);
                                springs[i]->node2->velocity = springs[i]->node2->velocity * ball->friction;
                            }
                        }
                    }
                    if (springs[i]->node2->TempSpringCount == 0) {
                        springs[i]->node2->integrate(timeStep);
                        springs[i]->node2->TempSpringCount = springs[i]->node2->SpringCount;
                        /** Ground collision **/
                        if (groundRender_visible) { 
                            if (getWorldPos(springs[i]->node2).y < ground->position.y) {
                                springs[i]->node2->position.y = ground->position.y - clothPos.y + 0.01;
                                springs[i]->node2->velocity = springs[i]->node2->velocity * ground->friction;
                            }
                        }
                        
                        /** Ball collision **/
                        if (ballRender_visible) {
                            Vec3 distVec = getWorldPos(springs[i]->node2) - ball->center;
                            double distLen = distVec.length();
                            double safeDist = ball->radius*1.1;
                            if (distLen < safeDist) {
                                distVec.normalize();
                                setWorldPos(springs[i]->node2, distVec*safeDist+ball->center);
                                springs[i]->node2->velocity = springs[i]->node2->velocity * ball->friction;
                            }
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
                        << "FPS: " << std::fixed << std::setprecision(2) << fps / iterationFreq << "\n"
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

    void setworldposandpin(Node* n, Vec3 pos)
    {
        setWorldPos(n, pos);
        n->isFixed = true;
    }

    void unpin(Node* n)
    {
        n->isFixed = false;
    }
    
	void collisionResponse(Ground* ground, Ball* ball)
	{
        for (int i = 0; i < nodes.size(); i++)
        {
            /** Ground collision **/
            // if (getWorldPos(nodes[i]).y < ground->position.y) {
            //     nodes[i]->position.y = ground->position.y - clothPos.y + 0.01;
            //     nodes[i]->velocity = nodes[i]->velocity * ground->friction;
            // }
            
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
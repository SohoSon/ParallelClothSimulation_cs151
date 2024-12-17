#pragma once

#include "Points.h"
#include "Vectors.h"

using namespace std;

class alignas(64) Spring
{
public:
    Node *node1;
    Node *node2;
	double restLen;
    double hookCoef;
    double dampCoef;
    
	Spring(Node *n1, Node *n2, double k)
	{
        node1 = n1;
        node2 = n2;
        node1->SpringCount++;
        node2->SpringCount++;
        node1->TempSpringCount++;
        node2->TempSpringCount++;
        Vec3 currSp = node2->position - node1->position;
        restLen = currSp.length();
        hookCoef = k;
        dampCoef = 5.0;
        //for node-based parallel
        node1->addSpring(this);
        node2->addSpring(this);
	}

	void applyInternalForce(double timeStep) // Compute spring internal force
	{
        double currLen = Vec3::dist(node1->position, node2->position);  // calculate the current length
        Vec3 fDir1 = (node2->position - node1->position)/currLen; // calculate the force direction
        Vec3 diffV1 = node2->velocity - node1->velocity; // calculate the velocity difference
        Vec3 f1 = fDir1 * ((currLen-restLen)*hookCoef + Vec3::dot(diffV1, fDir1)*dampCoef); // calculate the force
        node1->addForce(f1); // add force to node1
        node2->addForce(f1.minus()); // add force to node2
	}
    
    Vec3 computeForce(Node* node) {
        double currLen = Vec3::dist(node1->position, node2->position);  // calculate the current length
        if (currLen == 0) {
            std::cerr << "Warning: Zero length spring encountered!" << std::endl;
            return Vec3(0, 0, 0);
        }
        Vec3 fDir1 = (node2->position - node1->position)/currLen; // calculate the force direction
        Vec3 diffV1 = node2->velocity - node1->velocity; // calculate the velocity difference
        Vec3 f1 = fDir1 * ((currLen-restLen)*hookCoef + Vec3::dot(diffV1, fDir1)*dampCoef); // calculate the force

        // return the force corresponding to the node
        if (node == node1) {
            //cout << "node1 force: " << totalForce.x << ", " << totalForce.y << ", " << totalForce.z << endl;
            return f1;    // apply force to node1
        } else if (node == node2) {
            //cout << "node2 force: " << totalForce.x << ", " << totalForce.y << ", " << totalForce.z << endl;
            return f1.minus();   // apply negative force to node2
        } else {
            cout << "node not in spring" << endl;
            return Vec3(0.0, 0.0, 0.0); // if node not in spring, return zero force
        }
    }
};

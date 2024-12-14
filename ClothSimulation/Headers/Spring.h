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
        double currLen = Vec3::dist(node1->position, node2->position);  // 计算当前长度
        Vec3 fDir1 = (node2->position - node1->position)/currLen; // 计算力方向
        Vec3 diffV1 = node2->velocity - node1->velocity; // 计算速度差
        Vec3 f1 = fDir1 * ((currLen-restLen)*hookCoef + Vec3::dot(diffV1, fDir1)*dampCoef); // 计算力
        node1->addForce(f1); // 将力添加到节点1
        node2->addForce(f1.minus()); // 将力添加到节点2
	}
    
    Vec3 computeForce(Node* node) {
        double currLen = Vec3::dist(node1->position, node2->position);  // 计算当前长度
        if (currLen == 0) {
            std::cerr << "Warning: Zero length spring encountered!" << std::endl;
            return Vec3(0, 0, 0);
        }
        Vec3 fDir1 = (node2->position - node1->position)/currLen; // 计算力方向
        Vec3 diffV1 = node2->velocity - node1->velocity; // 计算速度差
        Vec3 f1 = fDir1 * ((currLen-restLen)*hookCoef + Vec3::dot(diffV1, fDir1)*dampCoef); // 计算力

        // 根据传入的节点返回相应的力
        if (node == node1) {
            //cout << "node1 force: " << totalForce.x << ", " << totalForce.y << ", " << totalForce.z << endl;
            return f1;    // 对 node1 施加正向力
        } else if (node == node2) {
            //cout << "node2 force: " << totalForce.x << ", " << totalForce.y << ", " << totalForce.z << endl;
            return f1.minus();   // 对 node2 施加反向力
        } else {
            cout << "node not in spring" << endl;
            return Vec3(0.0, 0.0, 0.0); // 如果节点不属于弹簧，返回零力
        }
    }
};

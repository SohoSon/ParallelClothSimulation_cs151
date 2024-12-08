#pragma once

#include <vector>
#include <omp.h>
#include <chrono>
#include "Spring.h"
#include "Rigid.h"
#include <iomanip> 
#include <unordered_set>
#include <deque>
#include <atomic>
#include <mutex>

class WorkStealingQueue {
private:
    // Task package that groups related springs with adaptive granularity
    struct SpringTask {
        static const int MIN_BATCH_SIZE = 64;
        static const int MAX_BATCH_SIZE = 256;
        
        std::vector<Spring*> springs;
        std::unordered_set<Node*> shared_nodes;
        double complexity_score;
        int batch_size;  // Dynamic batch size based on system load
        
        SpringTask() : batch_size(MIN_BATCH_SIZE) {
            springs.reserve(MAX_BATCH_SIZE);
        }
        
        // Calculate task complexity considering multiple factors
        void calculateComplexity() {
            complexity_score = 0.0;
            std::unordered_map<Node*, int> node_usage;
            
            // Count node usage frequency
            for (Spring* spring : springs) {
                node_usage[spring->node1]++;
                node_usage[spring->node2]++;
            }
            
            for (Spring* spring : springs) {
                // Base complexity from spring properties
                double spring_complexity = spring->hookCoef * 0.6 + spring->dampCoef * 0.4;
                
                // Node state factors
                if (!spring->node1->isFixed) {
                    spring_complexity *= (1.0 + 0.1 * node_usage[spring->node1]);
                }
                if (!spring->node2->isFixed) {
                    spring_complexity *= (1.0 + 0.1 * node_usage[spring->node2]);
                }
                
                // Memory locality factor
                spring_complexity *= (1.0 + 0.05 * shared_nodes.size());
                
                complexity_score += spring_complexity;
            }
            
            // Adjust for task size
            complexity_score *= (1.0 + springs.size() / static_cast<double>(MAX_BATCH_SIZE));
        }
        
        // Adaptive batch size adjustment based on execution history
        void adjustBatchSize(double execution_time) {
            static const double TARGET_TIME = 0.1;  // milliseconds
            if (execution_time > TARGET_TIME * 1.2) {
                batch_size = std::max(MIN_BATCH_SIZE, batch_size - 2);
            } else if (execution_time < TARGET_TIME * 0.8) {
                batch_size = std::min(MAX_BATCH_SIZE, batch_size + 2);
            }
        }
    };

    // Lock-free task queue implementation
    class LockFreeQueue {
    private:
        static const int QUEUE_SIZE = 1024;
        std::atomic<SpringTask*> tasks[QUEUE_SIZE];
        std::atomic<int> head{0};
        std::atomic<int> tail{0};
        
    public:
        LockFreeQueue() {
            for (int i = 0; i < QUEUE_SIZE; i++) {
                tasks[i].store(nullptr, std::memory_order_relaxed);
            }
        }
        
        bool push(SpringTask* task) {
            int curr_tail = tail.load(std::memory_order_relaxed);
            int next_tail = (curr_tail + 1) % QUEUE_SIZE;
            
            if (next_tail == head.load(std::memory_order_acquire)) {
                return false;  // Queue is full
            }
            
            tasks[curr_tail].store(task, std::memory_order_release);
            tail.store(next_tail, std::memory_order_release);
            return true;
        }
        
        SpringTask* pop() {
            int curr_head = head.load(std::memory_order_relaxed);
            if (curr_head == tail.load(std::memory_order_acquire)) {
                return nullptr;  // Queue is empty
            }
            
            SpringTask* task = tasks[curr_head].load(std::memory_order_acquire);
            head.store((curr_head + 1) % QUEUE_SIZE, std::memory_order_release);
            return task;
        }
        
        bool isEmpty() const {
            return head.load(std::memory_order_relaxed) == 
                   tail.load(std::memory_order_relaxed);
        }
    };

    LockFreeQueue task_queue;
    std::atomic<int> task_count{0};
    std::atomic<double> avg_execution_time{0.0};

public:
    bool addTask(SpringTask* task) {
        if (task_queue.push(task)) {
            task_count.fetch_add(1, std::memory_order_relaxed);
            return true;
        }
        return false;
    }
    
    SpringTask* stealTask() {
        SpringTask* task = task_queue.pop();
        if (task) {
            task_count.fetch_sub(1, std::memory_order_relaxed);
        }
        return task;
    }
    
    void updateExecutionStats(double execution_time) {
        double curr_avg = avg_execution_time.load(std::memory_order_relaxed);
        double new_avg = curr_avg * 0.9 + execution_time * 0.1;
        avg_execution_time.store(new_avg, std::memory_order_relaxed);
    }
};

class Cloth
{
public:
    const int nodesDensity = 10;
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
    int max_color;
    
    std::vector<Node*> nodes;
	std::vector<Spring*> springs;
	std::vector<Node*> faces;
    
    std::vector<std::vector<Spring*>> conflict_graph;
    std::vector<int> color_graph;

    std::vector<int> colors;
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
        conflict_graph = conflictGraph(springs.size());
        color_graph = colorGraph(conflict_graph);
        max_color = *std::max_element(color_graph.begin(), color_graph.end());
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

    std::vector<Spring*> springs_connected_to(Node* point) {
        std::vector<Spring*> connected_springs;
        for (int i = 0; i < springs.size(); i++) {
            if (springs[i]->node1 == point || springs[i]->node2 == point) {
                connected_springs.push_back(springs[i]);
            }
        }
        return connected_springs;
    }

    std::vector<std::vector<Spring*>> conflictGraph(int num_springs) {
        std::vector<std::vector<Spring*>> conflict_graph(num_springs);
        
        // Create a map of nodes to their connected springs for faster lookup
        std::unordered_map<Node*, std::vector<Spring*>> node_to_springs;
        for (int i = 0; i < springs.size(); i++) {
            node_to_springs[springs[i]->node1].push_back(springs[i]);
            node_to_springs[springs[i]->node2].push_back(springs[i]);
        }
        
        for (int spring_id = 0; spring_id < num_springs; ++spring_id) {
            Node* i = springs[spring_id]->node1;
            Node* j = springs[spring_id]->node2;
            
            // Add conflicts from node1's connected springs
            for (Spring* neighbor_spring : node_to_springs[i]) {
                if (neighbor_spring != springs[spring_id]) {
                    conflict_graph[spring_id].push_back(neighbor_spring);
                }
            }
            
            // Add conflicts from node2's connected springs
            for (Spring* neighbor_spring : node_to_springs[j]) {
                if (neighbor_spring != springs[spring_id] && 
                    std::find(conflict_graph[spring_id].begin(), 
                             conflict_graph[spring_id].end(), 
                             neighbor_spring) == conflict_graph[spring_id].end()) {
                    conflict_graph[spring_id].push_back(neighbor_spring);
                }
            }
        }
        
        return conflict_graph;
    }

    std::vector<int> colorGraph(std::vector<std::vector<Spring*>>& conflict_graph) {
        int num_springs = conflict_graph.size();
        std::vector<int> colors(num_springs, -1);
        
        // Calculate degrees for each spring
        std::vector<std::pair<int, int>> degrees; // (degree, spring_id)
        for (int i = 0; i < num_springs; i++) {
            degrees.push_back({conflict_graph[i].size(), i});
        }
        
        // Sort springs by degree in descending order
        std::sort(degrees.begin(), degrees.end(), 
            [](const std::pair<int, int>& a, const std::pair<int, int>& b) { 
                return a.first > b.first; 
            });
        
        // Color the graph
        for (const auto& deg_pair : degrees) {
            int spring_id = deg_pair.second;
            
            // Track used colors by neighbors
            std::vector<bool> used_colors(num_springs, false);
            for (Spring* neighbor_spring : conflict_graph[spring_id]) {
                int neighbor_id = std::find(springs.begin(), springs.end(), neighbor_spring) - springs.begin();
                if (neighbor_id < num_springs && colors[neighbor_id] != -1) {
                    used_colors[colors[neighbor_id]] = true;
                }
            }
            
            // Find first available color
            int color = 0;
            while (color < num_springs && used_colors[color]) {
                color++;
            }
            colors[spring_id] = color;
        }
        
        return colors;
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
                // Pre-fetch data for next iteration
                #pragma omp parallel for schedule(guided)
                for (int i = 0; i < springs.size(); i++) {
                    __builtin_prefetch(&springs[i+1], 0, 3);
                    // Process current spring
                }
                for (int color = 0; color <= max_color; ++color) {
                    #pragma omp parallel for num_threads(num_threads)
                    for (int i = 0; i < springs.size(); i++) {
                        if (color_graph[i] == color) {
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

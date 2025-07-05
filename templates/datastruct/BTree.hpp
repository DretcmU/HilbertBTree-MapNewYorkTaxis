#ifndef BTREE_HPP
#define BTREE_HPP

#include <vector>
#include <string>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <map>
#include <cmath>
#include <limits>
#include <algorithm>
#include <set>
#include <functional>

using namespace std;

struct Point {
    uint64_t key;                 // Índice de Hilbert posicion
    uint64_t indice_h;   
    double lat,lon;    // Coordenadas originales
    int cluster_id;

    Point() : key(0), coords(), cluster_id(-1) {}
    Point(uint64_t k, vector<double> c, int clus_id)
        : key(k), coords(c), cluster_id(clus_id) {}

    bool operator<(const Point& other) const {
        return key < other.key;
    } 

    bool operator==(const Point& other) const {
        return key == other.key;
    }
};

class BTreeNode {
public:
    vector<Point> keys;
    vector<BTreeNode*> children;
    bool leaf;
    int t;

    BTreeNode(int t, bool leaf);
    Point* search(uint64_t key);
    void insertNonFull(Point k);
    void splitChild(int i, BTreeNode* y);
    void traverse();
    ~BTreeNode();
};

class BTree {
public:
    BTreeNode* root;
    int t;

    BTree(int t);
    void insert(Point k);
    Point* search(uint64_t key);
    void traverse();
    ~BTree();
};

uint64_t hilbertIndexND(const vector<uint32_t>& coords, int bits);
vector<Point> leerDesdeBinario(const string& nombreArchivo);
std::vector<Point> buscarRangoHilbert(BTree& arbol, const vector<double>& coord_min, 
    const vector<double>& coord_max, int bits, const vector<double>&, const vector<double>&);
BTree* initHilbertTree(const string& pathBin, int orden);
void cargarMinMax(const string& filename, vector<double>& min_norm, vector<double>& max_norm);

#endif // BTREE_HPP

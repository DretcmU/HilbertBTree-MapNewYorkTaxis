#include "BTree.hpp"
#include <iostream>

BTreeNode::BTreeNode(int t, bool leaf) {
    this->t = t;
    this->leaf = leaf;
}

Point* BTreeNode::search(uint64_t key) {
    int i = 0;
    while (i < keys.size() && key > keys[i].key)
        i++;

    if (i < keys.size() && key == keys[i].key)
        return &keys[i];

    if (leaf)
        return nullptr;

    return children[i]->search(key);
}

void BTreeNode::insertNonFull(Point k) {
    int i = keys.size() - 1;

    if (leaf) {
        keys.push_back(k);
        while (i >= 0 && k < keys[i]) {
            keys[i + 1] = keys[i];
            i--;
        }
        keys[i + 1] = k;
    } else {
        while (i >= 0 && k < keys[i])
            i--;

        i++;
        if (children[i]->keys.size() == 2 * t - 1) {
            splitChild(i, children[i]);
            if (k.key > keys[i].key)
                i++;
        }
        children[i]->insertNonFull(k);
    }
}

void BTreeNode::splitChild(int i, BTreeNode* y) {
    BTreeNode* z = new BTreeNode(y->t, y->leaf);

    for (int j = 0; j < t - 1; j++)
        z->keys.push_back(y->keys[j + t]);

    if (!y->leaf) {
        for (int j = 0; j < t; j++)
            z->children.push_back(y->children[j + t]);
    }

    Point midKey = y->keys[t - 1];

    y->keys.resize(t - 1);
    if (!y->leaf)
        y->children.resize(t);

    children.insert(children.begin() + i + 1, z);
    keys.insert(keys.begin() + i, midKey);
}

void BTreeNode::traverse() {
    int i;
    for (i = 0; i < keys.size(); i++) {
        if (!leaf)
            children[i]->traverse();

        cout << "Key: " << keys[i].key << " | Coords: ";
        cout << keys[i].lat << ", "<<keys[i].lon;

        cout << " | cluster ID: "<<keys[i].cluster_id<<endl;
        // for (auto& d : keys[i].data) cout << d << " ";
        // cout << endl;
    }

    if (!leaf)
        children[i]->traverse();
}

BTree::BTree(int t) {
    root = nullptr;
    this->t = t;
}

void BTree::insert(Point k) {
    if (!root) {
        root = new BTreeNode(t, true);
        root->keys.push_back(k);
    } else {
        if (root->keys.size() == 2 * t - 1) {
            BTreeNode* s = new BTreeNode(t, false);
            s->children.push_back(root);
            s->splitChild(0, root);

            int i = (k.key > s->keys[0].key) ? 1 : 0;
            s->children[i]->insertNonFull(k);
            root = s;
        } else {
            root->insertNonFull(k);
        }
    }
}

Point* BTree::search(uint64_t key) {
    return root ? root->search(key) : nullptr;
}

void BTree::traverse() {
    if (root) root->traverse();
}

BTreeNode::~BTreeNode() {
    for (auto child : children)
        delete child;
}

BTree::~BTree() {
    delete root;
}


uint64_t hilbertIndexND(const vector<uint32_t>& coords, int bits) {
    int n = coords.size();
    vector<uint32_t> x(coords);
    uint64_t index = 0;
    uint32_t mask = 1U << (bits - 1);

    for (int i = 0; i < bits; ++i) {
        uint32_t bitword = 0;
        for (int j = 0; j < n; ++j) {
            uint32_t bit = (x[j] & mask) ? 1 : 0;
            bitword |= (bit << (n - j - 1));
        }

        uint32_t gray = bitword ^ (bitword >> 1);
        index = (index << n) | gray;

        for (int j = 0; j < n; ++j) {
            if ((bitword >> (n - j - 1)) & 1) {
                x[j] ^= mask;
            }
        }

        mask >>= 1;
    }

    return index;
}


std::vector<Point> buscarRangoHilbert(BTree& arbol, 
    const vector<double>& coord_min, const vector<double>& coord_max,
     int bits, const vector<double>& min_norm, const vector<double>& max_norm) {
    if (coord_min.size() != coord_max.size()) {
        cerr << "Error: las coordenadas deben tener la misma dimensión.\n";
        return std::vector<Point>{};
    }
    vector<double> hcoord_min(2),  hcoord_max(2);
    for(int i=0;i<2;i++){
        hcoord_min[i] = coord_min[i];
        hcoord_max[i] = coord_max[i];
    }

    double expand=0.01;
    for(int i=0;i<2;i++){
        hcoord_min[i] -= expand;
        hcoord_max[i] += expand;
    }

    // Normalización: convertir coordenadas reales a enteros entre 0 y 2^bits - 1
    auto normalizar = [bits](const vector<double>& coord, const vector<double>& min_vals, const vector<double>& max_vals) {
        vector<uint32_t> normalizada;
        for (size_t i = 0; i < coord.size(); ++i) {
            double norm = (coord[i] - min_vals[i]) / (max_vals[i] - min_vals[i]);
            norm = max(0.0, min(1.0, norm));  // clamp entre 0 y 1
            normalizada.push_back(static_cast<uint32_t>(norm * ((1 << bits) - 1)));
        }
        return normalizada;
    };

    vector<double> min_vals = min_norm;
    vector<double> max_vals = max_norm;

    vector<uint32_t> norm_min = normalizar(hcoord_min, min_vals, max_vals);
    vector<uint32_t> norm_max = normalizar(hcoord_max, min_vals, max_vals);

    uint64_t hilbert_min = hilbertIndexND(norm_min, bits);
    uint64_t hilbert_max = hilbertIndexND(norm_max, bits);

    if (hilbert_min > hilbert_max) std::swap(hilbert_min, hilbert_max);

    cout << "Buscando puntos con Hilbert entre " << hilbert_min << " y " << hilbert_max << endl;

    // Recorrer el árbol y verificar
    vector<Point> encontrados;

    // function<void(BTreeNode*)> buscarRec = [&](BTreeNode* nodo) {
    //     if (!nodo) return;

    //     int i = 0;
    //     for (; i < nodo->keys.size(); ++i) {
    //         if (!nodo->leaf)
    //             buscarRec(nodo->children[i]);

    //         if (nodo->keys[i].key >= hilbert_min && nodo->keys[i].key <= hilbert_max) {
    //             const auto& p = nodo->keys[i];
    //             bool dentro = true;
    //             //for (size_t d = 0; d < p.coords.size(); ++d) {
    //             for (size_t d = 0; d < 2; ++d) {
    //                 //cout<<p.coords[d] << " " <<coord_min[d] << " " << p.coords[d]<< " " <<coord_max[d]<<endl;
    //                 if (p.coords[d] < coord_min[d] || p.coords[d] > coord_max[d]) {
    //                     dentro = false;
    //                     break;
    //                 }
    //             }
    //             if (dentro)
    //                 encontrados.push_back(p);
    //         }
    //     }

    //     if (!nodo->leaf)
    //         buscarRec(nodo->children[i]);
    // };

    // buscarRec(arbol.root);
    
    function<void(BTreeNode*)> buscarRec = [&](BTreeNode* nodo) {
        if (!nodo) return;
    
        int i = 0;
    
        // Buscar la primera clave >= hilbert_min usando búsqueda secuencial (puedes usar binary search si hay muchos keys)
        while (i < nodo->keys.size() && nodo->keys[i].key < hilbert_min) {
            if (!nodo->leaf)
                buscarRec(nodo->children[i]);
            ++i;
        }
    
        // Desde aquí las claves están >= hilbert_min
        while (i < nodo->keys.size() && nodo->keys[i].key <= hilbert_max) {
            // Revisar hijo izquierdo
            if (!nodo->leaf)
                buscarRec(nodo->children[i]);
    
            // Revisar si el punto está dentro del rango original (por coordenadas reales)
            const auto& p = nodo->keys[i];

            if (!(p.lat < coord_min[0] || p.lat > coord_max[0] 
                || p.lon < coord_min[1] || p.lon > coord_max[1])) {
                encontrados.push_back(p);
            }            
    
            ++i;
        }
    
        // Después de hilbert_max, no vale la pena seguir
        if (!nodo->leaf && i < nodo->children.size())
            buscarRec(nodo->children[i]);
    };
    buscarRec(arbol.root);

    cout << "Se encontraron " << encontrados.size() << " puntos dentro del rango:\n";
    // for (const auto& p : encontrados) {
    //     cout << "Key: " << p.key << " | Coords: ";
    //     for (auto v : p.coords) cout << v << " ";
    //     cout << " | Cluster ID: " << p.cluster_id << endl;
    // }

    return encontrados;
}

vector<Point> leerDesdeBinario(const string& nombreArchivo) {
    vector<Point> datos;
    ifstream inBin(nombreArchivo, ios::binary);
    if (!inBin.is_open()) {
        cerr << "No se pudo abrir el archivo binario: " << nombreArchivo << endl;
        return datos;
    }

    uint64_t n_datos, n_columnas;
    inBin.read(reinterpret_cast<char*>(&n_datos), sizeof(uint64_t));
    inBin.read(reinterpret_cast<char*>(&n_columnas), sizeof(uint64_t));

    for (uint64_t i = 0; i < n_datos; ++i) {
        Point d;
        inBin.read(reinterpret_cast<char*>(&d.key), sizeof(uint64_t));
        inBin.read(reinterpret_cast<char*>(&d.indice_h), sizeof(uint64_t));
        inBin.read(reinterpret_cast<char*>(&d.lat), sizeof(double));
        inBin.read(reinterpret_cast<char*>(&d.lon), sizeof(double));
        inBin.read(reinterpret_cast<char*>(&d.cluster_id), sizeof(int));
        datos.push_back(d);
    }

    inBin.close();
    cout << "Archivo binario leído correctamente: " << nombreArchivo << "\n";
    return datos;
}


BTree* initHilbertTree(const string& pathBin = "../preprocessing/BinDatos.bin", int orden = 2) {
    BTree* tree = new BTree(orden);
    vector<Point> puntos = leerDesdeBinario(pathBin);


    for (const auto& punto : puntos) {
        tree->insert(punto);
    }

    cout << "Árbol Hilbert-BTree inicializado con " << puntos.size() << " puntos." << endl;
    return tree;
}

void cargarMinMax(const string& filename, vector<double>& min_norm, vector<double>& max_norm) {
    ifstream archivo(filename);
    if (!archivo.is_open()) {
        cerr << "No se pudo abrir el archivo para leer: " << filename << endl;
        return;
    }

    // min_norm.resize(2);
    // max_norm.resize(2);

    archivo >> min_norm[0] >> min_norm[1];
    archivo >> max_norm[0] >> max_norm[1];

    archivo.close();
}
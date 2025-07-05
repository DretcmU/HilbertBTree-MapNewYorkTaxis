#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <map>
#include <string>
#include <cmath>
#include <limits>
#include <algorithm>
#include <cstdint>
#include <set>
#include <chrono>

using namespace std;

// ===================
// ESTRUCTURA PRINCIPAL
// ===================

struct Dato {
    double lat;
    double lon;
    vector<uint32_t> columnas;
    uint64_t indice_h = 0;
    uint64_t key;
    int cluster_id = -1; // -1 = sin asignar, -2 = ruido
    bool visitado = false;
};

// ===================
// DISTANCIA EUCLIDEANA
// ===================
double distancia(const Dato& a, const Dato& b) {
    double suma = 0.0;
    for (size_t i = 0; i < a.columnas.size(); ++i) {
        double d = a.columnas[i] - b.columnas[i];
        suma += d * d;
    }
    return sqrt(suma);
}


struct KDNode {
    int index;
    KDNode* left = nullptr;
    KDNode* right = nullptr;
    int dimension;

    KDNode(int idx, int dim) : index(idx), dimension(dim) {}
};

class KDTree {
public:
    KDTree(const vector<Dato*>& datos) : datos(datos), dimensiones(datos[0]->columnas.size()) {
        indices.resize(datos.size());
        for (int i = 0; i < indices.size(); ++i) indices[i] = i;
        root = build(0, datos.size(), 0);
    }

    vector<int> pointsInSphere(const Dato& q, double eps) const {
        vector<int> result;
        search(root, q, eps, result);
        return result;
    }

private:
    const vector<Dato*>& datos;
    vector<int> indices;
    KDNode* root;
    int dimensiones;

    KDNode* build(int l, int r, int depth) {
        if (l >= r) return nullptr;
        int dim = depth % dimensiones;
        int m = (l + r) / 2;

        auto comp = [&](int a, int b) {
            return datos[a]->columnas[dim] < datos[b]->columnas[dim];
        };

        nth_element(indices.begin() + l, indices.begin() + m, indices.begin() + r, comp);

        KDNode* node = new KDNode(indices[m], dim);
        node->left = build(l, m, depth + 1);
        node->right = build(m + 1, r, depth + 1);
        return node;
    }

    void search(KDNode* node, const Dato& q, double eps, vector<int>& result) const {
        if (!node) return;
        Dato* p = datos[node->index];
        if (distancia(q, *p) <= eps)
            result.push_back(node->index);

        double delta = q.columnas[node->dimension] - p->columnas[node->dimension];

        if (delta <= eps) search(node->left, q, eps, result);
        if (delta >= -eps) search(node->right, q, eps, result);
    }
};

// ===================
// HILBERT ND
// ===================
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

void guardarMinMax(const string& filename, const vector<double>& minVals, const vector<double>& maxVals) {
    ofstream archivo(filename);
    if (!archivo.is_open()) {
        cerr << "No se pudo abrir el archivo para escribir: " << filename << endl;
        return;
    }

    archivo << minVals[0] << " " << minVals[1] << endl;
    archivo << maxVals[0] << " " << maxVals[1] << endl;

    archivo.close();
}

vector<Dato> leerColumnasDeCSVConHilbert(
    const string& filename,
    const vector<string>& columnasDeseadas,
    int bitsHilbert = 16
) {
    vector<Dato> datos;
    ifstream archivo(filename);
    string linea;

    if (!archivo.is_open()) {
        cerr << "No se pudo abrir el archivo: " << filename << endl;
        return datos;
    }

    // Leer encabezado
    getline(archivo, linea);
    istringstream ssEncabezado(linea);
    string columna;
    unordered_map<string, int> indiceColumna;
    int index = 0;

    while (getline(ssEncabezado, columna, ',')) {
        indiceColumna[columna] = index++;
    }

    // Verificar que todas las columnas deseadas existen
    vector<int> indicesSeleccionados;
    for (const auto& nombre : columnasDeseadas) {
        if (indiceColumna.find(nombre) != indiceColumna.end()) {
            indicesSeleccionados.push_back(indiceColumna[nombre]);
        } else {
            cerr << "Columna no encontrada: " << nombre << endl;
            return datos;
        }
    }

    // Leer todas las filas y filtrar válidas (columna 1 y 2 no deben ser 0.0 o vacías)
    vector<vector<string>> filasValidas;
    while (getline(archivo, linea)) {
        istringstream ss(linea);
        string valor;
        vector<string> filaCompleta;

        while (getline(ss, valor, ',')) {
            filaCompleta.push_back(valor);
        }

        // Validación: verificar que haya suficientes columnas
        if (filaCompleta.size() <= 2) continue;

        try {
            if (!filaCompleta[1].empty() && !filaCompleta[2].empty()) {
                double val1 = stod(filaCompleta[1]);
                double val2 = stod(filaCompleta[2]);

                if (val1 != 0.0 && val2 != 0.0) {
                    filasValidas.push_back(filaCompleta);
                }
            }
        } catch (...) {
            // Ignorar si no se puede convertir alguna columna
            continue;
        }
    }

    // Extraer valores numéricos para columnas deseadas
    vector<vector<double>> columnasRaw(columnasDeseadas.size());

    for (const auto& fila : filasValidas) {
        for (size_t i = 0; i < indicesSeleccionados.size(); ++i) {
            int idx = indicesSeleccionados[i];
            try {
                double val = stod(fila[idx]);
                columnasRaw[i].push_back(val);
            } catch (...) {
                columnasRaw[i].push_back(0.0);  // Si hay error, usa 0.0 (o podrías ignorar esa fila si prefieres)
            }
        }
    }

    // Calcular min y max para normalización
    vector<double> minVals, maxVals;
    for (const auto& col : columnasRaw) {
        double minV = *min_element(col.begin(), col.end());
        double maxV = *max_element(col.begin(), col.end());
        minVals.push_back(minV);
        maxVals.push_back(maxV);
    }

    int xdatoxd=0;
    // vector<double> min_norm = {39.9, -74.9};  // latitud mínima, longitud mínima
    // vector<double> max_norm = {41.1, -72.9};  // latitud máxima, longitud máxima
    // for(int i=0;i<2;i++){
    //     if(minVals[i]<min_norm[i]) minVals[i] = min_norm[i];
    //     if(maxVals[i]>max_norm[i]) maxVals[i] = max_norm[i];
    // }

    cout<<"Min: " <<minVals[0]<<" "<<minVals[1]<<endl;
    cout<<"Max: " <<maxVals[0]<<" "<<maxVals[1]<<endl;
    //cout<<"mira: ";cin>>xdatoxd;

    guardarMinMax("minmax.txt", minVals, maxVals);

    // Construir objetos Dato
    for (size_t i = 0; i < filasValidas.size(); ++i) {
        const auto& fila = filasValidas[i];
        if(minVals[0]>stod(fila[1]) || stod(fila[1])>maxVals[0] || 
        minVals[1]>stod(fila[2]) || stod(fila[2])>maxVals[1]){
            continue;
        }
        Dato dato;
        //vector<uint32_t> coords;
        vector<uint32_t> pos;

        for (size_t j = 0; j < indicesSeleccionados.size(); ++j) {
            int idx = indicesSeleccionados[j];
            double val = stod(fila[idx]);

            // Normalización
            uint32_t norm = static_cast<uint32_t>(
                ((val - minVals[j]) / (maxVals[j] - minVals[j])) * ((1 << bitsHilbert) - 1)
            );

            if (j < 2) {
                if (j==0)
                    dato.lat = val;
                else
                    dato.lon = val;

                pos.push_back(norm);
            } else {
                dato.columnas.push_back(norm);
            }
        }

        dato.key = hilbertIndexND(pos, bitsHilbert);
        dato.indice_h = hilbertIndexND(dato.columnas, bitsHilbert);
        datos.push_back(dato);
    }

    return datos;
}



// ============================
// EXPANSIÓN DEL CLUSTER (KD)
// ============================
void expandirClusterKD(vector<Dato*>& grupo, KDTree& tree, int index,
                       vector<int>& vecinos, int cluster_id, double eps, int minPts) {
    grupo[index]->cluster_id = cluster_id;

    size_t i = 0;
    while (i < vecinos.size()) {
        int idx = vecinos[i];
        Dato* actual = grupo[idx];
        if (!actual->visitado) {
            actual->visitado = true;
            auto nuevos_vecinos = tree.pointsInSphere(*actual, eps);
            if (nuevos_vecinos.size() >= minPts) {
                vecinos.insert(vecinos.end(), nuevos_vecinos.begin(), nuevos_vecinos.end());
            }
        }

        if (actual->cluster_id == -1 || actual->cluster_id == -2) {
            actual->cluster_id = cluster_id;
        }

        ++i;
    }
}

// ============================
// DBSCAN CON KD-TREE
// ============================
void dbscanKD(vector<Dato*>& grupo, double eps, int minPts, int& nextClusterId) {
    KDTree tree(grupo);

    for (int i = 0; i < grupo.size(); ++i) {
        Dato* punto = grupo[i];
        if (punto->visitado) continue;

        punto->visitado = true;
        auto vecinos = tree.pointsInSphere(*punto, eps);

        if (vecinos.size() < minPts) {
            punto->cluster_id = -2; // ruido
        } else {
            expandirClusterKD(grupo, tree, i, vecinos, nextClusterId, eps, minPts);
            nextClusterId++;
        }
    }
}

void guardarEnBinario(const vector<Dato>& datos, const string& nombreArchivo) {
    ofstream outBin(nombreArchivo, ios::binary);
    if (!outBin.is_open()) {
        cerr << "No se pudo crear el archivo binario de salida: " << nombreArchivo << endl;
        return;
    }

    // También abrimos un archivo .txt adicional
    string nombreTxt = nombreArchivo + ".txt";  // Ejemplo: "BinDatos.bin.txt"
    ofstream outTxt(nombreTxt);
    if (!outTxt.is_open()) {
        cerr << "No se pudo crear el archivo de texto: " << nombreTxt << endl;
        return;
    }

    // Encabezado binario
    uint64_t n_columnas = 5;
    uint64_t n_datos = datos.size();
    outBin.write(reinterpret_cast<const char*>(&n_datos), sizeof(uint64_t));
    outBin.write(reinterpret_cast<const char*>(&n_columnas), sizeof(uint64_t));

    // Datos
    for (const Dato& d : datos) {
        outBin.write(reinterpret_cast<const char*>(&d.key), sizeof(uint64_t));
        outBin.write(reinterpret_cast<const char*>(&d.indice_h), sizeof(uint64_t));
        outBin.write(reinterpret_cast<const char*>(&d.lat), sizeof(double));
        outBin.write(reinterpret_cast<const char*>(&d.lon), sizeof(double));
        outBin.write(reinterpret_cast<const char*>(&d.cluster_id), sizeof(int));

        // También escribir en .txt para revisión humana
        outTxt << d.key << " " << d.indice_h << " " << d.lat << " " << d.lon << " " << d.cluster_id << "\n";
    }

    outBin.close();
    outTxt.close();
    cout << "Archivo binario guardado como '" << nombreArchivo << "'\n";
    cout << "Archivo de texto guardado como '" << nombreTxt << "'\n";
}

int estimarMinPts(size_t cantidadDatos) {
    return std::max(4, static_cast<int>(log2(cantidadDatos))); // al menos 4
}

double estimarEpsilon1D(const vector<Dato*>& grupo, int k) {
    if (grupo.size() < 2) return 1.0;

    // Primero, normalizamos los valores de indice_h
    uint64_t min_h = grupo[0]->indice_h;
    uint64_t max_h = grupo[0]->indice_h;

    for (Dato* d : grupo) {
        min_h = std::min(min_h, d->indice_h);
        max_h = std::max(max_h, d->indice_h);
    }

    double rango = static_cast<double>(max_h - min_h);
    if (rango == 0.0) return 0.01;  // todos iguales, no hay separación

    vector<double> distancias;
    k = std::min(k, static_cast<int>(grupo.size() - 1));

    for (Dato* p : grupo) {
        vector<double> vecinos;
        double p_norm = (p->indice_h - min_h) / rango;
        for (Dato* q : grupo) {
            if (p == q) continue;
            double q_norm = (q->indice_h - min_h) / rango;
            vecinos.push_back(abs(p_norm - q_norm));
        }

        std::sort(vecinos.begin(), vecinos.end());

        if (vecinos.size() >= k)
            distancias.push_back(vecinos[k - 1]);
        else if (!vecinos.empty())
            distancias.push_back(vecinos.back());
    }

    if (distancias.empty()) return 0.01;

    std::sort(distancias.begin(), distancias.end());
    return distancias[distancias.size() * 0.9]; // valor en el percentil 90
}

int main() {
    string archivoCSV = "processed_data_subset_500k.csv";
    vector<string> columnas = {
        "pickup_latitude", "pickup_longitude",
        "passenger_count", "trip_distance", "fare_amount",
        "extra", "mta_tax", "tip_amount", "tolls_amount",
        "improvement_surcharge", "total_amount"
    };
    int bitsHilbert = 16;

    //auto [minVals, maxVals] = calcularMinMax(archivoCSV, columnas);
    //guardarMinMax("minmax.txt", minVals, maxVals);

    //vector<Dato> datos = leerYConstruirDatos(archivoCSV, columnas, minVals, maxVals, bitsHilbert);

    int agrupamiento_por_bits = 2; // numero de digitios que restaremos al (indice hilbert x cantidad de columnas)

    vector<Dato> datos = leerColumnasDeCSVConHilbert(archivoCSV, columnas, bitsHilbert);
    
    // Agrupar punteros por bits altos
    map<uint64_t, vector<Dato*>> gruposHilbert;

    for (Dato& dato : datos) {
        uint64_t grupo_id = dato.indice_h >> (bitsHilbert * (columnas.size()) - agrupamiento_por_bits);
        gruposHilbert[grupo_id].push_back(&dato);
    }

    // Aplicar DBSCAN a cada grupo
    int nextClusterId = 0, minPts=4;
    double eps = 0.2;
    auto inicio = std::chrono::high_resolution_clock::now();

    // APLICAR DBSCAN A CADA GRUPO
    for (auto it = gruposHilbert.begin(); it != gruposHilbert.end(); ++it) {
        uint64_t grupo_id = it->first;
        vector<Dato*>& grupo = it->second;
        
        //cout<<"A\n";
        int minPts = estimarMinPts(grupo.size());
        //cout<<"B\n";
        double eps = estimarEpsilon1D(grupo, minPts);        
        //cout << "Grupo " << grupo_id << " -> size: " << grupo.size()<< " | minPts: " << minPts << " | eps: " << eps << endl;
        dbscanKD(grupo, eps, minPts, nextClusterId);
    }

    // FIN DEL CRONÓMETRO
    auto fin = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duracion = fin - inicio;

    // Mostrar resumen
    map<int, int> conteoClusters;
    for (const auto& d : datos) {
        conteoClusters[d.cluster_id]++;
    }

    // Reasignar clusters de un solo punto como ruido (-2)
    for (auto& par : conteoClusters) {
        int id = par.first;
        int count = par.second;

        if (id >= 0 && count == 1) { // Solo los clusters válidos de un solo punto
            for (Dato& d : datos) {
                if (d.cluster_id == id) {
                    d.cluster_id = -2;
                    break; // solo uno por definición
                }
            }
        }
    }

    conteoClusters.clear();
    for (const auto& d : datos) {
        conteoClusters[d.cluster_id]++;
    }
    long con=0;
    for (auto it = conteoClusters.begin(); it != conteoClusters.end(); ++it) {
        int id = it->first;
        int count = it->second;
        cout << "Cluster " << id << ": " << count << " puntos\n";
        con++;
    }
    //cout<<conteoClusters.size()<<endl;
    cout<<con<<endl;
    ofstream out("numGrups.txt");
    out << con<< endl;
    out.close();
    // for (size_t i = 0; i < 150 && i < datos.size(); ++i) {
    //     cout << "Fila " << i << " cluster = " << datos[i].cluster_id << endl;
    // }

    // sort(datos.begin(), datos.end(), [](const Dato& a, const Dato& b) {
    //     return a.key < b.key;
    // });
    guardarEnBinario(datos, "BinDatos.bin");
    cout<<"size: "<<datos.size()<<endl;
    std::cout << "Tiempo total de DBSCAN + KDTree: " << duracion.count() << " segundos" << std::endl;
    
    return 0;
}
#include "utilserver/httplib.h"
#include "utilserver/json.hpp"
#include "datastruct/BTree.hpp"
#include <iostream>
#include <random>

using json = nlohmann::json;

int main() {
    BTree* arbol = initHilbertTree("templates/preprocessing/BinDatos.bin",2);
    vector<double> min_norm(2), max_norm(2);
    cargarMinMax("templates/preprocessing/minmax.txt", min_norm, max_norm);

    std::cout << ">> Servidor iniciando..." << std::endl;

    httplib::Server svr;

    svr.Options("/consultar", [](const httplib::Request& req, httplib::Response& res) {
        res.set_header("Access-Control-Allow-Origin", "*");
        res.set_header("Access-Control-Allow-Methods", "POST, OPTIONS");
        res.set_header("Access-Control-Allow-Headers", "Content-Type");
        res.status = 204; // No Content
    });

    svr.Post("/consultar", [arbol,min_norm,max_norm](const httplib::Request& req, httplib::Response& res) {
        res.set_header("Access-Control-Allow-Origin", "*");

        try {
            json recibido = json::parse(req.body);
            double minLat = recibido["minLat"];
            double maxLat = recibido["maxLat"];
            double minLon = recibido["minLon"];
            double maxLon = recibido["maxLon"];

            vector<double> min_coord = {minLat, minLon}; // lat, lon
            vector<double> max_coord = {maxLat, maxLon}; // lat, lon
            

            std::cout << "Recibido: [" << minLat << ", " << maxLat << "] x ["
                      << minLon << ", " << maxLon << "]\n";

            // std::random_device rd;
            // std::mt19937 gen(rd());
            // std::uniform_int_distribution<> count_dist(3, 7);
            // std::uniform_real_distribution<> lat_dist(minLat, maxLat);
            // std::uniform_real_distribution<> lon_dist(minLon, maxLon);
            // std::uniform_int_distribution<> color_dist(0, 2);
            
            int bits = 16;
            auto encontrados = buscarRangoHilbert(*arbol, min_coord, 
                max_coord, bits, min_norm, max_norm);

            json respuesta = json::array();
            for (const auto& p : encontrados) {
                respuesta.push_back({
                    {"lat", p.lat},
                    {"lon", p.lon},
                    {"color", p.cluster_id},
                    {"indice_h", p.indice_h}
                });
            }

            res.set_content(respuesta.dump(), "application/json");
        } catch (const std::exception& e) {
            res.status = 400;
            res.set_content(std::string("Error: ") + e.what(), "text/plain");
        }
    });

    std::cout << ">> Servidor escuchando en http://localhost:8080" << std::endl;
    svr.listen("0.0.0.0", 8080);
    return 0;
}

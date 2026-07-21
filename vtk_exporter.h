#ifndef VTK_EXPORTER_H
#define VTK_EXPORTER_H

#include <string>
#include <eigen3/Eigen/Core>

class VTKExporter {
private:
    std::string output_dir;
    int frame_counter = 0;

    int nx;      // Número de nodos en X
    int ny;      // Número de nodos en Y
    double dx;   // Espaciado entre nodos en X
    double dy;   // Espaciado entre nodos en Y
    int print_every; // Cada cuántos frames se reporta progreso en consola 0 = nunca

public:
    // Constructor: directorio de salida y geometría de la malla
    explicit VTKExporter(const std::string& output_directory,
                          int nx, int ny, double dx, double dy, int print_every);

    void write_frame(const Eigen::VectorXd& field, double time,
                      const std::string& field_name = "Temperatura");

    
    void reset_counter();
};

#endif


# Se usa una imagen completa con el compilador solo para construir el binario.
FROM ubuntu:24.04 AS builder

# g++-14 es necesario porque el código usa <print> de C++23, que g++-13
# (el que viene por defecto en Ubuntu 24.04) todavía no soporta.

RUN apt-get update && apt-get install -y --no-install-recommends \
        g++-14 \
        libeigen3-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build
#COPY Neumann_VTK.cpp rk4.h vtk_exporter.cpp vtk_exporter.h ./
COPY Dirichlet_VTK.cpp euler.h rk4.h vtk_exporter.cpp vtk_exporter.h ./

#RUN g++-14 -std=c++23 -O3 -march=native Neumann_VTK.cpp vtk_exporter.cpp -o neumann.x
RUN g++-14 -std=c++23 -O3 -march=native Dirichlet_VTK.cpp vtk_exporter.cpp -o dirichlet.x


FROM ubuntu:24.04 AS runtime

RUN apt-get update && apt-get install -y --no-install-recommends \
        libstdc++6 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
#COPY --from=builder /build/neumann.x .
COPY --from=builder /build/dirichlet.x .

#CMD ["./neumann.x"]
CMD ["./dirichlet.x"]


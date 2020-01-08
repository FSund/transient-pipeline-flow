#include "heattransfer/pipewall.hpp"

const PipeWall PipeWall::defaultPipeWall = PipeWall(
    std::vector<Layer>(
        {
            Layer(0.024, Material::steel), // 24 mm steel
            Layer(0.007, Material::coating), // 7 mm coating
            Layer(0.08, Material::concrete) // 80 mm concrete
        }
    )
);

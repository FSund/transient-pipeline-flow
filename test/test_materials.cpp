#include "debug.hpp"
#include <iostream>
#include "heattransfer/pipewall.hpp"
#include "heattransfer/material.hpp"
#include "heattransfer/ambientfluid.hpp"
#include "heattransfer/burialmedium.hpp"

using namespace std;

TEST_CASE("Material")
{
    SUBCASE("Constructor")
    {
        Material mat(1, 2, 3);
        CHECK(mat.conductivity() == 1);
        CHECK(mat.density() == 2);
        CHECK(mat.heatCapacity() == 3);
    }
}

TEST_CASE("AmbientFluid")
{
    SUBCASE("Constructor 1")
    {
        AmbientFluid mat(1, 2, 3, 4, 5);
        CHECK(mat.velocity() == 1);
        CHECK(mat.viscosity() == 2);
        CHECK(mat.conductivity() == 3);
        CHECK(mat.density() == 4);
        CHECK(mat.heatCapacity() == 5);
    }

    SUBCASE("Constructor 2")
    {
        AmbientFluid mat(1, 2, Material(3, 4, 5));
        CHECK(mat.velocity() == 1);
        CHECK(mat.viscosity() == 2);
        CHECK(mat.conductivity() == 3);
        CHECK(mat.density() == 4);
        CHECK(mat.heatCapacity() == 5);
    }
}

TEST_CASE("BurialMedium")
{
    SUBCASE("Constructor 1")
    {
        BurialMedium mat(1, 2, 3);
        CHECK(mat.conductivity() == 1);
        CHECK(mat.density() == 2);
        CHECK(mat.heatCapacity() == 3);
    }

    SUBCASE("Constructor 2")
    {
        BurialMedium mat(Material(1, 2, 3));
        CHECK(mat.conductivity() == 1);
        CHECK(mat.density() == 2);
        CHECK(mat.heatCapacity() == 3);
    }
}

TEST_CASE("PipeWall")
{
    SUBCASE("Layer constructor 1")
    {
        PipeWall::Layer layer;
        CHECK(layer.thickness() == -1);
        CHECK(layer.conductivity() == -1);
        CHECK(layer.density() == -1);
        CHECK(layer.heatCapacity() == -1);
    }

    SUBCASE("Layer constructor 2")
    {
        PipeWall::Layer layer(1, 2, 3, 4);
        CHECK(layer.thickness() == 1);
        CHECK(layer.conductivity() == 2);
        CHECK(layer.density() == 3);
        CHECK(layer.heatCapacity() == 4);
    }

    SUBCASE("Layer constructor 3")
    {
        Material material(2, 3, 4);
        PipeWall::Layer layer(1, material);
        CHECK(layer.thickness() == 1);
        CHECK(layer.conductivity() == 2);
        CHECK(layer.density() == 3);
        CHECK(layer.heatCapacity() == 4);
    }

    SUBCASE("Layer setters")
    {
        PipeWall::Layer layer;
        layer.thickness() = 1;
        CHECK(layer.thickness() == 1);

        layer.conductivity() = 2;
        CHECK(layer.conductivity() == 2);

        layer.density() = 3;
        CHECK(layer.density() == 3);

        layer.heatCapacity() = 4;
        CHECK(layer.heatCapacity() == 4);
    }

    SUBCASE("Pipewall constructor 1")
    {
        PipeWall pipe(5);
        CHECK(pipe.layers().size() == 5);
        CHECK(pipe.size() == 5);
    }

    SUBCASE("Pipewall constructor 2")
    {
        PipeWall pipe({PipeWall::Layer(1, 2, 3, 4), PipeWall::Layer(2, 3, 4, 5)});
        CHECK(pipe.layers().size() == 2);
        CHECK(pipe.size() == 2);
        CHECK(pipe.layer(0).thickness() == 1);
        CHECK(pipe.layer(0).conductivity() == 2);
        CHECK(pipe.layer(0).density() == 3);
        CHECK(pipe.layer(0).heatCapacity() == 4);

        CHECK(pipe.layer(1).thickness() == 2);
        CHECK(pipe.layer(1).conductivity() == 3);
        CHECK(pipe.layer(1).density() == 4);
        CHECK(pipe.layer(1).heatCapacity() == 5);
    }
}

#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/grid/Partitioner.h"


#include "atlas/util/Config.h"
#include "tests/AtlasTestEnvironment.h"

using namespace atlas::mesh;

namespace atlas {
namespace test {

// ------------------------------------------------------------------

void testFunction(std::string csGridName, std::string gaussGridName) {
    const Grid grid               = Grid(csGridName);
    const util::Config meshConfig = util::Config("partitioner", "cubedsphere") | util::Config("halo", 1);
    const MeshGenerator meshGen   = MeshGenerator("cubedsphere", meshConfig);
    const Mesh mesh               = meshGen.generate(grid);
    const functionspace::NodeColumns CSFunctionSpace = functionspace::NodeColumns(mesh);

    const StructuredGrid gaussGrid(gaussGridName);

    const grid::MatchingMeshPartitioner CSPartitioner(CSFunctionSpace.mesh(), atlas::option::type("cubedsphere"));
    const grid::Distribution gaussToCSDistribution(gaussGrid, CSPartitioner);
    const Mesh structuredGaussMesh = atlas::StructuredMeshGenerator().generate(gaussGrid, gaussToCSDistribution); // Problem occurs here on non-rubix cube ranks.
}

CASE("CS12 O12") {
    testFunction("CS-LFR-12", "O12");
}

CASE("CS12 F12") {
    testFunction("CS-LFR-12", "F12");
}

CASE("CS15 O15") {
    testFunction("CS-LFR-15", "O15");
}

CASE("CS15 F15") {
    testFunction("CS-LFR-15", "O15");
}

CASE("CS12 O8") {
    testFunction("CS-LFR-12", "O8");
}

CASE("CS12 F8") {
    testFunction("CS-LFR-12", "F8");
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
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

CASE("fails on non rubix cube") {
    const Grid grid = Grid("CS-LFR-15");
    const util::Config meshConfig =
        util::Config("partitioner", "cubedsphere") | util::Config("halo", 1);
    const MeshGenerator meshGen   = MeshGenerator("cubedsphere", meshConfig);
    const Mesh mesh = meshGen.generate(grid);
    const functionspace::NodeColumns CSFunctionSpace = functionspace::NodeColumns(mesh);

    const StructuredGrid gaussGrid("F15");

    atlas::StructuredMeshGenerator meshGenerator;
    const grid::MatchingMeshPartitioner CSPartitioner =
        grid::MatchingMeshPartitioner(CSFunctionSpace.mesh(), atlas::option::type("cubedsphere"));
    const grid::Distribution gaussToCSDistribution = CSPartitioner.partition(gaussGrid);
    const Mesh structuredGaussMesh = meshGenerator.generate(gaussGrid, gaussToCSDistribution);
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
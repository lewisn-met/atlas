/*
 * (C) Crown Copyright 2026 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "eckit/filesystem/PathName.h"

#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Config.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

CASE("test_gmsh_custom_empty_partition_shards") {
    EXPECT_EQ(mpi::size(), 2);

    // All grid points are assigned to partition 0, leaving rank 1 empty.
    Grid grid("F4");
    std::vector<int> partitions(grid.size(), 0);
    grid::Distribution distribution(mpi::size(), grid.size(), partitions.data());

    // Match OOPS' custom mapping: empty ranks get partition number comm.size() + 1,
    // with only nonempty partitions counted in nb_parts.
    const int part = mpi::rank() == 0 ? 0 : mpi::size() + 1;
    util::Config mesh_config(grid.meshgenerator());
    mesh_config.set("part", part);
    mesh_config.set("nb_parts", 1);
    mesh_config.set("mpi_comm", mpi::comm().name());
    
    Mesh mesh = StructuredMeshGenerator(mesh_config).generate(grid, distribution);
    functionspace::StructuredColumns functionspace(grid, distribution, option::halo(0));
    Field field = functionspace.createField<double>(option::name("rank"));

    auto values = array::make_view<double, 1>(field);
    for (idx_t node = 0; node < field.shape(0); ++node) {
        values(node) = mpi::rank();
    }

    // Check the mesh and function space are populated only on the nonempty rank.
    EXPECT_EQ(mesh.metadata().get<int>("part"), part);
    EXPECT_EQ(mesh.metadata().get<int>("nb_parts"), 1);
    if (mpi::rank() == 0) {
        EXPECT(mesh.cells().size() > 0);
        EXPECT_EQ(functionspace.sizeOwned(), grid.size());
    }
    else {
        EXPECT_EQ(mesh.nodes().size(), 0);
        EXPECT_EQ(mesh.cells().size(), 0);
        EXPECT_EQ(functionspace.sizeOwned(), 0);
    }
    mpi::comm().barrier();

    // Write one shard per rank (no gather).
    const std::string filename = "test_gmsh_custom_empty_partition.msh";
    output::Gmsh gmsh(filename, util::Config("gather", false));
    gmsh.write(mesh);
    gmsh.write(field);
    mpi::comm().barrier();

    // Verify each shard contains every expected Gmsh section
    // Note this checks all shards in the output directory, make sure to empty this between runs.
    // (all shards should contain both mesh and field data).
    if (mpi::rank() == 0) {
        std::vector<eckit::PathName> shards;
        eckit::PathName::match(filename + ".p*", shards);
        EXPECT(!shards.empty());

        for (const auto& shard : shards) {
            Log::info() << "Checking mesh and field data in " << shard << std::endl;
            std::ifstream file(shard.localPath());

            // Construct a string `problems` which stores the output errors
            // We expect none, so EXPECT_EQ compares to the empty string
            // and prints `problems` if any are found.
            std::string problems;
            if (!file.is_open()) {
                problems += " [" + shard.asString() + ": cannot open file]";
            }
            std::ostringstream contents;
            contents << file.rdbuf();
            const auto data = contents.str();
            for (const std::string section : {"MeshFormat", "Nodes", "Elements", "NodeData"}) {
                if (data.find("$" + section + "\n") == std::string::npos) {
                    problems += " [" + shard.asString() + ": missing $" + section + "]";
                }
                if (data.find("$End" + section + "\n") == std::string::npos) {
                    problems += " [" + shard.asString() + ": missing $End" + section + "]";
                }
            }
            EXPECT_EQ(problems, std::string());
        }
    }
}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}

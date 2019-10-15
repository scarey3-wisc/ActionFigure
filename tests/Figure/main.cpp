#include "FiniteElementMesh.h"

#include <Eigen/Dense>

template<class T>
struct LatticeMesh : public FiniteElementMesh<T>
{
	using Base = FiniteElementMesh<T>;

	// from AnimatedMesh
	using Base::m_meshElements;
	using Base::m_particleX;
	using Base::initializeUSD;
	using Base::initializeTopology;
	using Base::initializeParticles;
	using Vector3 = typename Base::Vector3;

	// from FiniteElementMesh
	using Base::m_particleV;
	using Base::initializeUndeformedConfiguration;
	using Base::m_stepEndTime;

	std::array<int, 3> m_cellSize; // dimensions in grid cells
	T m_gridDX;

	std::vector<Vector3> m_particleUndeformedX;
	std::vector<int> righHand;

	std::vector<int> leftHand;

	std::vector<std::array<int, 3>> m_activeCells; // Marks the "active" cells in the lattice
	std::map<std::array<int, 3>, int> m_activeNodes; // Maps the "active" nodes to their particle index

	LatticeMesh()
		:Base(1.e3, 5., 1., .03)

	{
		
	}

	bool withinEllipsoid(std::array<int, 3> test, std::array<int, 3> cent, float xRad, float yRad, float zRad) {
		float r = 0;
		if (xRad != 0)
			r += 1.0 * (test[0] - cent[0]) * (test[0] - cent[0]) / (xRad * xRad);
		if(yRad != 0)
			r += 1.0 * (test[1] - cent[1]) * (test[1] - cent[1]) / (yRad * yRad);
		if(zRad != 0)
			r += 1.0 * (test[2] - cent[2]) * (test[2] - cent[2]) / (zRad * zRad);
		return r <= 1;
	}
	void initialize()
	{
		initializeUSD("ActionFigure.usda");


		for (int cell_i = 0; cell_i < m_cellSize[0]; cell_i++)
			for (int cell_j = 0; cell_j < m_cellSize[1]; cell_j++)
				for (int cell_k = 0; cell_k < m_cellSize[2]; cell_k++) {

					bool inside = false;

					if (withinEllipsoid({ cell_i, cell_j, cell_k }, { 20, 20, 10 }, 8, 12, 4) &&
						cell_i >= 14 && cell_i <= 26 && cell_j >= 12 && cell_j <= 30 && cell_k >= 7 && cell_k <= 13)
						inside = true;
					if (withinEllipsoid({ cell_i, cell_j, cell_k }, { 17, 6, 10 }, 2.5, 0, 2.5) &&
						cell_j <= 14)
						inside = true;
					if (withinEllipsoid({ cell_i, cell_j, cell_k }, { 23, 6, 10 }, 2.5, 0, 2.5) &&
						cell_j <= 14)
						inside = true;
					if (withinEllipsoid({ cell_i, cell_j, cell_k }, { 20, 34, 10 }, 4.5, 4.5, 4.5))
						inside = true;
					if (withinEllipsoid({ cell_i, cell_j, cell_k }, { 13, 27, 10 }, 2.5, 2.5, 0) &&
						cell_k >= 10 && cell_k <= 18)
						inside = true;
					if (withinEllipsoid({ cell_i, cell_j, cell_k }, { 27, 27, 10 }, 2.5, 2.5, 0) &&
						cell_k >= 10 && cell_k <= 18)
						inside = true;
					if (withinEllipsoid({ cell_i, cell_j, cell_k }, { 13, 27, 9 }, 2.5, 2.5, 2.5))
						inside = true;
					if (withinEllipsoid({ cell_i, cell_j, cell_k }, { 27, 27, 9 }, 2.5, 2.5, 2.5))
						inside = true;
					if (withinEllipsoid({ cell_i, cell_j, cell_k }, { 13, 27, 21 }, 2.5, 2.5, 2.5))
						inside = true;
					if (withinEllipsoid({ cell_i, cell_j, cell_k }, { 27, 27, 21 }, 2.5, 2.5, 2.5))
						inside = true;

					if(inside)
						m_activeCells.push_back(std::array<int, 3>{cell_i, cell_j, cell_k});

				}

		std::cout << "Created a model including " << m_activeCells.size() << " lattice cells" << std::endl;

		// Create (uniquely numbered) particles at the node corners of active cells

		for (const auto& cell : m_activeCells) {
			std::array<int, 3> node;
			for (node[0] = cell[0]; node[0] <= cell[0] + 1; node[0]++)
				for (node[1] = cell[1]; node[1] <= cell[1] + 1; node[1]++)
					for (node[2] = cell[2]; node[2] <= cell[2] + 1; node[2]++) {
						auto search = m_activeNodes.find(node);
						if (search == m_activeNodes.end()) { // Particle not yet created at this lattice node location -> make one
							m_activeNodes.insert({ node, m_particleX.size() });
							
							m_particleX.emplace_back(m_gridDX * T(node[0]), m_gridDX * T(node[1]), m_gridDX * T(node[2]));

						}
					}
		}
		std::cout << "Model contains " << m_particleX.size() << " particles" << std::endl;

		//Find Particle Handles
		for (int e = 0; e < m_particleX.size(); e++) {
			Vector3 pos = m_particleX[e];
			Vector3 rhp = m_gridDX * Vector3(T(13 + 0.5), T(27 + 0.5), T(21 + 0.5));
			if ((pos - rhp).norm() < 3.5 * m_gridDX && pos(2) > m_gridDX * 23) {
				righHand.emplace_back(e);
			}
			if ((pos - rhp).norm() < 3.5 * m_gridDX && pos(0) > m_gridDX * 15) {
				righHand.emplace_back(e);
			}
			if ((pos - rhp).norm() < 3.5 * m_gridDX && pos(0) < m_gridDX * 12) {
				righHand.emplace_back(e);
			}
			if ((pos - rhp).norm() < 3.5 * m_gridDX && pos(1) > m_gridDX * 29) {
				righHand.emplace_back(e);
			}
			if ((pos - rhp).norm() < 3.5 * m_gridDX && pos(1) < m_gridDX * 26) {
				righHand.emplace_back(e);
			}
			Vector3 lhp = m_gridDX * Vector3(T(27 + 0.5), T(27 + 0.5), T(21 + 0.5));
			if ((pos - lhp).norm() < 3.5 * m_gridDX && pos(2) > m_gridDX * 23) {
				leftHand.emplace_back(e);
			}
			if ((pos - lhp).norm() < 3.5 * m_gridDX && pos(0) > m_gridDX * 29) {
				leftHand.emplace_back(e);
			}
			if ((pos - lhp).norm() < 3.5 * m_gridDX && pos(0) < m_gridDX * 26) {
				leftHand.emplace_back(e);
			}
			if ((pos - lhp).norm() < 3.5 * m_gridDX && pos(1) > m_gridDX * 29) {
				leftHand.emplace_back(e);
			}
			if ((pos - lhp).norm() < 3.5 * m_gridDX && pos(1) < m_gridDX * 26) {
				leftHand.emplace_back(e);
			}

		}


		// Make tetrahedra out of all active cells (6 tetrahedra per cell)

		for (const auto& cell : m_activeCells) {
			int vertexIndices[2][2][2];
			for (int i = 0; i <= 1; i++)
				for (int j = 0; j <= 1; j++)
					for (int k = 0; k <= 1; k++) {
						std::array<int, 3> node{ cell[0] + i, cell[1] + j, cell[2] + k };
						auto search = m_activeNodes.find(node);
						if (search != m_activeNodes.end())
							vertexIndices[i][j][k] = search->second;
						else
							throw std::logic_error("particle at cell vertex not found");
					}

			m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][0][0], vertexIndices[1][1][0], vertexIndices[1][1][1]});
			m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][0][0], vertexIndices[1][1][1], vertexIndices[1][0][1]});
			m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][0][1], vertexIndices[1][1][1], vertexIndices[0][0][1]});
			m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][1][1], vertexIndices[0][1][1], vertexIndices[0][0][1]});
			m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][1][1], vertexIndices[0][1][0], vertexIndices[0][1][1]});
			m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][1][0], vertexIndices[0][1][0], vertexIndices[1][1][1]});
		}

		// Perform the USD-specific initialization of topology & particles
		// (this will also create a boundary *surface* to visualuze

		initializeTopology();
		initializeParticles();

		// Check particle indexing in mesh

		for (const auto& element : m_meshElements)
			for (const auto vertex : element)
				if (vertex < 0 || vertex >= m_particleX.size())
					throw std::logic_error("mismatch between mesh vertex and particle array");

		// Also resize the velocities to match
		m_particleV.resize(m_particleX.size(), Vector3::Zero());

		// Initialize rest shape matrices and particle mass
		initializeUndeformedConfiguration();
		m_particleUndeformedX = m_particleX;

	}
	void clearConstrainedParticles(std::vector<Vector3>& x) override
	{
		for (const auto v : righHand)
			x[v] = Vector3::Zero();
		for (const auto v : leftHand)
			x[v] = Vector3::Zero();
	}

	void setBoundaryConditions() override
	{
		if (m_stepEndTime <= 5) {
			T curTime = m_stepEndTime;
			for (const auto v : righHand) {
				m_particleX[v] = m_particleUndeformedX[v] + curTime * Vector3(T(0), T(0), T(0.02));
			}
			for (const auto v : leftHand) {
				m_particleX[v] = m_particleUndeformedX[v];
			}
		}
		else if (m_stepEndTime <= 8) {
			for (const auto v : righHand) {
				m_particleX[v] = m_particleUndeformedX[v] + 5 * Vector3(T(0), T(0), T(0.02));
			}
			for (const auto v : leftHand) {
				m_particleX[v] = m_particleUndeformedX[v];
			}
		}
		else if (m_stepEndTime <= 13) {
			T curTime = m_stepEndTime - 8;
			for (const auto v : righHand) {
				m_particleX[v] = m_particleUndeformedX[v] + (5 - curTime) * Vector3(T(0), T(0), T(0.02));
			}
			for (const auto v : leftHand) {
				m_particleX[v] = m_particleUndeformedX[v];
			}
		}
		else if (m_stepEndTime <= 15) {
			for (const auto v : righHand) {
				m_particleX[v] = m_particleUndeformedX[v];
			}
			for (const auto v : leftHand) {
				m_particleX[v] = m_particleUndeformedX[v];
			}
		}
		else if (m_stepEndTime <= 20) {
			T curTime = m_stepEndTime - 15;
			for (const auto v : righHand) {
				m_particleX[v] = m_particleUndeformedX[v] + curTime * Vector3(T(0), T(0.03), T(0));
			}
			for (const auto v : leftHand) {
				m_particleX[v] = m_particleUndeformedX[v] + curTime * Vector3(T(0), T(-0.03), T(0));
			}
		}
		else if (m_stepEndTime <= 25) {
			for (const auto v : righHand) {
				m_particleX[v] = m_particleUndeformedX[v] + 5 * Vector3(T(0), T(0.03), T(0));
			}
			for (const auto v : leftHand) {
				m_particleX[v] = m_particleUndeformedX[v] + 5 * Vector3(T(0), T(-0.03), T(0));
			}
		}
		else if (m_stepEndTime <= 35) {
			T curTime = m_stepEndTime - 25;
			for (const auto v : righHand) {
				m_particleX[v] = m_particleUndeformedX[v] + (5 - curTime) * Vector3(T(0), T(0.03), T(0));
			}
			for (const auto v : leftHand) {
				m_particleX[v] = m_particleUndeformedX[v] + (5 - curTime) * Vector3(T(0), T(-0.03), T(0));
			}
		}
		else if (m_stepEndTime <= 40) {
			for (const auto v : righHand) {
				m_particleX[v] = m_particleUndeformedX[v] -5 * Vector3(T(0), T(0.03), T(0));
			}
			for (const auto v : leftHand) {
				m_particleX[v] = m_particleUndeformedX[v] -5 * Vector3(T(0), T(-0.03), T(0));
			}
		}
		else if (m_stepEndTime <= 45) {
			T curTime = m_stepEndTime - 40;
			for (const auto v : righHand) {
				m_particleX[v] = m_particleUndeformedX[v] + (curTime - 5) * Vector3(T(0), T(0.03), T(0));
			}
			for (const auto v : leftHand) {
				m_particleX[v] = m_particleUndeformedX[v] + (curTime - 5) * Vector3(T(0), T(-0.03), T(0));
			}
		}
		else if (m_stepEndTime <= 55) {
			for (const auto v : righHand) {
				m_particleX[v] = m_particleUndeformedX[v];
			}
			for (const auto v : leftHand) {
				m_particleX[v] = m_particleUndeformedX[v];
			}
		}
		else if (m_stepEndTime <= 60) {
			T curTime = m_stepEndTime - 55;
			for (const auto v : righHand) {
				m_particleX[v] = m_particleUndeformedX[v] + curTime * Vector3(T(0.02), T(0), T(0));
			}
			for (const auto v : leftHand) {
				m_particleX[v] = m_particleUndeformedX[v] + curTime * Vector3(T(-0.02), T(0), T(0));
			}
		}
		else if (m_stepEndTime <= 65) {
			for (const auto v : righHand) {
				m_particleX[v] = m_particleUndeformedX[v] + 5 * Vector3(T(0.02), T(0), T(0));
			}
			for (const auto v : leftHand) {
				m_particleX[v] = m_particleUndeformedX[v] + 5 * Vector3(T(-0.02), T(0), T(0));
			}
		}
		else if (m_stepEndTime <= 70) {
			T curTime = m_stepEndTime - 65;
			for (const auto v : righHand) {
				m_particleX[v] = m_particleUndeformedX[v] + (5 - curTime) * Vector3(T(0.02), T(0), T(0));
			}
			for (const auto v : leftHand) {
				m_particleX[v] = m_particleUndeformedX[v] + (5 - curTime) * Vector3(T(-0.02), T(0), T(0));
			}
		}
		else {
			for (const auto v : righHand) {
				m_particleX[v] = m_particleUndeformedX[v];
			}
			for (const auto v : leftHand) {
				m_particleX[v] = m_particleUndeformedX[v];
			}
		}
	}
};

int main(int argc, char* argv[])
{
	LatticeMesh<float> simulationMesh;
	simulationMesh.m_cellSize = { 40, 40, 40 };
	simulationMesh.m_gridDX = 0.025;
	simulationMesh.m_nFrames = 200;
	simulationMesh.m_subSteps = 1;
	simulationMesh.m_frameDt = 0.5;

	// Initialize the simulation example
	simulationMesh.initialize();

	// Output the initial shape of the mesh
	simulationMesh.writeFrame(0);

	// Perform the animation, output results at each frame
	for (int frame = 1; frame <= simulationMesh.m_nFrames; frame++) {
		simulationMesh.simulateFrame(frame);
		simulationMesh.writeFrame(frame);
	}

	// Write the entire timeline to USD
	simulationMesh.writeUSD();

	return 0;
}
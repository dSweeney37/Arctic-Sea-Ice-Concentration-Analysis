/*
 Author: Daniel T. Sweeney
 Class: CSE310 - 83657
 Date: 12//2019
 Description:
*/
#include "Project3.h"





int main() {
	cell** cellMatrix = new cell * [CELL_DIMENSIONS_X];
	float*** sicArray = new float** [CELL_DIMENSIONS_X];
	int** componentMatrix = new int* [3]{ 0 };
	int** degreeMatrix = new int* [3]{ 0 };
	short*** distanceGraph = new short** [CELL_DIMENSIONS_X * CELL_DIMENSIONS_Y];



	for (int x = 0; x < CELL_DIMENSIONS_X; x++) {
		cellMatrix[x] = new cell[CELL_DIMENSIONS_Y];
		sicArray[x] = new float* [CELL_DIMENSIONS_Y];
		
		if (x < 3) {
			componentMatrix[x] = new int[CELL_DIMENSIONS_X * CELL_DIMENSIONS_Y]();
			degreeMatrix[x] = new int[CELL_DIMENSIONS_X * CELL_DIMENSIONS_Y]();
		}

		for (int y = 0; y < CELL_DIMENSIONS_Y; y++) { 
			sicArray[x][y] = new float[(MAX_YEAR - MIN_YEAR + 1) * NUM_OF_WEEKS];
			distanceGraph[(x * CELL_DIMENSIONS_X) + y] = new short* [CELL_DIMENSIONS_X * CELL_DIMENSIONS_Y];
			

			for (int z = 0; z < CELL_DIMENSIONS_X * CELL_DIMENSIONS_Y; z++) {
				distanceGraph[(x * CELL_DIMENSIONS_X) + y][z] = new short[3];
				fill(&distanceGraph[(x * CELL_DIMENSIONS_X) + y][z][0], &distanceGraph[(x * CELL_DIMENSIONS_X) + y][z][3], INF);
			}
			fill(&distanceGraph[y][y][0], &distanceGraph[y][y][3], 0);
		}
	}

	LoadData(sicArray, cellMatrix);
	LoadEdges(sicArray, cellMatrix, degreeMatrix, distanceGraph);
	DFS(cellMatrix, componentMatrix);
	PrintDegreeData(degreeMatrix);
	PrintComponentData(componentMatrix);
	CalculateClusteringCoefficients(cellMatrix);
	FloydWarshall(distanceGraph);

	delete[] cellMatrix;
	delete[] componentMatrix;
	delete[] degreeMatrix;
	delete[] sicArray;
	delete[] distanceGraph;
	return 0;
}





void CalculateClusteringCoefficients(cell** pCellMatrix) {
	float clusterCoefficient[3] = { 0 };
	int meanDivisor = 0;



	for (int y = 0; y < CELL_DIMENSIONS_X; y++) {
		for (int x = 0; x < CELL_DIMENSIONS_X; x++) {
			if (pCellMatrix[x][y].isSIC == true) {
				int numOfEdges[3] = { 0 }, numOfCells[3] = { 0 };


				numOfCells[0] = pCellMatrix[x][y].edgeNum[0] * (pCellMatrix[x][y].edgeNum[0] - 1);
				numOfCells[1] = pCellMatrix[x][y].edgeNum[1] * (pCellMatrix[x][y].edgeNum[1] - 1);
				numOfCells[2] = pCellMatrix[x][y].edgeNum[2] * (pCellMatrix[x][y].edgeNum[2] - 1);

				for (int j = 0; j < CELL_DIMENSIONS_Y; j++) {
					for (int i = 0; i < CELL_DIMENSIONS_X; i++) {
						if (pCellMatrix[x][y].adjArr[i][j][0] == true) {
							numOfEdges[0] += CountNeighboursEdges(pCellMatrix, pCellMatrix[x][y], i, j, 0);
						}
						if (pCellMatrix[x][y].adjArr[i][j][1] == true) {
							numOfEdges[1] += CountNeighboursEdges(pCellMatrix, pCellMatrix[x][y], i, j, 1);
						}
						if (pCellMatrix[x][y].adjArr[i][j][2] == true) {
							numOfEdges[2] += CountNeighboursEdges(pCellMatrix, pCellMatrix[x][y], i, j, 2);
						}
					}
				}

				for (int i = 0; i < 3; i++) {
					if (numOfCells[i] != 0) { clusterCoefficient[i] += (float)(2 * numOfEdges[i]) / numOfCells[i]; }
				}

				meanDivisor++;
			}
		}
	}
	cout << "-------------Mean Clustering Coefficients-------------\n";
	for (int i = 0; i < 3; i++) {
		clusterCoefficient[i] /= meanDivisor;
		cout << " Threshold of " << arrThreshold[i] << ": " << clusterCoefficient[i] * 100 << "%" << endl;
	}
	cout << "\n\n\n";
}





void CalculateCorrelation(float*** pSICArray, cell** pCellMatrix, short*** pDistanceGraph, int x1, int y1, int x2, int y2) {
	bool calcSxx = false, calcSyy = false;
	float correlation, Sxy = 0;
	


	if (pCellMatrix[x1][y1].squareSum == 0) { calcSxx = true; }
	if (pCellMatrix[x2][y2].squareSum == 0) { calcSyy = true; }

	for (int z = 0; z < (MAX_YEAR - MIN_YEAR + 1) * NUM_OF_WEEKS; z++) {
		if (calcSxx == true) { pCellMatrix[x1][y1].squareSum += pow(pSICArray[x1][y1][z] - pCellMatrix[x1][y1].bar, 2); }
		if (calcSyy == true) { pCellMatrix[x2][y2].squareSum += pow(pSICArray[x2][y2][z] - pCellMatrix[x2][y2].bar, 2); }
		Sxy += (pSICArray[x1][y1][z] - pCellMatrix[x1][y1].bar) * (pSICArray[x2][y2][z] - pCellMatrix[x2][y2].bar);
	}

	correlation = abs(Sxy / sqrt(pCellMatrix[x1][y1].squareSum * pCellMatrix[x2][y2].squareSum));
	
	if (correlation >= stof(arrThreshold[2])) {
		pCellMatrix[x1][y1].edgeNum[2]++;
		pCellMatrix[x1][y1].adjArr[x2][y2][2] = true;
		pCellMatrix[x2][y2].edgeNum[2]++;
		pCellMatrix[x2][y2].adjArr[x1][y1][2] = true;
		pDistanceGraph[(y1 * CELL_DIMENSIONS_X) + x1][(y2 * CELL_DIMENSIONS_X) + x2][2] = 1;
		pDistanceGraph[(y2 * CELL_DIMENSIONS_X) + x2][(y1 * CELL_DIMENSIONS_X) + x1][2] = 1;
		
		if (correlation >= stof(arrThreshold[1])) {
			pCellMatrix[x1][y1].edgeNum[1]++;
			pCellMatrix[x1][y1].adjArr[x2][y2][1] = true;
			pCellMatrix[x2][y2].edgeNum[1]++;
			pCellMatrix[x2][y2].adjArr[x1][y1][1] = true;
			pDistanceGraph[(y1 * CELL_DIMENSIONS_X) + x1][(y2 * CELL_DIMENSIONS_X) + x2][1] = 1;
			pDistanceGraph[(y2 * CELL_DIMENSIONS_X) + x2][(y1 * CELL_DIMENSIONS_X) + x1][1] = 1;
			
			if (correlation >= stof(arrThreshold[0])) {
				pCellMatrix[x1][y1].edgeNum[0]++;
				pCellMatrix[x1][y1].adjArr[x2][y2][0] = true;
				pCellMatrix[x2][y2].edgeNum[0]++;
				pCellMatrix[x2][y2].adjArr[x1][y1][0] = true;
				pDistanceGraph[(y1 * CELL_DIMENSIONS_X) + x1][(y2 * CELL_DIMENSIONS_X) + x2][0] = 1;
				pDistanceGraph[(y2 * CELL_DIMENSIONS_X) + x2][(y1 * CELL_DIMENSIONS_X) + x1][0] = 1;
			}
		}
	}
}




  `	1Q37
int CountNeighboursEdges(cell** pCellMatrix, cell pCell, int X, int Y, int pSlot) {
	int x = (X + 1 + (Y * CELL_DIMENSIONS_X)) % CELL_DIMENSIONS_X;
	int y = floor((X + 1 + (Y * CELL_DIMENSIONS_X)) / CELL_DIMENSIONS_X);
	int count = 0;

	

	while (y < CELL_DIMENSIONS_Y) {
		while (x < CELL_DIMENSIONS_X) {
			if (pCell.adjArr[x][y][pSlot] == true) { if (pCellMatrix[X][Y].adjArr[x][y][pSlot] == true) { count++; } }
			++x;
		}
		x = 0;
		++y;
	}
	return count;
}





void DFS(cell** pCellMatrix, int** pComponentMatrix) {
	for (int y = 0; y < CELL_DIMENSIONS_Y; y++) {
		for (int x = 0; x < CELL_DIMENSIONS_X; x++) {
			int arrCount[3];
			fill(begin(arrCount), end(arrCount), -1);
			


			if (pCellMatrix[x][y].isSIC == true) {
				if (pCellMatrix[x][y].explored[0] == false) { arrCount[0] = DFSVisit(pCellMatrix, x, y, 1, 0); }
				if (pCellMatrix[x][y].explored[1] == false) { arrCount[1] = DFSVisit(pCellMatrix, x, y, 1, 1); }
				if (pCellMatrix[x][y].explored[2] == false) { arrCount[2] = DFSVisit(pCellMatrix, x, y, 1, 2); }

				if (arrCount[0] > -1) { pComponentMatrix[0][arrCount[0]]++; }
				if (arrCount[1] > -1) { pComponentMatrix[1][arrCount[1]]++; }
				if (arrCount[2] > -1) { pComponentMatrix[2][arrCount[2]]++; }
			}
		}	
	}
}





int DFSVisit(cell** pCellMatrix, int X, int Y, int pCount, int pSlot) {
	pCellMatrix[X][Y].explored[pSlot] = true;

	for (int y = 0; y < CELL_DIMENSIONS_Y; y++) {
		for (int x = 0; x < CELL_DIMENSIONS_X; x++) {
			if (pCellMatrix[X][Y].adjArr[x][y][pSlot] == true && pCellMatrix[x][y].explored[pSlot] == false) {
				pCount = DFSVisit(pCellMatrix, x, y, ++pCount, pSlot);
			}
		}
	}
	return pCount;
}





void FloydWarshall(short*** pDistanceGraph) {
	float tempArr[3][2] = { { 0 } };



	for (int i = 0; i < CELL_DIMENSIONS_X * CELL_DIMENSIONS_Y; i++) {
		for (int s = 0; s < CELL_DIMENSIONS_X * CELL_DIMENSIONS_Y; s++) {
			for (int d = s + 1; d < CELL_DIMENSIONS_X * CELL_DIMENSIONS_Y; d++) {

				for (int z = 0; z < 3; z++) {
					if (pDistanceGraph[s][i][z] + pDistanceGraph[i][d][z] < pDistanceGraph[s][d][z]) {
						pDistanceGraph[s][d][z] = pDistanceGraph[s][i][z] + pDistanceGraph[i][d][z];
					}
					if (i == CELL_DIMENSIONS_X * CELL_DIMENSIONS_Y - 1 && pDistanceGraph[s][d][z] != INF) {
						tempArr[z][0] += pDistanceGraph[s][d][z];
						++tempArr[z][1];
					}
				}
			}
		}
	}
	for (int j = 0; j < 3; j++) {
		cout << tempArr[j][0]/tempArr[j][1] << ", " << tempArr[j][0] << ", " << tempArr[j][1] <<  endl;
	}
}





void LoadData(float*** pSICArray, cell** pCellMatrix) {
	string fileName = "CS310_project_subregion/####/Beaufort_Sea_diffw##y####+landmask";
	


	for (int z = 0; z < (MAX_YEAR - MIN_YEAR + 1) * NUM_OF_WEEKS; z++) {
		string week = to_string(z % NUM_OF_WEEKS + 1);


		if (z % NUM_OF_WEEKS < 9) { week = "0" + week; }
		fileName.replace(47, 2, week); // Replaces the w## part of 'fileName' w/the file year.
		fileName.replace(24, 4, to_string((int)(MIN_YEAR + floor(z / NUM_OF_WEEKS)))); // Replaces the y#### part of 'fileName' w/the file year.
		fileName.replace(50, 4, to_string((int)(MIN_YEAR + floor(z / NUM_OF_WEEKS)))); // Replaces the y#### part of 'fileName' w/the file year.
		ifstream inFile(fileName, ios::in | ios::binary);

		for (int y = 0; y < CELL_DIMENSIONS_Y; y++) {
			for (int x = 0; x < CELL_DIMENSIONS_X; x++) {
				float dataIn = 0;


				inFile.read((char*)&dataIn, 4); // Read 4 bytes of data.
				pSICArray[x][y][z] = dataIn; // Loads each cell value into its corresponding pSlot in the 3D array.

				if (dataIn == 168) { pCellMatrix[x][y].isSIC = false; }
				else if (dataIn != 157) {
					pCellMatrix[x][y].bar += dataIn;
					pCellMatrix[x][y].counter++;
				}

				if (z == ((MAX_YEAR - MIN_YEAR + 1) * NUM_OF_WEEKS) - 1 && pCellMatrix[x][y].counter != 0) {
					pCellMatrix[x][y].bar /= pCellMatrix[x][y].counter;
				}
			}
		}
		inFile.close();
	}
}





void LoadEdges(float*** pSICArray, cell** pCellMatrix, int** pDegreeMatrix, short*** pDistanceGraph) {
	int pos = 0, x1, y1, x2, y2;
	


	while (pos < CELL_DIMENSIONS_X * CELL_DIMENSIONS_Y) {
		x1 = pos % CELL_DIMENSIONS_X;
		y1 = floor(pos / CELL_DIMENSIONS_X);
		
		if (pCellMatrix[x1][y1].isSIC == false) { ++pos; }
		
		else if (pos == CELL_DIMENSIONS_X * CELL_DIMENSIONS_Y - 1) {
			pDegreeMatrix[0][pCellMatrix[x1][y1].edgeNum[0]]++;
			pDegreeMatrix[1][pCellMatrix[x1][y1].edgeNum[1]]++;
			pDegreeMatrix[2][pCellMatrix[x1][y1].edgeNum[2]]++;
			++pos;
		}

		else {
			x2 = (pos + 1) % CELL_DIMENSIONS_X;
			y2 = floor((pos + 1) / CELL_DIMENSIONS_X);

			while (y2 < CELL_DIMENSIONS_Y) {
				while (x2 < CELL_DIMENSIONS_X) {
					if (pCellMatrix[x2][y2].isSIC == true && pCellMatrix[x2][y2].explored[0] == false) {
						CalculateCorrelation(pSICArray, pCellMatrix, pDistanceGraph, x1, y1, x2, y2); 
					}
					++x2;
				}
				x2 = 0;
				++y2;
			}

			pDegreeMatrix[0][pCellMatrix[x1][y1].edgeNum[0]]++;
			pDegreeMatrix[1][pCellMatrix[x1][y1].edgeNum[1]]++;
			pDegreeMatrix[2][pCellMatrix[x1][y1].edgeNum[2]]++;
			++pos;
		}
	}
}





void PrintComponentData(int** pComponentMatrix) {
	for (int i = 0; i < 3; i++) {
		int count = 0;


		cout << "-------------Connected Components for " << arrThreshold[i] << " Threshold-------------\n";

		for (int j = 0; j < CELL_DIMENSIONS_X * CELL_DIMENSIONS_Y; j++) {
			if (pComponentMatrix[i][j] != 0) {
				cout << " Connected Cells of Length " << j << ": " << pComponentMatrix[i][j] << endl;
				count += pComponentMatrix[i][j];
			}
		}
		cout << " Number of Components: " << count << endl << "\n\n\n";
	}
}





void PrintDegreeData(int** pDegreeMatrix) {
	for (int i = 0; i < 3; i++) {
		cout << "-------------Degree Distribution for " << arrThreshold[i] << " Threshold-------------\n";

		for (int j = 0; j < CELL_DIMENSIONS_X * CELL_DIMENSIONS_Y; j++) {
			if (pDegreeMatrix[i][j] != 0) {
				//string str(pDegreeMatrix[i][j], '*');
				//cout << " Degree " << j << ": " << str << endl;
				cout << " Degree " << j << ": " << pDegreeMatrix[i][j] << endl;
			}
		}
		cout << "\n\n\n";
	}
}
#define INF CELL_DIMENSIONS_X * CELL_DIMENSIONS_Y + 1
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
using namespace std;



const int CELL_DIMENSIONS_X = 63; // Large data set (304), Small data set (63)
const int CELL_DIMENSIONS_Y = 63; // Large data set (448), Small data set (63)
const int MIN_YEAR = 1990; // Large data set (1979 - MAX_YEAR), Small data set (1990 - MAX_YEAR)
const int MAX_YEAR = 2005;
const int NUM_OF_WEEKS = 52;
string arrThreshold[3]{ ".95", ".925", ".9" };



struct cell {
	//.95 values -> arr[x][y][0], .925 values -> arr[x][y][1], .9 values -> arr[x][y][2]
	bool adjArr[CELL_DIMENSIONS_X][CELL_DIMENSIONS_Y][3] = { { { false } } };
	bool explored[3] = { false };
	bool isSIC = true;
	float bar = 0;
	float squareSum = 0;
	int counter = 0;
	int edgeNum[3] = { 0 };
};





void CalculateClusteringCoefficients(cell** pCellMatrix);
void CalculateCorrelation(float*** pSICArray, cell** pCellMatrix, short*** pDistanceGraph, int x1, int y1, int x2, int y2);
int CountNeighboursEdges(cell** pCellMatrix, cell pCell, int X, int Y, int pSlot);
void DFS(cell** pCellMatrix, int** pComponentMatrix);
int DFSVisit(cell** pCellMatrix, int x, int y, int pCount, int pSlot);
void FloydWarshall(short*** pDistanceGraph);
void LoadData(float*** pSICArray, cell** pCellMatrix);
void LoadEdges(float*** pSICArray, cell** pCellMatrix, int** pDegreeMatrix, short*** pDistanceGraph);
void PrintComponentData(int** pComponentMatrix);
void PrintDegreeData(int** pDegreeMatrix);
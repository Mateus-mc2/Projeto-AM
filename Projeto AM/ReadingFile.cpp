#include "ReadingFile.h"

using namespace std;

//x,x,x,x,o,o,x,o,o,positive to [1,1,1,1,0,0,1,0,0]
int* lineToVector(string line){
	int count = 0;
	int vectorReturn[9];

	for (int i = 0; i < 18; i = i + 2)
	{
		if (line[i] == 'x'){
			vectorReturn[count] = 1;
		}
		else if (line[i] == 'o'){
			vectorReturn[count] = 0;
		}
		else{
			vectorReturn[count] = -1;
		}
		count++;
	}
	return vectorReturn;
}
//return matrix [num. lines][9]
int** matrixExamples()
{
	string strInput;
	int** matrizPrinc = new int*[9];
	int line = 0;

	ifstream inf("TestFile.txt");

	if (!inf)
	{
		cout << "//could not be opened for reading!\\" << endl;
		exit(1);
	}

	while (getline(inf, strInput))
	{
		matrizPrinc[line] = new int[9];
		for (int j = 0; j < 9; j++)
		{
			matrizPrinc[line][j] = lineToVector(strInput)[j];
			//cout << matrizPrinc[line][j] ;
		}
		//cout << "--------------------------" << endl;
		line++;
	}

	return matrizPrinc;
}

int** matrixDissimilarity(int** m1, int totalLines){
	int auxLine = 0;
	int** result = new int*[totalLines];
	int** m2 = m1;

	for (int i = 0; i < totalLines; i++)
	{
		result[i] = new int[totalLines];
		for (int j = 0; j < totalLines; j++)
		{
			result[i][j] = 0;
		}
	}

	while (auxLine < totalLines){
		for (int i = 0; i < totalLines; i++)
		{
			for (int j = 0; j < 9; j++)
			{
				if (m2[auxLine][j] != m1[i][j]){
					result[auxLine][i] += 1;
				}
			}
			cout << " " << result[auxLine][i];
		}
		cout << endl;
		auxLine++;
	}

	return result;
}
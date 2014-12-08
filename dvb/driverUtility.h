
#pragma once

#include <fstream>
#include <iostream>
using namespace std;

#define		WRITE_FILE_FOR_DRIVER	1

template <typename T>
void 	readArray(T* pArray, int nSize, char* strFileName)
{
	ifstream file(strFileName);

	if( file == NULL )
		cout << "failed to open: " << strFileName << endl;

	for (int i=0;i<nSize;i++)
		file >> pArray[i];

	file.close();
}

template <typename T>
void	writeArray(T* pArray, int nSize, char* strFileName)
{
	ofstream file(strFileName);

	if( file == NULL )
		cout << "failed to open: " << strFileName << endl;

	for (int i=0;i<nSize;i++)
		file << pArray[i];

	file.close();
}

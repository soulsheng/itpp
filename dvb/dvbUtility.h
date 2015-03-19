
#pragma once


#define		FILENAME_IT		"../data/random_3_6_16200.it"
#define		FILENAME_ALIST	"../data/dvbs2_r12.alist"

#define		COUNT_REPEAT_DEF	1	// repeat time 
#define		SIZE_PACKET		188

#define		N_BCH			31
#define		T_BCH			2
#define		K_BCH			21

#define		VAR_SIZE_CODE		16200
#define		CHECK_SIZE_CODE		8100//8073

#define		REMOVE_NOISE		0

#if REMOVE_NOISE
#define		EBNO			10
#else
#define		EBNO			2.2
#endif

#define		REMOVE_BCH			1

#define		MOD_TYPE_DEFAULT	2

#include <itpp/itcomm.h>
#include <iostream>
using namespace std;
using namespace itpp;

void convertBufferToVec( char* buffer, bvec& a );

void convertBufferToVec( double* buffer, vec& a );

void convertBufferToVec( double* buffer, cvec& a );

void convertVecToBuffer( char* buffer, bvec& a );

void convertVecToBuffer( double* buffer, vec& a );

void convertVecToBuffer( double* buffer, cvec& a );


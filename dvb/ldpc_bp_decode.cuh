
#pragma once
#include <itpp/itcomm.h>

class ldpc_gpu
{
protected:
bool syndrome_check_gpu( ) ;

void updateVariableNode_gpu( ) ;

void updateCheckNode_gpu( );

void initializeMVC_gpu( );

bool check_parity_cpu(char *LLR);

public:
int bp_decode(int *LLRin, int *LLRout,
	bool psc = true,			//!< check syndrom after each iteration
	int max_iters = 50 );		//!< Maximum number of iterations

	/*!
	   * LDPC����   *
		* \param		LLRin	�������룺������Ȼ��ֵ����
		* \param		bitout	�������������������Ϣλ�����ַ�����
		* \param		psc	    �������룺�Ƿ�ÿ�ε�����������żУ��
		* \param		max_iters �������룺����������
	*/
int bp_decode_once(int *LLRin, char *LLRout,
	bool psc = true,			//!< check syndrom after each iteration
	int max_iters = 50 );		//!< Maximum number of iterations
int bp_decode_once(itpp::vec& softbits, char *LLRout);		//!< Maximum number of iterations
int bp_decode_once(double* softbits, char *LLRout);		//!< Maximum number of iterations

	/*!
	   * LDPC��������ʼ��   *
		* \param	nvar 	�������룺�����ڵ���Ŀ�����볤��
		* \param	ncheck	�������룺У��ڵ���Ŀ����Ϣ����
		* \param	nmaxX1	�������룺У�����ÿ�С�1���ĸ��������ֵ
		* \param	nmaxX2  �������룺У�����ÿ�С�1���ĸ��������ֵ
		* \param	V 		�������룺У�����ÿ�С�1����������
		* \param	sumX1	�������룺У�����ÿ�С�1���ĸ���
		* \param	sumX2	�������룺У�����ÿ�С�1���ĸ���
		* \param	iind	�������룺У�����ÿ�С�1��������
		* \param	jind	�������룺У�����ÿ�С�1��������
		* \param	Dint1/2/3		�������룺ͬ������Ȼ��class LLR_calc_unit
		* \param	logexp_table	�������룺������Ȼ�Ȳ��ұ�
	*/
	bool	initialize( );

	~ldpc_gpu();

	int get_nvar() const { return nvar; }
	int get_ncheck() const { return ncheck; }
	int get_ninfo() const { return nvar - ncheck; }
	float get_rate();

private:
	bool	release();

private:
	int* d_synd ;

	int* d_sumX1 ;
	int* d_sumX2 ;
	int* d_mcv ;
	int* d_mvc ;
	int* d_iind ;
	int* d_jind ;
	int* d_V ;

	int* d_logexp_table ;
	
	//int *d_ml, *d_mr ;
	
	int* d_LLRin ;
	char* d_LLRout ;
	
	int* h_V, *h_sumX2;
	int* h_mcv, *h_mvc ;

private:
	int nvar, ncheck;
	int nmaxX1, nmaxX2; // max(sumX1) max(sumX2)
	short int Dint1, Dint2, Dint3;	//! Decoder (lookup-table) parameters
	//int max_cnd;	//! Maximum check node degree that the class can handle
	int QLLR_MAX;

	itpp::LDPC_Code ldpc;

};
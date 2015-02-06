#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "apsk.h"


void apsk16(char *apsk16_input, double *apsk16_output_real, double *apsk16_output_imaginary)
{
	double kl = sqrt(2.0f);
	double pl = sin((double)(PI / 12.0f));
	double ql = cos((double)(PI / 12.0f));

	double a = 3.15f;

	if(apsk16_input[0] == '0' && apsk16_input[1] == '0' && apsk16_input[2] == '0' && apsk16_input[3] == '0')
	{
		*apsk16_output_real = a / 2.0f * kl;
		*apsk16_output_imaginary = a / 2.0f * kl;
	}
	else if(apsk16_input[0] == '0' && apsk16_input[1] == '0' && apsk16_input[2] == '0' && apsk16_input[3] == '1')
	{
		*apsk16_output_real = a / 2.0f * kl;
		*apsk16_output_imaginary = -a / 2.0f * kl;
	}
	else if(apsk16_input[0] == '0' && apsk16_input[1] == '0' && apsk16_input[2] == '1' && apsk16_input[3] == '0')
	{
		*apsk16_output_real = -a / 2.0f * kl;
		*apsk16_output_imaginary = a / 2.0f * kl;
	}

	else if(apsk16_input[0] == '0' && apsk16_input[1] == '0' && apsk16_input[2] == '1' && apsk16_input[3] == '1')
	{
		*apsk16_output_real = -a / 2.0f * kl;
		*apsk16_output_imaginary = -a / 2.0f * kl;
	}
	else if(apsk16_input[0] == '0' && apsk16_input[1] == '1' && apsk16_input[2] == '0' && apsk16_input[3] == '0')
	{
		*apsk16_output_real = a * ql;
		*apsk16_output_imaginary = a * pl;
	}
	else if(apsk16_input[0] == '0' && apsk16_input[1] == '1' && apsk16_input[2] == '0' && apsk16_input[3] == '1')
	{
		*apsk16_output_real = a * ql;
		*apsk16_output_imaginary = -a * pl;
	}
	else if(apsk16_input[0] == '0' && apsk16_input[1] == '1' && apsk16_input[2] == '1' && apsk16_input[3] == '0')
	{
		*apsk16_output_real = -a * ql;
		*apsk16_output_imaginary = a * pl;
	}
	else if(apsk16_input[0] == '0' && apsk16_input[1] == '1' && apsk16_input[2] == '1' && apsk16_input[3] == '1')
	{
		*apsk16_output_real = -a * ql;
		*apsk16_output_imaginary = -a * pl;
	}
	
	else if(apsk16_input[0] == '1' && apsk16_input[1] == '0' && apsk16_input[2] == '0' && apsk16_input[3] == '0')
	{
		*apsk16_output_real = a * pl;
		*apsk16_output_imaginary = a * ql;
	}
	else if(apsk16_input[0] == '1' && apsk16_input[1] == '0' && apsk16_input[2] == '0' && apsk16_input[3] == '1')
	{
		*apsk16_output_real = a * pl;
		*apsk16_output_imaginary = -a * ql;
	}
	else if(apsk16_input[0] == '1' && apsk16_input[1] == '0' && apsk16_input[2] == '1' && apsk16_input[3] == '0')
	{
		*apsk16_output_real = -a * pl;
		*apsk16_output_imaginary = a * ql;
	}
	else if(apsk16_input[0] == '1' && apsk16_input[1] == '0' && apsk16_input[2] == '1' && apsk16_input[3] == '1')
	{
		*apsk16_output_real = -a * pl;
		*apsk16_output_imaginary = -a * ql;
	}
	else if(apsk16_input[0] == '1' && apsk16_input[1] == '1' && apsk16_input[2] == '0' && apsk16_input[3] == '0')
	{
		*apsk16_output_real = kl / 2.0f;
		*apsk16_output_imaginary = kl / 2.0f;
	}
	else if(apsk16_input[0] == '1' && apsk16_input[1] == '1' && apsk16_input[2] == '0' && apsk16_input[3] == '1')
	{
		*apsk16_output_real = kl / 2.0f;
		*apsk16_output_imaginary = -kl / 2.0f;
	}
	else if(apsk16_input[0] == '1' && apsk16_input[1] == '1' && apsk16_input[2] == '1' && apsk16_input[3] == '0')
	{
		*apsk16_output_real = -kl / 2.0f;
		*apsk16_output_imaginary = kl / 2.0f;
	}
	else if(apsk16_input[0] == '1' && apsk16_input[1] == '1' && apsk16_input[2] == '1' && apsk16_input[3] == '1')
	{
		*apsk16_output_real = -kl / 2.0f;
		*apsk16_output_imaginary = -kl / 2.0f;
	}
	
}

void encode16(char *intrlvrOut, double *apsk16_output, int bit_num)
{
	int symbol_num = bit_num / 4;

	char *apsk16_input = (char*)malloc(4 * sizeof(char));
	zeros(apsk16_input, 4);

	for(int m=0;m<symbol_num*2;m+=2)
	{
		for(int c=0;c<4;c++)
		{
			apsk16_input[c] = intrlvrOut[4*(m/2)+c];
		}
		apsk16(apsk16_input, &apsk16_output[m], &apsk16_output[m+1]);
	}
	free(apsk16_input);
}


void de_apsk16(double *apsk16_output_real, double *apsk16_output_imaginary, char *de_apsk16_output4)
{
	double z1 = *apsk16_output_real;
	double z2 = *apsk16_output_imaginary;

	double kl = sqrt(2.0f);
	double pl = sin((double)(PI / 12.0f));
	double ql = cos((double)(PI / 12.0f));
	double jl = asin((double)(z2 / 3.0f));
	double hl = acos((double)(z2 / 3.0f));

#if 1
	double z = sqrt(pow(z1,2)+pow(z2,2));

	double v1 = PI/6.0f;
	double v2 = PI/3.0f;
	double v3 = PI/2.0f;

	if((z-2.0f)<=EPSINON)
	{
		if(z1>=EPSINON)
			z1 = kl/2.0f;
		else z1 = -kl/2.0f;

		if(z2>=EPSINON)
			z2 = kl/2.0f;
		else z2 = -kl/2.0f;
	}

	else 
	{
		if (z1>=EPSINON)
		{
			if (jl>=EPSINON && (jl-v1)<EPSINON) 
			{
				z1=3.0f*ql;
				z2=3.0f*pl;
			}
			else if ((jl-v1)>=EPSINON && (jl-v2)<EPSINON)
			{
				z1=3.0f/2.0f*kl;
				z2=3.0f/2.0f*kl;
			}
			else if ((jl-v2)>=EPSINON && (jl-v3)<EPSINON)
			{
				z1=3.0f*pl;
				z2=3.0f*ql;
			}
			else if ((jl+v1)>=EPSINON && jl<EPSINON)
			{
				z1=3*ql;
				z2=-3*pl;
			}
			else if ((jl+v2)>=EPSINON && (jl+v1)<EPSINON)
			{
				z1=3/2.0f*kl;
				z2=-3/2.0f*kl;
			}
			else if ((jl+v3)>=EPSINON && (jl-v2)<EPSINON)
			{
				z1=3.0f*pl;
				z2=-3.0f*ql;
			}
		}
		
		else 
		{
			if (hl>=EPSINON && (hl-v1)<EPSINON) 
			{
				z1=-3.0f*pl;
				z2=3.0f*ql;
			}
			else if ((hl-v1)<=EPSINON && (hl-v2)<EPSINON)
			{
				z1=-3/2.0f*kl;
				z2=3/2.0f*kl;
			}
			else if ((hl-v2)>=EPSINON && (hl-v3)<EPSINON)
			{
				z1=-3.0f*ql;
				z2=3.0f*pl;
			}
			else if ((hl-v3)>=EPSINON && (hl-2.0f*v2)<EPSINON)
			{
				z1=-3.0f*ql;
				z2=-3.0f*pl;
			}
			else if ((hl-2*v2)>=EPSINON && (hl-5.0f*v1)<EPSINON)
			{
				z1=-3/2.0f*kl;
				z2=-3/2.0f*kl;
			}
			else if ((hl-5.0f*v1)>=EPSINON && (hl-PI)<EPSINON)
			{
				z1=-3.0f*pl;
				z2=-3.0f*ql;
			}
		}
	}

#endif
	
	double *z_real = (double*)malloc(16 * sizeof(double));
	double *z_imaginary = (double*)malloc(16 * sizeof(double));
	zeros(z_real, 16);
	zeros(z_imaginary, 16);

	z_real[0] = 3.0f/2.0f*kl;
	z_imaginary[0] = 3.0f/2.0f*kl;

	z_real[1] = 3.0f/2.0f*kl;
	z_imaginary[1] = -3.0f/2.0f*kl;

	z_real[2] = -3.0f/2.0f*kl;
	z_imaginary[2] = 3.0f/2.0f*kl;

	z_real[3] = -3.0f/2.0f*kl;
	z_imaginary[3] = -3.0f/2.0f*kl;

	z_real[4] = 3.0f*ql;
	z_imaginary[4] = 3.0f*pl;

	z_real[5] = 3.0f*ql;
	z_imaginary[5] = -3.0f*pl;

	z_real[6] = -3.0f*ql;
	z_imaginary[6] = 3.0f*pl;

	z_real[7] = -3.0f*ql;
	z_imaginary[7] = -3.0f*pl;

	z_real[8] = 3.0f*pl;
	z_imaginary[8] = 3.0f*ql;

	z_real[9] = 3.0f*pl;
	z_imaginary[9] = -3.0f*ql;

	z_real[10] = -3.0f*pl;
	z_imaginary[10] = 3.0f*ql;

	z_real[11] = -3.0f*pl;
	z_imaginary[11] = -3.0f*ql;

	z_real[12] = kl/2.0f;
	z_imaginary[12] = kl/2.0f;

	z_real[13] = kl/2.0f;
	z_imaginary[13] = -kl/2.0f;

	z_real[14] = -kl/2.0f;
	z_imaginary[14] = kl/2.0f;

	z_real[15] = -kl/2.0f;
	z_imaginary[15] = -kl/2.0f;


	if (fabs(z1 - z_real[0]) <= EPSINON && fabs(z2 - z_imaginary[0]) <= EPSINON) 
	{
		de_apsk16_output4[0] = '0';
		de_apsk16_output4[1] = '0';
		de_apsk16_output4[2] = '0';
		de_apsk16_output4[3] = '0';
	}
	else if (fabs(z1-z_real[1]) <= EPSINON && fabs(z2-z_imaginary[1]) <= EPSINON) 
	{
		de_apsk16_output4[0] = '0';
		de_apsk16_output4[1] = '0';
		de_apsk16_output4[2] = '0';
		de_apsk16_output4[3] = '1';
	}
	else if (fabs(z1-z_real[2]) <= EPSINON && fabs(z2-z_imaginary[2]) <= EPSINON) 
	{
		de_apsk16_output4[0] = '0';
		de_apsk16_output4[1] = '0';
		de_apsk16_output4[2] = '1';
		de_apsk16_output4[3] = '0';
	}
	else if (fabs(z1-z_real[3]) <= EPSINON && fabs(z2-z_imaginary[3]) <= EPSINON) 
	{
		de_apsk16_output4[0] = '0';
		de_apsk16_output4[1] = '0';
		de_apsk16_output4[2] = '1';
		de_apsk16_output4[3] = '1';
	}
	else if (fabs(z1-z_real[4]) <= EPSINON && fabs(z2-z_imaginary[4]) <= EPSINON) 
	{
		de_apsk16_output4[0] = '0';
		de_apsk16_output4[1] = '1';
		de_apsk16_output4[2] = '0';
		de_apsk16_output4[3] = '0';
	}
	else if (fabs(z1-z_real[5]) <= EPSINON && fabs(z2-z_imaginary[5]) <= EPSINON) 
	{
		de_apsk16_output4[0] = '0';
		de_apsk16_output4[1] = '1';
		de_apsk16_output4[2] = '0';
		de_apsk16_output4[3] = '1';
	}
	else if (fabs(z1-z_real[6]) <= EPSINON && fabs(z2-z_imaginary[6]) <= EPSINON) 
	{
		de_apsk16_output4[0] = '0';
		de_apsk16_output4[1] = '1';
		de_apsk16_output4[2] = '1';
		de_apsk16_output4[3] = '0';
	}
	else if (fabs(z1-z_real[7]) <= EPSINON && fabs(z2-z_imaginary[7]) <= EPSINON) 
	{
		de_apsk16_output4[0] = '0';
		de_apsk16_output4[1] = '1';
		de_apsk16_output4[2] = '1';
		de_apsk16_output4[3] = '1';
	}
	else if (fabs(z1-z_real[8]) <= EPSINON && fabs(z2-z_imaginary[8]) <= EPSINON) 
	{
		de_apsk16_output4[0] = '1';
		de_apsk16_output4[1] = '0';
		de_apsk16_output4[2] = '0';
		de_apsk16_output4[3] = '0';
	}
	else if (fabs(z1-z_real[9]) <= EPSINON && fabs(z2-z_imaginary[9]) <= EPSINON) 
	{
		de_apsk16_output4[0] = '1';
		de_apsk16_output4[1] = '0';
		de_apsk16_output4[2] = '0';
		de_apsk16_output4[3] = '1';
	}
	else if (fabs(z1-z_real[10]) <= EPSINON && fabs(z2-z_imaginary[10]) <= EPSINON) 
	{
		de_apsk16_output4[0] = '1';
		de_apsk16_output4[1] = '0';
		de_apsk16_output4[2] = '1';
		de_apsk16_output4[3] = '0';
	}
	else if (fabs(z1-z_real[11]) <= EPSINON && fabs(z2-z_imaginary[11]) <= EPSINON) 
	{
		de_apsk16_output4[0] = '1';
		de_apsk16_output4[1] = '0';
		de_apsk16_output4[2] = '1';
		de_apsk16_output4[3] = '1';
	}
	else if (fabs(z1-z_real[12]) <= EPSINON && fabs(z2-z_imaginary[12]) <= EPSINON) 
	{
		de_apsk16_output4[0] = '1';
		de_apsk16_output4[1] = '1';
		de_apsk16_output4[2] = '0';
		de_apsk16_output4[3] = '0';
	}
	else if (fabs(z1-z_real[13]) <= EPSINON && fabs(z2-z_imaginary[13]) <= EPSINON) 
	{
		de_apsk16_output4[0] = '1';
		de_apsk16_output4[1] = '1';
		de_apsk16_output4[2] = '0';
		de_apsk16_output4[3] = '1';
	}
	else if (fabs(z1-z_real[14]) <= EPSINON && fabs(z2-z_imaginary[14]) <= EPSINON) 
	{
		de_apsk16_output4[0] = '1';
		de_apsk16_output4[1] = '1';
		de_apsk16_output4[2] = '1';
		de_apsk16_output4[3] = '0';
	}
	else if (fabs(z1-z_real[15]) <= EPSINON && fabs(z2-z_imaginary[15]) <= EPSINON) 
	{
		de_apsk16_output4[0] = '1';
		de_apsk16_output4[1] = '1';
		de_apsk16_output4[2] = '1';
		de_apsk16_output4[3] = '1';
	}
	
	free(z_real);
	free(z_imaginary);

}

void decode16(double *apsk16_output, char *de_apsk16_output, int bit_num)
{
	int symbol_num = bit_num / 4;

	char *de_apsk16_output4 = (char*)malloc(4 * sizeof(char));
	zeros(de_apsk16_output4, 4);

	for(int i=0;i<symbol_num*2;i+=2)
	{
		zeros(de_apsk16_output4, 4);

		de_apsk16(&apsk16_output[i], &apsk16_output[i+1], de_apsk16_output4);

		for(int j=0;j<4;j++)
		{
			de_apsk16_output[4*(i/2)+j] = de_apsk16_output4[j];
		}
	}

	free(de_apsk16_output4);
}


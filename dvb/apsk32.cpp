#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "apsk.h"

void zeros(char *input,int length)
{
	int i;
	for(i=0;i<length;i++)
	{
		input[i] = '0';
	}
}

void zeros(double *input,int length)
{
	int i;
	for(i=0;i<length;i++)
	{
		input[i] = 0.0f;
	}
}

void apsk32(char *apsk32_input, double *apsk32_output_real, double *apsk32_output_imaginary)
{
	double kl = sqrt(2.0f);
	double pl = sin((double)(PI / 12.0f));
	double ql = cos((double)(PI / 12.0f));
	double hl = sin((double)(PI / 8.0f));
	double jl = cos((double)(PI / 8.0f));

	double a = 2.84f;
	double b = 5.27f;

	if(apsk32_input[0] == '0' && apsk32_input[1] == '0' && apsk32_input[2] == '0' && apsk32_input[3] == '0' && apsk32_input[4] == '0')
	{
		*apsk32_output_real = a * kl / 2.0f;
		*apsk32_output_imaginary = a * kl / 2.0f;
	}
	else if(apsk32_input[0] == '0' && apsk32_input[1] == '0' && apsk32_input[2] == '0' && apsk32_input[3] == '0' && apsk32_input[4] == '1')
	{
		*apsk32_output_real = a * pl;
		*apsk32_output_imaginary = a * ql;
	}
	else if(apsk32_input[0] == '0' && apsk32_input[1] == '0' && apsk32_input[2] == '0' && apsk32_input[3] == '1' && apsk32_input[4] == '0')
	{
		*apsk32_output_real = a * kl / 2.0f;
		*apsk32_output_imaginary = -a * kl / 2.0f;
	}

	else if(apsk32_input[0] == '0' && apsk32_input[1] == '0' && apsk32_input[2] == '0' && apsk32_input[3] == '1' && apsk32_input[4] == '1')
	{
		*apsk32_output_real = a * pl;
		*apsk32_output_imaginary = -a * ql;
	}
	else if(apsk32_input[0] == '0' && apsk32_input[1] == '0' && apsk32_input[2] == '1' && apsk32_input[3] == '0' && apsk32_input[4] == '0')
	{
		*apsk32_output_real = -a * kl / 2.0f;
		*apsk32_output_imaginary = a * kl / 2.0f;
	}
	else if(apsk32_input[0] == '0' && apsk32_input[1] == '0' && apsk32_input[2] == '1' && apsk32_input[3] == '0' && apsk32_input[4] == '1')
	{
		*apsk32_output_real = -a * pl;
		*apsk32_output_imaginary = a * ql;
	}
	else if(apsk32_input[0] == '0' && apsk32_input[1] == '0' && apsk32_input[2] == '1' && apsk32_input[3] == '1' && apsk32_input[4] == '0')
	{
		*apsk32_output_real = -a * kl / 2.0f;
		*apsk32_output_imaginary = -a * kl / 2.0f;
	}
	else if(apsk32_input[0] == '0' && apsk32_input[1] == '0' && apsk32_input[2] == '1' && apsk32_input[3] == '1' && apsk32_input[4] == '1')
	{
		*apsk32_output_real = -a * pl;
		*apsk32_output_imaginary = -a * ql;
	}
	else if(apsk32_input[0] == '0' && apsk32_input[1] == '1' && apsk32_input[2] == '0' && apsk32_input[3] == '0' && apsk32_input[4] == '0')
	{
		*apsk32_output_real = b * jl;
		*apsk32_output_imaginary = b * hl;
	}
	else if(apsk32_input[0] == '0' && apsk32_input[1] == '1' && apsk32_input[2] == '0' && apsk32_input[3] == '0' && apsk32_input[4] == '1')
	{
		*apsk32_output_real = b * hl;
		*apsk32_output_imaginary = b * jl;
	}
	else if(apsk32_input[0] == '0' && apsk32_input[1] == '1' && apsk32_input[2] == '0' && apsk32_input[3] == '1' && apsk32_input[4] == '0')
	{
		*apsk32_output_real = b * kl / 2.0f;
		*apsk32_output_imaginary = -b * kl / 2.0f;
	}
	else if(apsk32_input[0] == '0' && apsk32_input[1] == '1' && apsk32_input[2] == '0' && apsk32_input[3] == '1' && apsk32_input[4] == '1')
	{
		*apsk32_output_real = 0.0f;
		*apsk32_output_imaginary = -b;
	}
	else if(apsk32_input[0] == '0' && apsk32_input[1] == '1' && apsk32_input[2] == '1' && apsk32_input[3] == '0' && apsk32_input[4] == '0')
	{
		*apsk32_output_real = -b * kl / 2.0f;
		*apsk32_output_imaginary = b * kl / 2.0f;
	}
	else if(apsk32_input[0] == '0' && apsk32_input[1] == '1' && apsk32_input[2] == '1' && apsk32_input[3] == '0' && apsk32_input[4] == '1')
	{
		*apsk32_output_real = 0.0f;
		*apsk32_output_imaginary = b;
	}
	else if(apsk32_input[0] == '0' && apsk32_input[1] == '1' && apsk32_input[2] == '1' && apsk32_input[3] == '1' && apsk32_input[4] == '0')
	{
		*apsk32_output_real = -b * jl;
		*apsk32_output_imaginary = -b * hl;
	}
	else if(apsk32_input[0] == '0' && apsk32_input[1] == '1' && apsk32_input[2] == '1' && apsk32_input[3] == '1' && apsk32_input[4] == '1')
	{
		*apsk32_output_real = -b * hl;
		*apsk32_output_imaginary = -b * jl;
	}
	else if(apsk32_input[0] == '1' && apsk32_input[1] == '0' && apsk32_input[2] == '0' && apsk32_input[3] == '0' && apsk32_input[4] == '0')
	{
		*apsk32_output_real = a * ql;
		*apsk32_output_imaginary = a * pl;
	}
	else if(apsk32_input[0] == '1' && apsk32_input[1] == '0' && apsk32_input[2] == '0' && apsk32_input[3] == '0' && apsk32_input[4] == '1')
	{
		*apsk32_output_real = kl / 2.0f;
		*apsk32_output_imaginary = kl / 2.0f;
	}
	else if(apsk32_input[0] == '1' && apsk32_input[1] == '0' && apsk32_input[2] == '0' && apsk32_input[3] == '1' && apsk32_input[4] == '0')
	{
		*apsk32_output_real = a * ql;
		*apsk32_output_imaginary = -a * pl;
	}
	else if(apsk32_input[0] == '1' && apsk32_input[1] == '0' && apsk32_input[2] == '0' && apsk32_input[3] == '1' && apsk32_input[4] == '1')
	{
		*apsk32_output_real = kl / 2.0f;
		*apsk32_output_imaginary = -kl / 2.0f;
	}
	else if(apsk32_input[0] == '1' && apsk32_input[1] == '0' && apsk32_input[2] == '1' && apsk32_input[3] == '0' && apsk32_input[4] == '0')
	{
		*apsk32_output_real = -a * ql;
		*apsk32_output_imaginary = a * pl;
	}
	else if(apsk32_input[0] == '1' && apsk32_input[1] == '0' && apsk32_input[2] == '1' && apsk32_input[3] == '0' && apsk32_input[4] == '1')
	{
		*apsk32_output_real = -kl / 2.0f;
		*apsk32_output_imaginary = kl / 2.0f;
	}
	else if(apsk32_input[0] == '1' && apsk32_input[1] == '0' && apsk32_input[2] == '1' && apsk32_input[3] == '1' && apsk32_input[4] == '0')
	{
		*apsk32_output_real = -a * ql;
		*apsk32_output_imaginary = -a * pl;
	}
	else if(apsk32_input[0] == '1' && apsk32_input[1] == '0' && apsk32_input[2] == '1' && apsk32_input[3] == '1' && apsk32_input[4] == '1')
	{
		*apsk32_output_real = -kl / 2.0f;
		*apsk32_output_imaginary = -kl / 2.0f;
	}
	else if(apsk32_input[0] == '1' && apsk32_input[1] == '1' && apsk32_input[2] == '0' && apsk32_input[3] == '0' && apsk32_input[4] == '0')
	{
		*apsk32_output_real = b;
		*apsk32_output_imaginary = 0.0f;
	}
	else if(apsk32_input[0] == '1' && apsk32_input[1] == '1' && apsk32_input[2] == '0' && apsk32_input[3] == '0' && apsk32_input[4] == '1')
	{
		*apsk32_output_real = b * kl / 2.0f;
		*apsk32_output_imaginary = b * kl / 2.0f;
	}
	else if(apsk32_input[0] == '1' && apsk32_input[1] == '1' && apsk32_input[2] == '0' && apsk32_input[3] == '1' && apsk32_input[4] == '0')
	{
		*apsk32_output_real = b * jl;
		*apsk32_output_imaginary = -b * hl;
	}
	else if(apsk32_input[0] == '1' && apsk32_input[1] == '1' && apsk32_input[2] == '0' && apsk32_input[3] == '1' && apsk32_input[4] == '1')
	{
		*apsk32_output_real = b * hl;
		*apsk32_output_imaginary = -b * jl;
	}
	else if(apsk32_input[0] == '1' && apsk32_input[1] == '1' && apsk32_input[2] == '1' && apsk32_input[3] == '0' && apsk32_input[4] == '0')
	{
		*apsk32_output_real = -b * jl;
		*apsk32_output_imaginary = b * hl;
	}
	else if(apsk32_input[0] == '1' && apsk32_input[1] == '1' && apsk32_input[2] == '1' && apsk32_input[3] == '0' && apsk32_input[4] == '1')
	{
		*apsk32_output_real = -b * hl;
		*apsk32_output_imaginary = b * jl;
	}
	else if(apsk32_input[0] == '1' && apsk32_input[1] == '1' && apsk32_input[2] == '1' && apsk32_input[3] == '1' && apsk32_input[4] == '0')
	{
		*apsk32_output_real = -b;
		*apsk32_output_imaginary = 0.0f;
	}
	else if(apsk32_input[0] == '1' && apsk32_input[1] == '1' && apsk32_input[2] == '1' && apsk32_input[3] == '1' && apsk32_input[4] == '1')
	{
		*apsk32_output_real = -b * kl / 2.0f;
		*apsk32_output_imaginary = -b * kl / 2.0f;
	}
}

void encode32(char *intrlvrOut, double *apsk32_output, int bit_num)
{
	int symbol_num = bit_num / 5;

	char *apsk32_input = (char*)malloc(5 * sizeof(char));
	zeros(apsk32_input, 5);

	for(int m=0;m<symbol_num*2;m+=2)
	{
		for(int c=0;c<5;c++)
		{
			apsk32_input[c] = intrlvrOut[5*(m/2)+c];
		}
		apsk32(apsk32_input, &apsk32_output[m], &apsk32_output[m+1]);
	}
	free(apsk32_input);
}


void de_apsk32(double *apsk32_output_real, double *apsk32_output_imaginary, char *de_apsk32_output5)
{
	double kl = sqrt(2.0f);
	double pl = sin((double)(PI / 12.0f));
	double ql = cos((double)(PI / 12.0f));
	double hl = sin((double)(PI / 8.0f));
	double jl = cos((double)(PI / 8.0f));

	double a = 2.84f;
	double b = 5.27f;

	double x_real = *apsk32_output_real;
	double x_imaginary = *apsk32_output_imaginary;

	double *z_real = (double*)malloc(32 * sizeof(double));
	double *z_imaginary = (double*)malloc(32 * sizeof(double));
	zeros(z_real, 32);
	zeros(z_imaginary, 32);

	double *dl = (double*)malloc(32 * sizeof(double));
	zeros(dl, 32);

	z_real[0] = a*kl/2.0f;
	z_imaginary[0] = a*kl/2.0f;

	z_real[1] = a*pl;
	z_imaginary[1] = a*ql;

	z_real[2] = a*kl/2.0f;
	z_imaginary[2] = -a*kl/2.0f;

	z_real[3] = a*pl;
	z_imaginary[3] = -a*ql;

	z_real[4] = -a*kl/2.0f;
	z_imaginary[4] = a*kl/2.0f;

	z_real[5] = -a*pl;
	z_imaginary[5] = a*ql;

	z_real[6] = -a*kl/2.0f;
	z_imaginary[6] = -a*kl/2.0f;

	z_real[7] = -a*pl;
	z_imaginary[7] = -a*ql;

	z_real[8] = b*jl;
	z_imaginary[8] = b*hl;

	z_real[9] = b*hl;
	z_imaginary[9] = b*jl;

	z_real[10] = b*kl/2.0f;
	z_imaginary[10] = -b*kl/2.0f;

	z_real[11] = 0.0f;
	z_imaginary[11] = -b;

	z_real[12] = -b*kl/2.0f;
	z_imaginary[12] = b*kl/2.0f;

	z_real[13] = 0.0f;
	z_imaginary[13] =b;

	z_real[14] = -b*jl;
	z_imaginary[14] = -b*hl;

	z_real[15] = -b*hl;
	z_imaginary[15] = -b*jl;

	z_real[16] = a*ql;
	z_imaginary[16] = a*pl;

	z_real[17] = kl/2.0f;
	z_imaginary[17] = kl/2.0f;

	z_real[18] = a*ql;
	z_imaginary[18] = -a*pl;

	z_real[19] = kl/2.0f;
	z_imaginary[19] = -kl/2.0f;

	z_real[20] = -a*ql;
	z_imaginary[20] = a*pl;

	z_real[21] = -kl/2.0f;
	z_imaginary[21] = kl/2.0f;

	z_real[22] = -a*ql;
	z_imaginary[22] = -a*pl;

	z_real[23] = -kl/2.0f;
	z_imaginary[23] = -kl/2.0f;

	z_real[24] = b;
	z_imaginary[24] = 0.0f;

	z_real[25] = b*kl/2.0f;
	z_imaginary[25] = b*kl/2.0f;

	z_real[26] = b*jl;
	z_imaginary[26] = -b*hl;

	z_real[27] = b*hl;
	z_imaginary[27] = -b*jl;

	z_real[28] = -b*jl;
	z_imaginary[28] = b*hl;

	z_real[29] = -b*hl;
	z_imaginary[29] = b*jl;

	z_real[30] = -b;
	z_imaginary[30] = 0.0f;

	z_real[31] = -b*kl/2.0f;
	z_imaginary[31] = -b*kl/2.0f;

#if 0
	for(int nl=0;nl<32;nl++)
	{
		dl[nl] = sqrt(pow((x_real-z_real[nl]),2) + pow((x_imaginary-z_imaginary[nl]),2));
	}

	int dl1= dl[0];
	int nl1 = 0;
	for(int nl=0;nl<32;nl++)
	{
		if(dl[nl] < dl1)
		{
			dl1 = dl[nl];
			nl1 = nl;
		}
	}

	x_real = z_real[nl1];
	x_imaginary = z_imaginary[nl1];
#endif

	if (fabs(x_real - z_real[0]) <= EPSINON && fabs(x_imaginary - z_imaginary[0]) <= EPSINON) 
	{
		de_apsk32_output5[0] = '0';
		de_apsk32_output5[1] = '0';
		de_apsk32_output5[2] = '0';
		de_apsk32_output5[3] = '0';
		de_apsk32_output5[4] = '0';
	}
	else if (fabs(x_real-z_real[1]) <= EPSINON && fabs(x_imaginary-z_imaginary[1]) <= EPSINON) 
	{
		de_apsk32_output5[0] = '0';
		de_apsk32_output5[1] = '0';
		de_apsk32_output5[2] = '0';
		de_apsk32_output5[3] = '0';
		de_apsk32_output5[4] = '1';
	}
	else if (fabs(x_real-z_real[2]) <= EPSINON && fabs(x_imaginary-z_imaginary[2]) <= EPSINON) 
	{
		de_apsk32_output5[0] = '0';
		de_apsk32_output5[1] = '0';
		de_apsk32_output5[2] = '0';
		de_apsk32_output5[3] = '1';
		de_apsk32_output5[4] = '0';
	}
	else if (fabs(x_real-z_real[3]) <= EPSINON && fabs(x_imaginary-z_imaginary[3]) <= EPSINON) 
	{
		de_apsk32_output5[0] = '0';
		de_apsk32_output5[1] = '0';
		de_apsk32_output5[2] = '0';
		de_apsk32_output5[3] = '1';
		de_apsk32_output5[4] = '1';
	}
	else if (fabs(x_real-z_real[4]) <= EPSINON && fabs(x_imaginary-z_imaginary[4]) <= EPSINON) 
	{
		de_apsk32_output5[0] = '0';
		de_apsk32_output5[1] = '0';
		de_apsk32_output5[2] = '1';
		de_apsk32_output5[3] = '0';
		de_apsk32_output5[4] = '0';
	}
	else if (fabs(x_real-z_real[5]) <= EPSINON && fabs(x_imaginary-z_imaginary[5]) <= EPSINON) 
	{
		de_apsk32_output5[0] = '0';
		de_apsk32_output5[1] = '0';
		de_apsk32_output5[2] = '1';
		de_apsk32_output5[3] = '0';
		de_apsk32_output5[4] = '1';
	}
	else if (fabs(x_real-z_real[6]) <= EPSINON && fabs(x_imaginary-z_imaginary[6]) <= EPSINON) 
	{
		de_apsk32_output5[0] = '0';
		de_apsk32_output5[1] = '0';
		de_apsk32_output5[2] = '1';
		de_apsk32_output5[3] = '1';
		de_apsk32_output5[4] = '0';
	}
	else if (fabs(x_real-z_real[7]) <= EPSINON && fabs(x_imaginary-z_imaginary[7]) <= EPSINON) 
	{
		de_apsk32_output5[0] = '0';
		de_apsk32_output5[1] = '0';
		de_apsk32_output5[2] = '1';
		de_apsk32_output5[3] = '1';
		de_apsk32_output5[4] = '1';
	}
	else if (fabs(x_real-z_real[8]) <= EPSINON && fabs(x_imaginary-z_imaginary[8]) <= EPSINON) 
	{
		de_apsk32_output5[0] = '0';
		de_apsk32_output5[1] = '1';
		de_apsk32_output5[2] = '0';
		de_apsk32_output5[3] = '0';
		de_apsk32_output5[4] = '0';
	}
	else if (fabs(x_real-z_real[9]) <= EPSINON && fabs(x_imaginary-z_imaginary[9]) <= EPSINON) 
	{
		//de_apsk32_output5 = "01001";

		de_apsk32_output5[0] = '0';
		de_apsk32_output5[1] = '1';
		de_apsk32_output5[2] = '0';
		de_apsk32_output5[3] = '0';
		de_apsk32_output5[4] = '1';
	}
	else if (fabs(x_real-z_real[10]) <= EPSINON && fabs(x_imaginary-z_imaginary[10]) <= EPSINON) 
	{
		//de_apsk32_output5 = "01010";

		de_apsk32_output5[0] = '0';
		de_apsk32_output5[1] = '1';
		de_apsk32_output5[2] = '0';
		de_apsk32_output5[3] = '1';
		de_apsk32_output5[4] = '0';
	}
	else if (fabs(x_real-z_real[11]) <= EPSINON && fabs(x_imaginary-z_imaginary[11]) <= EPSINON) 
	{
		//de_apsk32_output5 = "01011";

		de_apsk32_output5[0] = '0';
		de_apsk32_output5[1] = '1';
		de_apsk32_output5[2] = '0';
		de_apsk32_output5[3] = '1';
		de_apsk32_output5[4] = '1';
	}
	else if (fabs(x_real-z_real[12]) <= EPSINON && fabs(x_imaginary-z_imaginary[12]) <= EPSINON) 
	{
		//de_apsk32_output5 = "01100";

		de_apsk32_output5[0] = '0';
		de_apsk32_output5[1] = '1';
		de_apsk32_output5[2] = '1';
		de_apsk32_output5[3] = '0';
		de_apsk32_output5[4] = '0';
	}
	else if (fabs(x_real-z_real[13]) <= EPSINON && fabs(x_imaginary-z_imaginary[13]) <= EPSINON) 
	{
		//de_apsk32_output5 = "01101";

		de_apsk32_output5[0] = '0';
		de_apsk32_output5[1] = '1';
		de_apsk32_output5[2] = '1';
		de_apsk32_output5[3] = '0';
		de_apsk32_output5[4] = '1';
	}
	else if (fabs(x_real-z_real[14]) <= EPSINON && fabs(x_imaginary-z_imaginary[14]) <= EPSINON) 
	{
		//de_apsk32_output5 = "01110";

		de_apsk32_output5[0] = '0';
		de_apsk32_output5[1] = '1';
		de_apsk32_output5[2] = '1';
		de_apsk32_output5[3] = '1';
		de_apsk32_output5[4] = '0';
	}
	else if (fabs(x_real-z_real[15]) <= EPSINON && fabs(x_imaginary-z_imaginary[15]) <= EPSINON) 
	{
		//de_apsk32_output5 = "01111";

		de_apsk32_output5[0] = '0';
		de_apsk32_output5[1] = '1';
		de_apsk32_output5[2] = '1';
		de_apsk32_output5[3] = '1';
		de_apsk32_output5[4] = '1';
	}
	else if (fabs(x_real-z_real[16]) <= EPSINON && fabs(x_imaginary-z_imaginary[16]) <= EPSINON) 
	{
		//de_apsk32_output5 = "10000";

		de_apsk32_output5[0] = '1';
		de_apsk32_output5[1] = '0';
		de_apsk32_output5[2] = '0';
		de_apsk32_output5[3] = '0';
		de_apsk32_output5[4] = '0';
	}
	else if (fabs(x_real-z_real[17]) <= EPSINON && fabs(x_imaginary-z_imaginary[17]) <= EPSINON) 
	{
		//de_apsk32_output5 = "10001";

		de_apsk32_output5[0] = '1';
		de_apsk32_output5[1] = '0';
		de_apsk32_output5[2] = '0';
		de_apsk32_output5[3] = '0';
		de_apsk32_output5[4] = '1';
	}
	else if (fabs(x_real-z_real[18]) <= EPSINON && fabs(x_imaginary-z_imaginary[18]) <= EPSINON) 
	{
		//de_apsk32_output5 = "10010";

		de_apsk32_output5[0] = '1';
		de_apsk32_output5[1] = '0';
		de_apsk32_output5[2] = '0';
		de_apsk32_output5[3] = '1';
		de_apsk32_output5[4] = '0';
	}
	else if (fabs(x_real-z_real[19]) <= EPSINON && fabs(x_imaginary-z_imaginary[19]) <= EPSINON) 
	{
		//de_apsk32_output5 = "10011";

		de_apsk32_output5[0] = '1';
		de_apsk32_output5[1] = '0';
		de_apsk32_output5[2] = '0';
		de_apsk32_output5[3] = '1';
		de_apsk32_output5[4] = '1';
	}
	else if (fabs(x_real-z_real[20]) <= EPSINON && fabs(x_imaginary-z_imaginary[20]) <= EPSINON) 
	{
		//de_apsk32_output5 = "10100";

		de_apsk32_output5[0] = '1';
		de_apsk32_output5[1] = '0';
		de_apsk32_output5[2] = '1';
		de_apsk32_output5[3] = '0';
		de_apsk32_output5[4] = '0';
	}
	else if (fabs(x_real-z_real[21]) <= EPSINON && fabs(x_imaginary-z_imaginary[21]) <= EPSINON) 
	{
		//de_apsk32_output5 = "10101";

		de_apsk32_output5[0] = '1';
		de_apsk32_output5[1] = '0';
		de_apsk32_output5[2] = '1';
		de_apsk32_output5[3] = '0';
		de_apsk32_output5[4] = '1';
	}
	else if (fabs(x_real-z_real[22]) <= EPSINON && fabs(x_imaginary-z_imaginary[22]) <= EPSINON) 
	{
		//de_apsk32_output5 = "10110";

		de_apsk32_output5[0] = '1';
		de_apsk32_output5[1] = '0';
		de_apsk32_output5[2] = '1';
		de_apsk32_output5[3] = '1';
		de_apsk32_output5[4] = '0';
	}
	else if (fabs(x_real-z_real[23]) <= EPSINON && fabs(x_imaginary-z_imaginary[23]) <= EPSINON) 
	{
		//de_apsk32_output5 = "10111";

		de_apsk32_output5[0] = '1';
		de_apsk32_output5[1] = '0';
		de_apsk32_output5[2] = '1';
		de_apsk32_output5[3] = '1';
		de_apsk32_output5[4] = '1';
	}
	else if (fabs(x_real-z_real[24]) <= EPSINON && fabs(x_imaginary-z_imaginary[24]) <= EPSINON) 
	{
		//de_apsk32_output5 = "11000";

		de_apsk32_output5[0] = '1';
		de_apsk32_output5[1] = '1';
		de_apsk32_output5[2] = '0';
		de_apsk32_output5[3] = '0';
		de_apsk32_output5[4] = '0';
	}
	else if (fabs(x_real-z_real[25]) <= EPSINON && fabs(x_imaginary-z_imaginary[25]) <= EPSINON) 
	{
		//de_apsk32_output5 = "11001";

		de_apsk32_output5[0] = '1';
		de_apsk32_output5[1] = '1';
		de_apsk32_output5[2] = '0';
		de_apsk32_output5[3] = '0';
		de_apsk32_output5[4] = '1';
	}
	else if (fabs(x_real-z_real[26]) <= EPSINON && fabs(x_imaginary-z_imaginary[26]) <= EPSINON) 
	{
		//de_apsk32_output5 = "11010";

		de_apsk32_output5[0] = '1';
		de_apsk32_output5[1] = '1';
		de_apsk32_output5[2] = '0';
		de_apsk32_output5[3] = '1';
		de_apsk32_output5[4] = '0';
	}
	else if (fabs(x_real-z_real[27]) <= EPSINON && fabs(x_imaginary-z_imaginary[27]) <= EPSINON) 
	{
		//de_apsk32_output5 = "11011";

		de_apsk32_output5[0] = '1';
		de_apsk32_output5[1] = '1';
		de_apsk32_output5[2] = '0';
		de_apsk32_output5[3] = '1';
		de_apsk32_output5[4] = '1';
	}
	else if (fabs(x_real-z_real[28]) <= EPSINON && fabs(x_imaginary-z_imaginary[28]) <= EPSINON) 
	{
		//de_apsk32_output5 = "11100";

		de_apsk32_output5[0] = '1';
		de_apsk32_output5[1] = '1';
		de_apsk32_output5[2] = '1';
		de_apsk32_output5[3] = '0';
		de_apsk32_output5[4] = '0';
	}
	else if (fabs(x_real-z_real[29]) <= EPSINON && fabs(x_imaginary-z_imaginary[29]) <= EPSINON) 
	{
		//de_apsk32_output5 = "11101";

		de_apsk32_output5[0] = '1';
		de_apsk32_output5[1] = '1';
		de_apsk32_output5[2] = '1';
		de_apsk32_output5[3] = '0';
		de_apsk32_output5[4] = '1';
	}
	else if (fabs(x_real-z_real[30]) <= EPSINON && fabs(x_imaginary-z_imaginary[30]) <= EPSINON) 
	{
		//de_apsk32_output5 = "11110";

		de_apsk32_output5[0] = '1';
		de_apsk32_output5[1] = '1';
		de_apsk32_output5[2] = '1';
		de_apsk32_output5[3] = '1';
		de_apsk32_output5[4] = '0';
	}
	else if (fabs(x_real-z_real[31]) <= EPSINON && fabs(x_imaginary-z_imaginary[31]) <= EPSINON) 
	{
		//de_apsk32_output5 = "11111";

		de_apsk32_output5[0] = '1';
		de_apsk32_output5[1] = '1';
		de_apsk32_output5[2] = '1';
		de_apsk32_output5[3] = '1';
		de_apsk32_output5[4] = '1';
	}

	free(z_real);
	free(z_imaginary);

	free(dl);

}

void decode32(double *apsk32_output, char *de_apsk32_output, int bit_num)
{
	int symbol_num = bit_num / 5;

	char *de_apsk32_output5 = (char*)malloc(5 * sizeof(char));
	zeros(de_apsk32_output5, 5);

	for(int i=0;i<symbol_num*2;i+=2)
	{
		zeros(de_apsk32_output5, 5);

		de_apsk32(&apsk32_output[i], &apsk32_output[i+1], de_apsk32_output5);

		for(int j=0;j<5;j++)
		{
			de_apsk32_output[5*(i/2)+j] = de_apsk32_output5[j];
		}
	}

	free(de_apsk32_output5);
}


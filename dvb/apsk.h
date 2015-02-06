#ifndef PI
#define PI 3.1415926
#endif

#ifndef EPSINON
#define EPSINON 0.0000001
#endif

void zeros(char *input,int length);

void zeros(double *input,int length);

void apsk32(char *apsk32_input, double *apsk32_output_real, double *apsk32_output_imaginary);

void encode32(char *intrlvrOut, double *apsk32_output, int bit_num);

void de_apsk32(double *apsk32_output_real, double *apsk32_output_imaginary, char *de_apsk32_output5);

void decode32(double *apsk32_output, char *de_apsk32_output, int bit_num);

void apsk16(char *apsk16_input, double *apsk16_output_real, double *apsk16_output_imaginary);

void encode16(char *intrlvrOut, double *apsk16_output,int bit_num);

void de_apsk16(double *apsk16_output_real, double *apsk16_output_imaginary, char *de_apsk16_output4);

void decode16(double *apsk16_output, char *de_apsk16_output, int bit_num);
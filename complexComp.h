#ifndef COMPLEXCOMP_H
#define COMPLEXCOMP_H

#define pi 3.14159265358

#include <stdio.h>
#include <math.h>

typedef struct complex{
    double real;
    double img;
    double mag;
    double arg; // radians
} complex;

#define iota(I) complex I = {.real = 0.0, .img = 1.0, .mag = 1.0, .arg = pi/2}

// Initialize the complex number
complex c_set(double x, double y, char i){
    complex z;

    if(i=='c'){
        z.real = x;
        z.img = y;
        z.mag = sqrt(x*x + y*y);

        double angle = atan2(y, x);  // Calculate initial angle
        // Adjust angle to the principal range (-π, π]
        if (angle < 0) {
            angle += 2 * 3.14159;
        }
        z.arg = angle;
    }
    if(i=='p'){
        z.mag = x;
        z.arg = y;
        z.real = x*cos(y);
        z.img = x*sin(y);
    }
    return z;
}
void c_to_gen(complex *x){
    complex z = *x;
    *x = c_set(z.real,z.img,'c');
}

// Get the real part of an complex number
double c_real(complex z){
    return z.real;
}

// Get the imaginary part of an complex number
double c_img(complex z){
    return z.img;
}

// Magnitude of the complex number
double c_mag(complex z){
    c_to_gen(&z);
    return z.mag;
}

// Argument of the complex number
double c_arg(complex z){
    c_to_gen(&z);
    return z.arg;
}

// Conjugate of the complex number
complex c_conj(complex z){
    complex z1 = c_set(z.real,-z.img,'c');
    return z1;
}

// Sum
complex c_sum(complex x, complex y){
    double real,img;
    real = x.real+y.real;
    img = x.img+y.img;
    complex z = c_set(real,img,'c');
    return z;
}

// Difference
complex c_diff(complex x, complex y){
    double real,img;
    real = x.real-y.real;
    img = x.img-y.img;
    complex z = c_set(real,img,'c');
    return z;
}

// Product
complex c_prod(complex x, complex y){
    double real,img;
    real = x.real*y.real - x.img*y.img;
    img = x.real*y.img + x.img*y.real;
    complex z = c_set(real,img,'c');
    return z;
}

// Division
complex c_div(complex x, complex y){
    complex z,result;
    z = c_prod(x,c_conj(y));
    result = c_set(z.real/(c_mag(y)*c_mag(y)),z.img/(c_mag(y)*c_mag(y)),'c');
    return result;
}

// Power of the complex number
complex c_pow(complex base, complex exponent){
    // Calculate the magnitude and argument of the base
    double base_mag = c_mag(base);
    double base_arg = c_arg(base);

    // Calculate the exponent's effect on magnitude and argument
    double result_mag = pow(base_mag, exponent.real) * exp(-exponent.img * base_arg);
    double result_arg = (exponent.real * base_arg) + (exponent.img * log(base_mag));

    complex result = c_set(result_mag,result_arg,'p');
    return result;
}

// Exponent of the complex number
complex c_exp(complex z){

    double result_real = exp(z.real) * cos(z.img);
    double result_img = exp(z.real) * sin(z.img);

    complex result = c_set(result_real,result_img,'c');
    return result;
}

//Compute complex sine
complex c_sin(complex z){

    double result_real = sin(z.real)*cosh(z.img);
    double result_img = cos(z.real)*sinh(z.img);

    complex result = c_set(result_real,result_img,'c');
    return result;
}

//Compute complex cosine
complex c_cos(complex z){

    double result_real = cos(z.real)*cosh(z.img);
    double result_img = -sin(z.real)*sinh(z.img);

    complex result = c_set(result_real,result_img,'c');
    return result;
}

//Compute complex tangent
complex c_tan(complex z){

    complex result = c_div(c_sin(z),c_cos(z));
    return result;
}

// Print the complex number
void c_print(complex z, char i){
    if(i=='c'){
        if(z.img<0){
            printf("%.3f - j%.3f",z.real, -z.img);
        }
        else{
            printf("%.3f + j%.3f", z.real, z.img);
        }
    }
    else if(i=='p'){
        printf("%.3f rad%.3f",z.mag,z.arg);
    }
}

#endif
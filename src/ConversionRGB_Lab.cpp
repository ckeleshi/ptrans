#include <math.h>
#include <iostream>
#include <vector>

std::vector<double> RGB_to_lab(std::vector<double> rgb){
    int unsigned i=0;
    for (double value : rgb){
        value /= 255;
        value = (value > 0.04045)? pow((value+0.055)/1.055,2.4) :   value /= 12.92;
        rgb[i] = value*100;
        i++;
    }

    std::vector<double> XYZ{0, 0, 0};
  /*  std::cout <<  "RGB_to_lab"
    << rgb[0] << std::endl
    << rgb[1] << std::endl
    << rgb[2] << std::endl;*/

    XYZ[0] = rgb[0]*0.4124 + rgb[1]*0.3576 + rgb[2] * 0.1805;
    XYZ[1] = rgb[0]*0.2126 + rgb[1]*0.7152 + rgb[2] * 0.0722;
    XYZ[2] = rgb[0]*0.0193 + rgb[1]*0.1192 + rgb[2] * 0.9505;

    XYZ[0] /= 95.047;
    XYZ[1] /= 100.0;
    XYZ[2] /= 108.883;


    i = 0;
    for(double value : XYZ){
        value = (value >0.008856)?  pow(value,0.33333333333) : (7.787 * value) + (16./116.);
        XYZ[i] = value;
        i++;
    }

    std::vector<double> LAB{0,0,0};
    LAB[0] = (116 * XYZ[1]) - 16.;
    LAB[1] = 500 * (XYZ[0]-XYZ[1]);
    LAB[2] = 200 * (XYZ[1] - XYZ[2]);

    return LAB;

}
std::vector<double> LAB_to_RGB(std::vector<double> LAB){

  std::vector<double> XYZ {0,0,0};
  XYZ[1] = (LAB[0]+16.)/116.;
  XYZ[0] = LAB[1]/500. + XYZ[1];
  XYZ[2] = XYZ[1] - LAB[2]/200.;

  double epsi = pow(0.008856,3);
  int unsigned i =0;
  for(double value : XYZ){
     value = (value < epsi) ? pow(value,3) : (value-16./116.)/7.787;
     XYZ[i] = value;
     i++;
  }

  XYZ[0] *= 95.047;
  XYZ[1] *= 100.;
  XYZ[2] *= 108.883;

  double r = 3.2406*XYZ[0] + (-1.5372*XYZ[1]) + (-0.4986*XYZ[2]);
  double g = -0.9689*XYZ[0] + 1.8757*XYZ[1] + 0.04152*XYZ[2];
  double b = 0.05571*XYZ[0] + (-0.2040*XYZ[1]) + 1.05699*XYZ[2];
  std::vector<double> RGB{r,g,b};
  /*std::cout << "LAB_to_RGB"
  << RGB[0] << std::endl
  << RGB[1] << std::endl
  << RGB[2] << std::endl;*/

  i =0;
  epsi = 0.0031308;
  for(double value : RGB){
    value = (value < epsi)? pow(value,(1/2.4))/1.055 - 0.055 : value*12.92;
    value*=255;
    value /=100;
    RGB[i] = roundf(value);
    i++;
  }

   return RGB;
}

int main(int argc, char const *argv[])
{
  double r,g,b;
  std::cout<<"Choisissez composantes"<<std::endl<<"R : ";
  std::cin>>r;
  std::cout<<"G : ";
  std::cin>>g;
  std::cout<<"B : ";
  std::cin>>b;
    std::vector<double> RGB{r,g,b};
    std::vector<double> LAB = RGB_to_lab(RGB);
    std::cout<<"conversion RGB to lab"<<std::endl
             <<LAB[0]<<std::endl
             <<LAB[1]<<std::endl
             <<LAB[2]<<std::endl;
    std::cout<<"===================================="<<std::endl;
    std::vector<double> temp = LAB_to_RGB(LAB);
    std::cout<<"conversion lab to RGB"<<std::endl
              << temp[0]<< std::endl
              << temp[1] << std::endl
              << temp[2] << std::endl;
    return 0;
}

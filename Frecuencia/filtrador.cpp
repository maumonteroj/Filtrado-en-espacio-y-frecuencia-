//File Management
#include <iostream>
#include <fstream>
#include <string>

//Time management
#include <ltiTimer.h>

//Images
#include <ltilib/ltiIOImage.h>  // To read/write images from files (jpg, png, and bmp)
#include <ltiImage.h>
#include <ltilib/ltiViewer2D.h> //Display image

//Kernels
#include <ltiKernel2D.h>
#include <ltiGaussKernels.h>//Gausian Kernel
#include <ltiOctagonalKernel.h>//Octagonal

//Fourier
#include <ltiFFT.h>//Fast fourier transform
#include <ltiIFFT.h>//Inverse Fast fourier transform
#include <ltiBoundaryExpansion.h>//Padding

//Matrix handling
#include <ltiMatrix.h>//Allows for element-wise multiplication

//Math
#include <ltiMath.h>
#include <ltiRound.h>

//namespace po = boost::program_options;
using namespace std;
using namespace lti;

int main(int ac, char* av[])
{
    
    viewer2D view;
    ioImage loader; //An image loader
    image img;    //An image container
    
    int numberOfImages = 12;
    int numberOfKernels = 14;    
    
    string fileNames[numberOfImages];
    fileNames[0] = "64x48.png";
    fileNames[1] = "230x140.png";
    fileNames[2] = "400x230.jpg";
    fileNames[3] = "570x330.png";
    fileNames[4] = "740x420.png";
    fileNames[5] = "920x520.jpg";
    fileNames[6] = "1080x610.png";
    fileNames[7] = "1245x705.png";
    fileNames[8] = "1400x800.png";
    fileNames[9] = "1600x900.png";
    fileNames[10] = "1750x1000.png";
    fileNames[11] = "1920x1080.png";
    
    
    kernel2D<float> kern [numberOfKernels];
    int kernelSizes [numberOfKernels];
    kernelSizes[0] = 3;
    kernelSizes[1] = 5;
    kernelSizes[2] = 8;
    kernelSizes[3] = 12;
    kernelSizes[4] = 19;
    kernelSizes[5] = 29;
    kernelSizes[6] = 45;
    kernelSizes[7] = 70;
    kernelSizes[8] = 109;
    kernelSizes[9] = 171;
    kernelSizes[10] = 267;
    kernelSizes[11] = 418;
    kernelSizes[12] = 654;
    kernelSizes[13] = 1023;
    
    int variance;
    
    for(int i = 0; i < 14; i++){
        variance = (kernelSizes[i]+2)*(kernelSizes[i]+2)/36;
        gaussKernel2D<float> gauss(kernelSizes[i], variance);
        kern[i].castFrom(gauss);
    }
    
    /*
    Padding
    */    
    
    //Sets parameters
    eBoundaryType padding = Zero;
    boundaryExpansion::parameters bepar;
    bepar.boundaryType = padding;
    bepar.topBorder = 0;
    bepar.leftBorder = 0;
    
    boundaryExpansion::parameters kernelParm;
    kernelParm.boundaryType = padding;
    
    //    
    channel chnlA, chnlB;
    
    //FFT parameters
    channel img_F_R, img_F_I, kern_F_R, kern_F_I, productA, productB, result_R, result_I, filtered;
    //Fast fourier Trasform method
    fft fft2D;
    ifft ifft2D;

    double times [numberOfImages*numberOfKernels];
    
    int time = 500000;//500ms+-20us
    double n;//Holds image size
    int iterations;
    matrix<float> mat;
    kernel2D<float> kernel;
    

    int top;
    int side;

    for(int i = 0; i < numberOfImages; i++){//Iterates over number of images
        cout << i << endl;
        loader.load(fileNames[i], img);
        for(int j = 0; j < numberOfKernels; j++){//Iterates over number of kernelsZ
            cout << "\t" << j << endl;
            kernel = kern[j];
            
            chnlB.castFrom(kernel);
            iterations = 0;
            timer chron;
            chron.start();
            
            while(chron.getTime() < time){
            
                n = iround(pow(2, ceil(log(2*lti::max((img.rows() + kernel.rows()),
                                                           (img.columns() + kernel.columns())))/
                                          log(2.0f))));
                
                top = (n-kernel.columns())/2;
                side = (n-kernel.rows())/2;
                         
                kernelParm.topBorder = top;
                kernelParm.leftBorder = side;
                kernelParm.bottomBorder = top;
                kernelParm.rightBorder = side;
                
                bepar.bottomBorder = n-img.rows();
                bepar.rightBorder = n-img.columns();

                //Creates Padder
                boundaryExpansion be(bepar);
                boundaryExpansion beKern(kernelParm);


                //Applies Padding
                chnlA.castFrom(img);
                be.apply(chnlA);
                beKern.apply(chnlB);
               
               
                fft2D.apply(chnlA, img_F_R, img_F_I);
                fft2D.apply(chnlB, kern_F_R, kern_F_I);
                
                //Applies the element by element matrix multiplication
                productA = mat.emultiply(img_F_R, kern_F_R);
                productB = mat.emultiply(img_F_I, kern_F_I);
                result_R = mat.subtract(productA, productB);
                
                productA = mat.emultiply(img_F_I, kern_F_R);
                productB = mat.emultiply(img_F_R, kern_F_I);
                result_I = productA + productB;

                ifft2D.apply(result_R, result_I, filtered);
                
                
                //view.show(filtered);
                //view.waitKeyPressed();

                
                //Inverts padding
                channel chnlOut(img.rows(), img.columns());
                chnlOut.fill(filtered);
                
                iterations++;
            }
            chron.stop();
            
            times[i*numberOfImages + j] = chron.getTime() / iterations;
        }
    } 
    
    

    //Saves time parameters
    ofstream myfile;
    myfile.open ("time.txt");
    
    myfile << "# name: time\n# type: matrix\n# rows: "<< numberOfKernels <<"\n# columns: "<<  numberOfImages << "\n";
    
    for(int i = 0; i < numberOfImages; i++){
        for(int j = 0; j < numberOfKernels; j++){
            myfile << times[i*numberOfImages + j] << " ";
        }
        myfile << endl;
    } 
    myfile.close();
    
   
        
    //Saves gauss kernel size
    myfile.open ("kernel.txt");
    
    myfile << "# name: kernel\n# type: matrix\n# rows: "<< 1 <<"\n# columns: "<<  numberOfKernels << "\n";
    
    for(int i = 0; i < numberOfKernels; i++){
        myfile << kernelSizes[i]*kernelSizes[i] << " ";
    } 
    myfile.close();

    //Saves filesSize
    myfile.open ("imageSize.txt");
    myfile << "# name: imageSize\n# type: matrix\n# rows: "<< 1 <<"\n# columns: "<<  numberOfImages << "\n";
    
    myfile << (64*64)/2 << " ";//64x64
    myfile << (230*140)/2 << " ";//230x140
    myfile << (400*230)/2 << " ";//400x230
    myfile << (570*330)/2 << " ";//570x330
    myfile << (740*420)/2 << " ";//740x420
    myfile << (920*520)/2 << " ";//920x520
    myfile << (1080*610)/2 << " ";//1080x610
    myfile << (1245*705)/2 << " ";//1245x705
    myfile << (1400*800)/2 << " ";//1400x800
    myfile << (1600*900)/2 << " ";//1600x900
    myfile << (1750*1000)/2 << " ";//1750x1000
    myfile << (1920*1080)/2 << " ";//1920x1080
    
    
    myfile.close();    

}









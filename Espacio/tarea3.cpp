//Manejo del documento
#include <iostream>
#include <fstream>
#include <string>

//Manejo del tiempo
#include <ltiTimer.h>

//Imagenes
#include <ltiIOImage.h>
#include <ltiViewer2D.h>
#include <ltiImage.h>
#include <ltiChannel.h>

//Kenels
#include <ltiGaussKernels.h>
#include <ltiKernel2D.h>
#include <ltiIOLTI.h>
#include <cstdlib>
#include <ltiOctagonalKernel.h>

//Convolucion
#include <ltiConvolution.h>

#include <ltiGuiServer.h>

using namespace std;
using namespace lti;

int main(int argc,char* argv[]) {
	
	timer chron;
	ioImage loader; //Cargar la imagen
	channel data,resultado; //Contenedor de la imagen en blanco y negro

    string fileNames[12];
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
   

    int kernelSizes [14];
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
    

	int numberOfImages = 12;
    int numberOfKernels = 14;    
    double times [numberOfImages*numberOfKernels];
    double timesoc [numberOfImages*numberOfKernels];

	int time = 500000;//500ms+-20us
    int iterations;

 for(int i = 0; i < numberOfImages; i++){//Iterates over number of images
	loader.load(fileNames[i],data);
	cout << i << endl;
	 for(int j = 0; j < numberOfKernels; j++){//Iterates over number of kernelsZ
		gaussKernel2D<channel::value_type> kernel(kernelSizes[j],((kernelSizes[j]+2)*(kernelSizes[j]+2)/36));
		octagonalKernel<channel::value_type > kerneloc(kernelSizes[j]);
		iterations = 0;
		chron.start();
		while(chron.getTime() < time){
			convolution filter;                        // convolution operator
			convolution::parameters param;             // parameters
			param.setKernel(kernel);                        // use the gauss kernel
			filter.setParameters(param); 
			filter.apply(data);
			iterations++;
		}
		chron.stop();
		times[i*numberOfImages + j] = chron.getTime() / iterations;
		
		iterations = 0;
		chron.start();
		while(chron.getTime() < time){
			convolution filter;                        // convolution operator
			convolution::parameters param;             // parameters
			param.setKernel(kerneloc);                        // use the gauss kernel
			filter.setParameters(param); 
			filter.apply(data);
			iterations++;
		}
		chron.stop();
		timesoc[i*numberOfImages + j] = chron.getTime() / iterations;
		
	}
 }
  
  //Guardar datos
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
    
//Guardar datos octogonal
   myfile.open ("timeoc.txt");
    
    myfile << "# name: time\n# type: matrix\n# rows: "<< numberOfKernels <<"\n# columns: "<<  numberOfImages << "\n";
    
    for(int i = 0; i < numberOfImages; i++){
        for(int j = 0; j < numberOfKernels; j++){
            myfile << timesoc[i*numberOfImages + j] << " ";
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

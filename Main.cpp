#include "MatrixMult.h"
//#include <chrono>
int main(int argc, char** argv) //D,d,DoseLim,SearchLim,RefPath,EvaPath
{
        if (argc!=7) 
        {
                std::cout << "Number of input arguments incorrect!\n";
                return 1;
        }

        double D=std::strtod(argv[1],nullptr); //Dose criterion, in %
        double d=std::strtod(argv[2],nullptr); //Distance criterion, in mm
        double DoseLim=std::strtod(argv[3],nullptr); //Lower dose threshold, in %
        double SearchLim=std::strtod(argv[4],nullptr); //Search distance, in mm
        Matrix ReadMatrix(std::string); //forward declare
        Matrix Gamma3D(const Matrix&, const Matrix&, const double, const double, const double, const double);
        Matrix Ref=ReadMatrix(argv[5]); //Reference dose distribution
        Matrix Eva=ReadMatrix(argv[6]); //Evaluated dose distribution
        //Eva.PrintMatrix(1);
       //auto start = std::chrono::high_resolution_clock::now();
        Matrix Gamma=Gamma3D(Ref,Eva,D,d,DoseLim,SearchLim);
        //auto finish = std::chrono::high_resolution_clock::now();
        //std::chrono::duration<double> elapsed = finish - start;
        Gamma.PrintMatrix("GammaDistribution.dat");
        //std::cout << "Elapsed time: " << elapsed.count() << " s\n";
        return 0;
}

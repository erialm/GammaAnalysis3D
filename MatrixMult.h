#include <iostream>
#include <fstream>
#include <string>
class Matrix
{
        private:
                size_t Rows, Columns, Slices;
                double* Mat;
                double XCorner, YCorner, ZCorner;
                double SizeX, SizeY, SizeZ;

        public:
                Matrix(size_t Rows, size_t Columns, size_t Slices=1, double Val=0); //Constructor
                Matrix(const Matrix&) = delete; //No copy constructor
                Matrix (Matrix&&);    //Move constructor

                ~Matrix() {if (Mat!=nullptr) delete[] Mat;}      //Destructor
                
                size_t GetRows() const noexcept {return Rows;} 
                size_t GetColumns() const noexcept {return Columns;}
                size_t GetSlices() const noexcept {return Slices;}
                double GetXPosition(size_t Column) const noexcept {return XCorner+Column*SizeX;}
                double GetYPosition(size_t Row) const noexcept {return YCorner+Row*SizeY;}
                double GetZPosition(size_t Slice) const noexcept {return ZCorner+Slice*SizeZ;}
                double GetSizeX() const noexcept {return SizeX;}
                double GetSizeY() const noexcept {return SizeY;}
                double GetSizeZ() const noexcept {return SizeZ;}

                void SetXCorner(double X) noexcept {XCorner = X;}
                void SetYCorner(double Y) noexcept {YCorner = Y;}
                void SetZCorner(double Z) noexcept {ZCorner = Z;}
                void SetSizeX(double X) noexcept {SizeX = X;}
                void SetSizeY(double Y) noexcept {SizeY = Y;}
                void SetSizeZ(double Z) noexcept {SizeZ = Z;}
                 

                double  operator() (size_t Row, size_t Column, size_t Slice) const {return Mat[Slice+Slices*(Column+Columns*Row)];}  
                double  operator() (size_t Index) const {return Mat[Index];}  
                double& operator() (size_t Row, size_t Column, size_t Slice) {return Mat[Slice+Slices*(Column+Columns*Row)];}       

                Matrix& operator=(const Matrix&) = delete; //No copy assignment
                Matrix operator*(const Matrix&); //Matrix multiplication
                Matrix& operator=(Matrix&&);

                void PrintMatrix(size_t) const noexcept;
                void PrintMatrix(std::string) const noexcept;
                Matrix Transpose() const noexcept;
                void SwapRows(size_t,size_t) noexcept;
                void GaussianElimination(Matrix&) noexcept;
                Matrix BackSubstitute(Matrix&) noexcept;
};




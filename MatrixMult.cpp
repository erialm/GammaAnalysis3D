#include "MatrixMult.h"
Matrix::Matrix(size_t NoRows, size_t NoColumns, size_t NoSlices, double Val)
        :Rows{NoRows}, Columns{NoColumns}, Slices{NoSlices}, Mat{new double[NoRows*NoColumns*NoSlices]}
{
        for (size_t i=0;i<Rows*Columns*Slices;++i) Mat[i]=Val;
}


Matrix::Matrix(Matrix&& Rhs)
{
        Rows=Rhs.Rows;
        Columns=Rhs.Columns;
        Slices=Rhs.Slices;
        Mat=Rhs.Mat;
        Rhs.Mat=nullptr;
        Rhs.Rows=0;
        Rhs.Columns=0;
        Rhs.Slices=0;
}

Matrix& Matrix::operator=(Matrix&& Rhs)
{
        if (this !=&Rhs)
        {
                delete[] Mat;
                Mat = Rhs.Mat;
                Rows=Rhs.Rows;
                Columns=Rhs.Columns;
                Slices=Rhs.Slices;

                Rhs.Mat=nullptr;
                Rhs.Rows=0;
                Rhs.Columns=0;
                Rhs.Slices=0;
        }
        return *this;
}

void Matrix::PrintMatrix(size_t Slice) const noexcept
{
        using std::cout;
        for (size_t i=0;i<Rows;++i)
        {
                for (size_t j=0; j<Columns; ++j) cout << "  " << (*this)(i,j,Slice);
                cout << '\n';
        }
        cout << '\n';
}

Matrix Matrix::operator*(const Matrix& Rhs)
{
        using std::cout;
        if (Rows!=Rhs.GetColumns() && Columns!=Rhs.GetRows())
        {
                cout << "Dimensions of the matrices are incorrect!\n";
                exit (1);
        }
        Matrix ReturnMatrix(Rows,Rhs.GetColumns());
        for (size_t i=0;i<ReturnMatrix.GetRows();++i)
        {
                for (size_t j=0;j<ReturnMatrix.GetColumns();++j) 
                {
                        for (size_t k=0;k<Columns;++k) ReturnMatrix(i,j,0)+=(*this)(i,k,0)*Rhs(k,j,0);
                }
        }
        return ReturnMatrix;
}

Matrix Matrix::Transpose() const noexcept
{
        Matrix TempMatrix(Columns,Rows);
        for (size_t i=0;i<Columns;i++)
        {
                for (size_t j=0;j<Rows;++j) TempMatrix(i,j,0)=(*this)(j,i,0);
        }
        return TempMatrix;
}

void Matrix::SwapRows(size_t Row1, size_t Row2) noexcept
{
       if (Row1==Row2) return;
       double TempElement;
       for (size_t i=0;i<Columns;++i)
       {
               TempElement=(*this)(Row1,i,0); 
               (*this)(Row1,i,0)=(*this)(Row2,i,0);
               (*this)(Row2,i,0)=TempElement;
       }       
}

void Matrix::GaussianElimination(Matrix& B) noexcept
{
        double MaxEl;
        size_t MaxRow;
        size_t BColumns=B.GetColumns();
        for (size_t i=0;i<Rows;++i)
        {
                MaxEl=(*this)(i,i,0);
                MaxRow=i;
                for (size_t j=i+1;j<Rows;++j)
                {
                        if ((*this)(j,i,0)>MaxEl)
                        {
                                MaxEl=(*this)(j,i,0);
                                MaxRow=j;
                        }
                }
                SwapRows(i,MaxRow);
                B.SwapRows(i,MaxRow);
                double FirstElement;
                for (size_t j=i+1;j<Rows;++j)
                {
                        for (size_t k=0;k<BColumns;++k) B(j,k,0)=B(j,k,0)-B(i,k,0)*(*this)(j,i,0)/MaxEl;
                        FirstElement=(*this)(j,i,0);
                        for (size_t k=0;k<Columns;++k)
                        {
        
                                (*this)(j,k,0)=(*this)(j,k,0)-(*this)(i,k,0)*FirstElement/MaxEl;
                        }
                }       
        }
}

Matrix Matrix::BackSubstitute(Matrix& B) noexcept
{
        Matrix x(B.GetRows(),B.GetColumns());

        for (int i=(B.GetRows()-1);i>=0;--i) 
        {
                for (size_t j=0;j<B.GetColumns();++j)
                {
                        x(i,j,0)=B(i,j,0)/(*this)(i,i,0);
                        for (int k=i-1;k>=0;--k)
                        {
                                B(k,j,0)=B(k,j,0)-x(i,j,0)*(*this)(k,i,0);
                        }
                }
        }
        return x;
}

/*Matrix Matrix::MLDivide(Matrix& B) noexcept
{
       (*this).GaussianElimination(B);
       Matrix x = (*this).BackSubstitute(B);
       return x;
}*/



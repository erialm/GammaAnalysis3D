#include "MatrixMult.h"
#include <cmath>
#include <cfloat>
#include <vector>
#include <algorithm>
#include <iterator>
#include <limits>
#define DIM 4
struct Point
{
        double X, Y, Z, Dose;
};


Matrix ReadMatrix(std::string Path) noexcept
{
        //Change Geant4 output to write corner coordinates along with number of voxels before dose data
        std::ifstream InputFile{Path,std::ios::binary};
        size_t Dim[DIM-1];
        double Coor[DIM-1];
        double VoxSize[DIM-1];
        double InputD;
        int InputI;
        for (size_t i=0;i<(DIM-1);++i)
        {
                InputFile.read(reinterpret_cast<char*>(&InputI),sizeof(int));
                InputFile.read(reinterpret_cast<char*>(&InputD),sizeof(double));
                Dim[i]=InputI;
                Coor[i]=InputD;
                InputFile.read(reinterpret_cast<char*>(&InputD),sizeof(double));
                VoxSize[i]=InputD;
        }
        Matrix Image(Dim[0],Dim[1],Dim[2]);
        Image.SetXCorner(Coor[0]);
        Image.SetYCorner(Coor[1]);
        Image.SetZCorner(Coor[2]);
        Image.SetSizeX(VoxSize[0]);
        Image.SetSizeY(VoxSize[1]);
        Image.SetSizeZ(VoxSize[2]);
        for (size_t i=0;i<Dim[2];++i)
        {
                for (size_t j=0;j<Dim[1];++j)
                {
                        for (size_t k=0;k<Dim[0];++k)
                        {
                                InputFile.read(reinterpret_cast<char*>(&InputD),sizeof(double));
                                Image(j,k,i)=InputD;
                        }
                }
        }
        return Image;
}

double GetMax(const Matrix& M) noexcept
{
        size_t Rows=M.GetRows();
        size_t Columns=M.GetColumns();
        size_t Slices=M.GetSlices();
        double Max=0;
        for (size_t i=0;i<Rows;++i)
        {
                for (size_t j=0;j<Columns;++j)
                {
                        for (size_t k=0;k<Slices;++k) 
                        {
                                if (M(i,j,k)>Max) Max=M(i,j,k);
                        }
                }
        }
        return Max;
}

double GetMin(const Matrix& M) noexcept
{
        size_t Rows=M.GetRows();
        size_t Columns=M.GetColumns();
        size_t Slices=M.GetSlices();
        double Min=DBL_MAX;
        for (size_t i=0;i<Rows;++i)
        {
                for (size_t j=0;j<Columns;++j)
                {
                        for (size_t k=0;k<Slices;++k) 
                        {
                                if (std::isnan(M(i,j,k))) continue;
                                if (M(i,j,k)<Min) Min=M(i,j,k);
                        }
                }
        }
        return Min;
}


std::vector<Point> CreateCube(const Matrix& Eva, const size_t Row, const size_t Column,const size_t Slice, const double D, const double d) noexcept
{
        std::vector<Point> Cube;
        Point EntryPoint;

        //P1
        EntryPoint.X=Eva.GetXPosition(Column)/d;
        EntryPoint.Y=Eva.GetYPosition(Row)/d;
        EntryPoint.Z=Eva.GetZPosition(Slice)/d;
        EntryPoint.Dose=Eva(Row,Column,Slice)/D;
        Cube.push_back(EntryPoint);

        //P2
        EntryPoint.X=Eva.GetXPosition(Column+1)/d;
        EntryPoint.Y=Eva.GetYPosition(Row)/d;
        EntryPoint.Z=Eva.GetZPosition(Slice)/d;
        EntryPoint.Dose=Eva(Row,Column+1,Slice)/D;
        Cube.push_back(EntryPoint);

        //P3
        EntryPoint.X=Eva.GetXPosition(Column+1)/d;
        EntryPoint.Y=Eva.GetYPosition(Row+1)/d;
        EntryPoint.Z=Eva.GetZPosition(Slice)/d;
        EntryPoint.Dose=Eva(Row+1,Column+1,Slice)/D;
        Cube.push_back(EntryPoint);
       
        //P4
        EntryPoint.X=Eva.GetXPosition(Column)/d;
        EntryPoint.Y=Eva.GetYPosition(Row+1)/d;
        EntryPoint.Z=Eva.GetZPosition(Slice+1)/d;
        EntryPoint.Dose=Eva(Row+1,Column,Slice+1)/D;
        Cube.push_back(EntryPoint);

        //P5
        EntryPoint.X=Eva.GetXPosition(Column+1)/d;
        EntryPoint.Y=Eva.GetYPosition(Row)/d;
        EntryPoint.Z=Eva.GetZPosition(Slice+1)/d;
        EntryPoint.Dose=Eva(Row,Column+1,Slice+1)/D;
        Cube.push_back(EntryPoint);
        
        //P6
        EntryPoint.X=Eva.GetXPosition(Column)/d;
        EntryPoint.Y=Eva.GetYPosition(Row)/d;
        EntryPoint.Z=Eva.GetZPosition(Slice+1)/d;
        EntryPoint.Dose=Eva(Row,Column,Slice+1)/D;
        Cube.push_back(EntryPoint);
        
        //P7
        EntryPoint.X=Eva.GetXPosition(Column)/d;
        EntryPoint.Y=Eva.GetYPosition(Row+1)/d;
        EntryPoint.Z=Eva.GetZPosition(Slice)/d;
        EntryPoint.Dose=Eva(Row+1,Column,Slice)/D;
        Cube.push_back(EntryPoint);
        
        //P8
        EntryPoint.X=Eva.GetXPosition(Column+1)/d;
        EntryPoint.Y=Eva.GetYPosition(Row+1)/d;
        EntryPoint.Z=Eva.GetZPosition(Slice+1)/d;
        EntryPoint.Dose=Eva(Row+1,Column+1,Slice+1)/D;
        Cube.push_back(EntryPoint);
        
        return Cube;
}

Matrix SetUpP(std::vector<Point> SimplexVertices, const double* const R) noexcept
{
        size_t Size=SimplexVertices.size();//SimplexVertices.size();
        Matrix P(DIM,1);
        P(0,0,0)=R[0]-SimplexVertices[Size-1].X;
        P(1,0,0)=R[1]-SimplexVertices[Size-1].Y;
        P(2,0,0)=R[2]-SimplexVertices[Size-1].Z;
        P(3,0,0)=R[3]-SimplexVertices[Size-1].Dose;

        return P;
}

Matrix SetUpV(std::vector<Point> SimplexVertices) noexcept
{
        size_t Size=SimplexVertices.size();
        Matrix V(DIM,Size-1);
        for (size_t VCol=0;VCol<(Size-1);++VCol)
        {
                V(0,VCol,0)=SimplexVertices[VCol].X-SimplexVertices[Size-1].X;
                V(1,VCol,0)=SimplexVertices[VCol].Y-SimplexVertices[Size-1].Y;
                V(2,VCol,0)=SimplexVertices[VCol].Z-SimplexVertices[Size-1].Z;
                V(3,VCol,0)=SimplexVertices[VCol].Dose-SimplexVertices[Size-1].Dose;
        }
        return V;
}


Matrix MLDivide(Matrix& A,Matrix& B) noexcept
{
       A.GaussianElimination(B);
       Matrix x = A.BackSubstitute(B);
       return x;
}

Matrix GetW(Matrix& VTV, Matrix& VT, const Matrix& P) noexcept
{
        Matrix TmpW=(MLDivide(VTV,VT))*P;
        Matrix W(TmpW.GetRows()+1,1);
        double Sum=0;
        for (size_t i=0;i<TmpW.GetRows();++i) 
        {
                W(i,0,0)=TmpW(i,0,0);
                Sum+=TmpW(i,0,0);
        }
        W(TmpW.GetRows(),0,0)=1-Sum;
        return W;
}

Point ComputeNearestPoint(const Matrix& W, const std::vector<Point>& SimplexVertices) noexcept
{
        Point NearestPoint;
        for (size_t i=0;i<W.GetRows();++i) 
        {
                NearestPoint.X+=W(i,0,0)*SimplexVertices[i].X;
                NearestPoint.Y+=W(i,0,0)*SimplexVertices[i].Y;
                NearestPoint.Z+=W(i,0,0)*SimplexVertices[i].Z;
                NearestPoint.Dose+=W(i,0,0)*SimplexVertices[i].Dose;
        }
        return NearestPoint;
}

double GetGammaFromLine(const std::vector<Point>& Line, const Matrix& W, const double * const R) noexcept
{
        using std::pow;
        bool PositiveW=true;
        for (size_t i=0;i<W.GetRows();++i)
        {
                if (W(i,0,0)<0)
                {
                        PositiveW=false;
                        break;
                }
        }
        if (PositiveW) 
        {
                Point NP=ComputeNearestPoint(W,Line);
                return pow(R[0]-NP.X,2)+pow(R[1]-NP.Y,2)+pow(R[2]-NP.Z,2)+pow(R[3]-NP.Dose,2);
        }
        else 
        {
                Point P1, P2;
                P1.X=Line[0].X;
                P1.Y=Line[0].Y;
                P1.Z=Line[0].Z;
                P1.Dose=Line[0].Dose;
                P2.X=Line[1].X;
                P2.Y=Line[1].Y;
                P2.Z=Line[1].Z;
                P2.Dose=Line[1].Dose;
                double GAMMA1=pow(R[0]-P1.X,2)+pow(R[1]-P1.Y,2)+pow(R[2]-P1.Z,2)+pow(R[3]-P1.Dose,2);
                double GAMMA2=pow(R[0]-P2.X,2)+pow(R[1]-P2.Y,2)+pow(R[3]-P2.Z,2)+pow(R[3]-P2.Dose,2);
                return GAMMA1<GAMMA2 ? GAMMA1 : GAMMA2;
        }

}

double CheckSimplex(const std::vector<Point>&,const double* const) noexcept;

double GetGammaFromSimplex(const std::vector<Point>& Simplex,const Matrix& W, const double * const R) noexcept
{

        bool PositiveW=true;
        for (size_t i=0;i<W.GetRows();++i)
        {
                if (W(i,0,0)<0)
                {
                        PositiveW=false;
                        break;
                }
        }
        double GAMMA;
        if (PositiveW) 
        {
                Point NP=ComputeNearestPoint(W,Simplex);
                GAMMA=pow(R[0]-NP.X,2)+pow(R[1]-NP.Y,2)+pow(R[2]-NP.Z,2)+pow(R[3]-NP.Dose,2);
        }
        else 
        {
                std::vector<Point> NonNegW;
                Point Entry;
                std::vector<double> Gamma;
                for (size_t i=0;i<W.GetRows();++i)
                {
                        if (W(i,0,0)<0)
                        {
                                for (size_t j=0;j<W.GetRows();++j)
                                {
                                        if (i!=j)
                                        {
                                                Entry.X=Simplex[j].X;
                                                Entry.Y=Simplex[j].Y;
                                                Entry.Z=Simplex[j].Z;
                                                Entry.Dose=Simplex[j].Dose;
                                                NonNegW.push_back(Entry);
                                        } 
                                }
                                GAMMA=CheckSimplex(NonNegW, R);
                                NonNegW.clear();
                        }
                } 
        }
        return GAMMA;
}

double CheckSimplex(const std::vector<Point>& Simplex, const double* const R) noexcept
{

        Matrix P=SetUpP(Simplex, R);
        Matrix V=SetUpV(Simplex);
        Matrix VT=V.Transpose();
        Matrix VTV=VT*V;
        Matrix W=GetW(VTV,VT,P);

        double GAMMA;
        if (Simplex.size()>2) GAMMA=GetGammaFromSimplex(Simplex,W,R);
        else GAMMA=GetGammaFromLine(Simplex,W,R);
        
        return GAMMA;
}

double Get3DGamma(const Matrix& Eva, const double D, const double d, const double* const R, const size_t Row, const size_t Column, const size_t Slice) noexcept
{
        using std::vector;
        vector<Point> Cube=CreateCube(Eva,Row,Column,Slice,D,d);
        vector<Point> Tetrahedon1, Tetrahedon2, Tetrahedon3, Tetrahedon4, Tetrahedon5;
        Matrix Gammas(5,1,1);
        
        Tetrahedon1.push_back(Cube[6]);
        Tetrahedon1.push_back(Cube[3]);
        Tetrahedon1.push_back(Cube[5]);
        Tetrahedon1.push_back(Cube[7]);

        Gammas(0,0,0)=CheckSimplex(Tetrahedon1,R);

        Tetrahedon2.push_back(Cube[6]);
        Tetrahedon2.push_back(Cube[7]);
        Tetrahedon2.push_back(Cube[5]);
        Tetrahedon2.push_back(Cube[1]);

        Gammas(1,0,0)=CheckSimplex(Tetrahedon2,R);


        Tetrahedon3.push_back(Cube[0]);
        Tetrahedon3.push_back(Cube[6]);
        Tetrahedon3.push_back(Cube[1]);
        Tetrahedon3.push_back(Cube[5]);

        Gammas(2,0,0)=CheckSimplex(Tetrahedon3,R);

        Tetrahedon4.push_back(Cube[5]);
        Tetrahedon4.push_back(Cube[1]);
        Tetrahedon4.push_back(Cube[4]);
        Tetrahedon4.push_back(Cube[7]);

        Gammas(3,0,0)=CheckSimplex(Tetrahedon4,R);

        Tetrahedon5.push_back(Cube[6]);
        Tetrahedon5.push_back(Cube[1]);
        Tetrahedon5.push_back(Cube[2]);
        Tetrahedon5.push_back(Cube[7]);

        Gammas(4,0,0)=CheckSimplex(Tetrahedon5,R);
        return GetMin(Gammas);
}


Matrix Gamma3D(const Matrix& Ref, const Matrix& Eva, double Dose, const double d, const double DoseLim, const double SearchLim) noexcept
{
        using std::sqrt;
        using std::vector;
        using std::abs;
        using std::min_element;
        using std::begin;
        using std::end;

        Matrix Gamma(Ref.GetRows(),Ref.GetColumns(),Ref.GetSlices(),std::numeric_limits<double>::quiet_NaN());
        Gamma.SetXCorner(Ref.GetXPosition(0));
        Gamma.SetYCorner(Ref.GetYPosition(0));
        Gamma.SetZCorner(Ref.GetZPosition(0));
        
        Gamma.SetSizeX(Ref.GetSizeX());
        Gamma.SetSizeY(Ref.GetSizeY());
        Gamma.SetSizeZ(Ref.GetSizeZ());
        vector<double> ReturnGamma;
        ReturnGamma.reserve(1000);

        const double MaxDose=GetMax(Ref);
        const double D=((100+Dose)/100)*MaxDose-MaxDose;
        const double DL=((100+DoseLim)/100)*MaxDose-MaxDose;
        const double SL=SearchLim/d;

        double* Rr = new double[DIM];
        double* Er = new double[DIM];
        size_t RefRows=Ref.GetRows();
        size_t RefColumns=Ref.GetColumns();
        size_t RefSlices=Ref.GetSlices();
        size_t EvaRows=Eva.GetRows();
        size_t EvaColumns=Eva.GetColumns();
        size_t EvaSlices=Eva.GetSlices();
        for (size_t RefRow=0; RefRow<RefRows;++RefRow)
        {
                Rr[1]=Ref.GetYPosition(RefRow)/d;
                for (size_t RefCol=0; RefCol<RefColumns;++RefCol)
                {
                        Rr[0]=Ref.GetXPosition(RefCol)/d;
                        for (size_t RefSli=0;RefSli<RefSlices;++RefSli)
                        {

                                if (Ref(RefRow,RefCol,RefSli)<DL) continue;
                                Rr[2]=Ref.GetZPosition(RefSli)/d;
                                Rr[3]=Ref(RefRow,RefCol,RefSli)/D;
                                for (size_t EvaRow=0; EvaRow<(EvaRows-1);++EvaRow) 
                                {
                                        Er[1]=Eva.GetYPosition(EvaRow)/d;
                                        if (abs(Rr[1]-Er[1])>SL) continue;
                                        for (size_t EvaCol=0;EvaCol<(EvaColumns-1);++EvaCol)
                                        {
                                                Er[0]=Eva.GetXPosition(EvaCol)/d;
                                                if (abs(Rr[0]-Er[0])>SL) continue;
                                                for (size_t EvaSli=0;EvaSli<(EvaSlices-1);++EvaSli)
                                                {
                                                        Er[2]=Eva.GetZPosition(EvaSli)/d;
                                                        if (abs(Rr[2]-Er[2])>SL) continue;
                                                        ReturnGamma.push_back(Get3DGamma(Eva,D,d,Rr,EvaRow,EvaCol,EvaSli));
                                                }
                                        }
                                }
                                Gamma(RefRow,RefCol,RefSli)=sqrt(*min_element(begin(ReturnGamma),end(ReturnGamma)));
                                ReturnGamma.clear();
                        }
                }
        }
        return Gamma;
}

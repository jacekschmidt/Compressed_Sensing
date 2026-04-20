#include <iostream>
#include "mvector.h"
#include "mmatrix.h"
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>

// for contour plots
void residualcontour(const std::string& filename, const MMatrix& A, const MVector& b, double x0, double x1, double y0, double y1, int nx, int ny)
{
    std::ofstream demofile;
    demofile.open(filename);
    if (!demofile)
    {
        std::cout<<"unable to open file";
        return;
    }
    for(int i=0; i<nx; i++)
    {
        for(int j=0; j<ny; j++)
        {
            MVector x={x0+(i*(x1-x0))/(nx-1),y0+(j*(y1-y0))/(ny-1)};
            demofile<<dot(b-(A*x),b-(A*x))<<" "; 
        }
        demofile<<"\n";
    }
    demofile.close();
}

// for 3dcontour plots
void residual3dcontour(const std::string& filename, const MMatrix& A, const MVector& b, double x0, double x1, double y0, double y1, double z0, double z1, int nx, int ny, int nz)
{
    std::ofstream demofile;
    demofile.open(filename);
    if (!demofile)
    {
        std::cout<<"unable to open file";
        return;
    }
    for(int i=0; i<nx; i++)
    {
        for(int j=0; j<ny; j++)
        {
            for(int k=0; k<nz; k++)
            {
                MVector x={x0+(i*(x1-x0))/(nx-1),y0+(j*(y1-y0))/(ny-1),z0+(k*(z1-z0))/(nz-1)};
                demofile<<dot(b-(A*x),b-(A*x))<<" ";
            } 
        }
    }
    demofile.close();
}

// SDLS
int SDLS(const MMatrix& A, const MVector& b, MVector& x,
int maxIterations, double tol)
{
    if (A.Cols()!=x.size() || A.Rows()!=b.size())
    {
        std::cout<<"incompatible arguments"<<std::endl;
        return 0;
    }
    double alpha;
    MVector r=Transpose(A)*(b-A*x);
    for(int iter=0; iter<=maxIterations; iter++)
    {
        std::cout<<x[1]<<" ";
        if (r.L2Norm() < tol) 
        {
            return iter;
        }
        MVector q=A*r;
        alpha=dot(r,r)/dot(q,q);
        x=x+alpha*r;
        r=r-alpha*(Transpose(A)*q);
    }
    return 0;
}


//IHT
int IHT(const MMatrix& A, const MVector& b, MVector&x, int k,
int maxIterations, double tol) 
{
    if (A.Cols()!=x.size() || A.Rows()!=b.size())
    {
        std::cout<<"incompatible arguments"<<std::endl;
        return 0;
    }
    double alpha;
    MVector r=Transpose(A)*(b-A*x);
    for(int iter=0; iter<=maxIterations; iter++)
    {
        MVector q=A*r;
        alpha=dot(r,r)/dot(q,q);
        x=x+alpha*r;
        x.Threshold(k);
        r=Transpose(A)*(b-A*x);
        if (r.L2Norm() < tol) 
        {
            return iter;
        }
    }
    return 0;
}

//NIHT
int NIHT(const MMatrix& A, const MVector& b, MVector&x, int k,
int maxIterations, double tol) 
{
    if (A.Cols()!=x.size() || A.Rows()!=b.size())
    {
        std::cout<<"incompatible arguments"<<std::endl;
        return 0;
    }
    double alpha;
    MVector r=Transpose(A)*(b-A*x);
    for(int iter=0; iter<=maxIterations; iter++)
    {
        r.Threshold(k);
        MVector q=A*r;
        alpha=dot(r,r)/dot(q,q);
        x=x+alpha*r;
        x.Threshold(k);
        r=Transpose(A)*(b-A*x);
        if (r.L2Norm() < tol) 
        {
            return iter;
        }
    }
    return 0;
}

// succesful recovery 
int successfulrecoveryIHT(const MMatrix& A, const MVector& x, MVector&x0, int k,
int maxIterations, double tol1,double tol2)
{
    MVector b=A*x;
    IHT(A, b, x0, k, maxIterations, tol1);
    if(dot(x-x0, x-x0)/dot(x,x)<tol2) return 1;
    else return 0;
}

// succesful recovery 
int successfulrecoveryNIHT(const MMatrix& A, const MVector& x, MVector&x0, int k,
int maxIterations, double tol1,double tol2)
{
    MVector b=A*x;
    NIHT(A, b, x0, k, maxIterations, tol1);
    if(dot(x-x0, x-x0)/dot(x,x)<tol2) return 1;
    else return 0;
}

/* phase transition */
void phasetransitionIHT(const std::string& filename, int maxn, int repitions,
int maxIterations, double tol1,double tol2 )
{
    std::ofstream demofile;
    demofile.open(filename);
    if (!demofile)
    {
        std::cout<<"unable to open file";
        return;
    }
    
    for(int m=1; m<=maxn; m++)
    {
        for(int k=1; k<=maxn; k++)
        {
            MMatrix A(m,maxn,0);
            MVector x0(maxn,0);
            MVector zeros(maxn,0);
            MVector x(maxn);
            double sum=0.0;
            for(int rep=0; rep<repitions; rep++)
            {
                A.initialise_normal();
                /*A.normalise_columns();*/
                x0=zeros;
                x.initialise_sparse_normal(k);
                sum+=successfulrecoveryIHT(A, x, x0, k, maxIterations, tol1, tol2);
            }
            demofile<<sum/repitions<<" ";
        }
    }
    demofile.close();
}

/* phase transition */
void phasetransitionNIHT(const std::string& filename, int maxn, int repitions,
int maxIterations, double tol1,double tol2 )
{
    std::ofstream demofile;
    demofile.open(filename);
    if (!demofile)
    {
        std::cout<<"unable to open file";
        return;
    }
    
    for(int m=1; m<=maxn; m++)
    {
        for(int k=1; k<=maxn; k++)
        {
            MMatrix A(m,maxn,0);
            MVector x0(maxn,0);
            MVector zeros(maxn,0);
            MVector x(maxn);
            double sum=0.0;
            for(int rep=0; rep<repitions; rep++)
            {
                A.initialise_normal();
                A.normalise_columns();
                x0=zeros;
                x.initialise_sparse_normal(k);
                sum+=successfulrecoveryNIHT(A, x, x0, k, maxIterations, tol1, tol2);
            }
            demofile<<sum/repitions<<" ";
        }
    }
    demofile.close();
}



int main()
{
    /*SLDS
    {
        MMatrix A(3,2,0.0);
        MVector x={0,0};
        A(0,0)=1; A(0,1)=2; A(1,0)=2; A(1,1)=1; A(2,0)=-1; A(2,1)=0;
        MVector b={10,-1,4};
        residualcontour("dat1.txt", A, b,-10, 10, -10, 10, 100, 100);
        b={10,-1,0};
        x={0,0};
        A(0,0)=1; A(0,1)=2; A(1,0)=2; A(1,1)=1; A(2,0)=-1; A(2,1)=0;
        residualcontour("dat2.txt", A, b,-10, 10, -10, 10, 100, 100);
        x={0,0};
        A(0,0)=1; A(0,1)=2; A(1,0)=2; A(1,1)=1; A(2,0)=1.8; A(2,1)=-2;
        residualcontour("dat3.txt", A, b,-10, 10, -10, 10, 100, 100);
        x={0,0};
        A(0,0)=1; A(0,1)=2; A(1,0)=2; A(1,1)=1; A(2,0)=-2; A(2,1)=-2;
        residualcontour("dat4.txt", A, b,-10, 10, -10, 10, 100, 100);
    }
    */

    /*Threshold
    {
        int k;
        MVector v={-9, 6, -7, 8, -5};

        k=3;

        
        MVector w(v.size(),0);
        auto itmin = std::min_element(v.begin(), v.end());
        int indexmin = std::distance(v.begin(), itmin);
        auto itmax = std::max_element(v.begin(), v.end());
        int indexmax = std::distance(v.begin(), itmax);

        for(int i=0;i<k;i++)
        {
            if(std::abs(*itmax)>std::abs(*itmin))
            {
                w[indexmax]=*itmax;
                v[indexmax]=0;
                itmax = std::max_element(v.begin(), v.end());
                indexmax = std::distance(v.begin(), itmax);
            }
            else
            {
                w[indexmin]=*itmin;
                v[indexmin]=0;
                itmin = std::min_element(v.begin(), v.end());
                indexmin = std::distance(v.begin(), itmin);
            }

            
        }
            
        v.Threshold(k);
        std::cout<<v<<std::endl;
    }
    */

    /*2x3 matrices
    {
        MMatrix A(2,3,0);
        A.initialise_normal();
        std::cout<<A<<std::endl;
        MVector x(3);              // creates a vector of size 10, initially all zeros
        x.initialise_sparse_normal(1); // fills 3 random entries with N(0,1) samples
        std::cout << x << "\n";
        MVector b=A*x;
        x={0,0,0};
        std::cout<<IHT(A,b,x,1,1000,1e-2)<<std::endl;
        std::cout<<x<<std::endl;
    }
    */

    /*projection illustration
    {
        MMatrix A(2,3,0);
        A(0,0)=1; A(0,1)=2;A(0,2)=3;A(1,0)=4;A(1,1)=5;A(1,2)=6;
        MVector x={0.5,0.2,1};
        MVector b=A*x;
        x={0.2,0,0};
        //std::cout<<IHT(A,b,x,1,1000,1)<<std::endl;
        //std::cout<<x<<std::endl;
        residual3dcontour("dat5.txt",A,b,0,1.5,0,1.5,0,1.5,100,100,100);
    }
    */

    /*projection illustration
    {
        MMatrix A(2,3,0);
        A(0,0)=0.6428; A(0,1)=1; A(0,2)=0.6428; A(1,0)=-0.7660; A(1,1)=0; A(1,2)=0.7660;
        MVector b={0.7428,0.7660};
        MVector x={0.6, 0.4, 0.3};
        std::cout<<IHT(A,b,x,1,100,0.15);
        residual3dcontour("dat7.txt",A,b,-0.5,1,-0.2,0.5,-0.2,1.2,100,100,100);
    }
    */

    /*successful recovery test
    {
        
        MMatrix A(2,3,0);
        MVector x0(3,0);
        MVector x(3);
        A.initialise_normal();
        std::cout<<A<<std::endl;
        x.initialise_sparse_normal(1);
        std::cout<<x<<std::endl;
        std::cout<<x0<<std::endl;
        std::cout<<successfulrecoveryIHT(A, x, x0, 1, 1000000, 0.1, 0.2)<<std::endl;
        std::cout<<x0<<std::endl;
    }
    */

    /*phase transition*/
    {
        phasetransitionIHT("phasetransition40_20IHT.txt", 40, 20, 100, 0.1, 0.2);
    }
    
    return 0;
}
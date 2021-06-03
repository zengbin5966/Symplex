#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

class Simplex{
    private:
        int rows, cols;
        //stores coefficients of all the variables
        std::vector <std::vector<double> > A;
        //stores constants of constraints
        std::vector<double> B;
        //stores the coefficients of the objective function
        std::vector<double> C;
        double maximum = 0;
        bool isUnbounded;

    public:
        Simplex(std::vector <std::vector<double> > matrix,std::vector<double> b ,std::vector<double> c){
            maximum = 0;
            isUnbounded = false;
            rows = matrix.size();
            cols = matrix[0].size();
            A.resize( rows , vector<double>( cols , 0 ) );
            B.resize(b.size());
            C.resize(c.size());
            for(int i= 0;i<rows;i++){             //pass A[][] values to the metrix
                for(int j= 0; j< cols;j++ ){
                    A[i][j] = matrix[i][j];
                }
            }
            for(int i=0; i< c.size() ;i++ ){      //pass c[] values to the B vector
                C[i] = c[i] ;
            }
            for(int i=0; i< b.size();i++ ){      //pass b[] values to the B vector
                B[i] = b[i];
            }
        }

    /**
     * @brief 单纯性法检验数函数，检验当前结果是否为最优解
     * 
     * @return true 
     * @return false 
     */
    bool checkOptimality(){
        bool isOptimality = true;
        for (int i = 0; i < C.size(); i++){
            if (C[i]>0){
                isOptimality = false;
            }
        }
        return isOptimality;
    }
    
    /**
     * @brief 系数矩阵初始化函数
     * 
     */
    void matrixInint(){
        if (checkOptimality()) return ;
        int index = 0;
        for (int i = 0; i < rows; i++){
            double value = A[i][index];
            if (value == 0) continue;
            for (int j = 0; j < cols; j++){
               A[i][j] = A[i][j]/value;
            } //将某一行化为第index个元素为0的行
            B[i] = B[i]/value;
            for (int i1 = 0; i1 < rows; i1++){
                if (i1 == i) continue;
                if (index >= rows) continue;
                double rate = A[i1][index]/ A[i][index]; //计算相减行变换时的系数
                for (int jl = 0; jl < cols; jl++){
                    A[i1][jl] -= rate*A[i][jl];}
                B[i1] -= rate*B[i];
            } //对其他行进行减操作
            double rate = C[index]/ A[i][index];
            for (int jl = 0; jl < cols; jl++){
                C[jl] -= rate*A[i][jl];
            }
            maximum -= rate*B[i];
            index += 1;
        }

    }
    
    /**
     * @brief put the (index_i)th variable into Base Variables with putting the (index_j)th variable out of Base Variables.
     * 
     * @param index_i 
     * @param row_index 
     * @return true 
     * @return false 
     */
    bool cVectorUtil(int index_i, int row_index){
        aPrint();
        double value = A[row_index][index_i];
        double rate;
        for (int j = 0; j < cols; j++){
            A[row_index][j] = A[row_index][j]/value;
        } 
        for (int i = 0; i < rows; i++){
            if (i==row_index) continue;
            rate = A[i][index_i]/value;
            for (int j = 0; j < cols; j++){
                A[i][j] -= rate*A[row_index][j];
            }
            B[i] -= rate*B[row_index];
        rate = C[index_i]/A[row_index][index_i];
        cout << rate << "," <<  C[index_i] << "," <<  A[i][index_i] << endl;
        for (int jl = 0; jl < cols; jl++){
            C[jl] -= rate*A[row_index][jl];
        }
        maximum -= rate*B[row_index];
        }
        
    }

    bool aPrint(){
        for (int i = 0; i < rows; i++){
            for (int j = 0; j < cols; j++){
                cout << A[i][j] << "\t";
            }
            cout << B[i] << "\t";
            cout << endl;
        }
        for (int i = 0; i < cols; i++){
            cout << C[i] << "\t";
        }
        cout << maximum << "\t";
        cout << endl << endl;
    }

    bool doOptimize(){
        int index_i;
        int row_index;
        while (checkOptimality()==false){
            double tmp_max = 0;
            for (int i = 0; i < cols; i++){
                if (C[i]>=tmp_max){
                    tmp_max = C[i];
                    index_i = i;
                }
            }
            tmp_max = 0;
            double tmp;
            for (int j = 0; j < rows; j++){
                if (B[j]/A[j][index_i]>tmp_max){
                    tmp_max = B[j]/A[j][index_i];
                    row_index = j;
                }
            }
            cVectorUtil(index_i,row_index);
            aPrint();
        }
    }

    double const getMaximum(){
        return maximum;
    }
};



int main()
{

    int colSizeA=6;  //should initialise columns size in A
    int rowSizeA = 3;  //should initialise columns row in A[][] vector

    double C[]= {6,5,4,0,0,0};  //should initialis the c arry here
    double B[]={180,300,240};  // should initialis the b array here



    double a[3][6] = {    //should intialis the A[][] array here
                   { 2,  1,  1, 1, 0, 0},
                   { 1,  3,  2, 0, 1, 0},
                   { 2,  1,  2, 0, 0, 1}
             };


        std::vector <std::vector<double> > vec2D(rowSizeA, std::vector<double>(colSizeA, 0));

        std::vector<double> b(rowSizeA,0);
        std::vector<double> c(colSizeA,0);

       for(int i=0;i<rowSizeA;i++){         //make a vector from given array
            for(int j=0; j<colSizeA;j++){
                vec2D[i][j] = a[i][j];
            }
       }

       for(int i=0;i<rowSizeA;i++){
            b[i] = B[i];
       }

        for(int i=0;i<colSizeA;i++){
            c[i] = C[i];
       }


      // hear the make the class parameters with A[m][n] vector b[] vector and c[] vector
      Simplex simplex(vec2D,b,c);
      simplex.matrixInint();
      simplex.doOptimize();
    return 0;
}